#!/usr/bin/env python3
"""Gene neighbourhood of a MGnifam's member proteins on their contigs (operon check).

Each member protein (a MGnify90 cluster representative) is placed on its longest contig, and the genes around it are
read in its own direction (offset -1 = directly upstream, +1 = directly downstream):
- directon window: the same-strand run around it, up to 5 genes on each side, stopped by a gene on the other strand.
  A side is complete when it holds 5 same-strand genes or the strand switch; a contig edge before that is incomplete.
- census window (synteny_census.py): up to 5 genes on each side, either strand, within 10 kb of the member gene.
  A side is complete when its 5th gene or 10 kb of contig lie on it.
A neighbour gene is named by its Pfam accessions, by "MGnifam_<family>" when its cluster representative is a member,
and by "cluster:<rep>" when it has neither. Every name is ranked by the member proteins whose window holds it (the
denominator is the proteins whose window is complete on both sides), with enrichment over all genes on the same
contigs. Each --target group (PFxxxxx, CLxxxx expanded with the clans file, "family") tests one hypothesis.

Usage: neighbourhood.py MSA --clusters-metadata mgy_clusters_metadata.parquet --metadata mgy_proteins_metadata.parquet
       --pfam mgy_proteins_pfam.parquet --clans Pfam-A.clans.tsv.gz --out partners.tsv.gz --proteins proteins.tsv.gz
stdout: counts per step with timings, the targets, the top partners. Needs pyarrow (no numpy). Read-only.
"""

import argparse
import gzip
import statistics
import sys
import time
from array import array
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from contextlib import nullcontext
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # the two helpers are copied next to this script
from architectures import read_clans  # noqa: E402
from classify_multi_hits import read_msa, region  # noqa: E402

WINDOW_GENES = 5
WINDOW_BP = 10_000
TOP = 25
DEFAULT_TARGETS = ("CL0177", "PF00005", "PF00528,family")  # solute-binding, ATPase, second permease
OCCURRENCE = ["protein_id", "contig_id", "contig_length", "start_position", "end_position", "strand"]
GENE = ["contig_id", "start_position", "end_position", "strand", "protein_id", "cluster_rep"]
PARTNER_FIELDS = ("partner", "name", "directon_proteins", "directon_complete", "directon_pct", "census_proteins",
                  "census_complete", "census_pct", "modal_offset", "adjacent_pct", "median_gap_bp", "enrichment")

_SCAN = {}  # set before the worker pool forks: path, key column, sorted ids (array, shared read-only), columns


def _scan_row_group(rg):
    """Rows of one row group whose key is in the ids; only the ids inside the group's key min / max are looked up,
    so a sorted file is skipped and filtered cheaply, and an unsorted one costs one lookup of every id."""
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    if "file" not in _SCAN:  # one open file per worker process
        _SCAN["file"] = pq.ParquetFile(_SCAN["path"])
    pf, ids, key = _SCAN["file"], _SCAN["ids"], _SCAN["key"]
    lo, hi = 0, len(ids)
    stats = pf.metadata.row_group(rg).column(pf.schema_arrow.get_field_index(key)).statistics
    if stats is not None and stats.has_min_max:
        lo, hi = bisect_left(ids, stats.min), bisect_right(ids, stats.max)
        if lo == hi:
            return None
    if (lo, hi) == (0, len(ids)):
        if "all" not in _SCAN:
            _SCAN["all"] = pa.array(ids.tolist(), type=pa.int64())
        values = _SCAN["all"]
    else:
        values = pa.array(ids[lo:hi].tolist(), type=pa.int64())
    table = pf.read_row_group(rg, columns=_SCAN["columns"])
    return table.filter(pc.is_in(table.column(key), value_set=values))


def scan(parquet, key, ids, columns, cpus=1, label=""):
    """pyarrow Table of the rows whose `key` is in `ids`, one row group per worker (a single dataset scan of the
    117 GB Pfam parquet ran out of 64 GB)."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    n = pq.ParquetFile(str(parquet)).num_row_groups
    _SCAN.clear()
    _SCAN.update(path=str(parquet), key=key, ids=array("q", sorted(ids)), columns=columns)
    tables, t0 = [], time.time()
    with Pool(cpus) if cpus > 1 else nullcontext() as pool:
        results = pool.imap_unordered(_scan_row_group, range(n), chunksize=8) if pool else map(_scan_row_group, range(n))
        for i, table in enumerate(results, 1):
            if table is not None and table.num_rows:
                tables.append(table)
            if i % 5000 == 0 or i == n:
                print(f"#   {label} {i:,}/{n:,} row groups, {sum(t.num_rows for t in tables):,} rows, "
                      f"{time.time() - t0:.0f}s", flush=True)
    _SCAN.clear()
    if not tables:
        return pa.table({c: pa.array([], type=pa.int64()) for c in columns})
    return pa.concat_tables(tables)


def members_of(msa):
    """(family, {protein_id}) from a full MSA <family>.sto.gz."""
    family, _, rows = read_msa(msa)
    return family, {int(region(name, aligned)[0]) for name, aligned in rows}


def longest_contigs(table):
    """{protein: (contig, contig_length, start, end, strand)}, the longest contig of each (ties: lowest contig id)."""
    best = {}
    for protein, contig, length, start, end, strand in zip(*(table.column(c).to_pylist() for c in OCCURRENCE)):
        old = best.get(protein)
        if old is None or (length, -contig) > (old[1], -old[0]):
            best[protein] = (contig, length, start, end, strand)
    return best


def contig_genes(table):
    """{contig: [(start, end, strand, protein, cluster_rep)]} sorted by position."""
    genes = defaultdict(list)
    for contig, *gene in zip(*(table.column(c).to_pylist() for c in GENE)):
        genes[contig].append(tuple(gene))
    for rows in genes.values():
        rows.sort()
    return dict(genes)


def bp_gap(a_start, a_end, b_start, b_end):
    """Intergenic gap in bp; 0 when the genes overlap."""
    return max(0, b_start - a_end, a_start - b_end)


def windows(genes, i, contig_length, n=WINDOW_GENES, bp=WINDOW_BP):
    """(directon, census, directon_complete, census_complete) around genes[i]; windows are [(offset, index)] with the
    offset read in the gene's own direction."""
    start, end, strand = genes[i][:3]
    direction = 1 if strand >= 0 else -1
    directon, census, d_complete, c_complete = [], [], True, True
    for side in (-1, 1):
        step = side * direction  # index step of this side
        for k in range(1, n + 1):
            j = i + step * k
            if not 0 <= j < len(genes):
                d_complete = False
                break
            if genes[j][2] != strand:
                break
            directon.append((side * k, j))
        for k in range(1, n + 1):
            j = i + step * k
            if not 0 <= j < len(genes):
                break
            if bp_gap(start, end, genes[j][0], genes[j][1]) <= bp:
                census.append((side * k, j))
        edge = start - 1 if step < 0 else contig_length - end
        if not 0 <= i + step * n < len(genes) and edge < bp:
            c_complete = False
    return directon, census, d_complete, c_complete


def gene_names(gene, pfams, members, family):
    """Partner names of one gene: its Pfam accessions, MGnifam_<family> when its cluster rep is a member, else its
    cluster."""
    names = list(pfams.get(gene[3], ()))
    if gene[4] in members:
        names.append(f"MGnifam_{family}")
    return names or [f"cluster:{gene[4]}"]


def expand(target, clan_members, family):
    """Names one --target group stands for: PFxxxxx as is, CLxxxx -> its Pfam families, "family" -> MGnifam_<family>."""
    names = set()
    for token in target.split(","):
        token = token.strip()
        if token == "family":
            names.add(f"MGnifam_{family}")
        elif token.startswith("CL"):
            if token not in clan_members:
                raise ValueError(f"clan {token} is not in the clans file")
            names |= clan_members[token]
        else:
            names.add(token)
    return names


def tally(anchors, genes_by_contig, pfams, members, family):
    """Per-partner counts over the anchors [(protein, contig, contig_length, start)], and one row per protein."""
    stats = defaultdict(lambda: {"d": 0, "dc": 0, "c": 0, "cc": 0, "offsets": Counter(), "gaps": []})
    window_hits, window_genes, per_protein = Counter(), 0, []
    denominators = Counter()
    for protein, contig, length, start in anchors:
        genes = genes_by_contig.get(contig, [])
        i = next((k for k, g in enumerate(genes) if g[3] == protein and g[0] == start), None)
        if i is None:
            denominators["anchor_missing"] += 1
            continue
        directon, census, d_complete, c_complete = windows(genes, i, length)
        denominators["anchors"] += 1
        denominators["directon_complete"] += d_complete
        denominators["census_complete"] += c_complete
        nearest = {}
        for offset, j in directon:
            names = gene_names(genes[j], pfams, members, family)
            window_genes += 1
            window_hits.update(names)
            for name in names:
                if name not in nearest or abs(offset) < abs(nearest[name][0]):
                    nearest[name] = (offset, bp_gap(genes[i][0], genes[i][1], genes[j][0], genes[j][1]))
        for name, (offset, gap) in nearest.items():
            s = stats[name]
            s["d"] += 1
            s["dc"] += d_complete
            s["offsets"][offset] += 1
            s["gaps"].append(gap)
        census_names = {name for _, j in census for name in gene_names(genes[j], pfams, members, family)}
        for name in census_names:
            stats[name]["c"] += 1
            stats[name]["cc"] += c_complete
        per_protein.append((protein, contig, length, genes[i][2], len(genes), int(d_complete), int(c_complete),
                            " ".join(f"{n}@{o:+d}" for n, (o, _) in sorted(nearest.items(), key=lambda x: x[1][0])),
                            " ".join(sorted(census_names)), nearest, d_complete))
    background, total = Counter(), 0
    placed = {(a[1], a[0], a[3]) for a in anchors}  # the member genes themselves are not background
    for contig in {a[1] for a in anchors}:
        for gene in genes_by_contig.get(contig, []):
            if (contig, gene[3], gene[0]) in placed:
                continue
            total += 1
            background.update(gene_names(gene, pfams, members, family))
    return stats, denominators, per_protein, (window_hits, window_genes, background, total)


def partner_rows(stats, denominators, frequencies, names):
    """Rows of PARTNER_FIELDS, most directon proteins first."""
    window_hits, window_genes, background, total = frequencies
    rows = []
    for partner, s in stats.items():
        offset, n_modal = s["offsets"].most_common(1)[0] if s["offsets"] else (0, 0)
        observed = window_hits[partner] / window_genes if window_genes else 0
        expected = background[partner] / total if total else 0
        rows.append((partner, names.get(partner, ""), s["d"], s["dc"],
                     100 * s["dc"] / denominators["directon_complete"] if denominators["directon_complete"] else 0,
                     s["c"], s["cc"],
                     100 * s["cc"] / denominators["census_complete"] if denominators["census_complete"] else 0,
                     offset, 100 * (s["offsets"][-1] + s["offsets"][1]) / s["d"] if s["d"] else 0,
                     statistics.median(s["gaps"]) if s["gaps"] else "", observed / expected if expected else ""))
    rows.sort(key=lambda r: (-r[3], -r[2], r[0]))
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("msa", help="full MSA <family>.sto.gz")
    parser.add_argument("--clusters-metadata", required=True, help="mgy_clusters_metadata.parquet (representatives)")
    parser.add_argument("--metadata", required=True, help="mgy_proteins_metadata.parquet (every gene)")
    parser.add_argument("--pfam", required=True, help="mgy_proteins_pfam.parquet")
    parser.add_argument("--clans", required=True, help="Pfam-A.clans.tsv.gz of the parquet's Pfam release")
    parser.add_argument("--target", action="append", help=f"hypothesis group, repeatable (default {DEFAULT_TARGETS})")
    parser.add_argument("--out", required=True, help="every partner (.tsv.gz)")
    parser.add_argument("--proteins", required=True, help="one row per member protein (.tsv.gz)")
    parser.add_argument("--cpus", type=int, default=1, help="row-group workers")
    args = parser.parse_args(argv)
    t0 = time.time()
    clan_of, names = read_clans(args.clans)
    clan_members = defaultdict(set)
    for acc, clan in clan_of.items():
        clan_members[clan].add(acc)
    family, members = members_of(args.msa)
    targets = [(t, expand(t, clan_members, family)) for t in (args.target or DEFAULT_TARGETS)]
    names[f"MGnifam_{family}"] = f"MGnifam {family}"
    print(f"# family {family}: {len(members):,} member proteins, {time.time() - t0:.0f}s", flush=True)

    occurrences = scan(args.clusters_metadata, "protein_id", members, OCCURRENCE, args.cpus, "occurrences")
    best = longest_contigs(occurrences)
    del occurrences
    print(f"# {len(best):,} member proteins on a contig ({len(members) - len(best):,} not found), "
          f"{len({b[0] for b in best.values()}):,} contigs, {time.time() - t0:.0f}s", flush=True)

    genes_by_contig = contig_genes(scan(args.metadata, "contig_id", {b[0] for b in best.values()}, GENE, args.cpus,
                                        "genes"))
    n_genes = sum(map(len, genes_by_contig.values()))
    print(f"# {n_genes:,} genes on {len(genes_by_contig):,} contigs, {time.time() - t0:.0f}s", flush=True)

    proteins = {g[3] for rows in genes_by_contig.values() for g in rows}
    hits = scan(args.pfam, "protein_id", proteins, ["protein_id", "pfam_accession"], args.cpus, "pfam")
    pfams = defaultdict(set)
    for protein, acc in zip(hits.column("protein_id").to_pylist(), hits.column("pfam_accession").to_pylist()):
        pfams[protein].add("PF%05d" % acc)
    pfams = {p: tuple(sorted(a)) for p, a in pfams.items()}
    del hits
    print(f"# Pfam domains on {len(pfams):,} of {len(proteins):,} genes, {time.time() - t0:.0f}s", flush=True)

    anchors = [(p, b[0], b[1], b[2]) for p, b in best.items()]
    stats, denominators, per_protein, frequencies = tally(anchors, genes_by_contig, pfams, members, family)
    rows = partner_rows(stats, denominators, frequencies, names)
    dc, cc = denominators["directon_complete"], denominators["census_complete"]
    print(f"# anchors {denominators['anchors']:,} (missing on their contig {denominators['anchor_missing']:,}); "
          f"directon window complete {dc:,}, census window complete {cc:,}")

    print("target\tdirecton_complete_proteins\tdirecton_pct\tany_window_proteins")
    complete = [p for p in per_protein if p[-1]]
    for label, group in targets:
        hit = sum(1 for p in complete if group & p[-2].keys())
        anywhere = sum(1 for p in per_protein if group & p[-2].keys())
        print(f"{label}\t{hit}\t{100 * hit / dc if dc else 0:.10g}\t{anywhere}")
    every = sum(1 for p in complete if all(group & p[-2].keys() for _, group in targets))
    print(f"all {len(targets)} targets\t{every}\t{100 * every / dc if dc else 0:.10g}\t")

    print("\t".join(PARTNER_FIELDS))
    for r in rows[:TOP]:
        print("\t".join(f"{v:.10g}" if isinstance(v, float) else str(v) for v in r))
    with gzip.open(args.out, "wt") as out:
        out.write("\t".join(PARTNER_FIELDS) + "\n")
        for r in rows:
            out.write("\t".join(f"{v:.10g}" if isinstance(v, float) else str(v) for v in r) + "\n")
    with gzip.open(args.proteins, "wt") as out:
        out.write("protein\tcontig\tcontig_length\tstrand\tgenes_on_contig\tdirecton_complete\tcensus_complete"
                  "\tdirecton\tcensus\n")
        for p in per_protein:
            out.write("\t".join(map(str, p[:9])) + "\n")
    print(f"# done {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
