#!/usr/bin/env python3
"""Full domain-architecture counts of MGnifams member proteins.

One architecture per distinct member protein: its MGnifam copies and Pfam domains, N- to C-terminal by start.
- MGnifam copies come from the full MSA by the copy rules of classify_multi_hits.py: a split copy is one token.
  Full / partial and curved / jagged ends do not change the architecture; each row counts them per MGnifam token.
- Pfam domains are the MGnify Pfam parquet hits (GA thresholds, every hit). Of overlapping hits from one Pfam clan
  (a family without clan is its own), the most significant (i-Evalue, then score) is kept, as Pfam does.

Usage: architectures.py MSA... --pfam mgy_proteins_pfam.parquet --clans Pfam-A.clans.tsv.gz --out arch.tsv.gz
stdout: per family, proteins, architectures, rank at which 95% of proteins are covered, proteins with no Pfam,
then the top rows. --out: every architecture row. Needs pyarrow (no numpy). Read-only.
"""

import argparse
import gzip
import sys
import time
from bisect import bisect_left
from collections import defaultdict
from contextlib import nullcontext
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # classify_multi_hits.py is copied next to this script
from classify_multi_hits import chain, ends, is_full, model_span, read_msa, region  # noqa: E402

COVER = 0.95
TOP = 15
FIELDS = ("family", "rank", "proteins", "pct", "cum_pct", "all_full_pct", "architecture", "names", "mgnifam_ends")


def mgnifam_copies(path):
    """(family, {protein_id: [(start, end, full, n_curved, c_curved), ...]}), one tuple per copy."""
    family, rf, rows = read_msa(path)
    model_length = sum(c != "." for c in rf)
    proteins = defaultdict(list)
    for name, aligned in rows:
        span = model_span(aligned, rf)
        if span is None:
            raise ValueError(f"{path}: row {name} has no residue in a match column")
        mgyp, start, end = region(name, aligned)
        proteins[int(mgyp)].append((start, end) + span)
    return family, {protein: [copy_tuple(c, model_length) for c in chain(hits, model_length)]
                    for protein, hits in proteins.items()}


def copy_tuple(copy, model_length):
    """(start, end, full, n_curved, c_curved) of one copy; its ends are those of its first and last model positions."""
    n_end, c_end = ends(min(r[2] for r in copy), max(r[3] for r in copy), model_length)
    return (min(r[0] for r in copy), max(r[1] for r in copy), is_full(copy, model_length),
            n_end == "curved", c_end == "curved")


def resolve(hits, clan_of):
    """Pfam clan rule on (acc, i_evalue, score, env_from, env_to) hits; returns the kept ones by start.
    shortcut: no nesting exceptions (Pfam allows listed nested families), upgrade if a case study shows one."""
    kept = []
    for hit in sorted(hits, key=lambda h: (h[1], -h[2])):
        clan = clan_of.get(hit[0], hit[0])
        if not any(clan_of.get(k[0], k[0]) == clan and hit[3] <= k[4] and k[3] <= hit[4] for k in kept):
            kept.append(hit)
    return sorted(kept, key=lambda h: (h[3], h[4]))


def architecture(family, copies, pfams):
    """(tokens N to C, [(full, n_curved, c_curved) per MGnifam token])."""
    tokens = sorted([(c[0], c[1], family, c[2:]) for c in copies] + [(h[3], h[4], h[0], None) for h in pfams],
                    key=lambda t: (t[0], t[1]))
    return tuple(t[2] for t in tokens), [t[3] for t in tokens if t[3] is not None]


def aggregate(family, copies, pfams, clan_of):
    """Architecture rows of one family, most proteins first (ties by architecture)."""
    rows = {}
    for protein, protein_copies in copies.items():
        key, mgnifam = architecture(family, protein_copies, resolve(pfams.get(protein, []), clan_of))
        row = rows.setdefault(key, {"architecture": key, "proteins": 0, "all_full": 0,
                                    "ends": [[0, 0, 0] for _ in mgnifam]})
        row["proteins"] += 1
        row["all_full"] += all(m[0] for m in mgnifam)
        for counts, flags in zip(row["ends"], mgnifam):
            for i, flag in enumerate(flags):
                counts[i] += flag
    return sorted(rows.values(), key=lambda r: (-r["proteins"], r["architecture"]))


def rank_at(counts, cover):
    """Smallest number of rows (sorted, most first) whose proteins reach `cover` of the total."""
    total, running = sum(counts), 0
    for rank, n in enumerate(counts, 1):
        running += n
        if running >= cover * total:
            return rank
    return len(counts)


PFAM_COLUMNS = ["protein_id", "pfam_accession", "i_evalue", "score", "env_from", "env_to"]
_SCAN = {}  # set before the worker pool forks: path, sorted ids, value set, protein_id column index


def _scan_row_group(rg):
    """Member rows of one row group; skipped by its protein_id min / max when no member falls in the range."""
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    if "file" not in _SCAN:  # one open file per worker process
        _SCAN["file"] = pq.ParquetFile(_SCAN["path"])
    pf, ids = _SCAN["file"], _SCAN["ids"]
    stats = pf.metadata.row_group(rg).column(_SCAN["column"]).statistics
    if stats is not None and stats.has_min_max:
        i = bisect_left(ids, stats.min)
        if i == len(ids) or ids[i] > stats.max:
            return []
    table = pf.read_row_group(rg, columns=PFAM_COLUMNS)
    table = table.filter(pc.is_in(table.column("protein_id"), value_set=_SCAN["values"]))
    return list(zip(*(table.column(c).to_pylist() for c in PFAM_COLUMNS)))


def pfam_hits(parquet, proteins, cpus=1):
    """{protein_id: [(PFxxxxx, i_evalue, score, env_from, env_to)]}, one row group at a time as pfam_recount.sh
    (a single dataset scan of the 117 GB parquet ran out of 64 GB)."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    pf = pq.ParquetFile(str(parquet))
    n, schema = pf.num_row_groups, pf.schema_arrow
    ids = sorted(proteins)
    _SCAN.clear()
    _SCAN.update(path=str(parquet), ids=ids, column=schema.get_field_index("protein_id"),
                 values=pa.array(ids, type=schema.field("protein_id").type))
    del pf
    hits, t0 = defaultdict(list), time.time()
    with Pool(cpus) if cpus > 1 else nullcontext() as pool:
        results = pool.imap_unordered(_scan_row_group, range(n), chunksize=8) if pool else map(_scan_row_group, range(n))
        for i, rows in enumerate(results, 1):
            for protein, acc, i_evalue, score, start, end in rows:
                hits[protein].append(("PF%05d" % acc, i_evalue, score, start, end))
            if i % 2000 == 0 or i == n:
                print(f"#   {i:,}/{n:,} row groups, {time.time() - t0:.0f}s", flush=True)
    for rows in hits.values():
        rows.sort(key=lambda h: (h[3], h[4], h[0]))  # workers return row groups in any order
    return dict(hits)


def read_clans(path):
    """({acc: clan acc}, {acc: Pfam id}) from Pfam-A.clans.tsv.gz."""
    clan_of, names = {}, {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if f[1]:
                clan_of[f[0]] = f[1]
            names[f[0]] = f[3]
    return clan_of, names


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("msa", nargs="+", help="full MSA files <family>.sto.gz")
    parser.add_argument("--pfam", required=True, help="mgy_proteins_pfam.parquet")
    parser.add_argument("--clans", required=True, help="Pfam-A.clans.tsv.gz of the parquet's Pfam release")
    parser.add_argument("--out", required=True, help="every architecture row (.tsv.gz)")
    parser.add_argument("--cpus", type=int, default=1, help="row-group workers")
    args = parser.parse_args(argv)
    t0 = time.time()
    families = [mgnifam_copies(path) for path in args.msa]
    proteins = {p for _, copies in families for p in copies}
    print(f"# {len(families)} families, {len(proteins)} distinct proteins, {time.time() - t0:.0f}s", flush=True)
    pfams = pfam_hits(args.pfam, proteins, args.cpus)
    print(f"# Pfam rows for {len(pfams)} proteins, {sum(map(len, pfams.values()))} hits, {time.time() - t0:.0f}s",
          flush=True)
    clan_of, names = read_clans(args.clans)
    print("family\tproteins\tarchitectures\trank_at_95\tno_pfam_proteins\ttop_pct")
    with gzip.open(args.out, "wt") as out:
        out.write("\t".join(FIELDS) + "\n")
        for family, copies in families:
            rows = aggregate(family, copies, pfams, clan_of)
            total, running = sum(r["proteins"] for r in rows), 0
            no_pfam = sum(r["proteins"] for r in rows if all(t == family for t in r["architecture"]))
            for rank, r in enumerate(rows, 1):
                running += r["proteins"]
                out.write("\t".join(map(str, (
                    family, rank, r["proteins"], f"{100 * r['proteins'] / total:.10g}", f"{100 * running / total:.10g}",
                    f"{100 * r['all_full'] / r['proteins']:.10g}", " ".join(r["architecture"]),
                    " ".join(names.get(t, t) for t in r["architecture"]),
                    ",".join("/".join(map(str, e)) for e in r["ends"])))) + "\n")
            print(f"{family}\t{total}\t{len(rows)}\t{rank_at([r['proteins'] for r in rows], COVER)}\t{no_pfam}\t"
                  + " | ".join(f"{' '.join(r['architecture'])} {r['proteins']} ({100 * r['proteins'] / total:.1f}%)"
                               for r in rows[:5]), flush=True)
    print(f"# done {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
