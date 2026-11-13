#!/usr/bin/env python3
"""Source contigs of MGnifam member proteins, for a BLASTx check of their UniProtKB hits.

Each member protein (a MGnify90 cluster representative) is placed on its longest contig (the rule of
neighbourhood.py), with the assembly (ERZ), contig name and pipeline version that the MGnify API or the
processed-contigs FASTA needs to return the contig sequence. With --genes, every gene call on those contigs is
written too, to see whether one UniProtKB protein spans several called genes (exons split by the gene caller).

Usage: contig_lookup.py MSA [MSA ...] --clusters-metadata mgy_clusters_metadata.parquet --out members.tsv
       [--metadata mgy_proteins_metadata.parquet --genes genes.tsv.gz] [--cpus N]
stdout: per family, members placed, contigs, assemblies, pipeline versions. Needs pyarrow. Read-only.
"""

import argparse
import gzip
import statistics
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # neighbourhood.py and its two helpers are copied here
from neighbourhood import OCCURRENCE, longest_contigs, members_of, scan  # noqa: E402

EXTRA = ["assembly_accession", "contig_name", "pipeline_version", "truncation", "study_accession", "biome"]
MEMBER_FIELDS = ("family", "protein", "contig_id", "assembly_accession", "contig_name", "pipeline_version",
                 "contig_length", "start", "end", "strand", "truncation", "occurrences", "study_accession", "biome")
GENE_COLUMNS = ["contig_id", "start_position", "end_position", "strand", "protein_id"]


def mgyp(protein_id):
    return f"MGYP{protein_id:012d}"


def place(table, family_of):
    """Member rows: the longest contig of each protein, with its assembly, name, version, truncation, study and biome."""
    best = longest_contigs(table)
    extra, seen = {}, Counter()
    cols = [table.column(c).to_pylist() for c in ["protein_id", "contig_id", "start_position"] + EXTRA]
    for protein, contig, start, *rest in zip(*cols):
        seen[protein] += 1
        extra.setdefault((protein, contig, start), rest)
    rows = []
    for protein, (contig, length, start, end, strand) in best.items():
        assembly, name, version, truncation, study, biome = extra[(protein, contig, start)]
        for family in family_of[protein]:
            rows.append((family, mgyp(protein), contig, assembly, name, version, length, start, end, strand,
                         truncation, seen[protein], study, biome))
    return sorted(rows, key=lambda r: (int(r[0]), r[1]))


def summary(rows, members):
    print("family\tmembers\tplaced\tcontigs\tassemblies\tpipeline_versions\tcontig_length_median")
    for family in sorted(members, key=int):
        fam = [r for r in rows if r[0] == family]
        versions = Counter(r[5] for r in fam)
        print(f"{family}\t{len(members[family])}\t{len(fam)}\t{len({r[2] for r in fam})}\t{len({r[3] for r in fam})}"
              f"\t{' '.join(f'{v}:{n}' for v, n in sorted(versions.items()))}"
              f"\t{statistics.median(r[6] for r in fam) if fam else 'NA'}")
    assemblies = {r[3]: r[5] for r in rows}
    print(f"# all: {len({r[1] for r in rows}):,} proteins placed, {len({r[2] for r in rows}):,} contigs, "
          f"{len(assemblies):,} assemblies; assemblies by pipeline version: "
          f"{dict(sorted(Counter(assemblies.values()).items()))}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("msa", nargs="+", help="full MSAs <family>.sto.gz")
    parser.add_argument("--clusters-metadata", required=True, help="mgy_clusters_metadata.parquet (representatives)")
    parser.add_argument("--out", required=True, help="one row per family member (.tsv)")
    parser.add_argument("--metadata", help="mgy_proteins_metadata.parquet (every gene), with --genes")
    parser.add_argument("--genes", help="every gene call on the members' contigs (.tsv.gz)")
    parser.add_argument("--cpus", type=int, default=1, help="row-group workers")
    args = parser.parse_args(argv)
    if bool(args.metadata) != bool(args.genes):
        parser.error("--metadata and --genes go together")
    t0 = time.time()
    members, family_of = {}, {}
    for msa in args.msa:
        family, proteins = members_of(msa)
        members[family] = proteins
        for p in proteins:
            family_of.setdefault(p, []).append(family)
    print(f"# {len(members)} families, {len(family_of):,} member proteins", flush=True)

    table = scan(args.clusters_metadata, "protein_id", set(family_of), OCCURRENCE + EXTRA, args.cpus, "occurrences")
    rows = place(table, family_of)
    with open(args.out, "w") as f:
        f.write("\t".join(MEMBER_FIELDS) + "\n")
        f.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    summary(rows, members)
    print(f"# members -> {args.out}, {time.time() - t0:.0f}s", flush=True)

    if args.genes:
        genes = scan(args.metadata, "contig_id", {r[2] for r in rows}, GENE_COLUMNS, args.cpus, "genes")
        cols = [genes.column(c).to_pylist() for c in GENE_COLUMNS]
        out = sorted(zip(*cols))
        with gzip.open(args.genes, "wt") as f:
            f.write("contig_id\tstart\tend\tstrand\tprotein\n")
            f.writelines(f"{c}\t{s}\t{e}\t{st}\t{mgyp(p)}\n" for c, s, e, st, p in out)
        print(f"# {len(out):,} gene calls on {len({o[0] for o in out}):,} contigs -> {args.genes}, "
              f"{time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
