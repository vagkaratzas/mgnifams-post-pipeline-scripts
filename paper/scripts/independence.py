#!/usr/bin/env python3
"""Are a family's neighbourhood partners found in independent samples?

For the members whose same-strand run (directon) holds a partner, counts distinct studies, assemblies and biomes of
their contigs, and the pairwise identity of their family regions in the full MSA (match columns where both rows have
a residue; the longest region per protein). One genome sampled many times would show few studies and ~100% identity.

Usage: independence.py PROTEINS.tsv.gz MEMBERS.tsv MSA.sto.gz [--partner PF01183 ...]
  PROTEINS: neighbourhood.py --proteins; MEMBERS: contig_lookup.py members.tsv (with study_accession, biome).
Stdlib only. Read-only.
"""

import argparse
import csv
import gzip
import statistics
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # classify_multi_hits.py is copied next to this script
sys.path.insert(0, str(Path(__file__).resolve().parent))
from classify_multi_hits import read_msa, region  # noqa: E402

MIN_SHARED = 30  # identity over fewer shared match columns is noise


def pid(row):
    """Integer protein id from 'MGYP000…', '123' or 123."""
    return int(str(row).replace("MGYP", ""))


def regions(msa):
    """{protein: match-column string of its longest region} from the full MSA."""
    _, rf, rows = read_msa(msa)
    best = {}
    for name, aligned in rows:
        mgyp, start, end = region(name, aligned)
        match = "".join(c for c, r in zip(aligned, rf) if r != ".")
        p = pid(mgyp)
        if p not in best or end - start > best[p][0]:
            best[p] = (end - start, match.upper())
    return {p: m for p, (_, m) in best.items()}


def identity(a, b):
    """Identical / shared residue columns, or None below MIN_SHARED."""
    shared = [(x, y) for x, y in zip(a, b) if x.isalpha() and y.isalpha()]
    return sum(x == y for x, y in shared) / len(shared) if len(shared) >= MIN_SHARED else None


def describe(label, proteins, meta, seqs):
    """One line: n, studies, assemblies, biomes, pairwise identity quartiles and max, pairs >= 95%."""
    rows = [meta[p] for p in proteins if p in meta]
    ids = [i for a, b in combinations(sorted(p for p in proteins if p in seqs), 2)
           if (i := identity(seqs[a], seqs[b])) is not None]
    q = statistics.quantiles(ids, n=4, method="inclusive") if len(ids) > 1 else [float("nan")] * 3
    print(f"{label}: n={len(proteins)} studies={len({r['study_accession'] for r in rows})} "
          f"assemblies={len({r['assembly_accession'] for r in rows})} "
          f"pairs={len(ids)} identity q1={q[0]:.10g} median={q[1]:.10g} q3={q[2]:.10g} "
          f"max={max(ids, default=float('nan')):.10g} pairs>=0.95={sum(i >= 0.95 for i in ids)}")
    print(f"  biomes: {dict(Counter(r['biome'] for r in rows).most_common())}")
    print(f"  top studies: {Counter(r['study_accession'] for r in rows).most_common(5)}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("proteins", help="neighbourhood.py --proteins output (.tsv.gz)")
    parser.add_argument("members", help="contig_lookup.py members.tsv")
    parser.add_argument("msa", help="full MSA of the family (.sto.gz)")
    parser.add_argument("--partner", action="append", default=[], help="partner name as in the directon column")
    args = parser.parse_args(argv)
    placed = list(csv.DictReader(gzip.open(args.proteins, "rt"), delimiter="\t"))
    meta = {pid(r["protein"]): r for r in csv.DictReader(open(args.members), delimiter="\t")}
    seqs = regions(args.msa)
    describe("all members", [pid(r["protein"]) for r in placed], meta, seqs)
    complete = [r for r in placed if r["directon_complete"] == "1"]
    describe("complete directon", [pid(r["protein"]) for r in complete], meta, seqs)
    for partner in args.partner:
        has = [r for r in placed if any(t.split("@")[0] == partner for t in r["directon"].split())]
        describe(f"{partner} in directon", [pid(r["protein"]) for r in has], meta, seqs)
        describe(f"{partner} in complete directon",
                 [pid(r["protein"]) for r in has if r["directon_complete"] == "1"], meta, seqs)


def selftest():
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        t = Path(tmp)
        core = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY"  # 40 match columns
        other = "W" * 10 + core[10:]  # 75% identical to core over 40 columns
        with gzip.open(t / "f.sto.gz", "wt") as h:
            h.write("# STOCKHOLM 1.0\n"
                    f"1/1-40 {core}\n2/1-40 {core}\n3/1-40 {other}\n3/50-60 {'-' * 30}KLMNPQRSTV\n"
                    f"#=GC RF {'x' * 40}\n//\n")
        with gzip.open(t / "p.tsv.gz", "wt") as h:
            h.write("protein\tcontig\tcontig_length\tstrand\tgenes_on_contig\tdirecton_complete\tcensus_complete\t"
                    "directon\tcensus\n1\t10\t5000\t1\t5\t1\t1\tPF1@+1 PF2@+2\tPF1\n"
                    "2\t20\t5000\t1\t5\t0\t0\tPF1@+1\tPF1\n3\t30\t900\t1\t1\t0\t0\t\t\n")
        (t / "m.tsv").write_text("family\tprotein\tassembly_accession\tstudy_accession\tbiome\n"
                                 "7\tMGYP000000000001\tERZ1\tMGYS1\troot:A\n7\tMGYP000000000002\tERZ2\tMGYS1\troot:A\n"
                                 "7\tMGYP000000000003\tERZ3\tMGYS2\troot:B\n")
        seqs = regions(t / "f.sto.gz")
        assert set(seqs) == {1, 2, 3} and seqs[3] == other, seqs  # longest region of protein 3 kept
        assert identity(seqs[1], seqs[2]) == 1.0 and identity(seqs[1], seqs[3]) == 0.75
        assert identity(core, "-" * 11 + core[11:20] + "-" * 20) is None  # 9 shared columns
        main([str(t / "p.tsv.gz"), str(t / "m.tsv"), str(t / "f.sto.gz"), "--partner", "PF1"])
    print("selftest ok")


if __name__ == "__main__":
    selftest() if sys.argv[1:] == ["--selftest"] else main()
