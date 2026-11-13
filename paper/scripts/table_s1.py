#!/usr/bin/env python3
"""Supplementary Table S1: annotation evidence of the MGnifams of one release, from the annotation_evidence log.

Input: the log of `annotation_evidence.py` (committed in data/); its
"== <release>" block gives every count. Two rows are derived: PDB and AlphaFoldDB = PDB + AlphaFoldDB - either,
PDB, not AlphaFoldDB = either - AlphaFoldDB. Stdlib only.

Usage (from the repository root):
  python paper/scripts/table_s1.py [--release 1.1] [--selftest]
"""
import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(HERE, "..", "data", "mgnifams_annotation_evidence_20730227.log")


def block(text, release):
    m = re.search(rf"^== {re.escape(release)}: (\d+) families\n(.*?)(?=^-- |^== |\Z)", text, re.M | re.S)
    if not m:
        sys.exit(f"no '== {release}' block in the log")
    return int(m.group(1)), m.group(2)


def num(pattern, text):
    m = re.search(pattern, text)
    if not m:
        sys.exit(f"pattern not in the log block: {pattern}")
    return int(m.group(1))


def rows(n, b):
    pfam = num(r"Pfam (\d+), profile only", b)
    profile_only = num(r"profile only (\d+) ", b)
    funfam = num(r"FunFam (\d+), FunFam only", b)
    funfam_only = num(r"FunFam only (\d+)", b)
    family = num(r"family resource (\d+) ", b)
    late = sum(int(x) for x in re.findall(r"late Pfam \([^)]*\): (\d+) of", b))
    pdb, afdb = num(r"PDB (\d+), AlphaFoldDB", b), num(r"AlphaFoldDB (\d+), either", b)
    either = num(r"either (\d+) ", b)
    afdb_not_pdb = num(r"AlphaFoldDB not PDB (\d+)", b)
    afdb_only = num(r"AlphaFoldDB only evidence (\d+)", b)
    none = num(r"no match at all (\d+) ", b)
    no_fam_pdb = num(r"no Pfam \(seq or profile\), FunFam or PDB (\d+) ", b)
    both = pdb + afdb - either
    assert afdb - both == afdb_not_pdb, (afdb, both, afdb_not_pdb)  # the log's own overlap agrees
    assert family == pfam + profile_only + funfam_only, family  # disjoint: Pfam, profile only, FunFam only
    assert none <= no_fam_pdb, (none, no_fam_pdb)
    out = [("Pfam match at the sequence level", pfam),
           ("Pfam match at the profile level only", profile_only),
           ("FunFam match", funfam),
           ("FunFam match, no Pfam match", funfam_only),
           ("Pfam or FunFam match (protein family resource)", family),
           ("Pfam sequence-level matches only to Pfam families added after the masking release", late),
           ("PDB match", pdb),
           ("AlphaFoldDB match", afdb),
           ("PDB or AlphaFoldDB match", either),
           ("PDB and AlphaFoldDB match", both),
           ("PDB match, no AlphaFoldDB match", either - afdb),
           ("AlphaFoldDB match is the only evidence", afdb_only),
           ("no Pfam, FunFam or PDB match", no_fam_pdb),
           ("no match to any resource", none)]
    return [(label, c, f"{100 * c / n:.1f}") for label, c in out]


def selftest():
    text = open(LOG).read()
    n, b = block(text, "1.1")
    r = {label: c for label, c, _ in rows(n, b)}
    assert n == 35459 and r["Pfam or FunFam match (protein family resource)"] == 17203, r
    assert r["Pfam sequence-level matches only to Pfam families added after the masking release"] == 636  # 434 + 202
    assert r["PDB and AlphaFoldDB match"] == 13693 and r["PDB match, no AlphaFoldDB match"] == 30, r
    n0, b0 = block(text, "1.0")
    assert {label: c for label, c, _ in rows(n0, b0)}["no match to any resource"] == 7266
    print("selftest ok")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--release", default="1.1")
    p.add_argument("--log", default=LOG)
    p.add_argument("--out", default="supplementary_tables/table-s1_annotation-evidence.tsv")
    p.add_argument("--selftest", action="store_true")
    a = p.parse_args()
    selftest()
    if a.selftest:
        return
    n, b = block(open(a.log).read(), a.release)
    out = a.out if a.release == "1.1" else a.out.replace(".tsv", f"_{a.release}.tsv")
    with open(out, "w") as fh:
        fh.write(f"evidence\tfamilies\t% of {n:,}\n")
        for label, c, pc in rows(n, b):
            fh.write(f"{label}\t{c}\t{pc}\n")
            print(f"{label:<85} {c:>7,} {pc:>5}")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
