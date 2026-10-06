#!/usr/bin/env python3
"""
Pfam version checks for the MGnifams manuscript.

infer       Which Pfam release annotated a set of MGnify proteins? Pfam accessions are
            release-specific: an annotation made with release X has no accession that
            X lacks, and contains hits to every family X added. For each candidate
            release it counts the accessions missing from that release.

    pfam_version_checks.py infer --input chunk_000001.csv.gz [more.csv.gz ...] \\
        --release 35.0=Pfam35.0/Pfam-A.clans.tsv.gz --release 36.0=Pfam36.0/Pfam-A.clans.tsv.gz

late-hits   MGnifams are built from regions left after masking MGnify's own Pfam
            annotation, but their representatives are searched against a newer Pfam.
            Count the families with a representative Pfam hit to an accession the
            masking release did not have (so it could not have been masked).

    pfam_version_checks.py late-hits --pfams mgnifam_pfams.csv \\
        --baseline Pfam36.0/Pfam-A.clans.tsv.gz --output late_pfam_hits.tsv

A release's family list is its Pfam-A.clans.tsv.gz (~300 KB) on
https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam<version>/ ; a Pfam-A.hmm(.gz)
also works. Input CSVs are MGnify protein tables whose `metadata` column holds a JSON
object with the Pfam hits under "p" as [accession, E-value, score, ...].
"""

import argparse
import csv
import gzip
import json
import sys
from collections import Counter

csv.field_size_limit(sys.maxsize)


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_accessions(path):
    """Accessions of one release: first column of a clans file, or ACC lines of an HMM file."""
    out = set()
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("ACC"):
                out.add(line.split()[1].split(".")[0])
            elif line.startswith("PF"):
                out.add(line.split("\t")[0].split(".")[0])
    return out


def run_infer(args):
    hits = Counter()
    for path in args.input:
        with open_text(path) as fh:
            for row in csv.DictReader(fh):
                for hit in json.loads(row["metadata"]).get("p", []):
                    hits[hit[0].split(".")[0]] += 1
    print("release\tfamilies\taccessions_missing\thits_missing\tdistinct_accessions\thits")
    for spec in args.release:
        name, path = spec.split("=", 1)
        known = load_accessions(path)
        missing = [a for a in hits if a not in known]
        print("\t".join(map(str, [name, len(known), len(missing), sum(hits[a] for a in missing),
                                  len(hits), sum(hits.values())])))


def run_late_hits(args):
    baseline = load_accessions(args.baseline)
    fams = {}
    with open_text(args.pfams) as fh:
        for row in csv.DictReader(fh):
            fams.setdefault(row["id"], []).append(row["pfam"].split(".")[0])
    n_newer = n_only = 0
    with open(args.output, "w") as out:
        out.write("id\tn_hits\tn_newer\tnewer_accessions\n")
        for fam, accs in fams.items():
            newer = sorted({a for a in accs if a not in baseline})
            n_newer += bool(newer)
            n_only += bool(newer) and all(a not in baseline for a in accs)
            out.write("%s\t%d\t%d\t%s\n" % (fam, len(accs), sum(a not in baseline for a in accs),
                                             ",".join(newer)))
    print("baseline families: %d" % len(baseline))
    print("families with a Pfam hit: %d" % len(fams))
    print("with a hit to a family newer than the baseline: %d" % n_newer)
    print("with only newer-family hits: %d" % n_only)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    i = sub.add_parser("infer", help="which Pfam release annotated these proteins")
    i.add_argument("--input", nargs="+", required=True)
    i.add_argument("--release", action="append", required=True, metavar="NAME=FAMILY_LIST")
    i.set_defaults(func=run_infer)
    d = sub.add_parser("late-hits", help="MGnifam Pfam hits newer than the masking release")
    d.add_argument("--pfams", required=True, help="mgnifam_pfams.csv (id,pfam,...)")
    d.add_argument("--baseline", required=True, help="family list of the masking release")
    d.add_argument("--output", required=True)
    d.set_defaults(func=run_late_hits)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
