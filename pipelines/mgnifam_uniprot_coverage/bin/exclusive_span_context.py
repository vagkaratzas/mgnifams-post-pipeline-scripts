#!/usr/bin/env python3
"""Where do MGnifam-exclusive residues sit relative to Pfam: extending a Pfam
alignment, or in a separate region of the sequence?

Reads one chunk's pair of domtbls (the same pair the coverage pipeline maps),
merges each sequence's MGnifam hits into blocks, and classifies every block:

  extension   the block overlaps or abuts a Pfam span; its exclusive residues
              lengthen a region Pfam already covers
  separate    the sequence has Pfam hits but none touches this block
  no_pfam     the sequence has no Pfam hit (the newly covered sequences)

Exclusive residues per class are binned by the block's exclusive length
(<30, 30-74, >=75 aa; 75 aa is the shortest MGnifam representative), so
short boundary differences can be told apart from domain-sized additions.
The three classes partition the MGnifam-exclusive residues of the coverage
report, which is the check that this reads the same hits.

Per-chunk output is one TSV; sum chunks with:
    awk -F'\\t' 'NR>1 && $1!="class"{r[$1"\\t"$2]+=$3; b[$1"\\t"$2]+=$4}
                 END{for(k in r) print k"\\t"r[k]"\\t"b[k]}' */*.context.tsv
"""
import argparse
import importlib.util
import os
from collections import defaultdict

_spec = importlib.util.spec_from_file_location(
    "cov", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "mgnifam_uniprot_coverage_stats_from_domtbl.py"))
cov = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cov)

BINS = ((30, "<30"), (75, "30-74"), (None, ">=75"))


def length_bin(n):
    for hi, name in BINS:
        if hi is None or n < hi:
            return name


def read_spans(path, env, evalue_scale=None, max_evalue=None):
    """{target: [(start, end), ...]} from a domtbl, with the same filters as
    the coverage map step (Pfam: as searched with --cut_ga; MGnifams: rescaled
    sequence E and c-Evalue <= max_evalue)."""
    lo, hi = (19, 20) if env else (17, 18)
    spans = defaultdict(list)
    with cov.open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split(None, 21)
            if len(f) < 21:
                continue
            if evalue_scale is not None and (
                    float(f[6]) * evalue_scale > max_evalue
                    or float(f[11]) * evalue_scale > max_evalue):
                continue
            spans[f[0]].append((int(f[lo]), int(f[hi])))
    return spans


def classify(mgnifam_spans, pfam_spans):
    """{(class, length_bin): [residues, blocks]} over all targets."""
    out = defaultdict(lambda: [0, 0])
    for target, ivs in mgnifam_spans.items():
        blocks, _ = cov.merge_intervals(ivs)
        mask, _ = cov.merge_intervals(pfam_spans.get(target, []))
        for s, e in blocks:
            _, nres = cov.subtract_intervals([(s, e)], mask)
            if not nres:
                continue
            if not mask:
                cls = "no_pfam"
            # touching = overlap or abut (merge_intervals joins abutting spans too)
            elif any(ms <= e + 1 and me >= s - 1 for ms, me in mask):
                cls = "extension"
            else:
                cls = "separate"
            key = (cls, length_bin(nres))
            out[key][0] += nres
            out[key][1] += 1
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mgnifams-domtbl", required=True)
    ap.add_argument("--pfam-domtbl", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--evalue-scale", type=float, required=True,
                    help="target_z / search_z, as in the coverage run")
    ap.add_argument("--max-evalue", type=float, default=0.001)
    ap.add_argument("--env", action="store_true",
                    help="envelope instead of alignment coordinates")
    a = ap.parse_args()

    res = classify(
        read_spans(a.mgnifams_domtbl, a.env, a.evalue_scale, a.max_evalue),
        read_spans(a.pfam_domtbl, a.env))
    with open(a.out, "w") as fh:
        fh.write("class\tlength_bin\texclusive_residues\tblocks\n")
        for (cls, b), (nres, nblk) in sorted(res.items()):
            fh.write("%s\t%s\t%d\t%d\n" % (cls, b, nres, nblk))


if __name__ == "__main__":
    main()
