#!/usr/bin/env python3
"""How much of each member protein do a family's regions cover?

Usage: python -I protein_cover.py ROWS.tsv.gz FAMILY OUT.tsv [--sample N] [--seed 1]   (or --selftest)
ROWS: classify_multi_hits.py --rows output (family, mgyp, seq_from, seq_to, ...). A random sample of N member
proteins (default 2000) is looked up in the MGnify Proteins API for its length and full_length flag. Cover =
residues in the union of the protein's family regions / protein length. Read-only GETs, 4 at a time.
"""
import csv
import gzip
import json
import random
import statistics
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

API = "https://www.ebi.ac.uk/metagenomics/proteins/api/v1/protein/"


def union_len(spans):
    """Residues in the union of closed 1-based spans."""
    n, end = 0, 0
    for a, b in sorted(spans):
        if b > end:
            n += b - max(a, end + 1) + 1
            end = b
    return n


def get(mgyp, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(API + mgyp, timeout=60) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
        except (urllib.error.URLError, TimeoutError):
            pass
        time.sleep(5 * (i + 1))
    return None


def summary(label, xs):
    if not xs:
        return f"{label}: 0 proteins"
    q = statistics.quantiles(xs, n=4) if len(xs) > 1 else [xs[0]] * 3
    return (f"{label}: {len(xs)} proteins; cover median {q[1]:.3f} (IQR {q[0]:.3f}–{q[2]:.3f}); "
            f">= 0.8: {sum(x >= 0.8 for x in xs)} ({100 * sum(x >= 0.8 for x in xs) / len(xs):.1f}%); "
            f">= 0.5: {sum(x >= 0.5 for x in xs)} ({100 * sum(x >= 0.5 for x in xs) / len(xs):.1f}%)")


def main(rows_path, family, out, sample=2000, seed=1):
    spans = defaultdict(list)
    with gzip.open(rows_path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["family"] == family:
                spans[f"MGYP{int(r['mgyp']):012d}"].append((int(r["seq_from"]), int(r["seq_to"])))
    picked = random.Random(seed).sample(sorted(spans), min(sample, len(spans)))
    with ThreadPoolExecutor(4) as ex:
        recs = list(ex.map(get, picked))
    res = []
    for m, rec in zip(picked, recs):
        seq = rec.get("sequence") if rec else None
        cov = union_len(spans[m]) / len(seq) if seq else ""
        res.append({"mgyp": m, "regions": len(spans[m]), "covered": union_len(spans[m]),
                    "length": len(seq) if seq else "", "full_length": rec.get("full_length", "") if rec else "",
                    "cover": f"{cov:.4f}" if seq else ""})
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(res[0]), delimiter="\t")
        w.writeheader()
        w.writerows(res)
    ok = [r for r in res if r["length"]]
    print(f"family {family}: {len(spans)} member proteins, {len(picked)} sampled (seed {seed}), "
          f"{len(ok)} with a sequence; median length {statistics.median(r['length'] for r in ok):g}")
    print(summary("all", [float(r["cover"]) for r in ok]))
    print(summary("full_length True", [float(r["cover"]) for r in ok if r["full_length"] is True]))
    print(summary("full_length False", [float(r["cover"]) for r in ok if r["full_length"] is False]))
    print(f"wrote {out}")


def selftest():
    assert union_len([(1, 10)]) == 10 and union_len([(1, 10), (5, 20)]) == 20
    assert union_len([(30, 40), (1, 10)]) == 21 and union_len([(1, 10), (2, 5)]) == 10
    assert union_len([(1, 10), (11, 12)]) == 12
    print("selftest ok")


if __name__ == "__main__":
    a = sys.argv[1:]
    if a == ["--selftest"]:
        selftest()
    else:
        main(a[0], a[1], a[2], sample=int(a[a.index("--sample") + 1]) if "--sample" in a else 2000,
             seed=int(a[a.index("--seed") + 1]) if "--seed" in a else 1)
