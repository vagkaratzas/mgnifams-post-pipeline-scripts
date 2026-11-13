#!/usr/bin/env python3
"""Per-family details from a classify_multi_hits.py --rows file: copies per protein (proteins with >= 2 copies)
and the model gap between consecutive pieces of split copies (negative = overlap on the model).
Usage: rows_stats.py multi_hits_rows_<job>.tsv.gz      Self-check: rows_stats.py --selftest"""

import gzip
import sys
from collections import Counter, defaultdict


def stats(lines):
    proteins = defaultdict(list)
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if f[0] != "family":
            proteins[(f[0], f[1])].append((int(f[2]), int(f[4]), int(f[5]), int(f[8]), f[10]))
    copies, gaps = defaultdict(Counter), defaultdict(list)
    for (family, _), hits in proteins.items():
        hits.sort()
        n = max(h[3] for h in hits)
        if n >= 2:
            copies[family][n] += 1
        if hits[0][4] == "split":
            gaps[family] += [b[1] - a[2] - 1 for a, b in zip(hits, hits[1:])]
    return copies, gaps


def quartiles(values):
    v = sorted(values)
    return v[len(v) // 4], v[len(v) // 2], v[3 * len(v) // 4]


def selftest():
    rows = ["family\tmgyp\tseq_from\tseq_to\thmm_from\thmm_to\tn_end\tc_end\tcopy\tcopy_full\tclass",
            "f\t1\t1\t60\t1\t60\tcurved\tjagged\t1\tTrue\tsplit", "f\t1\t70\t120\t71\t100\tjagged\tcurved\t1\tTrue\tsplit",
            "f\t2\t1\t100\t1\t100\tcurved\tcurved\t1\tTrue\ttandem", "f\t2\t110\t210\t1\t100\tcurved\tcurved\t2\tTrue\ttandem"]
    copies, gaps = stats(rows)
    assert copies["f"] == Counter({2: 1}) and gaps["f"] == [10], (copies, gaps)
    print("selftest ok")


def main():
    if sys.argv[1:] == ["--selftest"]:
        return selftest()
    copies, gaps = stats(gzip.open(sys.argv[1], "rt"))
    for family in sorted(set(copies) | set(gaps), key=lambda x: (len(x), x)):
        g = gaps.get(family)
        print(family, "copies_per_protein", dict(sorted(copies[family].items())),
              "split_model_gap_q1_med_q3", quartiles(g) if g else "-", "n", len(g or []))


if __name__ == "__main__":
    main()
