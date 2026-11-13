#!/usr/bin/env python3
"""Biome breadth vs family size: Spearman rho and median leaves per size decade.

Usage: python -I breadth_size.py LEAVES.tsv.gz CATALOGUE.csv.gz   (or --selftest)
LEAVES: id, n_leaves from biome_breadth.py --dump. CATALOGUE: paper/data/catalogue_1.1.csv.gz (id, full_size, ...).
Size = regions in the full alignment (Figure 3A). Families without biome leaves are counted, not ranked.
"""
import csv
import gzip
import statistics
import sys
from collections import defaultdict


def ranks(xs):
    """Average ranks (1-based), ties share the mean rank."""
    order = sorted(range(len(xs)), key=xs.__getitem__)
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def spearman(x, y):
    return statistics.correlation(ranks(x), ranks(y))  # Pearson on average ranks


def table(path, delim):
    with gzip.open(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def main(leaves_path, catalogue_path):
    leaves = {int(r["id"]): int(r["n_leaves"]) for r in table(leaves_path, "\t")}
    size = {int(r["id"]): int(r["full_size"]) for r in table(catalogue_path, ",")}
    ids = sorted(set(leaves) & set(size))
    print(f"families: catalogue {len(size)}, with leaves {len(leaves)}, both {len(ids)}, "
          f"catalogue without leaves {len(set(size) - set(leaves))}")
    x, y = [size[i] for i in ids], [leaves[i] for i in ids]
    print(f"Spearman rho (full_size, n_leaves): {spearman(x, y):.4f} (n = {len(ids)})")
    by = defaultdict(list)
    for s, n in zip(x, y):
        by[len(str(s)) - 1].append(n)  # size decade, left-closed as in Figure 3A
    print("size decade\tfamilies\tmedian leaves\tmin\tmax")
    for d in sorted(by):
        v = by[d]
        print(f"{10 ** d:,}–{10 ** (d + 1) - 1:,}\t{len(v)}\t{statistics.median(v):g}\t{min(v)}\t{max(v)}")


def selftest():
    assert ranks([10, 20, 20, 30]) == [1.0, 2.5, 2.5, 4.0]
    assert abs(spearman([1, 2, 3, 4], [1, 4, 9, 16]) - 1) < 1e-12
    assert abs(spearman([1, 2, 3, 4], [4, 3, 2, 1]) + 1) < 1e-12
    # scipy.stats.spearmanr([1,2,3,4,5],[2,1,4,3,5]).statistic == 0.8
    assert abs(spearman([1, 2, 3, 4, 5], [2, 1, 4, 3, 5]) - 0.8) < 1e-12
    print("selftest ok")


if __name__ == "__main__":
    if sys.argv[1:] == ["--selftest"]:
        selftest()
    else:
        main(*sys.argv[1:3])
