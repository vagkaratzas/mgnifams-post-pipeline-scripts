#!/usr/bin/env python3
"""Biome breadth per release DB, exact values, using post-pipeline biome_analysis.py.

Usage: python3 biome_breadth.py DB.sqlite3... [--base OLD.sqlite3] [--family 25] [--dump leaves.tsv.gz]
--base fills families missing from a delta DB (the not-updated ones) from the previous release DB.
--dump writes id<TAB>n_leaves per family (last DB given) for breadth_size.py.
biome_analysis.py (mgnifams-post-pipeline-scripts bin/) must sit next to this file.
The 1.0 DB must give 121 leaf paths, median 56, max 85, 24 one-leaf families.
"""
import argparse
import gzip
import sqlite3
import statistics
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from biome_analysis import analyse_biomes, iter_biome_rows  # noqa: E402


def rows(db, base):
    seen = set()
    for fid, blob in iter_biome_rows(db):
        seen.add(fid)
        yield fid, blob
    if not base:
        return
    with sqlite3.connect(base) as con:
        ids = [i for (i,) in con.execute("SELECT id FROM mgnifam") if i not in seen]
        print(f"filled from {base}: {len(ids)} families")
        for j in range(0, len(ids), 500):  # id is the primary key: lookups, no full scan of the big DB
            c = ids[j:j + 500]
            q = f"SELECT id, biome_blob FROM mgnifam WHERE biome_blob IS NOT NULL AND id IN ({','.join('?' * len(c))})"
            yield from con.execute(q, c)


def report(db, family, base=None, dump=None):
    print(f"== {db}" + (f" + missing families from {base}" if base else ""))
    res = analyse_biomes(rows(db, base))
    fams = res["family_data"]
    n = [d["n_leaves"] for d in fams.values()]
    if dump:
        with gzip.open(dump, "wt") as fh:
            fh.write("id\tn_leaves\n")
            fh.writelines(f"{k}\t{d['n_leaves']}\n" for k, d in sorted(fams.items(), key=lambda kv: int(kv[0])))
        print(f"wrote {dump}: {len(fams)} families")
    tops_all = Counter(t for d in fams.values() for t in d["top_levels"])
    leaf_tops = Counter(str(p).split(":")[1] if ":" in str(p) else "Unknown" for p in res["leaf_df"]["path"])
    print(f"families with parseable leaves: {len(fams)}")
    print(f"leaf paths (most specific lineage of >= 1 family): {len(res['leaf_df'])}")
    print(f"leaf paths per top-level category: {dict(sorted(leaf_tops.items()))}")
    print(f"families per top-level category: {dict(sorted(tops_all.items()))}")
    print(f"leaves per family: median {statistics.median(n):.10g}, mean {statistics.mean(n):.10g}, max {max(n)}")
    print(f"families at max: {sorted(k for k, d in fams.items() if d['n_leaves'] == max(n))[:20]}")
    print(f"one-leaf families: {len(res['exclusive'])}")
    for k, leaf in sorted(res["exclusive"].items()):
        print(f"  one-leaf\t{k}\t{leaf[0]}\t{leaf[2]}")
    # paper claims around the numbers: broadest families' categories and terms, top leaves, repeated end terms
    term = lambda p: str(p).rsplit(":", 1)[-1]
    terms = {term(p) for p in res["leaf_df"]["path"]}
    top = [d for d in fams.values() if d["n_leaves"] == max(n)]
    print(f"families at max: {len(top)}; spanning all {len(tops_all)} top levels: "
          f"{sum(len(d['top_levels']) == len(tops_all) for d in top)}; "
          f"covering all {len(terms)} distinct end terms: {sum({term(l[0]) for l in d['leaves']} == terms for d in top)}")
    print("most common leaf paths (families):")
    for _, r in res["leaf_df"].sort_values("count", ascending=False).head(15).iterrows():
        print(f"  leaf\t{r['count']}\t{r['path']}")
    ends = Counter(term(p) for p in res["leaf_df"]["path"])
    print(f"end terms on > 1 leaf path: {sorted(((c, t) for t, c in ends.items() if c > 1), reverse=True)}")
    d = fams.get(family) or fams.get(str(family))
    if d:
        print(f"family {family}: {d['n_leaves']} leaves, top levels {sorted(d['top_levels'])}; "
              f"families with more leaves: {sum(x > d['n_leaves'] for x in n)}, with at least as many: {sum(x >= d['n_leaves'] for x in n)}")
        print(f"families per leaf count, top 8: {sorted(Counter(n).items(), reverse=True)[:8]}")
    else:
        print(f"family {family}: no biome leaves")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("dbs", nargs="+")
    ap.add_argument("--family", type=int, default=25)
    ap.add_argument("--base", help="previous release DB for families missing from a delta DB")
    ap.add_argument("--dump", help="write id<TAB>n_leaves (.tsv.gz) for breadth_size.py")
    a = ap.parse_args()
    for db in a.dbs:
        report(db, a.family, a.base, a.dump)
