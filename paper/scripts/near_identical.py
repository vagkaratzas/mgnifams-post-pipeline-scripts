"""Near-identical family pairs (residue Jaccard >= 0.95) in two `mgnifam compute_clans` runs.

Pairs are matched by family id (ids are stable between releases). For each pair only in one run, its residue Jaccard
in the other run is binned (>= 0.5, 0.1-0.5, or not written: < 0.1 with containment < 0.5, or no shared protein).
--families (families.tsv.gz of the newer release) splits the pairs by update status: both updated, one, none.
Groups = connected components of the near-identical pairs; "removable" = families beyond one per group.

  python -I near_identical.py <old>_clan_edges.csv <new>_clan_edges.csv [--families families.tsv.gz] [--cut 0.95]
  python -I near_identical.py --selftest
"""
import argparse
import csv
import gzip
import sys
from collections import Counter, defaultdict


def read_edges(path):
    with open(path, newline="") as f:
        return {tuple(sorted((r["family_1"], r["family_2"]))): float(r["residue_jaccard"]) for r in csv.DictReader(f)}


def groups(pairs):
    adj = defaultdict(set)
    for a, b in pairs:
        adj[a].add(b)
        adj[b].add(a)
    seen, sizes = set(), []
    for s in adj:
        if s in seen:
            continue
        stack, n = [s], 0
        seen.add(s)
        while stack:
            u = stack.pop()
            n += 1
            for v in adj[u] - seen:
                seen.add(v)
                stack.append(v)
        sizes.append(n)
    return sizes


def bin_of(j):
    if j is None:
        return "not_written"
    return ">=0.5" if j >= 0.5 else ("0.1-0.5" if j >= 0.1 else "<0.1_by_containment")


def compare(old, new, cut, status=None):
    po = {p for p, j in old.items() if j >= cut}
    pn = {p for p, j in new.items() if j >= cut}
    out = {"old_pairs": len(po), "new_pairs": len(pn), "both": len(po & pn),
           "only_old": len(po - pn), "only_new": len(pn - po)}
    for name, ps in (("old", po), ("new", pn)):
        sizes = groups(ps)
        out[f"{name}_families"] = sum(sizes)
        out[f"{name}_groups"] = len(sizes)
        out[f"{name}_removable"] = sum(sizes) - len(sizes)
    out["only_new_by_old_jaccard"] = dict(Counter(bin_of(old.get(p)) for p in pn - po))
    out["only_old_by_new_jaccard"] = dict(Counter(bin_of(new.get(p)) for p in po - pn))
    if status:
        def cls(p):
            return {2: "both_updated", 1: "one_updated", 0: "none_updated"}[sum(status[x] == "updated" for x in p)]
        for name, ps in (("old", po), ("new", pn), ("both", po & pn), ("only_new", pn - po), ("only_old", po - pn)):
            out[f"{name}_by_status"] = dict(Counter(cls(p) for p in ps))
        # two not-updated families have the same members in both runs, so the same Jaccard
        same = [p for p in set(old) & set(new) if cls(p) == "none_updated"]
        out["none_updated_pairs_in_both_edge_files"] = len(same)
        out["none_updated_pairs_jaccard_differs"] = sum(old[p] != new[p] for p in same)
    return out


def report(out):
    for k, v in out.items():
        if isinstance(v, dict):
            print(k)
            for kk, vv in sorted(v.items()):
                print(f"  {kk}\t{vv}")
        else:
            print(f"{k}\t{v}")


def selftest():
    old = {("a", "b"): 0.97, ("c", "d"): 0.96, ("d", "e"): 0.99, ("f", "g"): 0.3}
    new = {("a", "b"): 0.98, ("c", "d"): 0.6, ("f", "g"): 0.95, ("h", "i"): 0.96}
    st = {x: "updated" for x in "abcdefg"} | {"h": "not_updated", "i": "not_updated"}
    out = compare(old, new, 0.95, st)
    assert (out["old_pairs"], out["new_pairs"], out["both"], out["only_old"], out["only_new"]) == (3, 3, 1, 2, 2), out
    assert (out["old_groups"], out["old_families"], out["old_removable"]) == (2, 5, 3), out
    assert out["only_new_by_old_jaccard"] == {"0.1-0.5": 1, "not_written": 1}, out
    assert out["only_old_by_new_jaccard"] == {">=0.5": 1, "not_written": 1}, out
    assert out["new_by_status"] == {"both_updated": 2, "none_updated": 1}, out
    assert out["none_updated_pairs_in_both_edge_files"] == 0
    print("selftest ok")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
        sys.exit()
    ap = argparse.ArgumentParser()
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("--families")
    ap.add_argument("--cut", type=float, default=0.95)
    a = ap.parse_args()
    status = None
    if a.families:
        with gzip.open(a.families, "rt", newline="") as f:
            status = {r["family_id"]: r["status"] for r in csv.DictReader(f, delimiter="\t")}
    report(compare(read_edges(a.old), read_edges(a.new), a.cut, status))
