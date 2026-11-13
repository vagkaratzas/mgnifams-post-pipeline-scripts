"""Clan statistics from `mgnifam compute_clans` outputs: sizes, within-clan residue Jaccard, clan diameter.

Every family pair inside a clan falls in one class: a clan edge (residue_jaccard >= 0.5), a written edge below 0.5
(>= 0.1, or below 0.1 with containment >= 0.5), or not written (Jaccard < 0.1 and containment < 0.5, or no shared
protein). Unwritten pairs have no exact value, so the per-clan minimum is reported as "< 0.1" for them.

  python -I clan_stats.py <chunk>_clan_edges.csv <chunk>_clans.csv [--clan_jaccard 0.5] [--catalogue catalogue_<rel>.csv.gz]

--catalogue (paper/data/catalogue_1.1.csv.gz: id, ..., any_match) adds annotation per clan: clans with no family matching
Pfam, FunFam, PDB or AlphaFoldDB, and families in no clan with no match.
  python -I clan_stats.py --selftest
"""
import argparse
import csv
import gzip
import statistics
import sys
from collections import Counter, defaultdict, deque


def load(edges_path, clans_path, clan_jaccard):
    clan_of, members, reps = {}, defaultdict(list), {}
    with open(clans_path, newline="") as f:
        for r in csv.DictReader(f):
            clan_of[r["family_id"]] = r["clan_id"]
            members[r["clan_id"]].append(r["family_id"])
            if r["is_representative"] == "True":
                reps[r["clan_id"]] = r["family_id"]
    within = defaultdict(dict)  # clan -> {(a, b): (jaccard, containment)}
    adj = defaultdict(set)
    n_edges = 0
    with open(edges_path, newline="") as f:
        for r in csv.DictReader(f):
            a, b = r["family_1"], r["family_2"]
            j, c = float(r["residue_jaccard"]), float(r["residue_containment"])
            ca, cb = clan_of.get(a), clan_of.get(b)
            if j >= clan_jaccard:
                n_edges += 1
                assert ca is not None and ca == cb, f"clan edge {a},{b} not inside one clan"
                adj[a].add(b)
                adj[b].add(a)
            if ca is not None and ca == cb:
                within[ca][(a, b)] = (j, c)
    return members, reps, within, adj, n_edges


def diameter(fams, adj):
    best = 0
    for s in fams:
        dist, q = {s: 0}, deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if v not in dist:
                    dist[v] = dist[u] + 1
                    q.append(v)
        assert len(dist) == len(fams), "clan is not connected by its clan edges"
        best = max(best, max(dist.values()))
    return best


def summary(members, reps, within, adj, n_edges, clan_jaccard):
    out = {"clans": len(members), "families_in_clans": sum(len(m) for m in members.values()), "clan_edges": n_edges}
    pairs = Counter()
    min_class, min_fams, clan_min, diam = Counter(), Counter(), [], Counter()
    all_j = []  # written within-clan Jaccards
    for clan, fams in members.items():
        n = len(fams)
        total = n * (n - 1) // 2
        w = within[clan]
        e = sum(1 for j, _ in w.values() if j >= clan_jaccard)
        lo = sum(1 for j, _ in w.values() if 0.1 <= j < clan_jaccard)
        cont = sum(1 for j, _ in w.values() if j < 0.1)
        pairs["clan_edge"] += e
        pairs["written_0.1_to_edge"] += lo
        pairs["written_by_containment_only"] += cont
        pairs["not_written"] += total - len(w)
        pairs["all"] += total
        all_j.extend(j for j, _ in w.values())
        if total == len(w) and cont == 0:
            m = min(j for j, _ in w.values())
            clan_min.append(m)
            cls = "all pairs >= edge" if m >= clan_jaccard else "min in [0.1, edge)"
        else:
            cls = "min < 0.1"
        min_class[cls] += 1
        min_fams[cls] += n
        diam[diameter(fams, adj)] += 1
    out["within_clan_pairs"] = dict(pairs)
    out["clans_by_min_jaccard"] = dict(min_class)
    out["families_by_clan_min_jaccard"] = dict(min_fams)
    out["min_jaccard_of_clans_with_exact_min"] = min(clan_min) if clan_min else None
    out["written_within_clan_jaccard"] = {
        "n": len(all_j), "min": min(all_j), "median": statistics.median(all_j), "max": max(all_j)} if all_j else None
    out["diameter_histogram"] = dict(sorted(diam.items()))
    out["representatives"] = len(reps)
    return out


def annotation(members, catalogue_path):
    with gzip.open(catalogue_path, "rt", newline="") as f:
        match = {r["id"]: r["any_match"] == "1" for r in csv.DictReader(f)}
    in_clan = {fam for fams in members.values() for fam in fams}
    assert in_clan <= match.keys(), "clan family missing from the catalogue"
    annotated = Counter(sum(match[x] for x in fams) for fams in members.values())
    lone = [x for x in match if x not in in_clan]
    return {"families": len(match), "families_no_match": sum(not v for v in match.values()),
            "clans_no_matched_family": annotated[0], "clans_with_matched_family": len(members) - annotated[0],
            "families_in_clans_with_no_match": sum(not match[x] for x in in_clan),
            "families_in_unmatched_clans": sum(len(f) for f in members.values() if not any(match[x] for x in f)),
            "families_in_no_clan": len(lone), "families_in_no_clan_no_match": sum(not match[x] for x in lone)}


def report(out):
    for k, v in out.items():
        if isinstance(v, dict):
            print(k)
            for kk, vv in v.items():
                print(f"  {kk}\t{vv:.10g}" if isinstance(vv, float) else f"  {kk}\t{vv}")
        else:
            print(f"{k}\t{v:.10g}" if isinstance(v, float) else f"{k}\t{v}")
    p = out["within_clan_pairs"]
    for k in ("clan_edge", "written_0.1_to_edge", "written_by_containment_only", "not_written"):
        print(f"share_{k}\t{100 * p[k] / p['all']:.10g}%")


def selftest():
    import os
    import tempfile
    d = tempfile.mkdtemp()
    e, c = os.path.join(d, "e.csv"), os.path.join(d, "c.csv")
    # clan A: a-b-c chain (a,c written at 0.2); clan B: d-e clique; f singleton with a weak edge to a
    with open(e, "w") as f:
        f.write("family_1,family_2,shared_proteins,protein_jaccard,residue_jaccard,residue_containment\n"
                "a,b,1,0.6,0.6,0.7\nb,c,1,0.6,0.55,0.7\na,c,1,0.2,0.2,0.4\nd,e,1,0.9,0.96,1\na,f,1,0.1,0.15,0.3\n")
    with open(c, "w") as f:
        f.write("clan_id,family_id,is_representative\nX,a,False\nX,b,True\nX,c,False\nY,d,True\nY,e,False\n")
    loaded = load(e, c, 0.5)
    out = summary(*loaded, 0.5)
    cat = os.path.join(d, "cat.csv.gz")
    with gzip.open(cat, "wt") as f:
        f.write("id,any_match\na,0\nb,1\nc,0\nd,0\ne,0\nf,0\ng,1\n")
    ann = annotation(loaded[0], cat)
    assert ann == {"families": 7, "families_no_match": 5, "clans_no_matched_family": 1, "clans_with_matched_family": 1,
                   "families_in_clans_with_no_match": 4, "families_in_unmatched_clans": 2, "families_in_no_clan": 2,
                   "families_in_no_clan_no_match": 1}, ann
    assert out["clans"] == 2 and out["families_in_clans"] == 5 and out["clan_edges"] == 3
    assert out["within_clan_pairs"] == {"clan_edge": 3, "written_0.1_to_edge": 1, "written_by_containment_only": 0,
                                        "not_written": 0, "all": 4}, out
    assert out["clans_by_min_jaccard"] == {"min in [0.1, edge)": 1, "all pairs >= edge": 1}, out
    assert out["min_jaccard_of_clans_with_exact_min"] == 0.2 and out["diameter_histogram"] == {1: 1, 2: 1}
    print("selftest ok")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
        sys.exit()
    ap = argparse.ArgumentParser()
    ap.add_argument("edges")
    ap.add_argument("clans")
    ap.add_argument("--clan_jaccard", type=float, default=0.5)
    ap.add_argument("--catalogue")
    a = ap.parse_args()
    loaded = load(a.edges, a.clans, a.clan_jaccard)
    out = summary(*loaded, a.clan_jaccard)
    if a.catalogue:
        out["annotation"] = annotation(loaded[0], a.catalogue)
    report(out)
