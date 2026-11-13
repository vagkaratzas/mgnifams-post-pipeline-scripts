#!/usr/bin/env python3
"""Pfam domains per family and per clan, from the pipeline's per-family query TSVs.

Input: query_results/<family>.tsv written by the pipeline's build_parquet_domain_queries.py, one line per distinct
member protein: `mgyp<TAB><TAB>[['PF%05d', i_evalue, score, hmm_from, hmm_to, env_from, env_to], ...]` (empty when
the protein has no Pfam hit). Overlapping hits of one Pfam clan are resolved as in architectures.resolve.
Per family: proteins, proteins with >= 1 Pfam, distinct Pfams on any protein and on >= 10% of proteins.
Clans (MGnifams clans.tsv.gz) roll up: a clan's Pfams are the union of its families' Pfams. Stdlib only. Read-only.

Usage: pfam_per_family.py QUERY_DIR... --clans Pfam-A.clans.tsv.gz --mgnifam-clans clans.tsv.gz --out x.tsv.gz
"""

import argparse
import gzip
import re
import statistics
import sys
import time
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # architectures.py (+ classify_multi_hits.py) next to this
from architectures import read_clans, resolve  # noqa: E402

CORE = 0.10
HIT = re.compile(r"\['(PF\d+)', ([^,]+), ([^,]+), \d+, \d+, (\d+), (\d+)\]")
FIELDS = ("family", "proteins", "with_pfam", "pct_with_pfam", "distinct_any", "distinct_core", "core_pfams")
_CLAN_OF = {}  # Pfam acc -> Pfam clan, set before the pool forks


def family_counts(path):
    """(family, proteins, proteins with >= 1 Pfam, {acc: proteins carrying it}) of one query TSV."""
    proteins = with_pfam = 0
    per_acc = Counter()
    with open(path) as handle:
        for line in handle:
            proteins += 1
            hits_text = line.rstrip("\n").split("\t")[2]
            if not hits_text:
                continue
            hits = [(a, float(e), float(s), int(f), int(t)) for a, e, s, f, t in HIT.findall(hits_text)]
            if len(hits) != hits_text.count("['PF"):
                raise ValueError(f"{path}: unparsed hit in {line[:200]!r}")
            with_pfam += 1
            per_acc.update({h[0] for h in resolve(hits, _CLAN_OF)})
    return Path(path).stem, proteins, with_pfam, per_acc


def core(per_acc, proteins):
    """Pfams on at least CORE of the family's proteins."""
    return {a for a, n in per_acc.items() if n >= CORE * proteins}


def summary(label, values):
    """One line: n, quartiles ('inclusive', as the other probes), mean, max, and the 0..5+ histogram."""
    values = sorted(values)
    if not values:
        print(f"{label}: n=0", flush=True)
        return
    q =statistics.quantiles(values, n=4, method="inclusive") if len(values) > 1 else [values[0]] * 3
    hist = Counter(min(v, 5) for v in values) if all(isinstance(v, int) for v in values) else None
    line = (f"{label}: n={len(values)} q1={q[0]:.10g} median={q[1]:.10g} q3={q[2]:.10g} "
            f"mean={statistics.fmean(values):.10g} max={values[-1]:.10g}")
    if hist:
        line += " | " + " ".join(f"{k if k < 5 else '5+'}:{hist[k]}" for k in range(6))
    print(line, flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("query_dirs", nargs="+", help="query_results dirs with <family>.tsv")
    parser.add_argument("--clans", required=True, help="Pfam-A.clans.tsv.gz of the parquet's Pfam release")
    parser.add_argument("--mgnifam-clans", required=True, help="MGnifams clans.tsv.gz (clan_id, family_id, ...)")
    parser.add_argument("--out", required=True, help="per-family rows (.tsv.gz)")
    parser.add_argument("--cpus", type=int, default=1)
    args = parser.parse_args(argv)
    t0 = time.time()

    paths = {}
    for d in args.query_dirs:
        for p in Path(d).glob("*.tsv"):
            if p.stem in paths:
                raise ValueError(f"family {p.stem} in two dirs: {paths[p.stem]} and {p}")
            paths[p.stem] = p
    print(f"# {len(paths)} family TSVs in {len(args.query_dirs)} dirs", flush=True)

    _CLAN_OF.update(read_clans(args.clans)[0])
    rows = {}
    with Pool(args.cpus) as pool:
        for i, (family, n, with_pfam, per_acc) in enumerate(
                pool.imap_unordered(family_counts, sorted(paths.values(), key=lambda p: -p.stat().st_size)), 1):
            rows[family] = (n, with_pfam, per_acc)
            if i % 2000 == 0 or i == len(paths):
                print(f"#   {i}/{len(paths)} families, {time.time() - t0:.0f}s", flush=True)

    empty = [f for f, r in rows.items() if r[0] == 0]
    if empty:
        print(f"# WARNING {len(empty)} families with no protein line, left out: {' '.join(sorted(empty)[:20])}")
    fams = {f: r for f, r in rows.items() if r[0]}
    with gzip.open(args.out, "wt") as out:
        out.write("\t".join(FIELDS) + "\n")
        for f in sorted(fams, key=int):
            n, w, per_acc = fams[f]
            c = core(per_acc, n)
            out.write("\t".join(map(str, (f, n, w, f"{100 * w / n:.10g}", len(per_acc), len(c),
                                          ",".join(f"{a}:{100 * per_acc[a] / n:.4g}" for a in
                                                   sorted(c, key=lambda a: -per_acc[a]))))) + "\n")

    print("\n## families")
    summary("pct_with_pfam", [100 * r[1] / r[0] for r in fams.values()])
    summary("distinct_any", [len(r[2]) for r in fams.values()])
    summary("distinct_core_10pct", [len(core(r[2], r[0])) for r in fams.values()])
    tot, tot_w = sum(r[0] for r in fams.values()), sum(r[1] for r in fams.values())
    print(f"family-protein pairs: {tot} with_pfam {tot_w} ({100 * tot_w / tot:.10g}%)")
    print(f"families with no Pfam on any protein: {sum(not r[2] for r in fams.values())}")
    print("most frequent core Pfams (families): " + " ".join(
        f"{a}:{n}" for a, n in Counter(a for r in fams.values() for a in core(r[2], r[0])).most_common(15)))

    clans = defaultdict(list)
    with gzip.open(args.mgnifam_clans, "rt") as handle:
        next(handle)
        for line in handle:
            clan, family = line.split("\t")[:2]
            clans[clan].append(family)
    missing = [f for members in clans.values() for f in members if f not in fams]
    print(f"\n## clans: {len(clans)} clans, {sum(map(len, clans.values()))} families, {len(missing)} without a TSV")
    in_clan = {f for members in clans.values() for f in members}
    units_any, units_core = [], []
    for members in clans.values():
        present = [fams[f] for f in members if f in fams]
        units_any.append(len(set().union(*(r[2] for r in present))))
        units_core.append(len(set().union(*(core(r[2], r[0]) for r in present))))
    summary("clan distinct_any", units_any)
    summary("clan distinct_core_10pct", units_core)
    solo = [r for f, r in fams.items() if f not in in_clan]
    summary("no-clan family distinct_any", [len(r[2]) for r in solo])
    summary("no-clan family distinct_core_10pct", [len(core(r[2], r[0])) for r in solo])
    summary("clans + no-clan families distinct_any", units_any + [len(r[2]) for r in solo])
    summary("clans + no-clan families distinct_core_10pct", units_core + [len(core(r[2], r[0])) for r in solo])
    print(f"# done {time.time() - t0:.0f}s")


def selftest():
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        q = Path(tmp, "q")
        q.mkdir()
        # family 7: 10 proteins; PF00001 on 2 (20%), PF00002 + same-clan overlapping PF00003 on 1 (one kept), 7 bare
        lines = ["1\t\t[['PF00001', 1e-10, 50.0, 1, 90, 5, 100]]",
                 "2\t\t[['PF00001', 2e-10, 40.0, 1, 90, 5, 100], ['PF00002', 1e-20, 80.0, 1, 50, 120, 170], "
                 "['PF00003', 1e-05, 20.0, 1, 40, 130, 160]]"] + [f"{i}\t\t" for i in range(3, 11)]
        Path(q, "7.tsv").write_text("\n".join(lines) + "\n")
        Path(q, "8.tsv").write_text("11\t\t\n")
        clans = Path(tmp, "pfam.tsv.gz")
        with gzip.open(clans, "wt") as h:
            h.write("PF00001\t\t\tA\t\nPF00002\tCL1\t\tB\t\nPF00003\tCL1\t\tC\t\n")
        mclans = Path(tmp, "clans.tsv.gz")
        with gzip.open(mclans, "wt") as h:
            h.write("clan_id\tfamily_id\tis_representative\nc1\t7\tTrue\nc1\t8\tFalse\n")
        _CLAN_OF.update(read_clans(clans)[0])
        fam, n, w, per_acc = family_counts(Path(q, "7.tsv"))
        assert (fam, n, w) == ("7", 10, 2), (fam, n, w)
        assert per_acc == Counter({"PF00001": 2, "PF00002": 1}), per_acc
        assert core(per_acc, n) == {"PF00001", "PF00002"}  # 1 of 10 is exactly 10%: kept (at least)
        out = Path(tmp, "o.tsv.gz")
        main([str(q), "--clans", str(clans), "--mgnifam-clans", str(mclans), "--out", str(out), "--cpus", "1"])
        rows = gzip.open(out, "rt").read().splitlines()
        assert rows[1].split("\t")[:6] == ["7", "10", "2", "20", "2", "2"], rows[1]
        assert rows[2].split("\t")[:6] == ["8", "1", "0", "0", "0", "0"], rows[2]
    print("selftest ok")


if __name__ == "__main__":
    selftest() if sys.argv[1:] == ["--selftest"] else main()
