"""Unannotated vs annotated families, release 1.0 and 1.1 (Results catalogue: unannotated-families and pLDDT
paragraphs, Figure 3).

Usage: python3 catalogue_probe.py --old 1.0/output --new 1.1/output --reann 1.1/output_annotate/annotate_not_updated \
           --ids 1.1/output/update_families/families.tsv.gz --prev 1.0/output_db/families.tsv.gz --pfam_median 1089

Evidence as annotation_evidence.py (a row in the table = a match). Two splits:
  unannotated = no Pfam (sequence or profile), FunFam, PDB or AlphaFoldDB match;
  no family match = no Pfam (sequence or profile) or FunFam match (the pLDDT paragraph's Mann-Whitney split).
Size = families.tsv.gz full_size (domain regions in the full alignment; --prev for 1.0); representative length,
pLDDT, pTM from table_data/mgnifam.csv (1.1: updated from --new, the 1,040 from --old, whose structures are 1.0's).
Prints the paragraph numbers, then the Figure 3 bins (size by order of magnitude, length per 100 aa, pLDDT per
5 units) as TSV. With --dump DIR, also writes catalogue_<release>.csv.gz per release (id, full_size, rep_length,
plddt, ptm, any_match, family_match; 1/0 flags as the two splits above): the input of figure3.py. The 1.0 block reproduces the release 1.0 draft figures (3,628 / 14,222 "sequences" there may be another
size unit; 157 / 228 aa; 64.2% / 12.4%; 60.2 / 76.5; 48.2 / 73.5; pTM 0.23 / 0.62; 77.1% / 1.0%). Read-only.
"""
import argparse
import csv
import gzip
import math
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from annotation_evidence import evidence, load, load_models  # noqa: E402

csv.field_size_limit(sys.maxsize)


def sizes(path):
    with gzip.open(path, "rt") as fh:  # column 9 = full_size, as pass2.sh reads it
        rows = csv.reader(fh, delimiter="\t")
        assert next(rows)[8] == "full_size", path
        return {r[0]: int(r[8]) for r in rows}


def scores(path, keep):
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            if r["id"] in keep:
                try:
                    out[r["id"]] = (int(float(r["rep_length"])), float(r["plddt"]), float(r["ptm"]))
                except (ValueError, KeyError):
                    pass
    return out


def q(v):
    q1, med, q3 = statistics.quantiles(v, n=4, method="inclusive")
    return f"{med:.10g} (IQR {q1:.10g}-{q3:.10g}, n = {len(v)})"  # exact: paper numbers are read from here


def mann_whitney(a, b):
    """U of a vs b with tie-corrected normal z; rank-biserial r = 2U/(n1 n2) - 1 (> 0: a larger)."""
    pooled = sorted([(x, 0) for x in a] + [(x, 1) for x in b])
    ranks, ties, i = [0.0] * len(pooled), 0.0, 0
    while i < len(pooled):
        j = i
        while j + 1 < len(pooled) and pooled[j + 1][0] == pooled[i][0]:
            j += 1
        for k in range(i, j + 1):
            ranks[k] = (i + j) / 2 + 1
        t = j - i + 1
        ties += t ** 3 - t
        i = j + 1
    n1, n2 = len(a), len(b)
    r1 = sum(r for r, (_, g) in zip(ranks, pooled) if g == 0)
    u = r1 - n1 * (n1 + 1) / 2
    n = n1 + n2
    sd = math.sqrt(n1 * n2 / 12 * ((n + 1) - ties / (n * (n - 1))))
    z = (u - n1 * n2 / 2) / sd
    return u, z, math.erfc(abs(z) / math.sqrt(2)), 2 * u / (n1 * n2) - 1


def bins(values, key):
    out = {}
    for v in values:
        b = key(v)
        out[b] = out.get(b, 0) + 1
    return out


def share(unann, ann, pred):
    u, a = sum(1 for x in unann if pred(x)), sum(1 for x in ann if pred(x))
    return f"{100 * u / (u + a):.1f}% ({u} of {u + a})" if u + a else "empty"


def dump(path, ids, size, sc, ann, fam):
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "full_size", "rep_length", "plddt", "ptm", "any_match", "family_match"])
        for i in sorted(ids, key=int):
            w.writerow([i, size[i], *sc.get(i, ("", "", "")), int(i in ann), int(i in fam)])


def block(name, ids, e, size, sc, pfam_median, out=None):
    none = ids - (e["family"] | e["structure"])
    ann = ids - none
    nofam = ids - e["family"]
    if out:
        dump(f"{out}/catalogue_{name}.csv.gz", ids, size, sc, ann, e["family"])
    print(f"== {name}: {len(ids)} families; unannotated {len(none)} ({100 * len(none) / len(ids):.1f}%), annotated {len(ann)}; "
          f"no family match {len(nofam)}, family match {len(e['family'])}; missing scores {len(ids - set(sc))}")
    su, sa = [size[i] for i in none], [size[i] for i in ann]
    lu, la = [sc[i][0] for i in none if i in sc], [sc[i][0] for i in ann if i in sc]
    print(f"size: unannotated {q(su)}; annotated {q(sa)}; all {q(su + sa)}, range {min(su + sa)}-{max(su + sa)}")
    print(f"  unannotated median / Pfam median ({pfam_median}): {statistics.median(su) / pfam_median:.2f}")
    print(f"  unannotated share: size < 100 {share(su, sa, lambda x: x < 100)}; > 100,000 {share(su, sa, lambda x: x > 100000)}")
    print(f"rep length: unannotated {q(lu)}; annotated {q(la)}; all {q(lu + la)}, range {min(lu + la)}-{max(lu + la)}")
    pf, pn = [sc[i][1] for i in e["family"] if i in sc], [sc[i][1] for i in nofam if i in sc]
    u, z, p, r = mann_whitney(pf, pn)
    print(f"pLDDT: no family match {q(pn)}; family match {q(pf)}; Mann-Whitney U {u:.0f}, z {z:.2f}, p {p:.3g}, "
          f"rank-biserial {r:.2f}")
    pu, pa = [sc[i][1] for i in none if i in sc], [sc[i][1] for i in ann if i in sc]
    tu, ta = [sc[i][2] for i in none if i in sc], [sc[i][2] for i in ann if i in sc]
    print(f"  unannotated {q(pu)}; annotated {q(pa)}; pTM unannotated {statistics.median(tu):.3g}, annotated {statistics.median(ta):.3g}")
    print(f"  unannotated share: pLDDT 30-35 {share(pu, pa, lambda x: 30 <= x < 35)}; > 85 {share(pu, pa, lambda x: x > 85)}")
    print(f"-- {name} Figure 3 bins (panel, bin, annotated, unannotated)")
    for panel, vu, va, key in (("size", su, sa, lambda x: 10 ** int(math.log10(x))), ("length", lu, la, lambda x: x // 100 * 100),
                               ("plddt", pu, pa, lambda x: int(x // 5 * 5))):
        bu, ba = bins(vu, key), bins(va, key)
        for b in sorted(set(bu) | set(ba)):
            print(f"{panel}\t{b}\t{ba.get(b, 0)}\t{bu.get(b, 0)}")


def selftest():
    u, z, p, r = mann_whitney([3, 4, 5], [1, 2, 3])
    assert u == 8.5 and round(r, 3) == 0.889, (u, r)
    assert bins([5, 50, 99, 100], lambda x: 10 ** int(math.log10(x))) == {1: 1, 10: 2, 100: 1}
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        dump(f"{d}/t.csv.gz", {"10", "2"}, {"10": 5, "2": 7}, {"2": (100, 50.5, 0.3)}, {"2"}, set())
        with gzip.open(f"{d}/t.csv.gz", "rt") as fh:
            assert fh.read().splitlines()[1:] == ["2,7,100,50.5,0.3,1,0", "10,5,,,,0,0"]


def main():
    selftest()
    p = argparse.ArgumentParser()
    for k in ("old", "new", "reann", "ids", "prev"):
        p.add_argument(f"--{k}", required=True)
    p.add_argument("--pfam_median", type=float, default=1089)
    p.add_argument("--dump", help="directory for catalogue_<release>.csv.gz (Figure 3 input)")
    a = p.parse_args()
    with gzip.open(a.ids, "rt") as fh:
        status = {r["family_id"]: r["status"] for r in csv.DictReader(fh, delimiter="\t")}
    ids = set(status)
    upd = {i for i, s in status.items() if s == "updated"}
    nup = ids - upd
    models = load_models(f"{a.old}/table_data")

    d10, _, f10, _ = load(f"{a.old}/table_data", ids)
    sc10 = scores(f"{a.old}/table_data/mgnifam.csv", ids)
    block("1.0", ids, evidence(ids, d10, models, f10), sizes(a.prev), sc10, a.pfam_median, a.dump)

    du, _, fu, _ = load(f"{a.new}/table_data", upd)
    dn, _, fn, _ = load(f"{a.reann}/table_data", nup)
    d11 = {k: du[k] | dn[k] for k in ("pfam", "funfam")}
    sc11 = {**scores(f"{a.new}/table_data/mgnifam.csv", upd), **{i: sc10[i] for i in nup if i in sc10}}
    block("1.1", ids, evidence(ids, d11, models, {**fu, **fn}), sizes(a.ids), sc11, a.pfam_median, a.dump)


if __name__ == "__main__":
    main()
