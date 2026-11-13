"""Transmembrane / disordered classes of a release (Methods "Classification of transmembrane
and disordered families"; Results catalogue paragraph 2).

Usage: python3 tm_classes.py --old 1.0/table_data/mgnifam.csv --new 1.1/table_data/mgnifam.csv \
           --ids 1.1/update_families/families.tsv.gz

Prints the classes for 1.0 (all rows of --old; must reproduce the published 5,186 / 372 / 989,
+510, 1,777) and for 1.1: updated families from --new, not-updated ones from --old (same
DeepTMHMM / s4pred versions). Also the 1.0 -> 1.1 transitions. Read-only, except
--lists-out DIR: writes the Table 2 category lists (membrane_a.txt, membrane_b.txt, disorder.txt;
one MGYF%010d per line, the form of the UniProtKB coverage pipeline's --lists_dir) into DIR/1.0 and
DIR/1.1. The 1.0 set must equal the lists of the published Table 2 run before 1.1 is used.
"""
import argparse
import csv
import gzip
import os
import re
import sys
from collections import Counter

csv.field_size_limit(sys.maxsize)  # the table carries blob columns
FIELDS = ("rep_length", "membrane_alpha_percent", "membrane_beta_percent", "inside_percent",
          "outside_percent", "coil_percent", "plddt", "ptm")


def classes(r, alpha_long=5.0, plddt_max=50.0):
    """Set of classes of one family row (floats). alpha_long / plddt_max are the sensitivity knobs."""
    L, a, b = r["rep_length"], r["membrane_alpha_percent"], r["membrane_beta_percent"]
    out = set()
    need = 12.0 if 75 <= L <= 200 else alpha_long if L > 200 else None
    if need is not None and a >= need and r["inside_percent"] >= 1 and r["outside_percent"] >= 1:
        out.add("alpha")
    if b >= 20 and L >= 100:
        out.add("beta")
    if r["coil_percent"] >= 80 and r["plddt"] <= plddt_max and r["ptm"] <= 0.5 and a <= 1:
        out.add("disordered")
    return out


def read(path):
    rows, bad = {}, 0
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            try:
                rows[r["id"]] = {k: float(r[k]) for k in FIELDS}
            except (ValueError, TypeError):
                bad += 1
    return rows, bad


def summary(name, rows):
    n = len(rows)
    base = {i: classes(r) for i, r in rows.items()}
    alpha = sum("alpha" in c for c in base.values())
    beta = sum("beta" in c for c in base.values())
    both = sum({"alpha", "beta"} <= c for c in base.values())
    dis = sum("disordered" in c for c in base.values())
    tm_dis = sum(("alpha" in c or "beta" in c) and "disordered" in c for c in base.values())
    add3 = sum("alpha" in classes(r, alpha_long=3.0) and "alpha" not in base[i] for i, r in rows.items())
    dis70 = sum("disordered" in classes(r, plddt_max=70.0) for r in rows.values())
    print(f"== {name}: {n} families")
    print(f"membrane-alpha {alpha}  membrane-beta {beta}  both {both}  transmembrane {alpha + beta - both} "
          f"({100 * (alpha + beta - both) / n:.1f}%)")
    print(f"disordered {dis} ({100 * dis / n:.1f}%)  transmembrane and disordered {tm_dis}")
    print(f"sensitivity: alpha >200 aa at 3% adds {add3}; disorder pLDDT <= 70 gives {dis70} ({100 * dis70 / n:.1f}%)")
    for k in FIELDS:
        v = [r[k] for r in rows.values()]
        print(f"  range {k}: {min(v)} - {max(v)}")
    print(f"  rep_length < 75: {sum(r['rep_length'] < 75 for r in rows.values())}")
    return base


def label(c):
    return "+".join(sorted(c)) or "none"


def mgyf(family_id):
    """'444', 'MGYF0000000444' -> 'MGYF0000000444'."""
    return f"MGYF{int(re.sub(r'^MGYF', '', str(family_id))):010d}"


def write_lists(outdir, per_family):
    """per_family: {id: set of classes} -> membrane_a / membrane_b / disorder .txt, sorted."""
    os.makedirs(outdir, exist_ok=True)
    for name, cls in (("membrane_a", "alpha"), ("membrane_b", "beta"), ("disorder", "disordered")):
        ids = sorted(mgyf(i) for i, c in per_family.items() if cls in c)
        with open(os.path.join(outdir, f"{name}.txt"), "w") as fh:
            fh.write("".join(f"{i}\n" for i in ids))
        print(f"  wrote {outdir}/{name}.txt: {len(ids)}")


def selftest():
    assert mgyf("444") == mgyf("MGYF0000000444") == "MGYF0000000444"
    r = dict(rep_length=180, membrane_alpha_percent=12, membrane_beta_percent=0, inside_percent=1,
             outside_percent=5, coil_percent=10, plddt=80, ptm=0.8)
    assert classes(r) == {"alpha"}
    assert classes({**r, "inside_percent": 0}) == set()
    assert classes({**r, "rep_length": 440, "membrane_alpha_percent": 4}) == set()
    assert classes({**r, "rep_length": 440, "membrane_alpha_percent": 4}, alpha_long=3) == {"alpha"}
    d = {**r, "membrane_alpha_percent": 0, "coil_percent": 85, "plddt": 60, "ptm": 0.3}
    assert classes(d) == set() and classes(d, plddt_max=70) == {"disordered"}
    assert classes({**r, "membrane_alpha_percent": 0, "membrane_beta_percent": 25, "rep_length": 99}) == set()


def main():
    selftest()
    p = argparse.ArgumentParser()
    p.add_argument("--old", required=True)
    p.add_argument("--new", required=True)
    p.add_argument("--ids", required=True)
    p.add_argument("--lists-out")
    a = p.parse_args()
    with gzip.open(a.ids, "rt") as fh:
        status = {r["family_id"]: r["status"] for r in csv.DictReader(fh, delimiter="\t")}
    old, bad_old = read(a.old)
    new, bad_new = read(a.new)
    print(f"rows: old {len(old)} (unparsed {bad_old}), new {len(new)} (unparsed {bad_new}), "
          f"families.tsv {len(status)} {dict(Counter(status.values()))}")
    upd = [i for i, s in status.items() if s == "updated"]
    print(f"updated missing from --new: {sum(i not in new for i in upd)}; "
          f"not updated missing from --old: {sum(i not in old for i, s in status.items() if s != 'updated')}")
    c10 = summary("1.0 (all of --old)", old)
    rel = {i: (new[i] if s == "updated" else old[i]) for i, s in status.items()
           if (new if s == "updated" else old).get(i)}
    c11 = summary("1.1 (updated from --new, not updated from --old)", rel)
    summary("1.1 updated only", {i: new[i] for i in upd if i in new})
    if a.lists_out:
        write_lists(os.path.join(a.lists_out, "1.0"), c10)
        write_lists(os.path.join(a.lists_out, "1.1"), c11)
    t = Counter((label(c10.get(i, set())), label(c11[i])) for i in c11 if i in c10)
    print("== transitions 1.0 -> 1.1 (changed only)")
    for (x, y), n in sorted(t.items(), key=lambda kv: -kv[1]):
        if x != y:
            print(f"  {x} -> {y}: {n}")


if __name__ == "__main__":
    main()
