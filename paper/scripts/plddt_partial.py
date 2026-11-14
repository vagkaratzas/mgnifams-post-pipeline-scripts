#!/usr/bin/env python3
"""Model confidence by representative completeness: median (IQR) mean pLDDT and median pTM of the representative's
ESMFold model, for families whose representative MGnify Proteins flags full-length vs partial, split by match status.

Inputs: rep_resolve.py --all output (family_id, ..., full_length True/False/"" when missing) and the catalogue_probe.py
--dump output (id, ..., plddt, ptm, any_match).

  python -I plddt_partial.py rep_full_length_<rel>.tsv catalogue_<rel>.csv.gz
  python -I plddt_partial.py --selftest
"""
import csv
import gzip
import statistics
import sys

LABEL = {"True": "full-length", "False": "partial", "": "not in MGnify Proteins"}
MATCH = {"1": "any match", "0": "no match"}


def summarise(full_length, catalogue_rows):
    groups = {}
    for r in catalogue_rows:
        fl = full_length[r["id"]]  # KeyError on a family missing from rep_resolve: the two files must cover one release
        for key in ((fl, "all"), (fl, r["any_match"])):
            groups.setdefault(key, []).append((float(r["plddt"]), float(r["ptm"])))
    out = []
    for (fl, m), v in sorted(groups.items()):
        p = [a for a, _ in v]
        q = statistics.quantiles(p, n=4) if len(p) > 1 else [p[0]] * 3
        out.append((LABEL[fl], MATCH.get(m, m), len(v), statistics.median(p), q[0], q[2], statistics.median(b for _, b in v)))
    return out


def main(rep_path, catalogue_path):
    with open(rep_path, newline="") as f:
        full_length = {r["family_id"]: r["full_length"] for r in csv.DictReader(f, delimiter="\t")}
    with gzip.open(catalogue_path, "rt", newline="") as f:
        rows = summarise(full_length, csv.DictReader(f))
    print("representative\tmatch\tfamilies\tmedian_plddt\tq1\tq3\tmedian_ptm")
    for lab, m, n, med, q1, q3, ptm in rows:
        print(f"{lab}\t{m}\t{n}\t{med:.1f}\t{q1:.1f}\t{q3:.1f}\t{ptm:.2f}")


def selftest():
    fl = {"1": "True", "2": "False", "3": "True", "4": ""}
    cat = [{"id": "1", "plddt": "40", "ptm": "0.2", "any_match": "0"},
           {"id": "2", "plddt": "50", "ptm": "0.3", "any_match": "0"},
           {"id": "3", "plddt": "80", "ptm": "0.7", "any_match": "1"},
           {"id": "4", "plddt": "60", "ptm": "0.5", "any_match": "1"}]
    got = {(lab, m): (n, med) for lab, m, n, med, *_ in summarise(fl, cat)}
    assert got[("full-length", "all")] == (2, 60.0), got
    assert got[("full-length", "no match")] == (1, 40.0), got
    assert got[("partial", "no match")] == (1, 50.0), got
    assert got[("not in MGnify Proteins", "any match")] == (1, 60.0), got
    print("selftest ok")


if __name__ == "__main__":
    if sys.argv[1:] == ["--selftest"]:
        selftest()
    elif len(sys.argv) == 3:
        main(*sys.argv[1:])
    else:
        sys.exit(__doc__)
