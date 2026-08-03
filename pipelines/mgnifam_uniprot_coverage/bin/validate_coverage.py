#!/usr/bin/env python3
"""Internal-consistency checks over the reduced coverage tables.

Every number this pipeline emits is destined for a paper, so the failure mode
that matters is not a crash -- it is a plausible-looking wrong number. These
checks exist to turn each such failure into a non-zero exit.

Input files are the per-view reduced tables named <id>_<pass>_<coords>.reduced.tsv,
which is where the view identity comes from.
"""

import argparse
import csv
import sys
from pathlib import Path


def load_views(paths):
    """{(id, pass, coords): {category: row}} keyed off the filename."""
    views = {}
    for path in paths:
        stem = Path(path).name.split(".")[0]
        view_id, pass_name, coords = stem.rsplit("_", 2)
        with open(path) as fh:
            rows = {r["category"]: r for r in csv.DictReader(fh, delimiter="\t")}
        views[(view_id, pass_name, coords)] = rows
    return views


def num(row, field):
    value = row.get(field, "")
    return float(value) if value not in ("", None) else 0.0


class Checker:
    def __init__(self):
        self.lines = []
        self.failed = 0

    def check(self, ok, label, detail=""):
        status = "PASS" if ok else "FAIL"
        if not ok:
            self.failed += 1
        self.lines.append("%-4s  %s%s" % (status, label, "  -- %s" % detail if detail else ""))

    def note(self, label, detail=""):
        self.lines.append("NOTE  %s%s" % (label, "  -- %s" % detail if detail else ""))


def check_exclusive_within_total(chk, views):
    """Exclusive coverage is a subset of total coverage, category by category."""
    for (view_id, pass_name, coords), rows in sorted(views.items()):
        if pass_name != "exclusive":
            continue
        total = views.get((view_id, "total", coords))
        if not total:
            continue
        for cat, row in sorted(rows.items()):
            base = total.get(cat)
            if not base:
                chk.check(False, "%s/%s exclusive category absent from total" % (view_id, coords), cat)
                continue
            chk.check(num(row, "n_residues") <= num(base, "n_residues"),
                      "%s/%s/%s exclusive residues <= total" % (view_id, coords, cat),
                      "%s vs %s" % (row["n_residues"], base["n_residues"]))
            chk.check(num(row, "n_targets") <= num(base, "n_targets"),
                      "%s/%s/%s exclusive targets <= total" % (view_id, coords, cat),
                      "%s vs %s" % (row["n_targets"], base["n_targets"]))


def check_any_dominates(chk, views):
    """'any' is the union of all categories, so no category may exceed it."""
    for (view_id, pass_name, coords), rows in sorted(views.items()):
        any_row = rows.get("any")
        if not any_row:
            chk.check(False, "%s/%s/%s has an 'any' row" % (view_id, pass_name, coords))
            continue
        for cat, row in sorted(rows.items()):
            if cat == "any":
                continue
            chk.check(num(row, "n_residues") <= num(any_row, "n_residues"),
                      "%s/%s/%s/%s residues <= any" % (view_id, pass_name, coords, cat))
            chk.check(num(row, "n_targets") <= num(any_row, "n_targets"),
                      "%s/%s/%s/%s targets <= any" % (view_id, pass_name, coords, cat))


def check_unannotated(chk, views):
    """Newly annotated sequences are a subset of hit sequences, and are a
    presence/absence quantity -- so identical under ali and env coordinates.
    A difference there means a real bug, not a modelling choice."""
    for (view_id, pass_name, coords), rows in sorted(views.items()):
        if pass_name != "exclusive":
            continue
        for cat, row in sorted(rows.items()):
            chk.check(num(row, "n_seqs_unannotated_before") <= num(row, "n_targets"),
                      "%s/%s/%s newly annotated <= targets" % (view_id, coords, cat))

    by_coords = {}
    for (view_id, pass_name, coords), rows in views.items():
        if pass_name != "exclusive":
            continue
        for cat, row in rows.items():
            by_coords.setdefault((view_id, cat), {})[coords] = row["n_seqs_unannotated_before"]
    for (view_id, cat), seen in sorted(by_coords.items()):
        if len(seen) > 1:
            chk.check(len(set(seen.values())) == 1,
                      "%s/%s newly annotated is coordinate-independent" % (view_id, cat),
                      str(seen))


def check_pooling(chk, views):
    """The pooled UniProtKB view must be the sum of its subsets: hmmsearch
    chunks partition the target database, so no sequence is in two subsets."""
    subsets = sorted({v for v, _p, _c in views if v != "uniprotkb"})
    if not subsets:
        return
    for (view_id, pass_name, coords), rows in sorted(views.items()):
        if view_id != "uniprotkb":
            continue
        for cat, row in sorted(rows.items()):
            for field in ("n_targets", "n_residues"):
                parts = sum(num(views[(s, pass_name, coords)].get(cat, {}), field)
                            for s in subsets if (s, pass_name, coords) in views)
                chk.check(abs(num(row, field) - parts) < 1e-6,
                          "uniprotkb/%s/%s/%s %s == sum of subsets" % (pass_name, coords, cat, field),
                          "%s vs %s" % (row[field], parts))


def check_reference(chk, views, reference):
    """Golden reconciliation against the annotation_percentages pipeline.

    This is the one check that exercises the whole DAG end to end: if the
    SwissProt ali numbers land on the independently computed CSV to the
    residue, every stage in between is doing what it claims.
    """
    with open(reference) as fh:
        ref = next(csv.DictReader(fh))

    want_res = int(ref["annotated_amino_acids_after"]) - int(ref["annotated_amino_acids_before"])
    want_seq = int(ref["annotated_sequences_after"]) - int(ref["annotated_sequences_before"])
    want_pp_res = float(ref["annotated_amino_acid_percentage_point_increase"])
    want_pp_seq = float(ref["annotated_sequence_percentage_point_increase"])
    want_pfam_res = int(ref["annotated_amino_acids_before"])

    excl = views.get(("swissprot", "exclusive", "ali"))
    pfam = views.get(("swissprot", "pfam", "ali"))
    if not excl or not pfam:
        chk.check(False, "reference reconciliation needs the swissprot ali views")
        return

    row = excl["any"]
    chk.check(int(row["n_residues"]) == want_res,
              "reference: exclusive residues", "%s vs %s" % (row["n_residues"], want_res))
    chk.check(int(row["n_seqs_unannotated_before"]) == want_seq,
              "reference: newly annotated sequences",
              "%s vs %s" % (row["n_seqs_unannotated_before"], want_seq))
    chk.check(abs(float(row["pct_db_residues"]) - want_pp_res) < 0.01,
              "reference: residue percentage-point gain",
              "%s vs %s" % (row["pct_db_residues"], want_pp_res))
    chk.check(abs(float(row["pp_gain_sequences"]) - want_pp_seq) < 0.01,
              "reference: sequence percentage-point gain",
              "%s vs %s" % (row["pp_gain_sequences"], want_pp_seq))
    chk.check(int(pfam["any"]["n_residues"]) == want_pfam_res,
              "reference: Pfam residues", "%s vs %s" % (pfam["any"]["n_residues"], want_pfam_res))

    # A partial chunk set reduces cleanly and just understates everything, so
    # compare the Pfam share of the database against the known value. Warn only:
    # a legitimately different release would move it too.
    got = float(pfam["any"]["pct_db_residues"])
    want = float(ref["annotated_amino_acid_percentage_before"])
    if abs(got - want) > 1.0:
        chk.note("Pfam covers %.2f%% of the database, reference says %.2f%% -- "
                 "the signature of a partial chunk set" % (got, want))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reduced", nargs="+", required=True,
                    help="<id>_<pass>_<coords>.reduced.tsv files")
    ap.add_argument("--reference", default=None,
                    help="annotation_percentage_increase.csv to reconcile against")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    views = load_views(args.reduced)
    chk = Checker()
    chk.note("views", ", ".join("/".join(k) for k in sorted(views)))

    check_exclusive_within_total(chk, views)
    check_any_dominates(chk, views)
    check_unannotated(chk, views)
    check_pooling(chk, views)
    if args.reference:
        check_reference(chk, views, args.reference)
    else:
        chk.note("no --reference given; skipping the golden reconciliation")

    body = "\n".join(chk.lines)
    Path(args.output).write_text(body + "\n")
    print(body)

    if chk.failed:
        sys.exit("%d coverage validation check(s) failed" % chk.failed)


if __name__ == "__main__":
    main()
