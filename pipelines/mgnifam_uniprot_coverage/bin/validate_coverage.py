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


def load_reference(path):
    """The *_before/_after totals of one annotation_percentage_increase.csv."""
    with open(path) as fh:
        ref = next(csv.DictReader(fh))
    return {
        "exclusive_residues": int(ref["annotated_amino_acids_after"])
                              - int(ref["annotated_amino_acids_before"]),
        "new_sequences": int(ref["annotated_sequences_after"])
                         - int(ref["annotated_sequences_before"]),
        "pp_residues": float(ref["annotated_amino_acid_percentage_point_increase"]),
        "pp_sequences": float(ref["annotated_sequence_percentage_point_increase"]),
        "pfam_residues": int(ref["annotated_amino_acids_before"]),
        "pfam_pct_residues": float(ref["annotated_amino_acid_percentage_before"]),
    }


def check_reference(chk, views, subset, ref):
    """Golden reconciliation against the annotation_percentages pipeline.

    This is the check that exercises the whole DAG end to end: if a subset's
    ali numbers land on the independently computed CSV to the residue, every
    stage in between is doing what it claims. Alignment coordinates, because
    that is what the reference pipeline uses.
    """
    excl = views.get((subset, "exclusive", "ali"))
    pfam = views.get((subset, "pfam", "ali"))
    if not excl or not pfam:
        chk.check(False, "reference/%s reconciliation needs the ali views" % subset)
        return

    row = excl["any"]
    for label, got, want in [
        ("exclusive residues", int(row["n_residues"]), ref["exclusive_residues"]),
        ("newly annotated sequences", int(row["n_seqs_unannotated_before"]), ref["new_sequences"]),
        ("Pfam residues", int(pfam["any"]["n_residues"]), ref["pfam_residues"]),
    ]:
        chk.check(got == want, "reference/%s: %s" % (subset, label), "%s vs %s" % (got, want))
    for label, got, want in [
        ("residue percentage-point gain", float(row["pct_db_residues"]), ref["pp_residues"]),
        ("sequence percentage-point gain", float(row["pp_gain_sequences"]), ref["pp_sequences"]),
    ]:
        chk.check(abs(got - want) < 0.01, "reference/%s: %s" % (subset, label),
                  "%s vs %s" % (got, want))

    # A partial chunk set reduces cleanly and just understates everything, so
    # compare the Pfam share of the database against the known value. Warn only:
    # a legitimately different release would move it too.
    got = float(pfam["any"]["pct_db_residues"])
    if abs(got - ref["pfam_pct_residues"]) > 1.0:
        chk.note("%s: Pfam covers %.2f%% of the database, reference says %.2f%% -- "
                 "the signature of a partial chunk set" % (subset, got, ref["pfam_pct_residues"]))


def check_pooled_against_references(chk, views, refs):
    """The pooled view against the sum of the per-subset references.

    Only meaningful when every subset brought a reference; otherwise the sum
    is of a partial set and would fail for the wrong reason.
    """
    subsets = sorted({v for v, _p, _c in views if v != "uniprotkb"})
    if not subsets or set(subsets) != set(refs):
        chk.note("pooled reconciliation skipped: references cover %s, views cover %s"
                 % (sorted(refs) or "nothing", subsets))
        return
    pooled = views.get(("uniprotkb", "exclusive", "ali"))
    if not pooled:
        chk.check(False, "pooled reconciliation needs the uniprotkb ali exclusive view")
        return
    for field, key in [("n_residues", "exclusive_residues"),
                       ("n_seqs_unannotated_before", "new_sequences")]:
        want = sum(refs[s][key] for s in subsets)
        got = int(pooled["any"][field])
        chk.check(got == want, "reference/uniprotkb: %s == sum of subset references" % field,
                  "%s vs %s" % (got, want))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reduced", nargs="+", required=True,
                    help="<id>_<pass>_<coords>.reduced.tsv files")
    ap.add_argument("--references", nargs="*", default=[], metavar="SUBSET=CSV",
                    help="per-subset annotation_percentage_increase.csv to "
                         "reconcile against, e.g. swissprot=sprot.csv")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    views = load_views(args.reduced)
    refs = {}
    for spec in args.references:
        subset, _, path = spec.partition("=")
        if not path:
            sys.exit("--references wants SUBSET=CSV, got %r" % spec)
        refs[subset] = load_reference(path)

    chk = Checker()
    chk.note("views", ", ".join("/".join(k) for k in sorted(views)))

    check_exclusive_within_total(chk, views)
    check_any_dominates(chk, views)
    check_unannotated(chk, views)
    check_pooling(chk, views)
    if refs:
        for subset in sorted(refs):
            check_reference(chk, views, subset, refs[subset])
        check_pooled_against_references(chk, views, refs)
    else:
        chk.note("no references given; skipping the golden reconciliation")

    body = "\n".join(chk.lines)
    Path(args.output).write_text(body + "\n")
    print(body)

    if chk.failed:
        sys.exit("%d coverage validation check(s) failed" % chk.failed)


if __name__ == "__main__":
    main()
