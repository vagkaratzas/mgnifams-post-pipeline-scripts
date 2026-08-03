#!/usr/bin/env python3
"""Render the coverage tables into a publication-ready report.

Produces the same four tables plus a filled prose block as both Markdown and a
self-contained HTML page. The prose block exists so the paper's numbers are
transcribed by a program rather than by hand -- every figure in it is traceable
to a cell in one of the tables above it.

Primary coordinates are alignment ('ali'), which is what the published figures
use; envelope ('env') is carried alongside as the sensitivity check.
"""

import argparse
import csv
import html
from pathlib import Path

PRIMARY_COORDS = "ali"
VIEW_LABELS = {"swissprot": "SwissProt", "trembl": "TrEMBL", "uniprotkb": "UniProtKB"}
# categories that make up the "transmembrane or disordered" claim, plus the
# complement used for the "excluding these families" sentence
ENRICH_CATS = ["membrane_a", "membrane_b", "tm_any", "disorder", "tm_or_disorder"]


def load_views(paths):
    views = {}
    for path in paths:
        stem = Path(path).name.split(".")[0]
        view_id, pass_name, coords = stem.rsplit("_", 2)
        with open(path) as fh:
            views[(view_id, pass_name, coords)] = {
                r["category"]: r for r in csv.DictReader(fh, delimiter="\t")
            }
    return views


def load_totals(samplesheet):
    """Whole-database denominators per view, pooled view included."""
    totals = {}
    with open(samplesheet) as fh:
        for row in csv.DictReader(fh):
            totals[row["subset"]] = (int(row["total_sequences"]), int(row["total_residues"]))
    totals["uniprotkb"] = (sum(s for s, _r in totals.values()),
                           sum(r for _s, r in totals.values()))
    return totals


def load_library_sizes(path):
    with open(path) as fh:
        return {r["category"]: int(r["n_families"]) for r in csv.DictReader(fh, delimiter="\t")}


def get(views, view_id, pass_name, coords, cat="any", field=None):
    row = views.get((view_id, pass_name, coords), {}).get(cat)
    if row is None:
        return None
    return row if field is None else float(row[field] or 0)


def fmt(value, places=0):
    if value is None:
        return "n/a"
    return "{:,.{p}f}".format(value, p=places)


def order_views(views):
    present = {v for v, _p, _c in views}
    return [v for v in ("swissprot", "trembl", "uniprotkb") if v in present]


# ------------------------------------------------------------------- tables

def table_overview(views, totals, coords):
    header = ["database", "sequences", "residues", "Pfam residues", "Pfam %",
              "MGnifam residues", "MGnifam %", "MGnifam-exclusive residues",
              "pp gain (residues)", "relative increase %",
              "newly annotated sequences", "pp gain (sequences)"]
    rows = []
    for view in order_views(views):
        seqs, res = totals.get(view, (None, None))
        pfam_res = get(views, view, "pfam", coords, field="n_residues")
        mg_res = get(views, view, "total", coords, field="n_residues")
        ex_res = get(views, view, "exclusive", coords, field="n_residues")
        ex_pp = get(views, view, "exclusive", coords, field="pct_db_residues")
        new_seq = get(views, view, "exclusive", coords, field="n_seqs_unannotated_before")
        new_pp = get(views, view, "exclusive", coords, field="pp_gain_sequences")
        rows.append([
            VIEW_LABELS.get(view, view), fmt(seqs), fmt(res),
            fmt(pfam_res), fmt(100.0 * pfam_res / res, 2) if pfam_res and res else "n/a",
            fmt(mg_res), fmt(100.0 * mg_res / res, 2) if mg_res and res else "n/a",
            fmt(ex_res), fmt(ex_pp, 4),
            fmt(100.0 * ex_res / pfam_res, 2) if ex_res and pfam_res else "n/a",
            fmt(new_seq), fmt(new_pp, 4),
        ])
    return header, rows


def table_by_category(views, coords):
    header = ["database", "category", "families with hits", "sequences",
              "exclusive residues", "% of exclusive", "% of database residues",
              "newly annotated sequences"]
    rows = []
    for view in order_views(views):
        cats = views.get((view, "exclusive", coords), {})
        for cat in sorted(cats, key=lambda c: (c == "any", c)):
            row = cats[cat]
            rows.append([
                VIEW_LABELS.get(view, view), cat, fmt(float(row["families_with_hits"])),
                fmt(float(row["n_targets"])), fmt(float(row["n_residues"])),
                fmt(float(row["pct_residues"]), 2), fmt(float(row["pct_db_residues"]), 4),
                fmt(float(row["n_seqs_unannotated_before"])),
            ])
    return header, rows


def table_enrichment(views, sizes, coords, view="uniprotkb"):
    """Share of exclusive residues against share of the library.

    tm_or_disorder is a real union of the family sets, so this row is not the
    sum of the membrane and disorder rows -- summing them would double count
    residues wherever two such families overlap on the same protein.
    """
    library = sizes.get("__library__", 0)
    header = ["category", "families in library", "% of library",
              "% of exclusive residues", "fold enrichment"]
    rows = []
    cats = views.get((view, "exclusive", coords), {})
    for cat in ENRICH_CATS + ["not_tm_disorder"]:
        row = cats.get(cat)
        n_fam = sizes.get(cat)
        if row is None or n_fam is None or not library:
            continue
        share_lib = 100.0 * n_fam / library
        share_res = float(row["pct_residues"])
        rows.append([cat, fmt(n_fam), fmt(share_lib, 2), fmt(share_res, 2),
                     fmt(share_res / share_lib, 2) if share_lib else "n/a"])
    return header, rows


def table_sensitivity(views):
    header = ["database", "metric", "ali", "env"]
    rows = []
    metrics = [("MGnifam-exclusive residues", "exclusive", "n_residues", 0),
               ("pp gain (residues)", "exclusive", "pct_db_residues", 4),
               ("newly annotated sequences", "exclusive", "n_seqs_unannotated_before", 0),
               ("pp gain (sequences)", "exclusive", "pp_gain_sequences", 4),
               ("Pfam residues", "pfam", "n_residues", 0)]
    for view in order_views(views):
        for label, pass_name, field, places in metrics:
            ali = get(views, view, pass_name, "ali", field=field)
            env = get(views, view, pass_name, "env", field=field)
            if ali is None and env is None:
                continue
            rows.append([VIEW_LABELS.get(view, view), label, fmt(ali, places), fmt(env, places)])
    return header, rows


# -------------------------------------------------------------------- prose

def prose(views, totals, sizes, coords):
    """The paper paragraph, with every slot filled from the tables above."""
    out = []
    library = sizes.get("__library__", 0)

    for view in order_views(views):
        seqs, res = totals.get(view, (None, None))
        pfam_res = get(views, view, "pfam", coords, field="n_residues")
        pfam_seq = get(views, view, "pfam", coords, field="n_targets")
        mg_res = get(views, view, "total", coords, field="n_residues")
        ex = views.get((view, "exclusive", coords), {}).get("any")
        if not (res and pfam_res and mg_res and ex):
            continue
        ex_res = float(ex["n_residues"])
        out.append(
            "Of the {seqs} entries ({res} aa) in {label}, Pfam covered {pfam} residues "
            "({pfam_pct}%). MGnifams covered {mg} residues ({mg_pct}%), of which {ex} "
            "-- {ex_share}% of the MGnifam total -- fell outside any Pfam envelope. Total "
            "annotated residue space therefore rose from {pfam_pct}% to {after_pct}%, a "
            "{pp} percentage-point gain and a {rel}% relative increase over Pfam alone. At "
            "the sequence level, {pfam_seq} entries ({pfam_seq_pct}%) carried at least one "
            "Pfam hit; MGnifams added {new_seq} entries with no prior Pfam annotation "
            "(+{new_pp} percentage points), and contributed at least one residue of novel "
            "coverage to {ex_seq} entries ({ex_seq_pct}%).".format(
                label=VIEW_LABELS.get(view, view),
                seqs=fmt(seqs), res=fmt(res),
                pfam=fmt(pfam_res), pfam_pct=fmt(100.0 * pfam_res / res, 2),
                mg=fmt(mg_res), mg_pct=fmt(100.0 * mg_res / res, 2),
                ex=fmt(ex_res), ex_share=fmt(100.0 * ex_res / mg_res, 1),
                after_pct=fmt(100.0 * (pfam_res + ex_res) / res, 2),
                pp=fmt(float(ex["pct_db_residues"]), 2),
                rel=fmt(100.0 * ex_res / pfam_res, 2),
                pfam_seq=fmt(pfam_seq), pfam_seq_pct=fmt(100.0 * pfam_seq / seqs, 2),
                new_seq=fmt(float(ex["n_seqs_unannotated_before"])),
                new_pp=fmt(float(ex["pp_gain_sequences"]), 2),
                ex_seq=fmt(float(ex["n_targets"])),
                ex_seq_pct=fmt(100.0 * float(ex["n_targets"]) / seqs, 2)))

    cats = views.get(("uniprotkb", "exclusive", coords), {})
    if cats and library:
        def share(cat):
            return float(cats[cat]["pct_residues"]) if cat in cats else None

        tmd = share("tm_or_disorder")
        lib_share = 100.0 * sizes.get("tm_or_disorder", 0) / library
        rest = cats.get("not_tm_disorder")
        if tmd is not None and rest is not None:
            out.append(
                "Of the MGnifam-exclusive residues across UniProtKB, {tm}% were assigned by "
                "families flagged as transmembrane ({ma}% membrane-alpha, {mb}% membrane-beta) "
                "and {dis}% by predicted-disordered families. Taken as a union these are {tmd}% "
                "of exclusive residues against {lib}% of the library, an enrichment of "
                "{fold}-fold. Excluding these families, MGnifams contribute {rest_res} residues "
                "beyond Pfam, a {rest_pp} percentage-point gain in UniProtKB residue "
                "coverage.".format(
                    tm=fmt(share("tm_any"), 1), ma=fmt(share("membrane_a"), 1),
                    mb=fmt(share("membrane_b"), 1), dis=fmt(share("disorder"), 1),
                    tmd=fmt(tmd, 1), lib=fmt(lib_share, 1),
                    fold=fmt(tmd / lib_share, 2) if lib_share else "n/a",
                    rest_res=fmt(float(rest["n_residues"])),
                    rest_pp=fmt(float(rest["pct_db_residues"]), 2)))
    return out


# ------------------------------------------------------------------ render

def md_table(header, rows):
    lines = ["| " + " | ".join(header) + " |",
             "|" + "|".join(["---"] * len(header)) + "|"]
    lines += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(lines)


def html_table(header, rows):
    head = "".join("<th>%s</th>" % html.escape(str(c)) for c in header)
    body = "".join("<tr>%s</tr>" % "".join("<td>%s</td>" % html.escape(str(c)) for c in r)
                   for r in rows)
    return "<table><thead><tr>%s</tr></thead><tbody>%s</tbody></table>" % (head, body)


CSS = """
body{font-family:system-ui,-apple-system,Segoe UI,Helvetica,Arial,sans-serif;
     max-width:1100px;margin:2rem auto;padding:0 1rem;line-height:1.5;color:#111}
table{border-collapse:collapse;margin:1rem 0;font-size:.9rem;width:100%;
      display:block;overflow-x:auto}
th,td{border:1px solid #ccc;padding:.3rem .5rem;text-align:right;white-space:nowrap}
th:first-child,td:first-child,th:nth-child(2),td:nth-child(2){text-align:left}
th{background:#f2f2f2}
pre{background:#f7f7f7;padding:1rem;white-space:pre-wrap;border-left:3px solid #999}
.note{color:#555;font-size:.9rem}
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reduced", nargs="+", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--library-sizes", required=True)
    ap.add_argument("--validation", default=None)
    ap.add_argument("--prefix", default="coverage_report")
    args = ap.parse_args()

    views = load_views(args.reduced)
    totals = load_totals(args.samplesheet)
    sizes = load_library_sizes(args.library_sizes)
    coords = PRIMARY_COORDS if any(c == PRIMARY_COORDS for _v, _p, c in views) else \
        sorted({c for _v, _p, c in views})[0]

    tables = [
        ("Table 1. Coverage by database (%s coordinates)" % coords, table_overview(views, totals, coords)),
        ("Table 2. MGnifam-exclusive coverage by category (%s coordinates)" % coords, table_by_category(views, coords)),
        ("Table 3. Category enrichment among MGnifam-exclusive residues, UniProtKB", table_enrichment(views, sizes, coords)),
        ("Table 4. Coordinate sensitivity (alignment vs envelope)", table_sensitivity(views)),
    ]
    paragraphs = prose(views, totals, sizes, coords)
    validation = Path(args.validation).read_text() if args.validation else ""

    md = ["# MGnifam x UniProt coverage", "",
          "Primary coordinates: **%s**. Every number below is generated; do not edit "
          "this file by hand." % coords, ""]
    for title, (header, rows) in tables:
        md += ["## " + title, "", md_table(header, rows), ""]
    md += ["## Generated prose", "",
           "Numbers transcribed by the pipeline, wording to be edited by a human.", ""]
    md += ["> " + p + "\n" for p in paragraphs]
    if validation:
        md += ["## Validation", "", "```", validation.strip(), "```", ""]
    Path(args.prefix + ".md").write_text("\n".join(md) + "\n")

    body = ["<h1>MGnifam &times; UniProt coverage</h1>",
            "<p class='note'>Primary coordinates: <b>%s</b>. Every number below is "
            "generated; do not edit this file by hand.</p>" % html.escape(coords)]
    for title, (header, rows) in tables:
        body += ["<h2>%s</h2>" % html.escape(title), html_table(header, rows)]
    body += ["<h2>Generated prose</h2>"]
    body += ["<p>%s</p>" % html.escape(p) for p in paragraphs]
    if validation:
        body += ["<h2>Validation</h2>", "<pre>%s</pre>" % html.escape(validation.strip())]
    Path(args.prefix + ".html").write_text(
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<title>MGnifam x UniProt coverage</title><style>%s</style></head>"
        "<body>%s</body></html>\n" % (CSS, "".join(body)))


if __name__ == "__main__":
    main()
