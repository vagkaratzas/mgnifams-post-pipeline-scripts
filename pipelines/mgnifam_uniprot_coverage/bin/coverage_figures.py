#!/usr/bin/env python3
"""Manuscript figures for the MGnifam x UniProt coverage tables.

Vector PDF via plotnine. Panels 3 and 4 plot the same quantity on both axes
(share vs share, ali vs env), so they are square with a 1:1 reference line;
panels 1 and 2 are flat bars.
"""

import argparse
import csv
from pathlib import Path

import pandas as pd
from plotnine import (
    aes, element_blank, element_line, element_text, geom_abline, geom_bar,
    geom_col, geom_point, geom_text, ggplot, labs, coord_flip, scale_fill_manual,
    scale_colour_manual, scale_x_continuous, scale_y_continuous, theme, theme_bw,
)

# ------------------------------------------------------------ parameters
TEXT_SIZE_PT = 8            # every text element: axes, ticks, legend, strips
LINE_WIDTH = 0.3            # axes, grid, reference lines (mm)
POINT_SIZE = 1.8            # markers
WIDTH_MM = 85               # single column; 180 only if a venue asks
HEIGHT_MM = 55              # flat, except the square same-quantity panels
SQUARE_MM = 85              # width == height for share-vs-share panels
MM_PER_INCH = 25.4

# Bang Wong, colorblind-safe, used in this exact order for <= 8 discrete series
WONG = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7",
        "#D55E00", "#F0E442", "#0072B2", "#000000"]

# Continuous colour, if ever needed here, is viridis-family only (default
# "viridis") mapped through equal-count quantile bands computed from the data.
VIRIDIS_MAP = "viridis"

PRIMARY_COORDS = "ali"
VIEW_ORDER = ["swissprot", "trembl", "uniprotkb"]
VIEW_LABELS = {"swissprot": "SwissProt", "trembl": "TrEMBL", "uniprotkb": "UniProtKB"}


def base_theme():
    return theme_bw() + theme(
        text=element_text(family="Times New Roman", size=TEXT_SIZE_PT),
        axis_title=element_text(size=TEXT_SIZE_PT),
        axis_text=element_text(size=TEXT_SIZE_PT),
        legend_title=element_text(size=TEXT_SIZE_PT),
        legend_text=element_text(size=TEXT_SIZE_PT),
        strip_text=element_text(size=TEXT_SIZE_PT),
        plot_title=element_blank(),
        plot_subtitle=element_blank(),
        panel_grid_minor=element_blank(),
        panel_grid_major=element_line(size=LINE_WIDTH, colour="#DDDDDD"),
        panel_border=element_line(size=LINE_WIDTH, colour="#333333"),
        axis_ticks=element_line(size=LINE_WIDTH),
        legend_key_size=TEXT_SIZE_PT,
    )


def save(plot, path, width_mm=WIDTH_MM, height_mm=HEIGHT_MM):
    plot.save(path, width=width_mm / MM_PER_INCH, height=height_mm / MM_PER_INCH,
              units="in", dpi=600, verbose=False, transparent=True)


def load_frames(paths):
    rows = []
    for path in paths:
        stem = Path(path).name.split(".")[0]
        view, pass_name, coords = stem.rsplit("_", 2)
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                r.update(view=view, pass_name=pass_name, coords=coords)
                rows.append(r)
    frame = pd.DataFrame(rows)
    for col in ("n_targets", "n_residues", "pct_residues", "pct_db_residues",
                "pp_gain_sequences", "n_seqs_unannotated_before", "families_with_hits"):
        if col in frame:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
    return frame


def present_views(frame):
    return [v for v in VIEW_ORDER if v in set(frame["view"])]


# ------------------------------------------------------------- figure 1

def fig_coverage_stack(frame, totals, coords, out):
    """Where every residue of each database sits: Pfam, MGnifam-only, nothing."""
    records = []
    for view in present_views(frame):
        total = totals.get(view)
        pfam = frame.query("view == @view and pass_name == 'pfam' and coords == @coords "
                           "and category == 'any'")["n_residues"]
        excl = frame.query("view == @view and pass_name == 'exclusive' and coords == @coords "
                           "and category == 'any'")["n_residues"]
        if not total or pfam.empty or excl.empty:
            continue
        pfam, excl = float(pfam.iloc[0]), float(excl.iloc[0])
        for label, value in [("Pfam", pfam), ("MGnifam-exclusive", excl),
                             ("unannotated", max(total - pfam - excl, 0))]:
            records.append({"view": VIEW_LABELS[view], "layer": label,
                            "pct": 100.0 * value / total})
    if not records:
        return
    data = pd.DataFrame(records)
    data["layer"] = pd.Categorical(
        data["layer"], categories=["Pfam", "MGnifam-exclusive", "unannotated"])
    data["view"] = pd.Categorical(
        data["view"], categories=[VIEW_LABELS[v] for v in present_views(frame)])

    plot = (ggplot(data, aes("view", "pct", fill="layer"))
            + geom_col(width=0.65)
            + scale_fill_manual(values=WONG[:3], name="")
            + scale_y_continuous(expand=(0, 0), limits=(0, 100.5))
            + labs(x="", y="residues (% of database)")
            + base_theme())
    save(plot, out)


# ------------------------------------------------------------- figure 2

def fig_category_share(frame, coords, out):
    """Which family categories the exclusive residues come from."""
    data = frame.query("pass_name == 'exclusive' and coords == @coords and category != 'any'").copy()
    keep = ["novel", "novel_structure_any", "membrane_a", "membrane_b",
            "disorder", "tm_or_disorder", "not_tm_disorder"]
    data = data[data["category"].isin(keep)]
    if data.empty:
        return
    data["view"] = pd.Categorical(data["view"].map(VIEW_LABELS),
                                  categories=[VIEW_LABELS[v] for v in present_views(frame)])
    data["category"] = pd.Categorical(data["category"], categories=list(reversed(keep)))

    plot = (ggplot(data, aes("category", "pct_residues", fill="view"))
            + geom_col(position="dodge", width=0.7)
            + coord_flip()
            + scale_fill_manual(values=WONG[:data["view"].nunique()], name="")
            + labs(x="", y="share of MGnifam-exclusive residues (%)")
            + base_theme())
    save(plot, out, height_mm=70)


# ------------------------------------------------------------- figure 3

def fig_enrichment(frame, sizes, coords, out, view="uniprotkb"):
    """Share of exclusive residues against share of the library, 1:1 reference.

    Same quantity on both axes, so the panel is square.
    """
    library = sizes.get("__library__", 0)
    data = frame.query("view == @view and pass_name == 'exclusive' and coords == @coords").copy()
    data = data[data["category"].isin(sizes.keys()) & (data["category"] != "any")]
    if data.empty or not library:
        return
    data["lib_share"] = data["category"].map(lambda c: 100.0 * sizes[c] / library)

    limit = max(data["lib_share"].max(), data["pct_residues"].max()) * 1.15
    plot = (ggplot(data, aes("lib_share", "pct_residues"))
            + geom_abline(slope=1, intercept=0, size=LINE_WIDTH, colour="#999999")
            + geom_point(size=POINT_SIZE, colour=WONG[6])
            + geom_text(aes(label="category"), size=TEXT_SIZE_PT, nudge_y=limit * 0.035,
                        family="Times New Roman")
            + scale_x_continuous(limits=(0, limit))
            + scale_y_continuous(limits=(0, limit))
            + labs(x="share of the MGnifam library (%)",
                   y="share of MGnifam-exclusive residues (%)")
            + base_theme())
    save(plot, out, width_mm=SQUARE_MM, height_mm=SQUARE_MM)


# ------------------------------------------------------------- figure 4

def fig_coord_sensitivity(frame, out):
    """Alignment against envelope coordinates, per category and view."""
    data = frame.query("pass_name == 'exclusive'")
    wide = data.pivot_table(index=["view", "category"], columns="coords",
                            values="pct_db_residues").reset_index()
    if "ali" not in wide or "env" not in wide:
        return
    wide = wide.dropna(subset=["ali", "env"])
    if wide.empty:
        return
    wide["view"] = pd.Categorical(wide["view"].map(VIEW_LABELS),
                                  categories=[VIEW_LABELS[v] for v in present_views(frame)])

    limit = max(wide["ali"].max(), wide["env"].max()) * 1.15
    plot = (ggplot(wide, aes("ali", "env", colour="view"))
            + geom_abline(slope=1, intercept=0, size=LINE_WIDTH, colour="#999999")
            + geom_point(size=POINT_SIZE)
            + scale_colour_manual(values=WONG[:wide["view"].nunique()], name="")
            + scale_x_continuous(limits=(0, limit))
            + scale_y_continuous(limits=(0, limit))
            + labs(x="alignment coordinates (pp of database residues)",
                   y="envelope coordinates (pp of database residues)")
            + base_theme())
    save(plot, out, width_mm=SQUARE_MM, height_mm=SQUARE_MM)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reduced", nargs="+", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--library-sizes", required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    frame = load_frames(args.reduced)
    coords = PRIMARY_COORDS if PRIMARY_COORDS in set(frame["coords"]) \
        else sorted(set(frame["coords"]))[0]

    with open(args.samplesheet) as fh:
        totals = {r["subset"]: int(r["total_residues"]) for r in csv.DictReader(fh)}
    totals["uniprotkb"] = sum(totals.values())

    with open(args.library_sizes) as fh:
        sizes = {r["category"]: int(r["n_families"])
                 for r in csv.DictReader(fh, delimiter="\t")}

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    fig_coverage_stack(frame, totals, coords, out / "fig1_coverage_by_database.pdf")
    fig_category_share(frame, coords, out / "fig2_exclusive_by_category.pdf")
    fig_enrichment(frame, sizes, coords, out / "fig3_category_enrichment.pdf")
    fig_coord_sensitivity(frame, out / "fig4_coordinate_sensitivity.pdf")


if __name__ == "__main__":
    main()
