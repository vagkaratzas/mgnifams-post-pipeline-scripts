#!/usr/bin/env python3
# Family-level distribution figures for the MGnifams catalogue.
# Every panel is a stacked barplot of annotated vs unannotated (novel) families.
#
# Built to the scientific-publication-plotter standard: plotnine only, vector PDF output,
# one Times text size for every text element, no plot titles (context belongs in the
# manuscript caption), Wong colorblind-safe palette for the two discrete categories.
#
# Modes (--mode, default all):
#   figure        -> <out>/figures/            main-text Figure X
#                    A: family size, log10 bins, full range
#                    B: representative sequence length, 100-aa bins, full range (75-2,000)
#                    C: mean pLDDT, 5-unit bins
#                    each panel as its own PDF, plus the combined, tagged A/B/C
#                    figure_X_family_metadata.pdf
#   supplementary -> <out>/supplementary_figures/
#                    size split small/medium/large, length split short/medium/long
#                    (the fine-grained view behind panels A and B), and pTM
#
# Size uses log10 bins in the main figure because its range (29-1,515,677) cannot be binned
# legibly on a linear axis; the linear small/medium/large split is supplementary.
#
# Also writes <out>/mgnifam_metadata_stats.txt: min/Q1/median/Q3/max per metric, overall
# and split by annotated vs novel -- the numbers quoted in the manuscript paragraph.
#
# Production example (full catalogue):
# python bin/plot_mgnifam_metadata_distributions.py \
#   --metadata mgnifam_codon.csv \
#   --novel-ids assets/mgnifams_v2_results/generate_families/novel/mgnifams_no_annotation_ids.txt \
#   --output-dir output
#
# Local test example (small dataset, adjusted supplementary thresholds):
# python bin/plot_mgnifam_metadata_distributions.py \
#   --metadata assets/mgnifams_v2_results/table_data/mgnifam.csv \
#   --novel-ids assets/mgnifams_v2_results/generate_families/novel/mgnifams_no_annotation_ids.txt \
#   --output-dir output/test \
#   --size-small-max 100 --size-medium-max 1000 \
#   --size-small-bin 20 --size-medium-bin 200 --size-large-bin 1000

import argparse
import functools
import operator
import os
import warnings

import numpy as np
import pandas as pd
import polars as pl
import matplotlib

matplotlib.use("Agg")  # headless: never open a GUI canvas while composing panels
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
from plotnine import (aes, element_blank, element_line, element_rect, element_text, geom_col,
                      geom_text, ggplot, labs, scale_fill_manual, scale_y_continuous, theme,
                      theme_minimal)

# ----------------------------------------------------------------------------------------
# PARAMETER BLOCK -- every tunable that controls how the figures look lives here.
# ----------------------------------------------------------------------------------------
TEXT_SIZE_PT = 8            # one size for EVERY text element: axes, ticks, legend, labels, tags
FONT_STACK = ("Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif")

WIDTH_MM = 180              # double column: ~20 binned categories + in-bar counts need the width
PANEL_HEIGHT_MM = 65        # per panel; the combined figure is n_panels x this
PREVIEW_DPI = 600           # raster preview at journal resolution; the PDF is the deliverable
BAR_WIDTH = 0.75            # fraction of the category slot
LINE_WIDTH_MM = 0.2         # grid lines / bar outlines: thinner than the data
LABEL_PAD_MM = 0.5          # clearance a label needs inside a segment, beyond its own glyphs
STACK_GAP_MM = 1.5          # gap between the labels stacked above a bar, so they read apart
HEADROOM_FRAC = 0.04        # minimum y-axis headroom; the layout grows this to fit above-bar labels

# Discrete color -- Bang Wong colorblind-safe palette, used in its fixed order.
WONG = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#000000"]
ANNOTATED, NOVEL = "Annotated", "Unannotated (novel)"
FILL_COLORS = {ANNOTATED: WONG[0], NOVEL: WONG[1]}  # Wong 1, Wong 2
LABEL_COLOR = "#000000"     # both Wong fills are light enough to carry black text
LEGEND_KEY_PT = 10          # legend glyphs only; legend text stays at TEXT_SIZE_PT

# Continuous color -- viridis family with equal-count (quantile) breaks through the full ramp.
# No panel here encodes a continuous magnitude by color (the fill is the categorical annotation
# status above), so these are declared and unused.
CONTINUOUS_CMAP = "viridis"

# Every number is placed by measuring it: a count sits inside its segment only when the rendered
# glyphs fit the segment in both directions (horizontally if it can, upright if only that fits),
# and anything that does not fit is stacked upright above the bar instead of being dropped or
# spilled across its neighbours. LAYOUT_PASSES re-runs that layout after the y-axis has grown to
# hold the stacks, since taller axes shrink every segment in mm and can push one more label out.
LAYOUT_PASSES = 4
# ----------------------------------------------------------------------------------------


def resolve_font():
    """First available family in FONT_STACK, refreshing matplotlib's cache before giving up."""
    available = {f.name for f in fm.fontManager.ttflist}
    if not available & set(FONT_STACK):
        fm.fontManager = fm._load_fontmanager(try_read_cache=False)
        available = {f.name for f in fm.fontManager.ttflist}

    for family in FONT_STACK:
        if family in available:
            if family != FONT_STACK[0]:
                print(f"Note: '{FONT_STACK[0]}' unavailable, using metric-compatible '{family}'.")
            return family

    raise RuntimeError(f"None of the serif fonts {FONT_STACK} are installed.")


FONT_FAMILY = resolve_font()

# the log10 axis labels are mathtext ($10^{1}$); without this they would render in matplotlib's
# default DejaVu Sans and embed a second family in the PDF
matplotlib.rcParams.update({
    "mathtext.fontset": "custom",
    "mathtext.rm": FONT_FAMILY,
    "mathtext.it": f"{FONT_FAMILY}:italic",
    "mathtext.bf": f"{FONT_FAMILY}:bold",
})


def load_metadata(csv_path, novel_ids_path):
    with open(novel_ids_path) as f:
        novel_ids = {line.strip() for line in f if line.strip()}

    df = pl.read_csv(
        csv_path,
        columns=["id", "full_size", "rep_length", "plddt", "ptm"],
        comment_prefix="#",
    ).rename({"full_size": "size", "rep_length": "length"})

    return df.with_columns(
        ~pl.col("id").cast(pl.Utf8).is_in(novel_ids).alias("annotated")
    )


def summarise(handle, label, values):
    if len(values) == 0:
        return
    q1, med, q3 = np.quantile(values, [0.25, 0.5, 0.75])
    line = (f"{label}: n={len(values):,}  min={values.min():,.4g}  Q1={q1:,.4g}  "
            f"median={med:,.4g}  Q3={q3:,.4g}  max={values.max():,.4g}  mean={values.mean():,.4g}")
    print(line)
    handle.write(line + "\n")


def write_stats(df, output_path):
    n_total = len(df)
    n_annotated = int(df["annotated"].sum())
    n_novel = n_total - n_annotated

    with open(output_path, "w") as fh:
        header = (f"Total families: {n_total:,}\n"
                  f"Annotated (>=1 of Pfam/FunFam/AlphaFoldDB/PDB): {n_annotated:,} "
                  f"({100 * n_annotated / n_total:.1f}%)\n"
                  f"Unannotated / novel: {n_novel:,} ({100 * n_novel / n_total:.1f}%)\n")
        print(header, end="")
        fh.write(header)

        for col, label in [("size", "Family size (full-alignment sequences)"),
                           ("length", "Representative length (aa / HMM match states)"),
                           ("plddt", "mean pLDDT"),
                           ("ptm", "pTM")]:
            fh.write("\n")
            print()
            summarise(fh, f"{label} [all]", df[col].drop_nulls().to_numpy())
            summarise(fh, f"{label} [annotated]",
                      df.filter(pl.col("annotated"))[col].drop_nulls().to_numpy())
            summarise(fh, f"{label} [novel]",
                      df.filter(~pl.col("annotated"))[col].drop_nulls().to_numpy())

    print(f"\nSaved: {output_path}")


def bin_counts(df, value_col, bin_size, log=False, right_closed=False,
               lower_bound=None, upper_bound=None):
    """Bin one column, count annotated / novel families per bin, return a tidy pandas frame.

    log=True bins by decade (10^e <= v < 10^(e+1)) instead of by a fixed width, so a range
    spanning five orders of magnitude fits on one legible linear axis.

    right_closed=True bins integers as (start, start+bin_size] instead of [start, start+bin_size),
    so a column whose maximum is a round number closes exactly on it -- length runs to
    "1,901-2,000", not into a phantom "2,000-2,099" bin.

    lower_bound / upper_bound clip the printed interval labels to the range the subset really
    covers, so the first length bin reads "75-100" rather than "1-100". They default to the
    subset's own min / max.
    """
    values = df[value_col].drop_nulls()
    if values.is_empty():
        return pd.DataFrame()

    lo = values.min() if lower_bound is None else lower_bound
    hi = values.max() if upper_bound is None else upper_bound

    if log:
        df = df.filter(pl.col(value_col) > 0).with_columns(
            pl.col(value_col).log10().floor().cast(pl.Int64).alias("start"))
        label_of = lambda s: f"$10^{{{s}}}$-$10^{{{s + 1}}}$"  # noqa: E731
    elif right_closed:
        df = df.with_columns((((pl.col(value_col) - 1) // bin_size) * bin_size + 1).alias("start"))
        label_of = lambda s: f"{max(s, lo):,}-{min(s + bin_size - 1, hi):,}"  # noqa: E731
    elif df.schema[value_col] in (pl.Float32, pl.Float64):
        df = df.with_columns(((pl.col(value_col) // bin_size) * bin_size).alias("start"))
        label_of = lambda s: f"{s:g}-{s + bin_size:g}"  # noqa: E731
    else:
        df = df.with_columns(((pl.col(value_col) // bin_size) * bin_size).alias("start"))
        label_of = lambda s: f"{max(s, lo):,}-{min(s + bin_size - 1, hi):,}"  # noqa: E731

    agg = df.group_by(["start", "annotated"]).agg(pl.len().alias("count"))
    starts = sorted(agg["start"].unique().to_list())

    rows = []
    for start in starts:
        sub = agg.filter(pl.col("start") == start)
        counts = {ANNOTATED: sub.filter(pl.col("annotated"))["count"].sum(),
                  NOVEL: sub.filter(~pl.col("annotated"))["count"].sum()}
        total = counts[ANNOTATED] + counts[NOVEL]
        # stack annotated at the bottom: its label sits at half its own height, the novel
        # label above it at half of its own
        rows.append(dict(bin=label_of(start), status=ANNOTATED, count=counts[ANNOTATED],
                         mid=counts[ANNOTATED] / 2, total=total))
        rows.append(dict(bin=label_of(start), status=NOVEL, count=counts[NOVEL],
                         mid=counts[ANNOTATED] + counts[NOVEL] / 2, total=total,
                         pct=100 * counts[NOVEL] / total if total else 0))

    tidy = pd.DataFrame(rows)
    tidy["bin"] = pd.Categorical(tidy["bin"], categories=[label_of(s) for s in starts], ordered=True)
    # ggplot stacks the reverse of the level order, so NOVEL first puts ANNOTATED at the bottom
    tidy["status"] = pd.Categorical(tidy["status"], categories=[NOVEL, ANNOTATED], ordered=True)
    return tidy


@functools.lru_cache(maxsize=None)
def text_size_mm(text):
    """Rendered width and height of `text` in mm, at TEXT_SIZE_PT in the figure font."""
    prop = FontProperties(family=FONT_FAMILY, size=TEXT_SIZE_PT)
    extents = TextPath((0, 0), text, size=TEXT_SIZE_PT, prop=prop).get_extents()
    return extents.width / 72 * 25.4, extents.height / 72 * 25.4


def panel_area_mm(tidy, x_label, height_mm):
    """Width and height of the drawn plotting area, measured from a throwaway render.

    Label placement is a question about millimetres (does '10,357' fit across a 5 mm bar?),
    so it needs the real panel geometry rather than a guess at how much of the figure the
    rotated tick labels and axis titles consume. Measured with the legend on, which is the
    shortest the panel ever gets, so a label that fits here fits in the legend-less panels too.
    """
    plot = (base_panel(tidy, x_label, legend=True)
            + theme(figure_size=(WIDTH_MM / 25.4, height_mm / 25.4)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig = plot.draw(show=False)
        fig.canvas.draw()  # the panel only shrinks to fit the tick labels once it is rendered
    box = fig.axes[0].get_position()
    plt.close(fig)
    return box.width * WIDTH_MM, box.height * height_mm


def layout_labels(tidy, plot_w_mm, plot_h_mm, y_top):
    """Place every count and percentage; return (label rows, y the tallest stack reaches).

    Per bar: each count goes inside its own segment when the text fits there (horizontal if it
    can, upright if only that orientation fits), otherwise it joins the upright stack above the
    bar. The novel percentage always sits on top of that stack -- horizontal when it fits the
    category slot, upright when a horizontal one would run into the neighbouring bar.
    """
    units_per_mm = y_top / plot_h_mm
    bar_mm = BAR_WIDTH * plot_w_mm / tidy["bin"].nunique()

    rows, y_reach = [], 0.0
    for bin_label, group in tidy.groupby("bin", observed=True):
        indexed = group.set_index("status")
        total = indexed["count"].sum()
        stack = []  # labels that did not fit inside, bottom-to-top

        for status in (ANNOTATED, NOVEL):
            segment = indexed.loc[status]
            if segment["count"] == 0:
                continue
            label = f"{segment['count']:,}"
            w, h = text_size_mm(label)
            segment_mm = segment["count"] / units_per_mm
            fits = lambda across, along: (across + LABEL_PAD_MM <= bar_mm  # noqa: E731
                                          and along + LABEL_PAD_MM <= segment_mm)
            if fits(w, h):
                angle = 0
            elif fits(h, w):
                angle = 90
            else:
                stack.append(label)
                continue
            rows.append(dict(bin=bin_label, y=segment["mid"], label=label,
                             angle=angle, va="center"))

        stack.append(f"{indexed.loc[NOVEL, 'pct']:.1f}%")
        y = total + STACK_GAP_MM * units_per_mm
        if len(stack) == 1 and text_size_mm(stack[0])[0] + LABEL_PAD_MM <= bar_mm / BAR_WIDTH:
            w, h = text_size_mm(stack[0])
            rows.append(dict(bin=bin_label, y=y, label=stack[0], angle=0, va="bottom"))
            y_reach = max(y_reach, y + h * units_per_mm)
            continue
        for label in stack:  # upright, so a stack is never wider than the bar it sits on
            w, _ = text_size_mm(label)
            rows.append(dict(bin=bin_label, y=y, label=label, angle=90, va="bottom"))
            y = y + (w + STACK_GAP_MM) * units_per_mm
            y_reach = max(y_reach, y)

    return pd.DataFrame(rows), y_reach


def base_panel(tidy, x_label, tag=None, legend=True, y_top=None):
    """The bars, scales and theme -- everything except the number labels.

    y_top pins the axis top outright rather than expanding by a fraction: the label positions
    are themselves data, so a multiplicative expansion would stretch the axis past the height
    the labels were placed against and squeeze every gap the layout just measured out.
    """
    limits = (0, y_top) if y_top else None
    return (ggplot()
            + geom_col(tidy, aes("bin", "count", fill="status"), width=BAR_WIDTH)
            + scale_fill_manual(values=FILL_COLORS, breaks=[ANNOTATED, NOVEL])
            + scale_y_continuous(limits=limits, expand=(0, 0, 0 if y_top else HEADROOM_FRAC, 0),
                                 labels=lambda breaks: [f"{b:,.0f}" for b in breaks])
            + labs(x=x_label, y="Number of families", tag=tag)
            + theme_minimal()
            + theme(
                text=element_text(family=FONT_FAMILY, size=TEXT_SIZE_PT),
                axis_text_x=element_text(rotation=45, ha="right"),
                plot_title=element_blank(),
                plot_tag=element_text(family=FONT_FAMILY, size=TEXT_SIZE_PT, weight="bold"),
                plot_background=element_rect(fill="none", color="none"),
                panel_background=element_rect(fill="white", color="none"),
                panel_grid_minor=element_blank(),
                panel_grid_major_x=element_blank(),
                panel_grid_major_y=element_line(size=LINE_WIDTH_MM, color="#DDDDDD"),
                legend_position="top" if legend else "none",
                legend_title=element_blank(),
                legend_key_size=LEGEND_KEY_PT,
            ))


def build_panel(tidy, x_label, tag=None, legend=True, height_mm=PANEL_HEIGHT_MM):
    """One stacked barplot, every number placed to fit. No title -- context goes in the caption.

    The y-axis is grown to whatever the tallest above-bar stack needs, and the layout is then
    recomputed on the taller axis (a taller axis means fewer mm per family, which can push one
    more count out of its segment), repeating until the height settles.
    """
    plot_w_mm, plot_h_mm = panel_area_mm(tidy, x_label, height_mm)
    y_max = tidy["total"].max()

    # Iterate to a fixed point, then lay the labels out once more on the axis they are actually
    # drawn against: a label placed for a shorter axis would have its gaps squeezed by the
    # rescale, which is exactly the crowding this layout exists to avoid.
    y_top = y_max * (1 + HEADROOM_FRAC)
    for _ in range(LAYOUT_PASSES):
        _, y_reach = layout_labels(tidy, plot_w_mm, plot_h_mm, y_top)
        settled = y_reach + STACK_GAP_MM * y_top / plot_h_mm
        if settled <= y_top:
            break
        y_top = settled
    labels, _ = layout_labels(tidy, plot_w_mm, plot_h_mm, y_top)

    layers = [geom_text(group, aes("bin", "y", label="label"), size=TEXT_SIZE_PT,
                        color=LABEL_COLOR, angle=angle, va=va, ha="center")
              for (angle, va), group in labels.groupby(["angle", "va"])]

    return functools.reduce(operator.add, layers,
                            base_panel(tidy, x_label, tag, legend, y_top))


def save_figure(make_plot, path_no_ext, height_mm):
    """Vector PDF (the deliverable) plus a journal-resolution raster preview.

    Takes a factory, not a plot: drawing a composition consumes its layout registry, so a
    second save of the same object raises. Each format gets a freshly built figure.

    The size is set through the theme rather than save(width=, height=): a composition of
    panels ignores those arguments and would otherwise fall back to matplotlib's default
    6.4 x 4.8 in canvas, shrinking every panel out of the geometry the labels were placed for.
    """
    size = theme(figure_size=(WIDTH_MM / 25.4, height_mm / 25.4))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for ext in ("pdf", "png"):
            (make_plot() + size).save(f"{path_no_ext}.{ext}", dpi=PREVIEW_DPI, verbose=False)
    print(f"Saved: {path_no_ext}.pdf (+ .png preview)")


def figure_panels(plddt_bin):
    """The three main-text panels: (binning args, x axis label)."""
    return [
        (dict(value_col="size", bin_size=None, log=True),
         "Family size (sequences in full alignment)"),
        (dict(value_col="length", bin_size=100, right_closed=True),
         "Representative sequence length (aa)"),
        (dict(value_col="plddt", bin_size=plddt_bin),
         "Mean pLDDT of representative structure"),
    ]


def split_panels(df, col, x_label, name, small_max, medium_max, bins, labels):
    """3-way small/medium/large split of a wide-ranged column, for the supplement.

    Cuts are half-open ([0, small_max), [small_max, medium_max), [medium_max, inf)) so that
    each cut falls on a bin edge -- otherwise the boundary value lands in a bin of its own.
    """
    small_bin, medium_bin, large_bin = bins
    lo, mid, hi = labels
    return [
        (df.filter(pl.col(col) < small_max), dict(value_col=col, bin_size=small_bin,
                                                  upper_bound=small_max - 1),
         x_label, f"{name}_{lo}"),
        (df.filter((pl.col(col) >= small_max) & (pl.col(col) < medium_max)),
         dict(value_col=col, bin_size=medium_bin, lower_bound=small_max,
              upper_bound=medium_max - 1),
         x_label, f"{name}_{mid}"),
        (df.filter(pl.col(col) >= medium_max), dict(value_col=col, bin_size=large_bin,
                                                    lower_bound=medium_max),
         x_label, f"{name}_{hi}"),
    ]


def main():
    parser = argparse.ArgumentParser(
        description="Stacked barplots (annotated vs unannotated) of MGnifam family size, "
                    "representative length, pLDDT and pTM, plus a manuscript stats report.")
    parser.add_argument("--metadata", required=True,
                        help="mgnifam table CSV (id,full_size,rep_length,...,plddt,ptm)")
    parser.add_argument("--novel-ids", required=True,
                        help="TXT of unannotated (novel) family ids, one per line")
    parser.add_argument("--output-dir", default="output",
                        help="Figures go to <dir>/figures and <dir>/supplementary_figures "
                             "(default: output)")
    parser.add_argument("--mode", choices=["all", "figure", "supplementary"], default="all",
                        help="Which figures to produce (default: all)")

    parser.add_argument("--size-small-max", type=int, default=1000)
    parser.add_argument("--size-medium-max", type=int, default=50000)
    parser.add_argument("--size-small-bin", type=int, default=50)
    parser.add_argument("--size-medium-bin", type=int, default=2500)
    parser.add_argument("--size-large-bin", type=int, default=200000)

    parser.add_argument("--length-short-max", type=int, default=300)
    parser.add_argument("--length-medium-max", type=int, default=600)
    parser.add_argument("--length-short-bin", type=int, default=20)
    parser.add_argument("--length-medium-bin", type=int, default=25)
    parser.add_argument("--length-long-bin", type=int, default=100)

    parser.add_argument("--plddt-bin", type=float, default=5.0)
    parser.add_argument("--ptm-bin", type=float, default=0.05)
    args = parser.parse_args()

    fig_dir = os.path.join(args.output_dir, "figures")
    supp_dir = os.path.join(args.output_dir, "supplementary_figures")
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(supp_dir, exist_ok=True)

    df = load_metadata(args.metadata, args.novel_ids)
    write_stats(df, os.path.join(args.output_dir, "mgnifam_metadata_stats.txt"))

    if args.mode in ("all", "figure"):
        specs = list(zip("ABC", figure_panels(args.plddt_bin)))
        tidy_by_tag = {tag: bin_counts(df, **bin_args) for tag, (bin_args, _) in specs}

        for tag, (bin_args, x_label) in specs:
            save_figure(functools.partial(build_panel, tidy_by_tag[tag], x_label),
                        os.path.join(fig_dir, f"figure_X{tag.lower()}_{bin_args['value_col']}"),
                        PANEL_HEIGHT_MM)

        # `/` stacks panels vertically; the legend is identical on all three, so the combined
        # figure carries it once, on A
        def combined():
            return functools.reduce(operator.truediv, [
                build_panel(tidy_by_tag[tag], x_label, tag=tag, legend=tag == "A")
                for tag, (_, x_label) in specs])

        save_figure(combined, os.path.join(fig_dir, "figure_X_family_metadata"),
                    len(specs) * PANEL_HEIGHT_MM)

    if args.mode in ("all", "supplementary"):
        panels = split_panels(
            df, "size", "Family size (sequences in full alignment)", "size",
            args.size_small_max, args.size_medium_max,
            (args.size_small_bin, args.size_medium_bin, args.size_large_bin),
            ("small", "medium", "large"))
        panels += split_panels(
            df, "length", "Representative sequence length (aa)", "length",
            args.length_short_max, args.length_medium_max,
            (args.length_short_bin, args.length_medium_bin, args.length_long_bin),
            ("short", "medium", "long"))
        panels.append((df, dict(value_col="ptm", bin_size=args.ptm_bin),
                       "pTM of representative structure", "ptm"))

        for subset, bin_args, x_label, name in panels:
            tidy = bin_counts(subset, **bin_args)
            if tidy.empty:
                print(f"No families in range for: {name} - skipping.")
                continue
            save_figure(functools.partial(build_panel, tidy, x_label),
                        os.path.join(supp_dir, name), PANEL_HEIGHT_MM)


if __name__ == "__main__":
    main()
