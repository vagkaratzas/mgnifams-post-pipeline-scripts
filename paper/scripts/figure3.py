#!/usr/bin/env python3
"""Figure 3 of the MGnifams paper: family size, representative length and model confidence, as stacked bars of
families with at least one Pfam (sequence or profile), FunFam, PDB or AlphaFoldDB match vs families with none.

Input: catalogue_<release>.csv.gz from `catalogue_probe.py --dump`: id, full_size,
rep_length, plddt, ptm, any_match, family_match. Migrated from mgnifams-post-pipeline-scripts
`bin/plot_mgnifam_metadata_distributions.py` (`4e26598`); built to its scientific-publication-plotter standard:
plotnine only, vector PDF, one Times text size, no plot titles, Wong palette.

Usage (from the repository root):
  python paper/scripts/figure3.py \
      --catalogue paper/data/catalogue_1.1.csv.gz --release 1.1 \
      --probe-log tmp/mgnifams_catalogue_probe_<job>.log
  --selftest runs the binning checks only.

Writes <out>/figure-3_family-metadata.{pdf,png} (panels A size by order of magnitude, B length in 100-aa bins,
C mean pLDDT in 5-unit bins), <supp-out>/figure-s1..s3_*.{pdf,png} (Supplementary Figures S1 linear size
splits, S2 linear length splits, S3 pTM) and <stats-out>/figure-3_stats_<release>.txt (exact quantiles and the
bin counts). With --probe-log, the counts
must equal the probe's "Figure 3 bins" block for the release (same data, the probe's own bin keys), or it exits.
"""
import argparse
import functools
import math
import operator
import os
import re
import sys
import warnings

import matplotlib
import numpy as np
import pandas as pd
import polars as pl

matplotlib.use("Agg")  # headless
import matplotlib.font_manager as fm  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.font_manager import FontProperties  # noqa: E402
from matplotlib.textpath import TextPath  # noqa: E402
from plotnine import (aes, element_blank, element_line, element_rect, element_text, geom_col,  # noqa: E402
                      geom_text, ggplot, labs, scale_fill_manual, scale_y_continuous, theme, theme_minimal)

# ----------------------------------------------------------------------------------------
# PARAMETER BLOCK -- every tunable that controls how the figure looks lives here.
# ----------------------------------------------------------------------------------------
TEXT_SIZE_PT = 8            # one size for EVERY text element: axes, ticks, legend, labels, tags
FONT_STACK = ("Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif")

WIDTH_MM = 180              # double column: ~20 binned categories + in-bar counts need the width
PANEL_HEIGHT_MM = 65        # per panel; the combined figure is n_panels x this
PREVIEW_DPI = 600           # raster preview at journal resolution; the PDF is the deliverable
BAR_WIDTH = 0.75            # fraction of the category slot
LINE_WIDTH_MM = 0.2         # grid lines: thinner than the data
LABEL_PAD_MM = 0.5          # clearance a label needs inside a segment, beyond its own glyphs
STACK_GAP_MM = 1.5          # gap between the labels stacked above a bar
HEADROOM_FRAC = 0.04        # minimum y-axis headroom; the layout grows this to fit above-bar labels

# Discrete color -- Bang Wong colorblind-safe palette, in its fixed order.
WONG = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#000000"]
MATCH, NO_MATCH = "Pfam, FunFam, PDB or AlphaFoldDB match", "No match"
FILL_COLORS = {MATCH: WONG[0], NO_MATCH: WONG[1]}
LABEL_COLOR = "#000000"     # both Wong fills are light enough to carry black text
# A count stacked above its bar is drawn in a dark shade of its segment's fill (>= 4.5:1 on white), so the stack
# reads by colour as well as by order (bottom to top, as in the bar); the percentage stays black
STACK_COLORS = {MATCH: "#8F6000", NO_MATCH: WONG[6]}
TICK_ROTATION = 45          # x tick labels rotate only when the widest one does not fit its category slot
LEGEND_KEY_PT = 10          # legend glyphs only; legend text stays at TEXT_SIZE_PT

# Continuous color -- viridis family with equal-count (quantile) breaks. No panel encodes a magnitude by color.
CONTINUOUS_CMAP = "viridis"

# Labels are placed by measuring them; LAYOUT_PASSES re-runs the layout after the y-axis grows to hold the
# above-bar stacks, since a taller axis shrinks every segment in mm.
LAYOUT_PASSES = 4
DASH = "–"             # en dash in bin ranges
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


def load(path):
    df = pl.read_csv(path, schema_overrides={"plddt": pl.Float64, "ptm": pl.Float64})
    return df.rename({"full_size": "size", "rep_length": "length"}).with_columns(
        (pl.col("any_match") == 1).alias("matched"))


def bin_starts(df, col, kind, width=None):
    """Add `start` (the bin key) to df.

    decade: [10^e, 10^(e+1)), from the digit count, so no float log10 sits on a bin edge.
    right:  integers in (start - 1, start + width - 1], so a maximum of 2,000 closes the 1,901-2,000 bin.
    left:   [start, start + width); floats get an epsilon so 0.3 / 0.05 = 5.999... lands in bin 6.
    """
    c = pl.col(col)
    if kind == "decade":
        start = 10 ** (c.cast(pl.Int64).cast(pl.Utf8).str.len_chars().cast(pl.Int64) - 1)
    elif kind == "right":
        start = ((c - 1) // width) * width + 1
    elif df.schema[col].is_float():
        start = (c / width + 1e-9).floor() * width
    else:
        start = (c // width) * width
    return df.filter(c.is_not_null()).with_columns(start.alias("start"))


def bin_label(start, kind, width, lo, hi, is_float):
    if is_float:  # decimals of the width, so 0.05 bins read 0.10-0.15 and 5.0 bins read 30-35
        d = max(0, -math.floor(math.log10(width) + 1e-9)) if width % 1 else 0
        return f"{start:.{d}f}{DASH}{start + width:.{d}f}"
    end = start * 10 - 1 if kind == "decade" else start + width - 1
    if kind == "decade" and end >= hi and start > lo:  # open top decade: short enough for a horizontal tick
        return f"≥ {start:,}"
    return f"{max(start, lo):,}{DASH}{min(end, hi):,}"


def bin_counts(df, col, kind="left", width=None, lower=None, upper=None):
    """Tidy frame: one row per bin and status, with count, label mid-height and bin total.

    lower / upper clip the printed range to what the subset covers ("75-100", not "1-100").
    """
    values = df[col].drop_nulls()
    if values.is_empty():
        return pd.DataFrame()
    lo = values.min() if lower is None else lower
    hi = values.max() if upper is None else upper
    is_float = df.schema[col].is_float()
    agg = bin_starts(df, col, kind, width).group_by(["start", "matched"]).agg(pl.len().alias("count"))
    starts = sorted(agg["start"].unique().to_list())
    rows = []
    for s in starts:
        sub = agg.filter(pl.col("start") == s)
        m = int(sub.filter(pl.col("matched"))["count"].sum())
        n = int(sub.filter(~pl.col("matched"))["count"].sum())
        label = bin_label(s, kind, width, lo, hi, is_float)
        # matched at the bottom: its label sits at half its height, the no-match label above it
        rows.append(dict(bin=label, start=s, status=MATCH, count=m, mid=m / 2, total=m + n))
        rows.append(dict(bin=label, start=s, status=NO_MATCH, count=n, mid=m + n / 2, total=m + n,
                         pct=100 * n / (m + n)))
    tidy = pd.DataFrame(rows)
    tidy["bin"] = pd.Categorical(tidy["bin"], categories=list(dict.fromkeys(tidy["bin"])), ordered=True)
    # ggplot stacks the reverse of the level order, so NO_MATCH first puts MATCH at the bottom
    tidy["status"] = pd.Categorical(tidy["status"], categories=[NO_MATCH, MATCH], ordered=True)
    return tidy


@functools.lru_cache(maxsize=None)
def text_size_mm(text):
    """Rendered width and height of `text` in mm, at TEXT_SIZE_PT in the figure font."""
    prop = FontProperties(family=FONT_FAMILY, size=TEXT_SIZE_PT)
    extents = TextPath((0, 0), text, size=TEXT_SIZE_PT, prop=prop).get_extents()
    return extents.width / 72 * 25.4, extents.height / 72 * 25.4


def panel_area_mm(tidy, x_label, height_mm):
    """Width and height of the drawn plotting area, measured from a throwaway render (legend on: the
    shortest the panel gets, so a label that fits here fits in the legend-less panels too)."""
    plot = base_panel(tidy, x_label, legend=True) + theme(figure_size=(WIDTH_MM / 25.4, height_mm / 25.4))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig = plot.draw(show=False)
        fig.canvas.draw()
    box = fig.axes[0].get_position()
    plt.close(fig)
    return box.width * WIDTH_MM, box.height * height_mm


def layout_labels(tidy, plot_w_mm, plot_h_mm, y_top):
    """Place every count and percentage; return (label rows, y the tallest stack reaches).

    Each count goes inside its segment when it fits (horizontal if it can, upright if only that fits),
    otherwise into the upright stack above the bar; the no-match percentage tops that stack.
    """
    units_per_mm = y_top / plot_h_mm
    bar_mm = BAR_WIDTH * plot_w_mm / tidy["bin"].nunique()
    rows, y_reach = [], 0.0
    for bin_label_, group in tidy.groupby("bin", observed=True):
        indexed = group.set_index("status")
        total = indexed["count"].sum()
        stack = []
        for status in (MATCH, NO_MATCH):
            segment = indexed.loc[status]
            if segment["count"] == 0:
                continue
            label = f"{segment['count']:,}"
            w, h = text_size_mm(label)
            segment_mm = segment["count"] / units_per_mm

            def fits(across, along):
                return across + LABEL_PAD_MM <= bar_mm and along + LABEL_PAD_MM <= segment_mm

            if fits(w, h):
                angle = 0
            elif fits(h, w):
                angle = 90
            else:
                stack.append((label, STACK_COLORS[status]))
                continue
            rows.append(dict(bin=bin_label_, y=segment["mid"], label=label, angle=angle, va="center",
                             color=LABEL_COLOR))

        stack.append((f"{indexed.loc[NO_MATCH, 'pct']:.1f}%", LABEL_COLOR))
        y = total + STACK_GAP_MM * units_per_mm
        if len(stack) == 1 and text_size_mm(stack[0][0])[0] + LABEL_PAD_MM <= bar_mm / BAR_WIDTH:
            _, h = text_size_mm(stack[0][0])
            rows.append(dict(bin=bin_label_, y=y, label=stack[0][0], angle=0, va="bottom", color=LABEL_COLOR))
            y_reach = max(y_reach, y + h * units_per_mm)
            continue
        for label, color in stack:  # upright, so a stack is never wider than its bar
            w, _ = text_size_mm(label)
            rows.append(dict(bin=bin_label_, y=y, label=label, angle=90, va="bottom", color=color))
            y = y + (w + STACK_GAP_MM) * units_per_mm
            y_reach = max(y_reach, y)
    return pd.DataFrame(rows), y_reach


def tick_rotation(tidy, plot_w_mm=WIDTH_MM * 0.9):
    """0 when every x tick label fits its category slot horizontally (the width is an upper bound, so the
    throwaway render that measures the panel may rotate where the final one would not, never the reverse)."""
    slot = plot_w_mm / tidy["bin"].nunique()
    widest = max(text_size_mm(str(b))[0] for b in tidy["bin"].cat.categories)
    return 0 if widest + 2 * LABEL_PAD_MM <= slot else TICK_ROTATION


def base_panel(tidy, x_label, tag=None, legend=True, y_top=None, rotation=None):
    """Bars, scales and theme. y_top pins the axis top: the label positions are data, so a multiplicative
    expansion would stretch the axis past the height the labels were placed against."""
    limits = (0, y_top) if y_top else None
    rotation = tick_rotation(tidy) if rotation is None else rotation
    return (ggplot()
            + geom_col(tidy, aes("bin", "count", fill="status"), width=BAR_WIDTH)
            + scale_fill_manual(values=FILL_COLORS, breaks=[MATCH, NO_MATCH])
            + scale_y_continuous(limits=limits, expand=(0, 0, 0 if y_top else HEADROOM_FRAC, 0),
                                 labels=lambda breaks: [f"{b:,.0f}" for b in breaks])
            + labs(x=x_label, y="Number of families", tag=tag)
            + theme_minimal()
            + theme(
                text=element_text(family=FONT_FAMILY, size=TEXT_SIZE_PT),
                axis_text_x=element_text(rotation=rotation, ha="right" if rotation else "center"),
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


def composed_areas(tidies, x_labels, tags):
    """Plotting-area sizes (mm) of the panels stacked into one figure: the composition aligns them, so they
    differ from each panel drawn alone (by 2026-11-06: 37 mm each in the figure vs 28-44 mm alone)."""
    height = len(tidies) * PANEL_HEIGHT_MM
    comp = functools.reduce(operator.truediv, [base_panel(t, x, tag=tag, legend=i == 0)
                                               for i, (tag, t, x) in enumerate(zip(tags, tidies, x_labels))])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig = (comp + theme(figure_size=(WIDTH_MM / 25.4, height / 25.4))).draw(show=False)
        fig.canvas.draw()
    areas = [(ax.get_position().width * WIDTH_MM, ax.get_position().height * height) for ax in fig.axes]
    plt.close(fig)
    assert len(areas) == len(tidies), areas
    return areas


def build_panel(tidy, x_label, tag=None, legend=True, height_mm=PANEL_HEIGHT_MM, area=None):
    """One stacked barplot with every number placed to fit; the y-axis grows until the stacks fit.
    area: the plotting area in mm when the panel sits in a composition (composed_areas)."""
    plot_w_mm, plot_h_mm = area or panel_area_mm(tidy, x_label, height_mm)
    y_top = tidy["total"].max() * (1 + HEADROOM_FRAC)
    for _ in range(LAYOUT_PASSES):
        _, y_reach = layout_labels(tidy, plot_w_mm, plot_h_mm, y_top)
        settled = y_reach + STACK_GAP_MM * y_top / plot_h_mm
        if settled <= y_top:
            break
        y_top = settled
    labels, _ = layout_labels(tidy, plot_w_mm, plot_h_mm, y_top)
    layers = [geom_text(group, aes("bin", "y", label="label"), size=TEXT_SIZE_PT, color=color,
                        angle=angle, va=va, ha="center", family=FONT_FAMILY)
              for (angle, va, color), group in labels.groupby(["angle", "va", "color"])]
    return functools.reduce(operator.add, layers, base_panel(tidy, x_label, tag, legend, y_top))


def save_stack(tidies, x_labels, path_no_ext):
    """Panels stacked into one figure, tagged A, B, ... when there are several; only the top one carries the
    legend (the same on all)."""
    tags = list("ABCDEF"[:len(tidies)]) if len(tidies) > 1 else [None]
    areas = composed_areas(tidies, x_labels, tags)

    def combined():  # `/` stacks panels
        return functools.reduce(operator.truediv, [
            build_panel(t, x, tag=tag, legend=i == 0, area=area)
            for i, (tag, t, x, area) in enumerate(zip(tags, tidies, x_labels, areas))])

    save_figure(combined, path_no_ext, len(tidies) * PANEL_HEIGHT_MM)


def save_figure(make_plot, path_no_ext, height_mm):
    """Vector PDF plus a raster preview. Takes a factory: drawing a composition consumes it."""
    size = theme(figure_size=(WIDTH_MM / 25.4, height_mm / 25.4))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for ext in ("pdf", "png"):
            (make_plot() + size).save(f"{path_no_ext}.{ext}", dpi=PREVIEW_DPI, verbose=False)
    print(f"Saved: {path_no_ext}.pdf (+ .png)")


MAIN_PANELS = [  # (column, binning, x-axis label)
    ("size", dict(kind="decade"), "Family size (domain regions in the full alignment)"),
    ("length", dict(kind="right", width=100), "Representative length (aa)"),
    ("plddt", dict(kind="left", width=5.0), "Mean pLDDT of the representative's ESMFold model"),
]


SUPP_FIGURES = {"size": "figure-s1_family-size-linear", "length": "figure-s2_representative-length-linear",
                "ptm": "figure-s3_ptm"}  # Supplementary Figures S1-S3, cited in the Figure 3 caption


def supp_panels(df):
    """Linear small / medium / large splits of size and length, and pTM. Cuts sit on bin edges and close on the
    side the main panel's bins close: size left-closed (< 1,000), length right-closed (<= 300, as 201-300)."""
    out = []
    for col, x_label, kind, cuts, widths, names in (
            ("size", MAIN_PANELS[0][2], "left", (1000, 50000), (50, 2500, 200000), ("small", "medium", "large")),
            ("length", MAIN_PANELS[1][2], "right", (300, 600), (20, 25, 100), ("short", "medium", "long"))):
        (a, b), (w1, w2, w3), (n1, n2, n3) = cuts, widths, names
        c = pl.col(col)
        sh = 0 if kind == "left" else 1  # right-closed: the cut value belongs to the lower split
        out += [(df.filter(c < a + sh), dict(kind=kind, width=w1, upper=a - 1 + sh), x_label, f"{col}_{n1}"),
                (df.filter((c >= a + sh) & (c < b + sh)), dict(kind=kind, width=w2, lower=a + sh, upper=b - 1 + sh),
                 x_label, f"{col}_{n2}"),
                (df.filter(c >= b + sh), dict(kind=kind, width=w3, lower=b + sh), x_label, f"{col}_{n3}")]
    out.append((df, dict(width=0.05), "pTM of the representative's ESMFold model", "ptm"))
    return out


def quartiles(v):
    q1, med, q3 = np.quantile(v, [0.25, 0.5, 0.75])  # linear = statistics.quantiles(method="inclusive")
    return f"median {med:.10g} (IQR {q1:.10g}-{q3:.10g}), range {v.min():.10g}-{v.max():.10g}, n = {len(v)}"


def write_stats(df, release, tidies, path):
    n, m = len(df), int(df["matched"].sum())
    lines = [f"Release {release}: {n} families; match {m} ({100 * m / n:.1f}%), no match {n - m} ({100 * (n - m) / n:.1f}%)"]
    for col in ("size", "length", "plddt", "ptm"):
        for name, sub in (("all", df), ("match", df.filter(pl.col("matched"))), ("no match", df.filter(~pl.col("matched")))):
            lines.append(f"{col} [{name}]: {quartiles(sub[col].drop_nulls().to_numpy())}")
    lines.append("bins (panel, bin, match, no match, no-match %)")
    for (col, _, _), tidy in zip(MAIN_PANELS, tidies):
        for b, g in tidy.groupby("bin", observed=True):
            c = g.set_index("status")["count"]
            lines.append(f"{col}\t{b}\t{c[MATCH]}\t{c[NO_MATCH]}\t{100 * c[NO_MATCH] / (c[MATCH] + c[NO_MATCH]):.1f}")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines[:13]), f"\nSaved: {path}")


def probe_keys(df):
    """The probe's own bins (catalogue_probe.py block()): {(panel, key): (match, no match)}."""
    keys = {"size": lambda x: 10 ** int(math.log10(x)), "length": lambda x: x // 100 * 100,
            "plddt": lambda x: int(x // 5 * 5)}
    out = {}
    for col, key in keys.items():
        for v, matched in df.select(col, "matched").drop_nulls().iter_rows():
            k = (col, key(v))
            a, b = out.get(k, (0, 0))
            out[k] = (a + matched, b + (not matched))
    return out


def check_probe(df, release, log_path, tidies):
    """Exit unless the data gives the probe log's Figure 3 bins for this release, and every panel holds every
    family with a value."""
    want, on = {}, False
    for line in open(log_path):
        if line.startswith("-- "):
            on = line.startswith(f"-- {release} Figure 3 bins")
        elif line.startswith("=="):
            on = False
        elif on and line.strip():
            panel, key, a, b = line.split("\t")
            want[(panel, int(key))] = (int(a), int(b))
    if not want:
        sys.exit(f"No '{release} Figure 3 bins' block in {log_path}")
    got = probe_keys(df)
    if got != want:
        diff = sorted(k for k in set(got) | set(want) if got.get(k) != want.get(k))
        sys.exit(f"Bins differ from the probe log at {diff[:10]}: got {[got.get(k) for k in diff[:10]]}, "
                 f"log {[want.get(k) for k in diff[:10]]}")
    for (col, _, _), tidy in zip(MAIN_PANELS, tidies):
        assert tidy["count"].sum() == df[col].drop_nulls().len(), col
    print(f"Probe check: {len(want)} bins equal to {log_path}")


def selftest():
    df = pl.DataFrame({"size": [6, 9, 10, 99, 100, 1000, 1972419], "length": [75, 100, 101, 200, 1900, 1901, 2000],
                       "ptm": [0.0, 0.05, 0.3, 0.299, 0.15, 0.35, 0.999], "matched": [True] * 7})
    s = bin_starts(df, "size", "decade")["start"].to_list()
    assert s == [1, 1, 10, 10, 100, 1000, 1000000], s
    s = bin_starts(df, "length", "right", 100)["start"].to_list()
    assert s == [1, 1, 101, 101, 1801, 1901, 1901], s
    s = [round(x, 6) for x in bin_starts(df, "ptm", "left", 0.05)["start"].to_list()]
    assert s == [0.0, 0.05, 0.3, 0.25, 0.15, 0.35, 0.95], s
    labels = list(bin_counts(df, "size", "decade")["bin"].cat.categories)
    assert labels == ["6–9", "10–99", "100–999", "1,000–9,999", "≥ 1,000,000"], labels
    labels = list(bin_counts(df, "length", "right", 100)["bin"].cat.categories)
    assert labels[0] == "75–100" and labels[-1] == "1,901–2,000", labels
    assert list(bin_counts(df, "ptm", "left", 0.05)["bin"].cat.categories)[2] == "0.15–0.20"
    assert bin_label(30.0, "left", 5.0, 0, 100, True) == "30–35"
    supp = {name: (sub, kw) for sub, kw, _, name in supp_panels(df)}
    sub, kw = supp["length_long"]
    labels = list(bin_counts(sub, "length", **kw)["bin"].cat.categories)
    assert labels == ["1,801–1,900", "1,901–2,000"], labels  # 2,000 closes the last bin (S2C)
    assert supp["length_short"][0]["length"].to_list() == [75, 100, 101, 200], supp["length_short"][0]
    print("selftest ok")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--catalogue", help="catalogue_<release>.csv.gz from catalogue_probe.py --dump")
    p.add_argument("--release", default="1.1")
    p.add_argument("--probe-log", help="catalogue_probe log whose Figure 3 bins the data must reproduce")
    p.add_argument("--out", default="figures", help="main figure directory (default: figures)")
    p.add_argument("--supp-out", default="supplementary_figures", help="Supplementary Figures S1-S3 (default: supplementary_figures)")
    p.add_argument("--stats-out", default="tmp", help="figure-3_stats_<release>.txt (default: tmp)")
    p.add_argument("--selftest", action="store_true")
    a = p.parse_args()
    selftest()
    if a.selftest:
        return
    for d in (a.out, a.supp_out, a.stats_out):
        os.makedirs(d, exist_ok=True)

    df = load(a.catalogue)
    tidies = [bin_counts(df, col, **kw) for col, kw, _ in MAIN_PANELS]
    if a.probe_log:
        check_probe(df, a.release, a.probe_log, tidies)
    write_stats(df, a.release, tidies, os.path.join(a.stats_out, f"figure-3_stats_{a.release}.txt"))

    suffix = "" if a.release == "1.1" else f"_{a.release}"
    save_stack(tidies, [x for _, _, x in MAIN_PANELS], os.path.join(a.out, f"figure-3_family-metadata{suffix}"))

    groups = {}
    for sub, kw, x_label, name in supp_panels(df):
        tidy = bin_counts(sub, name.split("_")[0], **kw)
        if tidy.empty:
            print(f"No families in range for {name}; skipped.")
            continue
        groups.setdefault(name.split("_")[0], []).append((tidy, x_label))
    for col, panels in groups.items():
        save_stack([t for t, _ in panels], [x for _, x in panels], os.path.join(a.supp_out, f"{SUPP_FIGURES[col]}{suffix}"))


if __name__ == "__main__":
    main()
