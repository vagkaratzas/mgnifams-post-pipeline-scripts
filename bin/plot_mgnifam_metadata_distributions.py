#!/usr/bin/env python3
# Family-level distribution figures for the MGnifams catalogue.
# Every panel is a stacked barplot of annotated vs unannotated (novel) families.
#
# Modes (--mode, default all):
#   figure        -> <out>/figures/            main-text Figure X
#                    A: family size, log10 bins, full range
#                    B: representative sequence length, 100-aa bins, full range
#                    C: mean pLDDT, 5-unit bins
#                    saved as the three standalone panels plus the combined
#                    figure_X_family_metadata.{png,pdf} at 600 dpi
#   supplementary -> <out>/supplementary_figures/
#                    size split small/medium/large, length split short/medium/long
#                    (the fine-grained view behind panels A and B), and pTM
#
# Size uses log10 bins in the main figure because its range (29-1,515,677) cannot be
# binned legibly on a linear axis; the linear small/medium/large split is supplementary.
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
import os

import numpy as np
import polars as pl
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

ANNOTATED_COLOR = "#faad39"
NOVEL_COLOR = "#747b87"
FIGURE_DPI = 600


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


def bin_counts(df, value_col, bin_size, lower_bound=None, upper_bound=None, log=False):
    """Bin one column and count annotated / novel families per bin.

    log=True bins by decade (10^e <= v < 10^(e+1)) instead of by a fixed width, so a
    range spanning five orders of magnitude fits on one legible linear axis.

    lower_bound/upper_bound clip the printed interval labels to the range the subset
    really covers, so the first bin of the full-range length panel reads "75-99", not "1-99".
    """
    if lower_bound is None:
        lower_bound = df[value_col].min() - 1

    if log:
        df = df.filter(pl.col(value_col) > 0).with_columns(
            pl.col(value_col).log10().floor().cast(pl.Int64).alias("bin_start"))
    else:
        df = df.with_columns(((pl.col(value_col) // bin_size) * bin_size).alias("bin_start"))

    agg = df.group_by(["bin_start", "annotated"]).agg(pl.len().alias("count"))
    bins = sorted(agg["bin_start"].unique().to_list())

    if log:
        labels = [f"$10^{{{b}}}$-$10^{{{b + 1}}}$" for b in bins]
    elif df.schema[value_col] in (pl.Float32, pl.Float64):
        # continuous values: label the bin edges, e.g. "65-70"
        labels = [f"{b:g}-{b + bin_size:g}" for b in bins]
    else:
        # integer counts: label the closed interval the bin really holds, clipped to the
        # group's own range so boundary values never produce a phantom trailing bin
        hi = (lambda b: b + bin_size - 1) if upper_bound is None else \
             (lambda b: min(b + bin_size - 1, upper_bound))
        labels = [f"{max(b, lower_bound + 1):,}-{hi(b):,}" for b in bins]

    annotated, novel = [], []
    for b in bins:
        sub = agg.filter(pl.col("bin_start") == b)
        annotated.append(sub.filter(pl.col("annotated"))["count"].sum())
        novel.append(sub.filter(~pl.col("annotated"))["count"].sum())

    return labels, annotated, novel


def draw_panel(ax, labels, annotated, novel, x_label, title, legend=True):
    totals = [a + u for a, u in zip(annotated, novel)]
    x = np.arange(len(labels))

    ax.bar(x, annotated, 0.6, label="Annotated", color=ANNOTATED_COLOR)
    ax.bar(x, novel, 0.6, bottom=annotated, label="Unannotated (novel)", color=NOVEL_COLOR)

    # a segment thinner than this cannot hold its own count label without overlapping its
    # neighbour; a bar this thin is also too small to read a percentage off, and a bin holding
    # a handful of families would otherwise shout a meaningless "100.0% novel"
    min_labelled = 0.025 * max(totals)
    for i, (ann, unann, total) in enumerate(zip(annotated, novel, totals)):
        if ann >= min_labelled:
            ax.text(x[i], ann / 2, f"{ann:,}", ha="center", va="center",
                    fontsize=8, color="black", fontweight="bold")
        if unann >= min_labelled:
            ax.text(x[i], ann + unann / 2, f"{unann:,}", ha="center", va="center",
                    fontsize=8, color="white", fontweight="bold")
        if total >= min_labelled:
            ax.text(x[i], total + 0.05, f"{100 * unann / total:.1f}%",
                    ha="center", va="bottom", fontsize=8, color="black")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_xlabel(x_label, fontsize=11)
    ax.set_ylabel("Number of families", fontsize=11)
    if title:
        ax.set_title(title, fontsize=12)
    ax.set_ylim(0, max(totals) * 1.12)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    if legend:
        ax.legend(fontsize=9)
    ax.margins(x=0.02)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save_panel(df, output_path, value_col, bin_size, x_label, title,
               lower_bound=0, upper_bound=None, log=False):
    if df.is_empty():
        print(f"No families in range for: {title} - skipping.")
        return None

    binned = bin_counts(df, value_col, bin_size, lower_bound, upper_bound, log)
    labels = binned[0]

    fig, ax = plt.subplots(figsize=(max(10, len(labels) * 0.7), 6))
    draw_panel(ax, *binned, x_label, f"{title}\n(% above each bar = share of unannotated families)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=FIGURE_DPI)
    plt.close(fig)
    print(f"Saved: {output_path}")
    return binned


def figure_panels(df, plddt_bin):
    """The three main-text panels: (binning args, x label, standalone title)."""
    return [
        (dict(value_col="size", bin_size=None, log=True),
         "Family size (sequences in full alignment)",
         "Family size distribution (log$_{10}$ bins)"),
        (dict(value_col="length", bin_size=100),
         "Representative sequence length (aa)",
         "Representative sequence length distribution (100-aa bins)"),
        (dict(value_col="plddt", bin_size=plddt_bin),
         "Mean pLDDT of representative structure",
         f"Structural confidence distribution (pLDDT, {plddt_bin:g}-unit bins)"),
    ]


def build_figure(df, panels, output_prefix):
    """Combined publication figure: the panels stacked vertically, lettered A/B/C."""
    fig, axes = plt.subplots(len(panels), 1, figsize=(11, 5 * len(panels)))

    for ax, letter, (bin_args, x_label, _) in zip(axes, "ABCDEFG", panels):
        labels, annotated, novel = bin_counts(df, **bin_args)
        draw_panel(ax, labels, annotated, novel, x_label, title=None, legend=letter == "A")
        ax.text(-0.06, 1.02, letter, transform=ax.transAxes,
                fontsize=16, fontweight="bold", va="bottom", ha="right")

    fig.tight_layout(h_pad=3)
    for ext in ("png", "pdf"):
        path = f"{output_prefix}.{ext}"
        fig.savefig(path, dpi=FIGURE_DPI)
        print(f"Saved: {path}")
    plt.close(fig)


def split_panels(df, col, x_label, name, small_max, medium_max, bins, labels):
    """3-way small/medium/large split of a wide-ranged column, for the supplement.

    Cuts are half-open ([0, small_max), [small_max, medium_max), [medium_max, inf)) so that
    each cut falls on a bin edge -- otherwise the boundary value lands in a bin of its own.
    """
    small_bin, medium_bin, large_bin = bins
    lo, mid, hi = labels
    return [
        (df.filter(pl.col(col) < small_max), col, small_bin, x_label,
         f"Family {name} distribution - {lo} (< {small_max:,}, binned by {small_bin:,})",
         f"{name}_{lo}.png", None, small_max - 1),
        (df.filter((pl.col(col) >= small_max) & (pl.col(col) < medium_max)), col, medium_bin, x_label,
         f"Family {name} distribution - {mid} ({small_max:,}-{medium_max - 1:,}, binned by {medium_bin:,})",
         f"{name}_{mid}.png", small_max - 1, medium_max - 1),
        (df.filter(pl.col(col) >= medium_max), col, large_bin, x_label,
         f"Family {name} distribution - {hi} (>= {medium_max:,}, binned by {large_bin:,})",
         f"{name}_{hi}.png", medium_max - 1, None),
    ]


def main():
    parser = argparse.ArgumentParser(
        description="Stacked barplots (annotated vs unannotated) of MGnifam family size, "
                    "representative length, pLDDT and pTM, plus a manuscript stats report.")
    parser.add_argument("--metadata", required=True,
                        help="mgnifam metadata CSV (id,full_size,rep_length,...,plddt,ptm)")
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
        panels = figure_panels(df, args.plddt_bin)
        for letter, (bin_args, x_label, title) in zip("abc", panels):
            save_panel(df, os.path.join(fig_dir, f"figure_X{letter}_{bin_args['value_col']}.png"),
                       x_label=x_label, title=title, **bin_args)
        build_figure(df, panels, os.path.join(fig_dir, "figure_X_family_metadata"))

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
        panels.append((df, "ptm", args.ptm_bin, "pTM of representative structure",
                       f"Structural confidence distribution (pTM, binned by {args.ptm_bin:g})",
                       "ptm.png", 0, None))

        for subset, col, bin_size, x_label, title, name, lower, upper in panels:
            save_panel(subset, os.path.join(supp_dir, name), col, bin_size, x_label, title,
                       lower_bound=lower, upper_bound=upper)


if __name__ == "__main__":
    main()
