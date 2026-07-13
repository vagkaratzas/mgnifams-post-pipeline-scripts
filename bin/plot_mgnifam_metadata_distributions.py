#!/usr/bin/env python3
# Family-level distribution figures for the MGnifams catalogue (paper Figure X).
# Every panel is a stacked barplot of annotated vs unannotated (novel) families.
#
# Modes (--mode, default all):
#   size   : family size (full alignment sequences) -> 3 PNGs (small / medium / large)
#   length : representative sequence length (aa, == HMM match states) -> 3 PNGs (short / medium / long)
#   plddt  : mean pLDDT of the representative structure -> 1 PNG
#   ptm    : pTM of the representative structure -> 1 PNG
#
# Size and length are split into 3 PNGs each because their ranges (29-1,515,677 and
# 75-2,000) cannot be binned legibly in a single linear axis.
#
# Also writes <prefix>_stats.txt: min/Q1/median/Q3/max per metric, overall and split by
# annotated vs novel -- the numbers quoted in the manuscript paragraph.
#
# Production example (full catalogue):
# python bin/plot_mgnifam_metadata_distributions.py \
#   --metadata mgnifam_codon.csv \
#   --novel-ids assets/mgnifams_v2_results/generate_families/novel/mgnifams_no_annotation_ids.txt \
#   --output-prefix output/mgnifam_metadata
#
# Local test example (small dataset, adjusted thresholds):
# python bin/plot_mgnifam_metadata_distributions.py \
#   --metadata assets/mgnifams_v2_results/table_data/mgnifam.csv \
#   --novel-ids assets/mgnifams_v2_results/generate_families/novel/mgnifams_no_annotation_ids.txt \
#   --output-prefix output/test_mgnifam_metadata \
#   --size-small-max 100 --size-medium-max 1000 \
#   --size-small-bin 20 --size-medium-bin 200 --size-large-bin 1000

import argparse
import numpy as np
import polars as pl
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

ANNOTATED_COLOR = "#faad39"
NOVEL_COLOR = "#747b87"


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


def plot_group(df, value_col, bin_size, x_label, title, output_path, lower_bound=0, upper_bound=None):
    if df.is_empty():
        print(f"No families in range for: {title} - skipping.")
        return

    is_float = df.schema[value_col] in (pl.Float32, pl.Float64)
    df = df.with_columns(((pl.col(value_col) // bin_size) * bin_size).alias("bin_start"))
    agg = df.group_by(["bin_start", "annotated"]).agg(pl.len().alias("count"))

    bins = sorted(agg["bin_start"].unique().to_list())
    if is_float:
        # continuous values: label the bin edges, e.g. "65-70"
        bin_labels = [f"{b:g}-{b + bin_size:g}" for b in bins]
    else:
        # integer counts: label the closed interval the bin really holds, clipped to the
        # group's own range so boundary values never produce a phantom trailing bin
        hi = (lambda b: b + bin_size - 1) if upper_bound is None else \
             (lambda b: min(b + bin_size - 1, upper_bound))
        bin_labels = [f"{max(b, lower_bound + 1):,}-{hi(b):,}" for b in bins]

    annotated_counts, novel = [], []
    for b in bins:
        sub = agg.filter(pl.col("bin_start") == b)
        annotated_counts.append(sub.filter(pl.col("annotated"))["count"].sum())
        novel.append(sub.filter(~pl.col("annotated"))["count"].sum())

    totals = [a + u for a, u in zip(annotated_counts, novel)]
    x = np.arange(len(bins))
    bar_width = 0.6

    _, ax = plt.subplots(figsize=(max(10, len(bins) * 0.7), 6))
    ax.bar(x, annotated_counts, bar_width, label="Annotated", color=ANNOTATED_COLOR)
    ax.bar(x, novel, bar_width, bottom=annotated_counts, label="Unannotated (novel)", color=NOVEL_COLOR)

    # a segment thinner than this cannot hold its own count label without overlapping its neighbour
    min_labelled = 0.025 * max(totals)
    for i, (ann, unann, total) in enumerate(zip(annotated_counts, novel, totals)):
        if ann >= min_labelled:
            ax.text(x[i], ann / 2, f"{ann:,}", ha="center", va="center",
                    fontsize=8, color="black", fontweight="bold")
        if unann >= min_labelled:
            ax.text(x[i], ann + unann / 2, f"{unann:,}", ha="center", va="center",
                    fontsize=8, color="white", fontweight="bold")
        pct = f"{100 * unann / total:.1f}%" if total > 0 else "0%"
        ax.text(x[i], total + 0.05, pct, ha="center", va="bottom", fontsize=8, color="black")

    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels, rotation=45, ha="right", fontsize=8)
    ax.set_xlabel(x_label, fontsize=11)
    ax.set_ylabel("Number of families", fontsize=11)
    ax.set_title(f"{title}\n(% above each bar = share of unannotated families)", fontsize=13)
    ax.set_ylim(0, max(totals) * 1.12)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.legend()
    ax.margins(x=0.02)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Saved: {output_path}")


def split_groups(df, col, x_label, name, prefix, small_max, medium_max, bins, labels):
    """3-way small/medium/large split of a wide-ranged column.

    Cuts are half-open ([0, small_max), [small_max, medium_max), [medium_max, inf)) so that
    each cut falls on a bin edge -- otherwise the boundary value lands in a bin of its own.
    """
    small_bin, medium_bin, large_bin = bins
    lo, mid, hi = labels
    return [
        (df.filter(pl.col(col) < small_max), col, small_bin, x_label,
         f"Family {name} distribution - {lo} (< {small_max:,}, binned by {small_bin:,})",
         f"{prefix}_{name}_{lo}.png", 0, small_max - 1),
        (df.filter((pl.col(col) >= small_max) & (pl.col(col) < medium_max)), col, medium_bin, x_label,
         f"Family {name} distribution - {mid} ({small_max:,}-{medium_max - 1:,}, binned by {medium_bin:,})",
         f"{prefix}_{name}_{mid}.png", small_max - 1, medium_max - 1),
        (df.filter(pl.col(col) >= medium_max), col, large_bin, x_label,
         f"Family {name} distribution - {hi} (>= {medium_max:,}, binned by {large_bin:,})",
         f"{prefix}_{name}_{hi}.png", medium_max - 1, None),
    ]


def main():
    parser = argparse.ArgumentParser(
        description="Stacked barplots (annotated vs unannotated) of MGnifam family size, "
                    "representative length, pLDDT and pTM, plus a manuscript stats report.")
    parser.add_argument("--metadata", required=True,
                        help="mgnifam table CSV (id,full_size,rep_length,...,plddt,ptm)")
    parser.add_argument("--novel-ids", required=True,
                        help="TXT of unannotated (novel) family ids, one per line")
    parser.add_argument("--output-prefix", default="output/mgnifam_metadata",
                        help="Output path prefix (default: output/mgnifam_metadata)")
    parser.add_argument("--mode", choices=["all", "size", "length", "plddt", "ptm"], default="all",
                        help="Which panels to produce (default: all)")

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

    df = load_metadata(args.metadata, args.novel_ids)
    write_stats(df, f"{args.output_prefix}_stats.txt")

    groups = []
    if args.mode in ("all", "size"):
        groups += split_groups(
            df, "size", "Family size (number of sequences in full alignment)", "size",
            args.output_prefix, args.size_small_max, args.size_medium_max,
            (args.size_small_bin, args.size_medium_bin, args.size_large_bin),
            ("small", "medium", "large"))
    if args.mode in ("all", "length"):
        groups += split_groups(
            df, "length", "Representative sequence length (aa)", "length",
            args.output_prefix, args.length_short_max, args.length_medium_max,
            (args.length_short_bin, args.length_medium_bin, args.length_long_bin),
            ("short", "medium", "long"))
    if args.mode in ("all", "plddt"):
        groups.append((df, "plddt", args.plddt_bin, "Mean pLDDT of representative structure",
                       f"Family structural confidence distribution (pLDDT, binned by {args.plddt_bin:g})",
                       f"{args.output_prefix}_plddt.png", 0, None))
    if args.mode in ("all", "ptm"):
        groups.append((df, "ptm", args.ptm_bin, "pTM of representative structure",
                       f"Family structural confidence distribution (pTM, binned by {args.ptm_bin:g})",
                       f"{args.output_prefix}_ptm.png", 0, None))

    for subset, col, bin_size, x_label, title, path, lower, upper in groups:
        plot_group(subset, col, bin_size, x_label, title, path,
                   lower_bound=lower, upper_bound=upper)


if __name__ == "__main__":
    main()
