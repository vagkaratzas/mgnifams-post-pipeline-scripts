#!/usr/bin/env python3
# Produces 3 PNGs split by family length or size (--mode length|size):
#   short  : value <= short_max           (bin short_bin)
#   medium : short_max < value <= med_max (bin med_bin)
#   long   : value > med_max              (bin long_bin)
#
# Production example:
# python bin/plot_family_length_distribution.py \
#   --metadata <path/to/metadata_mqc.csv> \
#   --domtbl <path/to/funfams.domtbl.gz> <path/to/pfam.domtbl.gz> \
#   --output-prefix family_length \
#   --mode length \
#   --short-max 300 --med-max 1000 \
#   --short-bin 10 --med-bin 50 --long-bin 100
#
# Local test example (small dataset, adjusted thresholds):
# python bin/plot_family_length_distribution.py \
#   --metadata assets/mgnifams_v2_results/generate_families/metadata_mqc.csv \
#   --domtbl assets/mgnifams_v2_results/annotation/reps/funfams/reps_fasta.domtbl.gz \
#            assets/mgnifams_v2_results/annotation/reps/pfam/reps_fasta.domtbl.gz \
#   --output-prefix output/family_length \
#   --mode length \
#   --short-max 150 --med-max 500 \
#   --short-bin 10 --med-bin 50 --long-bin 100

import argparse
import gzip
import io
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def parse_annotated_ids(*domtbl_gz_paths):
    annotated = set()
    for path in domtbl_gz_paths:
        with gzip.open(path, "rt") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                annotated.add(line.split()[0])
    return annotated


def load_metadata(csv_path):
    with open(csv_path) as f:
        lines = [line for line in f if not line.startswith("#")]
    df = pl.read_csv(io.StringIO("".join(lines)))
    return df.with_columns([
        pl.col("Family Id").cast(pl.Utf8).alias("family_id"),
        pl.col("HMM consensus").str.len_chars().alias("hmm_length"),
        pl.col("Size").cast(pl.Int64).alias("size"),
    ]).select(["family_id", "hmm_length", "size"])


def plot_group(df, value_col, bin_size, x_label, title, output_path, lower_bound=0):
    if df.is_empty():
        print(f"No families in range for: {title} — skipping.")
        return

    df = df.with_columns(
        ((pl.col(value_col) // bin_size) * bin_size).alias("bin_start")
    )
    agg = (
        df.group_by(["bin_start", "annotated"])
        .agg(pl.len().alias("count"))
        .sort(["bin_start", "annotated"])
    )

    bins = sorted(agg["bin_start"].unique().to_list())
    bin_labels = [f"{max(b + 1, lower_bound + 1)}-{b + bin_size}" for b in bins]

    unannotated, annotated_counts = [], []
    for b in bins:
        sub = agg.filter(pl.col("bin_start") == b)
        annotated_counts.append(sub.filter(pl.col("annotated"))["count"].sum())
        unannotated.append(sub.filter(~pl.col("annotated"))["count"].sum())

    totals = [a + u for a, u in zip(annotated_counts, unannotated)]
    x = np.arange(len(bins))
    bar_width = 0.6

    fig, ax = plt.subplots(figsize=(max(10, len(bins) * 0.7), 6))
    ax.bar(x, annotated_counts, bar_width, label="Annotated", color="#faad39")
    ax.bar(x, unannotated, bar_width, bottom=annotated_counts, label="Unannotated", color="#747b87")

    for i, (ann, unann, total) in enumerate(zip(annotated_counts, unannotated, totals)):
        if ann > 0:
            ax.text(x[i], ann / 2, str(ann), ha="center", va="center",
                    fontsize=8, color="black", fontweight="bold")
        if unann > 0:
            ax.text(x[i], ann + unann / 2, str(unann), ha="center", va="center",
                    fontsize=8, color="white", fontweight="bold")
        pct = f"{100 * ann / total:.1f}%" if total > 0 else "0%"
        ax.text(x[i], total + 0.05, pct, ha="center", va="bottom",
                fontsize=8, color="black")

    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels, rotation=45, ha="right", fontsize=8)
    ax.set_xlabel(x_label, fontsize=11)
    ax.set_ylabel("Number of families", fontsize=11)
    ax.set_title(title, fontsize=13)
    ax.set_ylim(0, max(totals) * 1.1)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.legend()
    ax.margins(x=0.02)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Produce 3 stacked bar PNGs of family distribution by annotation status."
    )
    parser.add_argument("--metadata", required=True, help="metadata_mqc.csv")
    parser.add_argument("--domtbl", required=True, nargs="+", help="One or more .domtbl.gz annotation files")
    parser.add_argument("--mode", choices=["length", "size"], default="length",
                        help="Split and bin families by HMM consensus length (aa) or family size (# sequences) (default: length)")
    parser.add_argument("--short-max", type=int, default=300,
                        help="Upper bound for short-family PNG (default: 300)")
    parser.add_argument("--med-max", type=int, default=1000,
                        help="Upper bound for medium-family PNG (default: 1000)")
    parser.add_argument("--short-bin", type=int, default=10,
                        help="Bin size for short families (default: 10)")
    parser.add_argument("--med-bin", type=int, default=50,
                        help="Bin size for medium families (default: 50)")
    parser.add_argument("--long-bin", type=int, default=100,
                        help="Bin size for long families (default: 100)")
    parser.add_argument("--output-prefix", default="family_length",
                        help="Output path prefix; produces <prefix>_short.png, _medium.png, _long.png")
    args = parser.parse_args()

    value_col = "hmm_length" if args.mode == "length" else "size"
    x_label = "HMM consensus length (aa)" if args.mode == "length" else "Family size (number of sequences)"
    mode_label = "length" if args.mode == "length" else "size"

    annotated_ids = parse_annotated_ids(*args.domtbl)
    df = load_metadata(args.metadata)
    df = df.with_columns(
        pl.col("family_id").is_in(list(annotated_ids)).alias("annotated")
    )

    groups = [
        (
            df.filter(pl.col(value_col) <= args.short_max),
            args.short_bin,
            f"Family {mode_label} distribution — short (1–{args.short_max}, binned by {args.short_bin})",
            f"{args.output_prefix}_short.png",
            0,
        ),
        (
            df.filter((pl.col(value_col) > args.short_max) & (pl.col(value_col) <= args.med_max)),
            args.med_bin,
            f"Family {mode_label} distribution — medium ({args.short_max + 1}–{args.med_max}, binned by {args.med_bin})",
            f"{args.output_prefix}_medium.png",
            args.short_max,
        ),
        (
            df.filter(pl.col(value_col) > args.med_max),
            args.long_bin,
            f"Family {mode_label} distribution — long ({args.med_max + 1}+, binned by {args.long_bin})",
            f"{args.output_prefix}_long.png",
            args.med_max,
        ),
    ]

    for subset, bin_size, title, output_path, lower_bound in groups:
        plot_group(subset, value_col, bin_size, x_label, title, output_path, lower_bound=lower_bound)


if __name__ == "__main__":
    main()
