#!/usr/bin/env python3

import argparse
import io
import os
import sqlite3
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_OUTPUT_DIR = Path("output")
OUTPUT_FIGURE_NAME = "mgnifams_biome_distribution.png"
OUTPUT_REPORT_NAME = "mgnifams_biome_distribution.txt"

TOP_LEVEL_COLOURS = {
    "Host-associated": "#4C8BE0",
    "Environmental": "#E07C3A",
    "Engineered": "#4CAF76",
    "Mixed": "#D94F4F",
    "Unknown": "#AAAAAA",
}


def parse_biome_blob(blob):
    if isinstance(blob, bytes):
        blob = blob.decode("utf-8")
    if not blob or not blob.strip():
        return pd.DataFrame()

    try:
        df = pd.read_csv(io.StringIO(blob))
    except Exception:
        return pd.DataFrame()

    required_columns = {"ids", "labels", "parents", "counts"}
    if required_columns.difference(df.columns):
        return pd.DataFrame()

    df["counts"] = pd.to_numeric(df["counts"], errors="coerce").fillna(0).astype(int)
    df["depth"] = df["ids"].apply(lambda value: len(str(value).split(":")) - 1)
    return df


def get_leaf_biomes(df):
    """Return leaf rows, meaning nodes that are never another node's parent."""
    if df.empty:
        return pd.DataFrame()
    parent_ids = set(df["parents"].dropna())
    return df[~df["ids"].isin(parent_ids)].copy()


def get_top_level_label(leaf_id):
    parts = str(leaf_id).split(":")
    return parts[1] if len(parts) > 1 else "Unknown"


def fetch_biome_rows(db_path):
    with sqlite3.connect(db_path) as con:
        cur = con.cursor()
        cur.execute("SELECT id, biome_blob FROM mgnifam WHERE biome_blob IS NOT NULL")
        return cur.fetchall()


def analyse_biomes(rows):
    family_data = {}

    for fam_id, blob in rows:
        df = parse_biome_blob(blob)
        if df.empty:
            continue

        leaves_df = get_leaf_biomes(df)
        if leaves_df.empty:
            continue

        leaf_entries = []
        top_levels = set()
        for _, row in leaves_df.iterrows():
            top = get_top_level_label(row["ids"])
            count = int(row["counts"])
            leaf_entries.append((row["labels"], count, top))
            top_levels.add(top)

        family_data[fam_id] = {
            "leaves": leaf_entries,
            "top_levels": top_levels,
            "n_leaves": len(leaf_entries),
            "total_count": sum(count for _, count, _ in leaf_entries),
        }

    leaf_agg = defaultdict(lambda: {"family_ids": set(), "top_levels": set()})
    for fam_id, fam in family_data.items():
        for label, count, top in fam["leaves"]:
            leaf_agg[label]["family_ids"].add(fam_id)
            leaf_agg[label]["top_levels"].add(top)

    leaf_rows = []
    for label, info in leaf_agg.items():
        tops = info["top_levels"]
        colour_key = "Mixed" if len(tops) > 1 else next(iter(tops), "Unknown")
        leaf_rows.append(
            {
                "label": label,
                "count": len(info["family_ids"]),
                "colour": TOP_LEVEL_COLOURS.get(
                    colour_key, TOP_LEVEL_COLOURS["Unknown"]
                ),
                "top_level": colour_key,
            }
        )

    leaf_df = pd.DataFrame(leaf_rows, columns=["label", "count", "colour", "top_level"])
    if not leaf_df.empty:
        leaf_df = leaf_df.sort_values("count", ascending=True)

    exclusive = {
        fam_id: data["leaves"][0]
        for fam_id, data in family_data.items()
        if data["n_leaves"] == 1
    }
    exclusive_leaf_biomes = {
        label: sorted(info["family_ids"])
        for label, info in leaf_agg.items()
        if len(info["family_ids"]) == 1
    }

    return {
        "family_data": family_data,
        "total_families": len(family_data),
        "leaf_df": leaf_df,
        "exclusive": exclusive,
        "exclusive_leaf_biomes": exclusive_leaf_biomes,
    }


def write_leaf_distribution_plot(leaf_df, output_figure):
    if leaf_df.empty:
        return False

    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(12, max(8, len(leaf_df) * 0.28)))
    ax.barh(
        leaf_df["label"],
        leaf_df["count"],
        color=leaf_df["colour"],
        edgecolor="white",
        linewidth=0.4,
    )
    ax.set_xlabel("Number of families", fontsize=11)
    ax.set_title(
        "MGnifams leaf-level biome distribution", fontsize=13, fontweight="bold"
    )
    ax.spines[["top", "right"]].set_visible(False)

    legend_handles = [
        mpatches.Patch(color=colour, label=label)
        for label, colour in TOP_LEVEL_COLOURS.items()
        if label != "Unknown"
    ]
    ax.legend(handles=legend_handles, loc="lower right", framealpha=0.9, fontsize=9)

    output_figure.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_figure, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


def format_report(result, output_figure, figure_written):
    family_data = result["family_data"]
    total_families = result["total_families"]
    leaf_df = result["leaf_df"]
    exclusive = result["exclusive"]
    exclusive_leaf_biomes = result["exclusive_leaf_biomes"]
    lines = []

    lines.append(f"Successfully parsed: {total_families} families")
    lines.append("")

    if leaf_df.empty:
        lines.append("No parseable leaf-level biome data found.")
        return "\n".join(lines) + "\n"

    lines.append("=== Leaf-level biome distribution (top 30) ===")
    for _, row in leaf_df.sort_values("count", ascending=False).head(30).iterrows():
        lines.append(f"  {row['label']:<45} {row['count']:>8}  [{row['top_level']}]")

    if figure_written:
        lines.append("")
        lines.append(f"Figure saved -> {output_figure}")

    exclusive_by_biome = defaultdict(list)
    for fam_id, (label, count, top) in exclusive.items():
        exclusive_by_biome[label].append(fam_id)

    percent_exclusive = 100 * len(exclusive) / total_families
    lines.append("")
    lines.append("=== Niche / exclusive families (single leaf biome) ===")
    lines.append(
        f"  Total: {len(exclusive)} / {total_families} "
        f"({percent_exclusive:.1f}%)\n"
    )
    lines.append(f"  {'Biome':<45} {'Families':>10}")
    lines.append(f"  {'-' * 45} {'-' * 10}")
    for biome, fams in sorted(exclusive_by_biome.items(), key=lambda item: -len(item[1])):
        lines.append(f"  {biome:<45} {len(fams):>10}")

    lines.append("")
    lines.append("=== Exclusive leaf biomes (present in one family only) ===")
    lines.append(f"  Total: {len(exclusive_leaf_biomes)} leaf biomes\n")
    lines.append(f"  {'Biome':<45} Family ID")
    lines.append(f"  {'-' * 45} {'-' * 20}")
    for biome, fams in sorted(exclusive_leaf_biomes.items()):
        lines.append(f"  {biome:<45} {', '.join(map(str, fams))}")

    breadth = sorted(family_data.items(), key=lambda item: -item[1]["n_leaves"])
    lines.append("")
    lines.append("=== Broadest families (most leaf biomes) ===")
    lines.append(f"  {'Family ID':<20} {'# leaf biomes':>15}  Top-level categories")
    lines.append(f"  {'-' * 20} {'-' * 15}  {'-' * 30}")
    for fam_id, data in breadth[:25]:
        tops = ", ".join(sorted(data["top_levels"]))
        lines.append(f"  {fam_id:<20} {data['n_leaves']:>15}  {tops}")

    narrowest = sorted(
        family_data.items(), key=lambda item: (item[1]["n_leaves"], str(item[0]))
    )
    lines.append("")
    lines.append("=== Narrowest families (least leaf biomes) ===")
    lines.append(f"  {'Family ID':<20} {'# leaf biomes':>15}  Top-level categories")
    lines.append(f"  {'-' * 20} {'-' * 15}  {'-' * 30}")
    for fam_id, data in narrowest[:25]:
        tops = ", ".join(sorted(data["top_levels"]))
        lines.append(f"  {fam_id:<20} {data['n_leaves']:>15}  {tops}")

    n_leaves_dist = [data["n_leaves"] for data in family_data.values()]
    lines.append("")
    lines.append("  Leaf biome count distribution:")
    lines.append(f"    Median : {np.median(n_leaves_dist):.0f}")
    lines.append(f"    Mean   : {np.mean(n_leaves_dist):.1f}")
    lines.append(f"    Max    : {np.max(n_leaves_dist)}")
    lines.append(f"    > 5 biomes : {sum(1 for n in n_leaves_dist if n > 5)} families")
    lines.append(f"    > 10 biomes: {sum(1 for n in n_leaves_dist if n > 10)} families")
    return "\n".join(lines) + "\n"


def write_text_report(output_report, text):
    output_report.parent.mkdir(parents=True, exist_ok=True)
    output_report.write_text(text)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Analyse MGnifam biome_blob data and plot leaf-level biomes."
    )
    parser.add_argument(
        "db_path",
        type=Path,
        help="SQLite database path.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for output files. Defaults to {DEFAULT_OUTPUT_DIR}.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Print the biome report without writing a PNG.",
    )
    return parser.parse_args(argv)


def output_paths(args):
    return (
        args.output_dir / OUTPUT_FIGURE_NAME,
        args.output_dir / OUTPUT_REPORT_NAME,
    )


def main():
    args = parse_args()
    output_figure, output_report = output_paths(args)
    rows = fetch_biome_rows(args.db_path)
    print(f"Fetched {len(rows)} families with biome data.")

    result = analyse_biomes(rows)
    figure_written = False
    if not args.no_plot:
        figure_written = write_leaf_distribution_plot(result["leaf_df"], output_figure)
    report = format_report(result, output_figure, figure_written)
    print(report, end="")
    write_text_report(output_report, report)
    print(f"\nReport saved -> {output_report}")


if __name__ == "__main__":
    main()
