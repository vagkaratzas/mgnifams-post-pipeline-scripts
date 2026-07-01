#!/usr/bin/env python3
"""Analyse MGnifam biome blobs without collapsing duplicate leaf labels.

The production SQLite database can contain many families, and each `biome_blob`
is itself a small CSV-encoded tree. The analysis hot path therefore uses
streaming SQLite reads and stdlib CSV parsing; pandas is kept at the boundary for
the legacy parser helper and the final plot/report table.
"""

import argparse
import csv
import io
import logging
import os
import sqlite3
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_OUTPUT_DIR = Path("output")
OUTPUT_FIGURE_NAME = "mgnifams_biome_distribution.png"
OUTPUT_REPORT_NAME = "mgnifams_biome_distribution.txt"
DEFAULT_LOG_EVERY = 1000
LOGGER = logging.getLogger(__name__)

TOP_LEVEL_COLOURS = {
    "Host-associated": "#4C8BE0",
    "Environmental": "#E07C3A",
    "Engineered": "#4CAF76",
    "Mixed": "#D94F4F",
    "Unknown": "#AAAAAA",
}


def parse_biome_blob(blob):
    """Return a DataFrame view of a biome blob for legacy callers and tests."""
    records = parse_biome_records(blob)
    if not records:
        return pd.DataFrame()
    return pd.DataFrame(records)


def parse_biome_records(blob):
    """Parse a CSV biome blob into normalized records.

    Invalid or incomplete blobs are treated as absent biome data because a
    single malformed family should not stop a full database run.
    """
    if isinstance(blob, bytes):
        blob = blob.decode("utf-8")
    if not blob or not blob.strip():
        return []

    try:
        rows = list(csv.DictReader(io.StringIO(blob)))
    except Exception:
        return []

    required_columns = {"ids", "labels", "parents", "counts"}
    if not rows or required_columns.difference(rows[0]):
        return []

    records = []
    for row in rows:
        try:
            count = int(float(row.get("counts") or 0))
        except ValueError:
            count = 0
        leaf_id = str(row.get("ids") or "")
        records.append(
            {
                "ids": leaf_id,
                "labels": row.get("labels") or "",
                "parents": row.get("parents") or "",
                "counts": count,
                "depth": len(leaf_id.split(":")) - 1,
            }
        )
    return records


def get_leaf_biomes(df):
    """Return leaf rows, meaning tree nodes never referenced as a parent."""
    if df.empty:
        return pd.DataFrame()
    parent_ids = set(df["parents"].dropna())
    return df[~df["ids"].isin(parent_ids)].copy()


def get_leaf_biome_records(records):
    """Return leaf records from the blob's `ids`/`parents` tree encoding."""
    if not records:
        return []
    parent_ids = {record["parents"] for record in records if record["parents"]}
    return [record for record in records if record["ids"] not in parent_ids]


def get_top_level_label(leaf_id):
    parts = str(leaf_id).split(":")
    return parts[1] if len(parts) > 1 else "Unknown"


def count_biome_rows(db_path):
    """Return the progress denominator for families with biome data."""
    with sqlite3.connect(db_path) as con:
        cur = con.cursor()
        cur.execute("SELECT COUNT(*) FROM mgnifam WHERE biome_blob IS NOT NULL")
        return cur.fetchone()[0]


def iter_biome_rows(db_path):
    """Yield biome rows without materializing the full SQLite result set."""
    with sqlite3.connect(db_path) as con:
        cur = con.cursor()
        cur.execute("SELECT id, biome_blob FROM mgnifam WHERE biome_blob IS NOT NULL")
        yield from cur


def fetch_biome_rows(db_path):
    """Return all biome rows; kept for callers that need an in-memory list."""
    with sqlite3.connect(db_path) as con:
        cur = con.cursor()
        cur.execute("SELECT id, biome_blob FROM mgnifam WHERE biome_blob IS NOT NULL")
        return cur.fetchall()


def analyse_biomes(rows, total_rows=None, log_every=0):
    """Aggregate per-family biome leaves by full biome path.

    The full `root:...` path is the leaf identity. The terminal label is kept
    only as display metadata and to report labels that occur on multiple paths,
    such as `Sediment` under several aquatic branches.
    """
    family_data = {}
    processed = 0

    if total_rows is None:
        try:
            total_rows = len(rows)
        except TypeError:
            total_rows = None

    for fam_id, blob in rows:
        processed += 1
        if log_every and processed % log_every == 0:
            if total_rows is None:
                LOGGER.info("Processed %s families with biome data.", processed)
            else:
                LOGGER.info(
                    "Processed %s/%s families with biome data.",
                    processed,
                    total_rows,
                )

        records = parse_biome_records(blob)
        if not records:
            continue

        leaf_records = get_leaf_biome_records(records)
        if not leaf_records:
            continue

        leaf_entries = []
        top_levels = set()
        for row in leaf_records:
            leaf_path = row["ids"]
            top = get_top_level_label(leaf_path)
            count = int(row["counts"])
            leaf_entries.append((leaf_path, row["labels"], count, top))
            top_levels.add(top)

        family_data[fam_id] = {
            "leaves": leaf_entries,
            "top_levels": top_levels,
            "n_leaves": len(leaf_entries),
            "total_count": sum(count for _, _, count, _ in leaf_entries),
        }

    if log_every:
        if total_rows is None:
            LOGGER.info("Finished processing %s families with biome data.", processed)
        else:
            LOGGER.info(
                "Finished processing %s/%s families with biome data.",
                processed,
                total_rows,
            )

    # Full paths stay separate here; grouping by only `leaf_label` would turn
    # unrelated branches with the same terminal label into false `Mixed` leaves.
    leaf_agg = defaultdict(
        lambda: {"family_ids": set(), "top_levels": set(), "leaf_label": ""}
    )
    for fam_id, fam in family_data.items():
        for leaf_path, leaf_label, count, top in fam["leaves"]:
            leaf_agg[leaf_path]["family_ids"].add(fam_id)
            leaf_agg[leaf_path]["top_levels"].add(top)
            leaf_agg[leaf_path]["leaf_label"] = leaf_label

    leaf_rows = []
    for leaf_path, info in leaf_agg.items():
        tops = info["top_levels"]
        colour_key = "Mixed" if len(tops) > 1 else next(iter(tops), "Unknown")
        leaf_rows.append(
            {
                "label": leaf_path,
                "leaf_label": info["leaf_label"],
                "path": leaf_path,
                "count": len(info["family_ids"]),
                "colour": TOP_LEVEL_COLOURS.get(
                    colour_key, TOP_LEVEL_COLOURS["Unknown"]
                ),
                "top_level": colour_key,
            }
        )

    leaf_df = pd.DataFrame(
        leaf_rows,
        columns=["label", "leaf_label", "path", "count", "colour", "top_level"],
    )
    if not leaf_df.empty:
        leaf_df = leaf_df.sort_values("count", ascending=True)

    exclusive = {
        fam_id: data["leaves"][0]
        for fam_id, data in family_data.items()
        if data["n_leaves"] == 1
    }
    exclusive_leaf_biomes = {
        leaf_path: sorted(info["family_ids"])
        for leaf_path, info in leaf_agg.items()
        if len(info["family_ids"]) == 1
    }
    leaf_labels = defaultdict(lambda: {"paths": [], "family_ids": set()})
    for leaf_path, info in leaf_agg.items():
        label = info["leaf_label"]
        if not label:
            continue
        leaf_labels[label]["paths"].append(leaf_path)
        leaf_labels[label]["family_ids"].update(info["family_ids"])
    duplicate_leaf_labels = {
        label: {
            "paths": sorted(info["paths"]),
            "family_count": len(info["family_ids"]),
        }
        for label, info in leaf_labels.items()
        if len(info["paths"]) > 1
    }

    return {
        "family_data": family_data,
        "total_families": len(family_data),
        "processed_families": processed,
        "leaf_df": leaf_df,
        "exclusive": exclusive,
        "exclusive_leaf_biomes": exclusive_leaf_biomes,
        "duplicate_leaf_labels": duplicate_leaf_labels,
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
    bars = ax.barh(
        leaf_df["label"],
        leaf_df["count"],
        color=leaf_df["colour"],
        edgecolor="white",
        linewidth=0.4,
    )
    label_bar_values(ax, bars)
    ax.set_xlim(right=max(leaf_df["count"]) * 1.12)
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


def label_bar_values(ax, bars):
    """Draw the family count just beyond each horizontal bar."""
    for bar in bars:
        width = bar.get_width()
        ax.text(
            width + 0.2,
            bar.get_y() + bar.get_height() / 2,
            f"{int(width)}",
            va="center",
            ha="left",
            fontsize=8,
            color="#333333",
        )


def format_report(result, output_figure, figure_written):
    family_data = result["family_data"]
    total_families = result["total_families"]
    leaf_df = result["leaf_df"]
    total_leaf_paths = len(leaf_df)
    exclusive = result["exclusive"]
    exclusive_leaf_biomes = result["exclusive_leaf_biomes"]
    duplicate_leaf_labels = result["duplicate_leaf_labels"]
    lines = []

    lines.append(f"Successfully parsed: {total_families} families")
    lines.append("")

    if leaf_df.empty:
        lines.append("No parseable leaf-level biome data found.")
        return "\n".join(lines) + "\n"

    lines.append("=== Leaf-path biome distribution (top 30) ===")
    for _, row in leaf_df.sort_values("count", ascending=False).head(30).iterrows():
        lines.append(f"  {row['label']:<45} {row['count']:>8}  [{row['top_level']}]")

    if figure_written:
        lines.append("")
        lines.append(f"Figure saved -> {output_figure}")

    exclusive_by_biome = defaultdict(list)
    for fam_id, (leaf_path, leaf_label, count, top) in exclusive.items():
        exclusive_by_biome[leaf_path].append(fam_id)

    percent_exclusive = 100 * len(exclusive) / total_families
    lines.append("")
    lines.append("=== Niche / exclusive families (single leaf biome) ===")
    lines.append(
        f"  Total: {len(exclusive)} / {total_families} "
        f"({percent_exclusive:.1f}%)\n"
    )
    lines.append(f"  {'Biome path':<45} {'Families':>10}")
    lines.append(f"  {'-' * 45} {'-' * 10}")
    for biome, fams in sorted(exclusive_by_biome.items(), key=lambda item: -len(item[1])):
        lines.append(f"  {biome:<45} {len(fams):>10}")

    lines.append("")
    lines.append("=== Exclusive leaf paths (present in one family only) ===")
    lines.append(f"  Total: {len(exclusive_leaf_biomes)} leaf paths\n")
    lines.append(f"  {'Biome path':<45} Family ID")
    lines.append(f"  {'-' * 45} {'-' * 20}")
    for biome, fams in sorted(exclusive_leaf_biomes.items()):
        lines.append(f"  {biome:<45} {', '.join(map(str, fams))}")

    duplicate_path_count = sum(
        len(info["paths"]) for info in duplicate_leaf_labels.values()
    )
    lines.append("")
    lines.append("=== Duplicate terminal leaf labels ===")
    lines.append(
        f"  Total: {len(duplicate_leaf_labels)} labels across "
        f"{duplicate_path_count} leaf paths\n"
    )
    lines.append(f"  {'Leaf label':<45} {'Paths':>8} {'Families':>10}")
    lines.append(f"  {'-' * 45} {'-' * 8} {'-' * 10}")
    for label, info in sorted(
        duplicate_leaf_labels.items(),
        key=lambda item: (-len(item[1]["paths"]), item[0]),
    )[:30]:
        lines.append(
            f"  {label:<45} {len(info['paths']):>8} {info['family_count']:>10}"
        )

    breadth = sorted(family_data.items(), key=lambda item: -item[1]["n_leaves"])
    lines.append("")
    lines.append("=== Broadest families (most leaf paths) ===")
    lines.append(
        f"  {'Family ID':<20} {'# leaf paths':>15} "
        f"{'Available leaf paths':>20} {'% available':>12}  Top-level categories"
    )
    lines.append(
        f"  {'-' * 20} {'-' * 15} {'-' * 20} {'-' * 12}  {'-' * 30}"
    )
    for fam_id, data in breadth[:25]:
        tops = ", ".join(sorted(data["top_levels"]))
        percent_available = 100 * data["n_leaves"] / total_leaf_paths
        lines.append(
            f"  {fam_id:<20} {data['n_leaves']:>15} "
            f"{total_leaf_paths:>20} {percent_available:>11.1f}%  {tops}"
        )

    narrowest = sorted(
        family_data.items(), key=lambda item: (item[1]["n_leaves"], str(item[0]))
    )
    lines.append("")
    lines.append("=== Narrowest families (least leaf paths) ===")
    lines.append(
        f"  {'Family ID':<20} {'# leaf paths':>15} "
        f"{'Available leaf paths':>20} {'% available':>12}  Top-level categories"
    )
    lines.append(
        f"  {'-' * 20} {'-' * 15} {'-' * 20} {'-' * 12}  {'-' * 30}"
    )
    for fam_id, data in narrowest[:25]:
        tops = ", ".join(sorted(data["top_levels"]))
        percent_available = 100 * data["n_leaves"] / total_leaf_paths
        lines.append(
            f"  {fam_id:<20} {data['n_leaves']:>15} "
            f"{total_leaf_paths:>20} {percent_available:>11.1f}%  {tops}"
        )

    n_leaves_dist = [data["n_leaves"] for data in family_data.values()]
    lines.append("")
    lines.append("  Leaf path count distribution:")
    lines.append(f"    Median : {np.median(n_leaves_dist):.0f}")
    lines.append(f"    Mean   : {np.mean(n_leaves_dist):.1f}")
    lines.append(f"    Max    : {np.max(n_leaves_dist)}")
    lines.append(f"    > 5 paths : {sum(1 for n in n_leaves_dist if n > 5)} families")
    lines.append(f"    > 10 paths: {sum(1 for n in n_leaves_dist if n > 10)} families")
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
    parser.add_argument(
        "--log-every",
        type=int,
        default=DEFAULT_LOG_EVERY,
        help=(
            "Log progress every N families while analysing biome blobs. "
            f"Defaults to {DEFAULT_LOG_EVERY}; use 0 to disable."
        ),
    )
    return parser.parse_args(argv)


def output_paths(args):
    return (
        args.output_dir / OUTPUT_FIGURE_NAME,
        args.output_dir / OUTPUT_REPORT_NAME,
    )


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    args = parse_args()
    output_figure, output_report = output_paths(args)
    total_rows = count_biome_rows(args.db_path)
    rows = iter_biome_rows(args.db_path)
    LOGGER.info("Found %s families with biome data.", total_rows)

    result = analyse_biomes(rows, total_rows=total_rows, log_every=args.log_every)
    figure_written = False
    if not args.no_plot:
        figure_written = write_leaf_distribution_plot(result["leaf_df"], output_figure)
    report = format_report(result, output_figure, figure_written)
    print(report, end="")
    write_text_report(output_report, report)
    print(f"\nReport saved -> {output_report}")


if __name__ == "__main__":
    main()
