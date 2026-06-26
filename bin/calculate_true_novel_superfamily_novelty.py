#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Set


DEFAULT_NOVEL_FILTERED_MGNIFAMS = Path(
    "input/mgnifams_l100_plddt70_novel.csv"
)
DEFAULT_SUPERFAMILY_STATISTICS = Path(
    "assets/mgnifams_v2_results/generate_families/network/superfamily_statistics.csv"
)
DEFAULT_TRUE_NOVEL_SUPERFAMILIES = Path(
    "assets/mgnifams_v2_results/generate_families/novel/true_novel_superfamilies.txt"
)
DEFAULT_OUTPUT = Path(
    "assets/mgnifams_v2_results/generate_families/network/"
    "true_novel_superfamily_novelty_scores.csv"
)


def extract_integer_id(value: str) -> int:
    match = re.search(r"(\d+)", value)
    if not match:
        raise ValueError(f"Could not extract an integer ID from {value!r}")
    return int(match.group(1))


def read_novel_filtered_mgnifams(path: Path) -> Set[int]:
    novel_filtered_mgnifams: Set[int] = set()
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if header is None:
            return novel_filtered_mgnifams

        for row in reader:
            if not row:
                continue
            novel_filtered_mgnifams.add(extract_integer_id(row[0].strip()))

    return novel_filtered_mgnifams


def read_true_novel_superfamilies(path: Path) -> Set[str]:
    with open(path) as handle:
        return {line.strip() for line in handle if line.strip()}


def split_family_ids(value: str) -> List[str]:
    return [family_id.strip() for family_id in value.split(";") if family_id.strip()]


def count_novel_filtered_family_ids(
    family_ids: List[str],
    novel_filtered_mgnifams: Set[int],
) -> int:
    novel_count = 0
    for family_id in family_ids:
        try:
            if int(family_id) in novel_filtered_mgnifams:
                novel_count += 1
        except ValueError:
            continue
    return novel_count


def calculate_novelty_score(
    family_ids: List[str],
    novel_filtered_mgnifams: Set[int],
) -> float:
    if not family_ids:
        return 0.0

    novel_count = count_novel_filtered_family_ids(family_ids, novel_filtered_mgnifams)
    return 100.0 * novel_count / len(family_ids)


def true_novel_superfamily_novelty_scores(
    novel_filtered_mgnifams_csv: Path,
    superfamily_statistics_csv: Path,
    true_novel_superfamilies_txt: Path,
) -> List[Dict[str, str]]:
    novel_filtered_mgnifams = read_novel_filtered_mgnifams(novel_filtered_mgnifams_csv)
    true_novel_superfamilies = read_true_novel_superfamilies(true_novel_superfamilies_txt)

    with open(superfamily_statistics_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        required_columns = {"Cluster Id", "Family Size", "Family Ids"}
        missing_columns = required_columns.difference(reader.fieldnames or [])
        if missing_columns:
            raise ValueError(
                "Missing required columns in "
                f"{superfamily_statistics_csv}: {', '.join(sorted(missing_columns))}"
            )

        scored_rows: List[Dict[str, str]] = []
        for row in reader:
            if row["Cluster Id"] not in true_novel_superfamilies:
                continue

            family_ids = split_family_ids(row["Family Ids"])
            novel_filtered_count = count_novel_filtered_family_ids(
                family_ids,
                novel_filtered_mgnifams,
            )
            if novel_filtered_count == 0:
                continue

            scored_row = dict(row)
            scored_row["Novelty Score"] = (
                f"{100.0 * novel_filtered_count / len(family_ids):.6f}"
            )
            scored_rows.append(scored_row)

    scored_rows.sort(
        key=lambda row: (
            float(row["Novelty Score"]),
            int(row["Family Size"]),
        ),
        reverse=True,
    )
    return scored_rows


def write_scored_superfamilies(rows: List[Dict[str, str]], output_csv: Path) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        with open(output_csv, "w", newline="") as handle:
            handle.write("")
        return

    fieldnames = list(rows[0].keys())
    with open(output_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Filter true-novel superfamilies and add a novelty score based on the "
            "fraction of member family IDs present in the novel-filtered MGnifams set."
        )
    )
    parser.add_argument(
        "novel_filtered_mgnifams_csv",
        nargs="?",
        type=Path,
        default=DEFAULT_NOVEL_FILTERED_MGNIFAMS,
        help=f"Novel-filtered MGnifams CSV. Defaults to {DEFAULT_NOVEL_FILTERED_MGNIFAMS}.",
    )
    parser.add_argument(
        "superfamily_statistics_csv",
        nargs="?",
        type=Path,
        default=DEFAULT_SUPERFAMILY_STATISTICS,
        help=f"Superfamily statistics CSV. Defaults to {DEFAULT_SUPERFAMILY_STATISTICS}.",
    )
    parser.add_argument(
        "true_novel_superfamilies_txt",
        nargs="?",
        type=Path,
        default=DEFAULT_TRUE_NOVEL_SUPERFAMILIES,
        help=f"True-novel superfamily TXT. Defaults to {DEFAULT_TRUE_NOVEL_SUPERFAMILIES}.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output CSV path. Defaults to {DEFAULT_OUTPUT}.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    rows = true_novel_superfamily_novelty_scores(
        args.novel_filtered_mgnifams_csv,
        args.superfamily_statistics_csv,
        args.true_novel_superfamilies_txt,
    )
    write_scored_superfamilies(rows, args.output)

    print(f"True-novel superfamilies scored: {len(rows)} -> {args.output}")


if __name__ == "__main__":
    main()
