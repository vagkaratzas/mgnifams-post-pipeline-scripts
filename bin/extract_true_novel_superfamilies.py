#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
from typing import List, NamedTuple, Set


DEFAULT_NOVEL_IDS = Path(
    "assets/mgnifams_v2_results/generate_families/novel/novel_ids.txt"
)
DEFAULT_SUPERFAMILY_STATISTICS = Path(
    "assets/mgnifams_v2_results/generate_families/network/superfamily_statistics.csv"
)
DEFAULT_OUTPUT = Path(
    "assets/mgnifams_v2_results/generate_families/novel/"
    "true_novel_superfamilies.txt"
)


class Superfamily(NamedTuple):
    cluster_id: str
    family_ids: List[str]


def read_novel_ids(path: Path) -> Set[str]:
    with open(path) as handle:
        return {line.strip() for line in handle if line.strip()}


def split_family_ids(value: str) -> List[str]:
    return [family_id.strip() for family_id in value.split(";") if family_id.strip()]


def true_novel_superfamilies(
    superfamily_statistics_csv: Path,
    novel_ids_txt: Path,
) -> List[Superfamily]:
    novel_ids = read_novel_ids(novel_ids_txt)

    with open(superfamily_statistics_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        required_columns = {"Cluster Id", "Family Ids"}
        missing_columns = required_columns.difference(reader.fieldnames or [])
        if missing_columns:
            raise ValueError(
                "Missing required columns in "
                f"{superfamily_statistics_csv}: {', '.join(sorted(missing_columns))}"
            )

        superfamilies = []
        for row in reader:
            family_ids = split_family_ids(row["Family Ids"])
            if family_ids and set(family_ids).issubset(novel_ids):
                superfamilies.append(Superfamily(row["Cluster Id"], family_ids))

    return superfamilies


def write_superfamily_ids(superfamilies: List[Superfamily], output_txt: Path) -> None:
    output_txt.parent.mkdir(parents=True, exist_ok=True)
    with open(output_txt, "w") as handle:
        for superfamily in superfamilies:
            handle.write(f"{superfamily.cluster_id}\n")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Write superfamily ids where every member family is present in "
            "the MGnifam novel family id list."
        )
    )
    parser.add_argument(
        "novel_ids_txt",
        nargs="?",
        type=Path,
        default=DEFAULT_NOVEL_IDS,
        help=f"Novel family ids TXT. Defaults to {DEFAULT_NOVEL_IDS}.",
    )
    parser.add_argument(
        "superfamily_statistics_csv",
        nargs="?",
        type=Path,
        default=DEFAULT_SUPERFAMILY_STATISTICS,
        help=(
            "Superfamily statistics CSV containing Cluster Id and Family Ids. "
            f"Defaults to {DEFAULT_SUPERFAMILY_STATISTICS}."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output TXT path. Defaults to {DEFAULT_OUTPUT}.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    superfamilies = true_novel_superfamilies(
        args.superfamily_statistics_csv,
        args.novel_ids_txt,
    )
    write_superfamily_ids(superfamilies, args.output)

    print(f"True-novel superfamilies: {len(superfamilies)} -> {args.output}")


if __name__ == "__main__":
    main()

