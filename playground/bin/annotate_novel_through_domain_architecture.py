#!/usr/bin/env python3

import argparse
import csv
import json
import re
import sqlite3
from pathlib import Path
from typing import Dict, Iterable, List, NamedTuple, Set, Tuple


DEFAULT_SQLITE = Path("assets/mgnifams_v2_results/mgnifams_test.sqlite3")
DEFAULT_NOVEL_IDS = Path(
    "assets/mgnifams_v2_results/generate_families/novel/"
    "mgnifams_no_annotation_ids.txt"
)
DEFAULT_CLAN_MEMBERSHIP = Path(
    "assets/mgnifams_v2_results/generate_families/network/clan_membership.csv"
)
DEFAULT_FAMILY_PFAMS_OUTPUT = Path(
    "assets/mgnifams_v2_results/generate_families/novel/"
    "mgnifams_no_annotation_annotate_novel_through_domain_architecture.csv"
)
DEFAULT_CLAN_OUTPUT = Path(
    "assets/mgnifams_v2_results/generate_families/network/"
    "transiently_annotated_clans.csv"
)
DEFAULT_TRUE_NOVEL_WITHOUT_PFAMS = Path(
    "assets/mgnifams_v2_results/generate_families/novel/"
    "true_novel_without_annotate_novel_through_domain_architecture.txt"
)

PFAM_RE = re.compile(r"^PF\d{5}$")


class ClanAnnotation(NamedTuple):
    clan_id: str
    annotated: bool
    pfams: List[str]


def family_sort_key(family_id: str) -> Tuple[int, object]:
    try:
        return (0, int(family_id))
    except ValueError:
        return (1, family_id)


def normalise_family_id(value: str) -> str:
    match = re.search(r"(\d+)", value)
    if not match:
        raise ValueError(f"Could not extract an integer family ID from {value!r}")
    return str(int(match.group(1)))


def read_family_ids(path: Path) -> Set[str]:
    with open(path) as handle:
        return {normalise_family_id(line.strip()) for line in handle if line.strip()}


def split_family_ids(value: str) -> List[str]:
    return [
        normalise_family_id(family_id.strip())
        for family_id in value.split(";")
        if family_id.strip()
    ]


def parse_domain_blob(domain_blob) -> object:
    if domain_blob is None:
        return {}
    if isinstance(domain_blob, bytes):
        domain_blob = domain_blob.decode()
    return json.loads(domain_blob)


def extract_pfams_from_architecture(value: object) -> Set[str]:
    pfams: Set[str] = set()
    if isinstance(value, dict):
        domain_id = value.get("id")
        if isinstance(domain_id, str) and PFAM_RE.match(domain_id):
            pfams.add(domain_id)
        for child in value.values():
            pfams.update(extract_pfams_from_architecture(child))
    elif isinstance(value, list):
        for item in value:
            pfams.update(extract_pfams_from_architecture(item))
    return pfams


def novel_family_pfams(sqlite_db: Path, novel_ids_txt: Path) -> Dict[str, List[str]]:
    novel_ids = read_family_ids(novel_ids_txt)
    if not novel_ids:
        return {}

    family_pfams: Dict[str, List[str]] = {}
    with sqlite3.connect(sqlite_db) as connection:
        connection.execute("CREATE TEMP TABLE selected_mgnifam_ids (id INTEGER PRIMARY KEY)")
        connection.executemany(
            "INSERT INTO selected_mgnifam_ids (id) VALUES (?)",
            [(int(family_id),) for family_id in novel_ids],
        )
        rows = connection.execute(
            """
            SELECT CAST(m.id AS TEXT), m.domain_blob
            FROM mgnifam AS m
            JOIN selected_mgnifam_ids AS selected ON selected.id = m.id
            ORDER BY m.id
            """
        )
        for family_id, domain_blob in rows:
            architecture = parse_domain_blob(domain_blob)
            family_pfams[family_id] = sorted(extract_pfams_from_architecture(architecture))

    return family_pfams


def transiently_annotated_clans(
    clan_membership_csv: Path,
    novel_ids_txt: Path,
    family_pfams: Dict[str, List[str]],
) -> List[ClanAnnotation]:
    novel_ids = read_family_ids(novel_ids_txt)
    clan_rows: List[ClanAnnotation] = []

    with open(clan_membership_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        required_columns = {"Cluster Id", "Family Ids"}
        missing_columns = required_columns.difference(reader.fieldnames or [])
        if missing_columns:
            raise ValueError(
                "Missing required columns in "
                f"{clan_membership_csv}: {', '.join(sorted(missing_columns))}"
            )

        for row in reader:
            family_ids = split_family_ids(row["Family Ids"])
            if not family_ids or not set(family_ids).issubset(novel_ids):
                continue

            pfams = sorted(
                {
                    pfam
                    for family_id in family_ids
                    for pfam in family_pfams.get(family_id, [])
                }
            )
            clan_rows.append(ClanAnnotation(row["Cluster Id"], bool(pfams), pfams))

    return clan_rows


def write_family_pfams(family_pfams: Dict[str, List[str]], output_csv: Path) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["MGnifam id", "Pfam annotations"])
        for family_id in sorted(family_pfams, key=family_sort_key):
            writer.writerow([family_id, ";".join(family_pfams[family_id])])


def write_clan_annotations(
    clan_rows: Iterable[ClanAnnotation],
    output_csv: Path,
) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Clan Id", "Annotated", "Pfam annotations"])
        for clan_id, annotated, pfams in clan_rows:
            writer.writerow([clan_id, str(annotated), ";".join(pfams)])


def write_true_novel_without_pfams(
    clan_rows: Iterable[ClanAnnotation],
    output_txt: Path,
) -> None:
    output_txt.parent.mkdir(parents=True, exist_ok=True)
    with open(output_txt, "w") as handle:
        for clan_id, annotated, _ in clan_rows:
            if not annotated:
                handle.write(f"{clan_id}\n")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Extract Pfam IDs from MGnifam domain architecture blobs for novel "
            "families, then report true-novel clans with transient Pfam annotations."
        )
    )
    parser.add_argument(
        "--sqlite",
        type=Path,
        default=DEFAULT_SQLITE,
        help=f"MGnifam SQLite database. Defaults to {DEFAULT_SQLITE}.",
    )
    parser.add_argument(
        "--novel-ids",
        type=Path,
        default=DEFAULT_NOVEL_IDS,
        help=f"TXT file with one novel MGnifam id per line. Defaults to {DEFAULT_NOVEL_IDS}.",
    )
    parser.add_argument(
        "--clan-membership",
        type=Path,
        default=DEFAULT_CLAN_MEMBERSHIP,
        help=f"Clan membership CSV. Defaults to {DEFAULT_CLAN_MEMBERSHIP}.",
    )
    parser.add_argument(
        "--family-pfams-output",
        type=Path,
        default=DEFAULT_FAMILY_PFAMS_OUTPUT,
        help=f"Per-family Pfam CSV output. Defaults to {DEFAULT_FAMILY_PFAMS_OUTPUT}.",
    )
    parser.add_argument(
        "--clan-output",
        type=Path,
        default=DEFAULT_CLAN_OUTPUT,
        help=f"Per-clan transient annotation CSV output. Defaults to {DEFAULT_CLAN_OUTPUT}.",
    )
    parser.add_argument(
        "--true-novel-output",
        type=Path,
        default=DEFAULT_TRUE_NOVEL_WITHOUT_PFAMS,
        help=(
            "TXT output for true-novel clans/singletons with no recovered Pfams. "
            f"Defaults to {DEFAULT_TRUE_NOVEL_WITHOUT_PFAMS}."
        ),
    )
    return parser.parse_args()


def main():
    args = parse_args()
    family_pfams = novel_family_pfams(args.sqlite, args.novel_ids)
    clan_rows = transiently_annotated_clans(
        args.clan_membership,
        args.novel_ids,
        family_pfams,
    )

    write_family_pfams(family_pfams, args.family_pfams_output)
    write_clan_annotations(clan_rows, args.clan_output)
    write_true_novel_without_pfams(clan_rows, args.true_novel_output)

    annotated_clans = sum(1 for row in clan_rows if row.annotated)
    print(f"Novel MGnifams parsed: {len(family_pfams)} -> {args.family_pfams_output}")
    print(f"True-novel clans/singletons: {len(clan_rows)} -> {args.clan_output}")
    print(f"Transiently annotated clans/singletons: {annotated_clans}")
    print(f"Unannotated true-novel clans/singletons -> {args.true_novel_output}")


if __name__ == "__main__":
    main()
