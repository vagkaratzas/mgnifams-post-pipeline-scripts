#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
from typing import Iterable, Set, TextIO


def metadata_rows(handle: TextIO) -> csv.DictReader:
    return csv.DictReader(line for line in handle if not line.startswith("#"))


def read_family_ids(path: Path) -> Set[str]:
    family_ids = set()
    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            value = line.strip()
            if not value:
                continue
            parts = value.split("_", 1)
            if len(parts) != 2 or not parts[1]:
                raise ValueError(
                    f"Expected an underscore-delimited id at {path}:{line_number}: {value!r}"
                )
            family_ids.add(parts[1])
    return family_ids


def wrap_sequence(sequence: str, width: int = 80) -> Iterable[str]:
    for start in range(0, len(sequence), width):
        yield sequence[start : start + width]


def write_fasta(metadata_csv: Path, selected_ids_txt: Path, output_fasta: Path) -> int:
    selected_ids = read_family_ids(selected_ids_txt)
    count = 0

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(metadata_csv, newline="") as input_handle, open(output_fasta, "w") as output_handle:
        reader = metadata_rows(input_handle)
        required_columns = {"Family Id", "Sequence"}
        missing_columns = required_columns.difference(reader.fieldnames or [])
        if missing_columns:
            raise ValueError(
                f"Missing required columns in {metadata_csv}: {', '.join(sorted(missing_columns))}"
            )

        for row in reader:
            family_id = row["Family Id"].strip()
            if family_id not in selected_ids:
                continue

            sequence = row["Sequence"].strip()
            output_handle.write(f">{family_id}\n")
            output_handle.write("\n".join(wrap_sequence(sequence)))
            output_handle.write("\n")
            count += 1

    return count


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Write FASTA records from metadata_mqc.csv for family ids listed in "
            "an underscore-delimited clan/singleton id file."
        )
    )
    parser.add_argument(
        "metadata_csv",
        type=Path,
        help="Input metadata_mqc.csv containing Family Id and Sequence columns.",
    )
    parser.add_argument(
        "selected_ids_txt",
        type=Path,
        help="TXT file with ids like Singleton_243; the part after '_' is matched to Family Id.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output FASTA path.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    count = write_fasta(args.metadata_csv, args.selected_ids_txt, args.output)
    print(f"FASTA records written: {count} -> {args.output}")


if __name__ == "__main__":
    main()
