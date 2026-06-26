#!/usr/bin/env python3
import argparse
import csv
import gzip
import json
import sys
from pathlib import Path


def raise_csv_field_limit():
    """Protein metadata fields can exceed Python's default 131072-byte csv limit."""
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit = int(limit / 10)


raise_csv_field_limit()


STAT_FIELDS = [
    "label",
    "annotation_keys",
    "total_sequences",
    "annotated_sequences",
    "annotated_sequence_percentage",
    "total_amino_acids",
    "annotated_amino_acids",
    "annotated_amino_acid_percentage",
]


def open_text(path, mode="rt", newline=None):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode, newline=newline)
    return path.open(mode, newline=newline)


def percentage(numerator, denominator):
    if denominator == 0:
        return 0.0
    return round((numerator / denominator) * 100, 6)


def add_covered_range(covered, start, end, sequence_length):
    start = int(start)
    end = int(end)
    if start > end:
        start, end = end, start

    start = max(1, start)
    end = min(sequence_length, end)
    if start > end:
        return

    covered.update(range(start, end + 1))


def covered_positions(metadata, annotation_keys, sequence_length):
    covered = set()
    for key in annotation_keys:
        for annotation in metadata.get(key, []):
            if key == "p" and len(annotation) >= 7:
                add_covered_range(covered, annotation[5], annotation[6], sequence_length)
            elif key == "m" and len(annotation) >= 5:
                add_covered_range(covered, annotation[3], annotation[4], sequence_length)
            else:
                raise ValueError(f"Unsupported annotation key or malformed annotation: {key}")
    return covered


def calculate_stats(input_csv, annotation_keys, label):
    total_sequences = 0
    annotated_sequences = 0
    total_amino_acids = 0
    annotated_amino_acids = 0

    with open_text(input_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            total_sequences += 1
            sequence_length = len(row.get("sequence", ""))
            total_amino_acids += sequence_length

            raw_metadata = row.get("metadata") or "{}"
            metadata = json.loads(raw_metadata)
            covered = covered_positions(metadata, annotation_keys, sequence_length)

            if covered:
                annotated_sequences += 1
                annotated_amino_acids += len(covered)

    return {
        "label": label,
        "annotation_keys": ",".join(annotation_keys),
        "total_sequences": total_sequences,
        "annotated_sequences": annotated_sequences,
        "annotated_sequence_percentage": percentage(annotated_sequences, total_sequences),
        "total_amino_acids": total_amino_acids,
        "annotated_amino_acids": annotated_amino_acids,
        "annotated_amino_acid_percentage": percentage(
            annotated_amino_acids,
            total_amino_acids,
        ),
    }


def write_stats(output_csv, stats):
    with open_text(output_csv, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=STAT_FIELDS)
        writer.writeheader()
        writer.writerow(
            {
                key: f"{value:.6f}" if isinstance(value, float) else value
                for key, value in stats.items()
            }
        )


def parse_annotation_keys(value):
    keys = [key.strip() for key in value.split(",") if key.strip()]
    unsupported = sorted(set(keys) - {"p", "m"})
    if unsupported:
        raise argparse.ArgumentTypeError(
            f"Unsupported annotation key(s): {', '.join(unsupported)}"
        )
    if not keys:
        raise argparse.ArgumentTypeError("At least one annotation key is required")
    return keys


def main():
    parser = argparse.ArgumentParser(
        description="Calculate sequence and amino-acid annotation percentages from MGnify protein CSV metadata."
    )
    parser.add_argument("--input", required=True, help="Input protein CSV or CSV.GZ")
    parser.add_argument(
        "--annotation-keys",
        required=True,
        type=parse_annotation_keys,
        help="Comma-separated metadata annotation keys to count, for example p or p,m",
    )
    parser.add_argument("--label", required=True, help="Label to write in the stats row")
    parser.add_argument("--output", required=True, help="Output stats CSV or CSV.GZ")
    args = parser.parse_args()

    stats = calculate_stats(args.input, args.annotation_keys, args.label)
    write_stats(args.output, stats)


if __name__ == "__main__":
    main()
