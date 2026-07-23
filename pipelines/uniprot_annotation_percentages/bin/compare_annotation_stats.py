#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


COMPARISON_FIELDS = [
    "before_label",
    "after_label",
    "total_sequences_before",
    "total_sequences_after",
    "annotated_sequences_before",
    "annotated_sequences_after",
    "annotated_sequence_percentage_before",
    "annotated_sequence_percentage_after",
    "annotated_sequence_percentage_point_increase",
    "annotated_sequence_relative_increase_percentage",
    "total_amino_acids_before",
    "total_amino_acids_after",
    "annotated_amino_acids_before",
    "annotated_amino_acids_after",
    "annotated_amino_acid_percentage_before",
    "annotated_amino_acid_percentage_after",
    "annotated_amino_acid_percentage_point_increase",
    "annotated_amino_acid_relative_increase_percentage",
]


def read_stats(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        try:
            return next(reader)
        except StopIteration as exc:
            raise ValueError(f"No stats rows found in {path}") from exc


def to_int(row, key):
    return int(float(row[key]))


def to_float(row, key):
    return float(row[key])


def relative_increase(before, after):
    if before == 0:
        return 0.0
    return round(((after - before) / before) * 100, 6)


def compare_stats(before_csv, after_csv):
    before = read_stats(before_csv)
    after = read_stats(after_csv)

    before_sequence_pct = to_float(before, "annotated_sequence_percentage")
    after_sequence_pct = to_float(after, "annotated_sequence_percentage")
    before_residue_pct = to_float(before, "annotated_amino_acid_percentage")
    after_residue_pct = to_float(after, "annotated_amino_acid_percentage")

    return {
        "before_label": before["label"],
        "after_label": after["label"],
        "total_sequences_before": to_int(before, "total_sequences"),
        "total_sequences_after": to_int(after, "total_sequences"),
        "annotated_sequences_before": to_int(before, "annotated_sequences"),
        "annotated_sequences_after": to_int(after, "annotated_sequences"),
        "annotated_sequence_percentage_before": before_sequence_pct,
        "annotated_sequence_percentage_after": after_sequence_pct,
        "annotated_sequence_percentage_point_increase": round(
            after_sequence_pct - before_sequence_pct,
            6,
        ),
        "annotated_sequence_relative_increase_percentage": relative_increase(
            before_sequence_pct,
            after_sequence_pct,
        ),
        "total_amino_acids_before": to_int(before, "total_amino_acids"),
        "total_amino_acids_after": to_int(after, "total_amino_acids"),
        "annotated_amino_acids_before": to_int(before, "annotated_amino_acids"),
        "annotated_amino_acids_after": to_int(after, "annotated_amino_acids"),
        "annotated_amino_acid_percentage_before": before_residue_pct,
        "annotated_amino_acid_percentage_after": after_residue_pct,
        "annotated_amino_acid_percentage_point_increase": round(
            after_residue_pct - before_residue_pct,
            6,
        ),
        "annotated_amino_acid_relative_increase_percentage": relative_increase(
            before_residue_pct,
            after_residue_pct,
        ),
    }


def write_comparison(output_csv, comparison):
    with Path(output_csv).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COMPARISON_FIELDS)
        writer.writeheader()
        writer.writerow(
            {
                key: f"{value:.6f}" if isinstance(value, float) else value
                for key, value in comparison.items()
            }
        )


def main():
    parser = argparse.ArgumentParser(
        description="Compare Pfam-only and Pfam+MGnifam annotation statistics."
    )
    parser.add_argument("--before", required=True, help="Initial stats CSV")
    parser.add_argument("--after", required=True, help="Updated stats CSV")
    parser.add_argument("--output", required=True, help="Output comparison CSV")
    args = parser.parse_args()

    comparison = compare_stats(args.before, args.after)
    write_comparison(args.output, comparison)


if __name__ == "__main__":
    main()
