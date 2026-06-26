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


def open_text(path, mode="rt", newline=None):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode, newline=newline)
    return path.open(mode, newline=newline)


def parse_domtblout(domtblout, domain_evalue_threshold=0.001, annotations=None):
    if annotations is None:
        annotations = {}
    with open_text(domtblout) as file:
        for line in file:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.strip().split()
            if len(parts) < 22:
                continue

            target_name = parts[0]
            query_name = parts[3]
            domain_i_evalue = float(parts[12])
            domain_score = float(parts[13])
            ali_from = int(parts[17])
            ali_to = int(parts[18])

            if (
                domain_evalue_threshold is not None
                and domain_i_evalue > domain_evalue_threshold
            ):
                continue

            annotations.setdefault(target_name, []).append(
                [query_name, domain_i_evalue, domain_score, ali_from, ali_to]
            )

    return annotations


def parse_domtblouts(domtblouts, domain_evalue_threshold=0.001):
    annotations = {}
    for domtblout in domtblouts:
        parse_domtblout(
            domtblout,
            domain_evalue_threshold=domain_evalue_threshold,
            annotations=annotations,
        )
    return annotations


def update_csv_with_annotations(csv_file, annotation_data, output_file):
    updated_rows = []

    with open_text(csv_file, newline="") as file:
        reader = csv.DictReader(file)
        fieldnames = reader.fieldnames
        for row in reader:
            mgyp = row["mgyp"]
            metadata = json.loads(row.get("metadata") or "{}")
            metadata.pop("m", None)

            if mgyp in annotation_data:
                metadata["m"] = annotation_data[mgyp]

            row["metadata"] = json.dumps(metadata, separators=(",", ":"))
            updated_rows.append(row)

    with open_text(output_file, "wt", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(updated_rows)


def main():
    parser = argparse.ArgumentParser(
        description="Append MGnifam hmmsearch domtblout annotations to a protein CSV metadata column."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="One or more hmmsearch --domtblout files followed by input CSV and output CSV",
    )
    parser.add_argument(
        "--domain-evalue-threshold",
        type=float,
        default=0.001,
        help="Maximum independent domain E-value to append",
    )
    args = parser.parse_args()
    if len(args.inputs) < 3:
        parser.error("expected at least one domtblout, an input CSV, and an output CSV")

    domtblouts = args.inputs[:-2]
    csv_file = args.inputs[-2]
    output_file = args.inputs[-1]
    annotation_data = parse_domtblouts(
        domtblouts,
        domain_evalue_threshold=args.domain_evalue_threshold,
    )
    update_csv_with_annotations(csv_file, annotation_data, output_file)

    print(f"Updated CSV file saved as {output_file}")


if __name__ == "__main__":
    main()
