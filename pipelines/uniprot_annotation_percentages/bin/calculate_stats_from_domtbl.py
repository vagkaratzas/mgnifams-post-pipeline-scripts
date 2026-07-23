#!/usr/bin/env python3
"""Sequence- and residue-level annotation stats computed directly from hmmsearch domtblouts.

Coverage math (add_covered_range / integer position sets / percentage / STAT_FIELDS) is salvaged
from annotation_percentages/bin/calculate_annotation_stats.py; the domtblout column layout is
salvaged from annotation_percentages/bin/append_mgnifams_annot.py.
"""
import argparse
import gzip
from pathlib import Path


# STAT_FIELDS kept identical to calculate_annotation_stats.py so compare_annotation_stats.py
# consumes these rows unchanged.
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


def open_text(path, mode="rt"):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return path.open(mode)


def percentage(numerator, denominator):
    if denominator == 0:
        return 0.0
    return round((numerator / denominator) * 100, 6)


def read_fasta_lengths(fasta):
    """{seq_id: length}; id = header first token, matching hmmsearch domtblout target_name."""
    lengths = {}
    current = None
    with open_text(fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                current = line[1:].split()[0]
                lengths[current] = 0
            elif current is not None:
                lengths[current] += len(line.strip())
    return lengths


def parse_domtbl_ranges(domtbl, domain_evalue_threshold=None, ranges=None):
    """{target_name: [(ali_from, ali_to), ...]} across one domtblout."""
    if ranges is None:
        ranges = {}
    with open_text(domtbl) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 19:
                continue
            target_name = parts[0]
            domain_i_evalue = float(parts[12])
            ali_from = int(parts[17])
            ali_to = int(parts[18])
            if (
                domain_evalue_threshold is not None
                and domain_i_evalue > domain_evalue_threshold
            ):
                continue
            ranges.setdefault(target_name, []).append((ali_from, ali_to))
    return ranges


def parse_domtbls(domtbls, domain_evalue_threshold=None):
    ranges = {}
    for domtbl in domtbls:
        parse_domtbl_ranges(domtbl, domain_evalue_threshold, ranges)
    return ranges


def covered_positions(target_ranges, sequence_length):
    covered = set()
    for start, end in target_ranges:
        if start > end:
            start, end = end, start
        start = max(1, start)
        end = min(sequence_length, end)
        if start <= end:
            covered.update(range(start, end + 1))
    return covered


def calculate_stats(fasta_lengths, ranges, label, subtract_ranges=None):
    total_sequences = len(fasta_lengths)
    total_amino_acids = sum(fasta_lengths.values())
    annotated_sequences = 0
    annotated_amino_acids = 0
    exclusive_sequences = 0
    exclusive_amino_acids = 0

    for seq_id, seq_len in fasta_lengths.items():
        covered = covered_positions(ranges.get(seq_id, []), seq_len)
        if covered:
            annotated_sequences += 1
            annotated_amino_acids += len(covered)
        if subtract_ranges is not None:
            baseline = covered_positions(subtract_ranges.get(seq_id, []), seq_len)
            exclusive = covered - baseline
            if exclusive:
                exclusive_sequences += 1
                exclusive_amino_acids += len(exclusive)

    stats = {
        "label": label,
        "annotation_keys": label,
        "total_sequences": total_sequences,
        "annotated_sequences": annotated_sequences,
        "annotated_sequence_percentage": percentage(annotated_sequences, total_sequences),
        "total_amino_acids": total_amino_acids,
        "annotated_amino_acids": annotated_amino_acids,
        "annotated_amino_acid_percentage": percentage(annotated_amino_acids, total_amino_acids),
    }
    if subtract_ranges is not None:
        # Extra columns (e.g. Mgnifams-exclusive space); compare_annotation_stats.py ignores them.
        stats["exclusive_annotated_sequences"] = exclusive_sequences
        stats["exclusive_annotated_sequence_percentage"] = percentage(
            exclusive_sequences, total_sequences
        )
        stats["exclusive_annotated_amino_acids"] = exclusive_amino_acids
        stats["exclusive_annotated_amino_acid_percentage"] = percentage(
            exclusive_amino_acids, total_amino_acids
        )
    return stats


def write_stats(output, stats):
    import csv

    with open_text(output, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(stats.keys()))
        writer.writeheader()
        writer.writerow(
            {
                key: f"{value:.6f}" if isinstance(value, float) else value
                for key, value in stats.items()
            }
        )


def main():
    parser = argparse.ArgumentParser(
        description="Annotation stats (sequence + residue) from hmmsearch domtblouts and a FASTA."
    )
    parser.add_argument("--fasta", required=True, help="Sequence FASTA or FASTA.GZ (defines lengths)")
    parser.add_argument("--label", required=True, help="Stats row label")
    parser.add_argument("--output", required=True, help="Output stats CSV or CSV.GZ")
    parser.add_argument(
        "--domtbl", nargs="+", required=True, help="One or more hmmsearch --domtblout files (unioned)"
    )
    parser.add_argument(
        "--subtract",
        nargs="*",
        default=None,
        help="Domtblouts whose coverage is subtracted to also report exclusive residues",
    )
    parser.add_argument(
        "--domain-evalue-threshold",
        type=float,
        default=None,
        help="Max independent domain E-value; omit for --cut_ga searches (already filtered)",
    )
    args = parser.parse_args()

    fasta_lengths = read_fasta_lengths(args.fasta)
    ranges = parse_domtbls(args.domtbl, args.domain_evalue_threshold)
    subtract_ranges = (
        parse_domtbls(args.subtract, args.domain_evalue_threshold)
        if args.subtract
        else None
    )
    stats = calculate_stats(fasta_lengths, ranges, args.label, subtract_ranges)
    write_stats(args.output, stats)


if __name__ == "__main__":
    main()
