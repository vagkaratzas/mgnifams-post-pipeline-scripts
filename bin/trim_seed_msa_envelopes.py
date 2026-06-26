#!/usr/bin/env python3

import argparse
import gzip
import re
from pathlib import Path


HEADER_RE = re.compile(r"^(?P<name>[^/\s]+)(?:/(?P<start>\d+)-(?P<end>\d+))?(?P<rest>\s.*)?$")
GAP_CHARS = {"-", "."}


def read_rf(rf_path):
    rf = "".join(rf_path.read_text().split())
    invalid = set(rf) - {".", "x"}
    if invalid:
        invalid_symbols = "".join(sorted(invalid))
        raise ValueError(f"{rf_path} contains invalid RF symbol(s): {invalid_symbols}")
    if not rf:
        raise ValueError(f"{rf_path} is empty")
    if "x" not in rf:
        raise ValueError(f"{rf_path} does not contain any HMM match columns")
    return rf


def terminal_trim_bounds(rf):
    left = len(rf) - len(rf.lstrip("."))
    right = len(rf) - len(rf.rstrip("."))
    return left, len(rf) - right


def read_fasta_gz(msa_path):
    records = []
    header = None
    sequence_lines = []

    with gzip.open(msa_path, "rt") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence_lines)))
                header = line[1:]
                sequence_lines = []
            else:
                sequence_lines.append(line.strip())

    if header is not None:
        records.append((header, "".join(sequence_lines)))
    return records


def ungapped_length(sequence):
    return sum(1 for char in sequence if char not in GAP_CHARS)


def parse_header(header):
    match = HEADER_RE.match(header)
    if not match:
        raise ValueError(f"Unsupported FASTA header format: >{header}")

    name = match.group("name")
    start = match.group("start")
    end = match.group("end")
    rest = match.group("rest") or ""

    if start is None:
        return name, None, None, rest
    return name, int(start), int(end), rest


def format_header(name, start, end, rest, force_coordinates):
    if not force_coordinates:
        return f"{name}{rest}"
    return f"{name}/{start}-{end}{rest}"


def trim_record(header, sequence, keep_start, keep_end, force_coordinates):
    name, start, end, rest = parse_header(header)
    had_coordinates = start is not None
    clipped_left = sequence[:keep_start]
    clipped_right = sequence[keep_end:]
    trimmed_sequence = sequence[keep_start:keep_end]

    left_residues = ungapped_length(clipped_left)
    right_residues = ungapped_length(clipped_right)

    if start is None:
        start = 1
        end = ungapped_length(sequence)

    new_start = start + left_residues
    new_end = end - right_residues

    new_header = format_header(
        name, new_start, new_end, rest, had_coordinates or force_coordinates
    )
    return new_header, trimmed_sequence


def write_fasta_gz(records, output_path):
    with gzip.open(output_path, "wt") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n{sequence}\n")


def trim_msa(msa_path, rf_path, output_path):
    rf = read_rf(rf_path)
    keep_start, keep_end = terminal_trim_bounds(rf)
    records = read_fasta_gz(msa_path)
    if not records:
        raise ValueError(f"{msa_path} does not contain any FASTA records")

    force_coordinates = keep_start > 0 or keep_end < len(rf)
    trimmed_records = []
    for header, sequence in records:
        if len(sequence) != len(rf):
            raise ValueError(
                f"{msa_path} record >{header} has alignment length {len(sequence)} "
                f"but {rf_path} has length {len(rf)}"
            )
        trimmed_records.append(
            trim_record(header, sequence, keep_start, keep_end, force_coordinates)
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_fasta_gz(trimmed_records, output_path)


def seed_msa_basename(msa_path):
    name = msa_path.name
    if name.endswith(".fas.gz"):
        return name[: -len(".fas.gz")]
    return msa_path.stem


def trim_seed_msa_dir(seed_msa_dir, rf_dir, output_dir):
    seed_msa_dir = Path(seed_msa_dir)
    rf_dir = Path(rf_dir)
    output_dir = Path(output_dir)

    msa_paths = sorted(seed_msa_dir.glob("*.fas.gz"))
    if not msa_paths:
        raise ValueError(f"No .fas.gz seed MSA files found in {seed_msa_dir}")

    for msa_path in msa_paths:
        basename = seed_msa_basename(msa_path)
        rf_path = rf_dir / f"{basename}.txt"
        if not rf_path.exists():
            raise FileNotFoundError(f"Missing RF file for {msa_path.name}: {rf_path}")
        trim_msa(msa_path, rf_path, output_dir / msa_path.name)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Trim terminal envelope columns from gzipped seed MSAs using matching RF "
            "dot/x files, and recalculate sequence coordinates."
        )
    )
    parser.add_argument("seed_msa_dir", help="Directory containing seed MSA .fas.gz files")
    parser.add_argument("rf_dir", help="Directory containing matching RF .txt files")
    parser.add_argument("output_dir", help="Directory for filtered seed MSA .fas.gz files")
    return parser.parse_args()


def main():
    args = parse_args()
    trim_seed_msa_dir(args.seed_msa_dir, args.rf_dir, args.output_dir)


if __name__ == "__main__":
    main()
