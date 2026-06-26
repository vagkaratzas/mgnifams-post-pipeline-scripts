#!/usr/bin/env python3

import argparse
import csv
import gzip
import os
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


def extract_fasta(csv_file, fasta_file):
    """
    Extracts a FASTA file from a CSV file where the protein names are in the first column and the sequences are in the second.

    Args:
        csv_file (str): The path to the input CSV file.
        fasta_file (str): The path to the output FASTA file.
    """
    # Ensure the CSV file exists
    if not os.path.exists(csv_file):
        raise FileNotFoundError(f"The file {csv_file} does not exist.")

    with open_text(csv_file, newline='') as infile, open(fasta_file, 'w', encoding='utf-8') as outfile:
        reader = csv.DictReader(infile)
        
        for row in reader:
            protein_id = row['mgyp']
            sequence = row['sequence']
            
            # Skip if sequence is empty
            if not sequence:
                continue

            # Write to the FASTA file in the correct format
            outfile.write(f">{protein_id}\n")
            outfile.write(f"{sequence}\n")


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Extract FASTA file from a CSV file.")
    parser.add_argument("csv_file", help="Path to the input CSV file")
    parser.add_argument("fasta_file", help="Path to the output FASTA file")
    # Example call: ./bin/extract_fasta_from_proteins_csv.py assets/test_data/sequence_explorer_protein_test_10001.csv output/test.fasta

    # Parse the command-line arguments
    args = parser.parse_args()

    # Run the extraction
    extract_fasta(args.csv_file, args.fasta_file)


if __name__ == "__main__":
    main()
