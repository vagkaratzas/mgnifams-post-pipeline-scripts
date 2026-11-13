#!/usr/bin/env python3

import os
import sys
from Bio import SeqIO

def extract_first_sequence(msa_folder, output_fasta):
    with open(output_fasta, "w") as out_file:
        for root, _, files in os.walk(msa_folder):
            counter = 0
            for file in files:
                if file.endswith(".fas"):  # Assuming .fas is the extension for FASTA files
                    counter += 1
                    if counter % 500 == 0:
                        print(counter)
                    file_path = os.path.join(root, file)
                    basename = os.path.splitext(os.path.basename(file_path))[0]  # Get basename without extension
                    
                    # Check if file is empty
                    if os.path.getsize(file_path) == 0:
                        print(f"Empty file: {basename}")
                        continue  # Skip to the next file

                    with open(file_path, "r") as msa_file:
                        records = SeqIO.parse(msa_file, "fasta")
                        first_record = next(records)  # Get the first record
                        sequence = ''.join(filter(str.isalpha, str(first_record.seq))).upper()
                        # Write the basename and first sequence ID in FASTA format
                        out_file.write(f">{basename} {first_record.id}\n{sequence}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: extract_first_seq.py <msa_folder> <output_fasta>")
        sys.exit(1)

    # Local: ./bin/extract_reps_fasta_from_msa.py assets/test_data/test_output/families/msa output/mgnifams_reps.fasta

    msa_folder = sys.argv[1]
    output_fasta = sys.argv[2]

    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)

    extract_first_sequence(msa_folder, output_fasta)
