#!/usr/bin/env python3

import sys
import pandas as pd

def replace_file_basenames_with_id(input_csv, output_csv):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(input_csv)

    # List of columns to update with id-based basenames
    file_columns = [
        'cif_file', 'seed_msa_file', 'msa_file', 'hmm_file',
        'rf_file', 'biomes_file', 'domain_architecture_file'
    ]

    # Iterate through each row and update the file column basenames with the id
    for index, row in df.iterrows():
        id_value = str(row['id'])  # Convert the id to string for filename replacement
        for col in file_columns:
            extension = row[col].split('.')[-1]  # Get the file extension
            df.at[index, col] = f"{id_value}.{extension}"  # Replace with id + original extension

    # Save the modified DataFrame to a new CSV file
    df.to_csv(output_csv, index=False)
    print(f"Updated CSV has been written to {output_csv}")

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: python map_mgnifam_filepaths.py <mgnifam.csv> <outfile>")
        sys.exit(1)

    # Local: ./bin/map_mgnifam_filepaths.py output/redo/mgnifam_mapped.csv output/redo/mgnifam_mapped_correct.csv

    mgnifam_file = sys.argv[1]
    outfile = sys.argv[2]

    replace_file_basenames_with_id(mgnifam_file, outfile)
