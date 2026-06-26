#!/usr/bin/env python3

import sys
import os
import pandas as pd

def translate_ids(csv_file, mapping_file, id_column, outfile):
    # Create a dictionary from the mapping file (column 1 -> column 2), ensuring they are integers
    mapping_df = pd.read_csv(mapping_file, header=None, names=['old_id', 'new_id'], dtype=int)
    id_mapping = dict(zip(mapping_df['old_id'], mapping_df['new_id']))

    df = pd.read_csv(csv_file, dtype=str)  # Set all columns to string initially
    df[id_column] = pd.to_numeric(df[id_column], downcast='integer')  # Convert the id_column to integer

    # Translate the id_column using the mapping
    df[id_column] = df[id_column].map(id_mapping).fillna(df[id_column]).astype(int)

    # Write the translated DataFrame to a new CSV file
    df.to_csv(outfile, index=False)
    print(f"Translated data has been written to {outfile}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python translate_ids.py <csv_file> <mapping_file> <id_column> <outfile>")
        sys.exit(1)

    # Local: ./bin/map_ids.py output/redo/mgnifam_pfams_rem.csv output/redundancy_check/mapping_file.csv mgnifam_id output/redo/mgnifam_pfams_mapped.csv

    csv_file = sys.argv[1]
    mapping_file = sys.argv[2]
    id_column = sys.argv[3]
    outfile = sys.argv[4]

    # Create output directory if it doesn't exist
    directory_path = os.path.dirname(outfile)
    os.makedirs(directory_path, exist_ok=True)

    # Translate IDs in the CSV file based on the mapping file
    translate_ids(csv_file, mapping_file, id_column, outfile)
