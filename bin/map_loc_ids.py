#!/usr/bin/env python3

import sys
import os
import pandas as pd

def translate_ids(csv_file, mapping_file, id_column_position, outfile):
    # Create a dictionary from the mapping file (column 1 -> column 2), ensuring they are integers
    mapping_df = pd.read_csv(mapping_file, header=None, names=['old_id', 'new_id'], dtype=int)
    id_mapping = dict(zip(mapping_df['old_id'], mapping_df['new_id']))

    df = pd.read_csv(csv_file, dtype=str)  # Set all columns to string initially
    df.iloc[:, id_column_position] = pd.to_numeric(df.iloc[:, id_column_position], downcast='integer')  # Convert the id_column to integer

    # Translate the id_column using the mapping
    df.iloc[:, id_column_position] = df.iloc[:, id_column_position].map(id_mapping).fillna(df.iloc[:, id_column_position]).astype(int)

    # Write the translated DataFrame to a new CSV file
    df.to_csv(outfile, index=False)
    print(f"Translated data has been written to {outfile}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python translate_ids.py <csv_file> <mapping_file> <id_column_position> <outfile>")
        sys.exit(1)

    # Local: ./bin/map_loc_ids.py output/redo/pdb_scores_rem.csv output/redundancy_check/mapping_file.csv 0 output/redo/pdb_scores_mapped.csv

    csv_file = sys.argv[1]
    mapping_file = sys.argv[2]
    id_column_position = int(sys.argv[3])  # Convert the passed column position to an integer
    outfile = sys.argv[4]

    # Create output directory if it doesn't exist
    directory_path = os.path.dirname(outfile)
    os.makedirs(directory_path, exist_ok=True)

    # Translate IDs in the CSV file based on the mapping file
    translate_ids(csv_file, mapping_file, id_column_position, outfile)
