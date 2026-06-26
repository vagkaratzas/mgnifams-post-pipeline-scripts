#!/usr/bin/env python3

import sys
import os
import pandas as pd

def remove_redundant_ids(csv_file, redundant_ids_file, id_column_position, outfile):
    # Read the list of redundant IDs from the text file
    with open(redundant_ids_file, 'r') as f:
        redundant_ids = set(int(line.strip()) for line in f if line.strip())

    # Read the CSV file into a DataFrame, forcing all columns as string initially
    df = pd.read_csv(csv_file, dtype=str)

    # Convert the specified column (by position) to integer
    df.iloc[:, id_column_position] = pd.to_numeric(df.iloc[:, id_column_position], downcast='integer')

    # Remove rows where the column at the given position matches any of the redundant IDs
    filtered_df = df[~df.iloc[:, id_column_position].isin(redundant_ids)]

    # Write the filtered DataFrame to a new CSV file
    filtered_df.to_csv(outfile, index=False)
    print(f"Filtered data has been written to {outfile}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python remove_ids_from_csv.py <csv_file> <redundant_ids_file> <id_column_position> <outfile>")
        sys.exit(1)

    # Local: ./bin/remove_loc_ids_from_csv.py assets/test_data/structures/pdb_scores.csv output/redundancy_check/redundant_fam_ids.txt 0 output/redo/pdb_scores_rem.csv

    csv_file = sys.argv[1]
    redundant_ids_file = sys.argv[2]
    id_column_position = int(sys.argv[3])  # Convert the passed column position to an integer
    outfile = sys.argv[4]

    directory_path = os.path.dirname(outfile)
    os.makedirs(directory_path, exist_ok=True)

    remove_redundant_ids(csv_file, redundant_ids_file, id_column_position, outfile)
