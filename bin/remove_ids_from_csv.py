#!/usr/bin/env python3

import sys
import os
import pandas as pd

def remove_redundant_ids(csv_file, redundant_ids_file, id_column, outfile):
    # Read the list of redundant IDs from the text file
    with open(redundant_ids_file, 'r') as f:
        redundant_ids = set(int(line.strip()) for line in f if line.strip())

    # Read the CSV file into a DataFrame, forcing all columns as string except for the id_column
    df = pd.read_csv(csv_file, dtype=str)  # Set all columns to string initially
    df[id_column] = pd.to_numeric(df[id_column], downcast='integer')  # Convert the id_column to integer

    # Remove rows where the id_column matches any of the redundant IDs
    filtered_df = df[~df[id_column].isin(redundant_ids)]

    # Write the filtered DataFrame to a new CSV file
    filtered_df.to_csv(outfile, index=False)
    print(f"Filtered data has been written to {outfile}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python remove_ids_from_csv.py <csv_file> <redundant_ids_file> <id_column> <outfile>")
        sys.exit(1)

    # Local: ./bin/remove_ids_from_csv.py assets/test_data/tables/mgnifam_pfams.csv output/redundancy_check/redundant_fam_ids.txt mgnifam_id output/redo/mgnifam_pfams_rem.csv

    csv_file = sys.argv[1]
    redundant_ids_file = sys.argv[2]
    id_column = sys.argv[3]
    outfile = sys.argv[4]

    directory_path = os.path.dirname(outfile)
    os.makedirs(directory_path, exist_ok=True)

    remove_redundant_ids(csv_file, redundant_ids_file, id_column, outfile)
