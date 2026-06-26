#!/usr/bin/env python3

import sys
import os
import pandas as pd

def map_file_basenames(folder, mapping_file):
    # Create a dictionary from the mapping file (column 1 -> column 2)
    mapping_df = pd.read_csv(mapping_file, header=None, names=['old_id', 'new_id'], dtype=str)
    id_mapping = dict(zip(mapping_df['old_id'], mapping_df['new_id']))

    # Iterate over all files in the folder
    for filename in os.listdir(folder):
        basename, ext = os.path.splitext(filename)  # Get the file's basename and extension

        # Check if the basename is in the mapping
        if basename in id_mapping:
            new_basename = id_mapping[basename]  # Get the new basename from the mapping
            new_filename = new_basename + ext  # Keep the original extension

            old_filepath = os.path.join(folder, filename)
            new_filepath = os.path.join(folder, new_filename)

            # Rename the file
            os.rename(old_filepath, new_filepath)
            print(f"Renamed: {old_filepath} -> {new_filepath}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python map_file_ids.py <folder> <mapping_file>")
        sys.exit(1)

    # Local: ./bin/map_file_ids.py assets/test_data/structures/cif output/redundancy_check/mapping_file.csv

    folder = sys.argv[1]
    mapping_file = sys.argv[2]

    # Map and rename file basenames based on the mapping file
    map_file_basenames(folder, mapping_file)
