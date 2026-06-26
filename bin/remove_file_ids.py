#!/usr/bin/env python3

import sys
import os
import pandas as pd

def remove_files_with_redundant_ids(folder, redundant_ids_file):
    # Read the list of redundant IDs from the text file
    with open(redundant_ids_file, 'r') as f:
        redundant_ids = set(line.strip() for line in f if line.strip())

    # Iterate over all files in the folder
    for filename in os.listdir(folder):
        basename, _ = os.path.splitext(filename)  # Get the file's basename without the extension

        if basename in redundant_ids:
            file_path = os.path.join(folder, filename)
            os.remove(file_path)
            print(f"Removed file: {file_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python remove_file_ids.py <folder> <redundant_ids_file>")
        sys.exit(1)

    # Local: ./bin/remove_file_ids.py assets/test_data/structures/cif output/redundancy_check/redundant_fam_ids.txt

    folder = sys.argv[1]
    redundant_ids_file = sys.argv[2]

    # Remove files whose basenames match redundant IDs
    remove_files_with_redundant_ids(folder, redundant_ids_file)
