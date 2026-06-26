#!/usr/bin/env python3
import argparse
import os
import glob
import re

def check_versions_file(subfolder_path):
    versions_file_path = os.path.join(subfolder_path, 'versions.yml')
    if os.path.isfile(versions_file_path):
        try:
            if "PREDICT_STRUCTURES:ESMFOLD" in open(versions_file_path).read() or \
                "PREDICT_STRUCTURES:ESMFOLD_CPU" in open(versions_file_path).read():
                return True
        except Exception as e:
            print(f"An error occurred while reading {versions_file_path}: {e}")
    return False

def extract_scores(second_level_subfolder, output_file):
    # Find all files matching the pattern *_scores.txt in the second_level_subfolder
    score_files = glob.glob(os.path.join(second_level_subfolder, '*_scores.txt'))

    # Define a regex pattern to capture the relevant parts from the line
    pattern = re.compile(
        r"Predicted structure for (\d+) with length (\d+), pLDDT ([\d.]+), pTM ([\d.]+)"
    )

    # Open the output_file in append mode
    with open(output_file, 'a') as out_f:
        # Iterate through each file
        for score_file in score_files:
            # Read the content of the file
            with open(score_file, 'r') as in_f:
                for line in in_f:
                    # Search for the pattern in the line
                    match = pattern.search(line)
                    if match:
                        # Extract the name, length, pLDDT, and pTM
                        name, length, plddt, ptm = match.groups()
                        # Format the extracted data
                        formatted_data = f"{name},{length},{plddt},{ptm}\n"
                        # Append the formatted data to the output file
                        out_f.write(formatted_data)

def list_valid_subfolders(folder_path, output_file):
    valid_subfolders = []
    try:
        # Get first level subfolders
        first_level_subfolders = [f.name for f in os.scandir(folder_path) if f.is_dir()]
        for subfolder in first_level_subfolders:
            subfolder_path = os.path.join(folder_path, subfolder)
            # Get second level subfolders
            second_level_subfolders = [f.path for f in os.scandir(subfolder_path) if f.is_dir()]
            for second_level_subfolder in second_level_subfolders:
                if check_versions_file(second_level_subfolder):
                    valid_subfolders.append(second_level_subfolder)
                    extract_scores(second_level_subfolder, output_file)
    except Exception as e:
        print(f"An error occurred: {e}")
    return valid_subfolders

def main():
    parser = argparse.ArgumentParser(description="List valid second level subfolder paths of a given folder based on specific criteria.")
    parser.add_argument("folder_path", type=str, help="Path to the folder")
    parser.add_argument("output_file", type=str, help="Path to the output file")
    # Example call local: ./bin/extract_pdb_scores.py /home/vangelis/Desktop/Projects/mgnifams/work output/pdb_scores.csv
    
    args = parser.parse_args()
    folder_path = args.folder_path
    output_file = args.output_file

    if os.path.isdir(folder_path):
        valid_subfolders = list_valid_subfolders(folder_path, output_file)
        print("Valid second level subfolders containing the specific versions.yml:")
        for subfolder in valid_subfolders:
            print(subfolder)
    else:
        print("The provided path is not a valid directory.")

if __name__ == "__main__":
    main()
