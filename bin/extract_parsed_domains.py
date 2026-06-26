#!/usr/bin/env python3
import os
import argparse
import shutil
import glob

def check_versions_file(subfolder_path):
    versions_file_path = os.path.join(subfolder_path, '.command.run')
    if os.path.isfile(versions_file_path):
        try:
            if "EXPORT_DB_PARSE_DOMAINS" in open(versions_file_path).read():
                return True
        except Exception as e:
            print(f"An error occurred while reading {versions_file_path}: {e}")
    return False

def copy_files(subfolder_path, output_dir):
    # Define the path to the 'domain_results' subfolder
    domain_results_path = os.path.join(subfolder_path, 'domain_results')
    
    # Ensure the output directory exists, create it if it doesn't
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Find all files in the domain_results subfolder
    files = glob.glob(os.path.join(domain_results_path, '*'))

    # Iterate over each file and copy it to the output directory
    for file in files:
        if os.path.isfile(file):  # Ensure it's a file (not a directory)
            shutil.copy(file, output_dir)

def list_valid_subfolders(folder_path, output_dir):
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
                    copy_files(second_level_subfolder, output_dir)
    except Exception as e:
        print(f"An error occurred: {e}")
    return valid_subfolders

def main():
    parser = argparse.ArgumentParser(description="List valid second level subfolder paths of a given folder based on specific criteria.")
    parser.add_argument("folder_path", type=str, help="Path to the folder")
    parser.add_argument("output_dir", type=str, help="Path to the output directory where files will be copied")
    # Example call local: ./bin/extract_parsed_domains.py /home/vangelis/Desktop/Projects/mgnifams/work output/domain_results
    # local needs change to EXPORT_DB:PARSE_DOMAINS
    
    args = parser.parse_args()
    folder_path = args.folder_path
    output_dir = args.output_dir

    if os.path.isdir(folder_path):
        valid_subfolders = list_valid_subfolders(folder_path, output_dir)
        print("Valid second level subfolders containing the specific .command.run:")
        for subfolder in valid_subfolders:
            print(subfolder)
    else:
        print("The provided path is not a valid directory.")

if __name__ == "__main__":
    main()
