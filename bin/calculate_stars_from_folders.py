#!/usr/bin/env python3

import argparse
import pandas as pd
import csv
import os

def get_mgnifam_data(output_path):
    # Initialize lists and sets to hold data
    mgnifam_data = []
    
    # Define the file paths
    family_metadata_path = os.path.join(output_path, 'families', 'family_metadata.csv')
    pdb_scores_path = os.path.join(output_path, 'structures', 'pdb_scores.csv')

    # Read family_metadata.csv
    family_data = {}
    with open(family_metadata_path, mode='r') as family_file:
        reader = csv.reader(family_file)
        for row in reader:
            id = int(row[0])
            family_size = int(row[1])
            family_data[id] = family_size

    # Read pdb_scores.csv and match IDs
    with open(pdb_scores_path, mode='r') as pdb_file:
        reader = csv.reader(pdb_file)
        for row in reader:
            id = int(row[0])
            length = int(row[1])
            plddt = float(row[2])

            # Check if the ID exists in family_data
            if id in family_data:
                # Create a tuple in the format (id, plddt, family_size, rep_length)
                # rep_length is the length in this case
                family_size = family_data[id]
                mgnifam_data.append((id, plddt, family_size, length))

    return mgnifam_data

def get_mgnifam_pfams_ids(output_path):
    mgnifam_pfams_ids = set()  # Initialize an empty set

    # Define the file path
    annotated_models_path = os.path.join(output_path, 'hh', 'annotated_models.txt')

    # Read the annotated_models.txt file
    with open(annotated_models_path, mode='r') as file:
        for line in file:
            # Strip whitespace and convert the line to an integer
            mgnifam_id = int(line.strip())
            mgnifam_pfams_ids.add(mgnifam_id)  # Add the ID to the set

    return mgnifam_pfams_ids

def get_mgnifam_folds_ids(output_path):
    mgnifam_folds_ids = set()  # Initialize an empty set

    # Define the file path
    annotated_structures_path = os.path.join(output_path, 'structures', 'foldseek', 'annotated_structures.txt')

    # Read the annotated_models.txt file
    with open(annotated_structures_path, mode='r') as file:
        for line in file:
            # Strip whitespace and convert the line to an integer
            mgnifam_id = int(line.strip())
            mgnifam_folds_ids.add(mgnifam_id)  # Add the ID to the set

    return mgnifam_folds_ids

def calculate_quality_stars(plddt, family_size, rep_length):
    quality_stars = 0

    # Add stars based on plddt
    if (plddt is not None) and (plddt >= 70):
        quality_stars += 2
    elif (plddt is not None) and (plddt >= 50):
        quality_stars += 1

    # Add stars based on family_size
    if (family_size is not None) and (family_size >= 100):
        quality_stars += 1

    # Add stars based on rep_length
    if (rep_length is not None) and (rep_length >= 100):
        quality_stars += 1

    return quality_stars

def calculate_novelty_stars(mgnifam_id, mgnifam_pfams_ids, mgnifam_folds_ids):
    novelty_stars = 0

    # Add stars if mgnifam_id is not in mgnifam_pfams_ids
    if mgnifam_id not in mgnifam_pfams_ids:
        novelty_stars += 2

    # Add stars if mgnifam_id is not in mgnifam_folds_ids
    if mgnifam_id not in mgnifam_folds_ids:
        novelty_stars += 1

    return novelty_stars

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Calculate quality and novelty stars from the MGnifams output folders.")
    parser.add_argument("output_path", help="Path to the MGnifams output folder")

    # Parse command-line arguments
    args = parser.parse_args()

    get_mgnifam_data(args.output_path)

    # # Get the mgnifam data and the sets of mgnifam_ids from the other tables
    mgnifam_data      = get_mgnifam_data(args.output_path)
    mgnifam_pfams_ids = get_mgnifam_pfams_ids(args.output_path)
    mgnifam_folds_ids = get_mgnifam_folds_ids(args.output_path)

    # Prepare the data for the DataFrame
    ids = [row[0] for row in mgnifam_data]
    plddts = [row[1] for row in mgnifam_data]
    family_sizes = [row[2] for row in mgnifam_data]
    rep_lengths = [row[3] for row in mgnifam_data]
    quality_stars = [calculate_quality_stars(plddt, family_size, rep_length) for plddt, family_size, rep_length in zip(plddts, family_sizes, rep_lengths)]
    novelty_stars = [calculate_novelty_stars(mgnifam_id, mgnifam_pfams_ids, mgnifam_folds_ids) for mgnifam_id in ids]

    # Create the DataFrame
    df = pd.DataFrame({
        'MGnify Family ID': ids,
        'Quality Stars': quality_stars,
        'Novelty Stars': novelty_stars
    })

    # Print the DataFrame
    print(df)
    df.to_csv('output/stars.csv', index=False)

if __name__ == "__main__":
    main()
