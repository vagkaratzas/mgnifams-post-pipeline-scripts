#!/usr/bin/env python3

import sys
import pandas as pd

def calculate_jaccard_index(set1, set2):
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    jaccard_index = len(intersection) / len(union)

    return jaccard_index

def create_aa_set(df):
    aa_set = set()

    # Iterate through each row in the DataFrame
    for index, row in df.iterrows():
        protein = row['protein']
        region  = row['region']
        
        # Split the region_str into individual regions
        regions = region.split('-')
        # Extract start and end numbers from the first and last region
        start_num = int(regions[0])
        end_num   = int(regions[-1])
        
        # Generate all combinations and add them to the set
        for num in range(start_num, end_num + 1):
            aa_set.add(f"{protein}:{num}")
            
    return aa_set

def calculate_aa_jaccard_index(fam1, fam2, fam_proteins, fam_metadata):
    fam1_rep = fam_metadata.loc[fam_metadata['id'] == fam1, 'protein_rep'].values[0]
    if (fam_proteins.loc[fam_proteins['mgnifam_id'] == fam2, 'protein'].isin([fam1_rep]).any()):
        fam1_subset = fam_proteins[(fam_proteins['mgnifam_id'] == fam1) & (fam_proteins['protein'] == fam1_rep)]
        fam2_subset = fam_proteins[(fam_proteins['mgnifam_id'] == fam2) & (fam_proteins['protein'] == fam1_rep)]
        rep_length_value = str(int(fam_metadata.loc[fam_metadata['id'] == fam1, 'rep_length'].values[0]))
    else:
        fam2_rep = fam_metadata.loc[fam_metadata['id'] == fam2, 'protein_rep'].values[0]
        if (fam_proteins.loc[fam_proteins['mgnifam_id'] == fam1, 'protein'].isin([fam2_rep]).any()):
            fam1_subset = fam_proteins[(fam_proteins['mgnifam_id'] == fam1) & (fam_proteins['protein'] == fam2_rep)]
            fam2_subset = fam_proteins[(fam_proteins['mgnifam_id'] == fam2) & (fam_proteins['protein'] == fam2_rep)]
            rep_length_value = str(int(fam_metadata.loc[fam_metadata['id'] == fam2, 'rep_length'].values[0]))
        else:
            return 0
    
    
    # When and if whole protein, replace with its total length
    fam1_subset.loc[fam1_subset['region'] == "-", 'region'] = f'1-{rep_length_value}'
    fam2_subset.loc[fam2_subset['region'] == "-", 'region'] = f'1-{rep_length_value}'
    aa_set1 = create_aa_set(fam1_subset)
    aa_set2 = create_aa_set(fam2_subset)
    aa_jaccard_index = calculate_jaccard_index(aa_set1, aa_set2)

    return aa_jaccard_index

def main(fam_proteins_file, fam_metadata_file, log_file, output_file):
    fam_proteins = pd.read_csv(fam_proteins_file)
    fam_metadata = pd.read_csv(fam_metadata_file)
    fam_metadata = fam_metadata[['id', 'protein_rep', 'rep_region', 'rep_length']]

    # Extract unique mgnifam_id values
    unique_mgnifam_ids = fam_proteins['mgnifam_id'].unique()

    # Convert to list
    unique_mgnifam_ids_list = list(unique_mgnifam_ids)

    # Iterate through unique mgnifam_ids with a nested loop to avoid redundant comparisons
    for i in range(len(unique_mgnifam_ids_list)):
        fam1 = unique_mgnifam_ids_list[i]
        set_fam1 = set(fam_proteins[fam_proteins['mgnifam_id'] == fam1]['protein'])

        for j in range(i + 1, len(unique_mgnifam_ids_list)):
            fam2 = unique_mgnifam_ids_list[j]
            set_fam2 = set(fam_proteins[fam_proteins['mgnifam_id'] == fam2]['protein'])

            with open(log_file, 'a') as file:
                file.write(f"\n{fam1} {fam2} ")

            # Calculate Jaccard index
            jaccard_index = calculate_jaccard_index(set_fam1, set_fam2)
            if (jaccard_index >= 0.5):
                with open(log_file, 'a') as file:
                    file.write(f"jaccard_index = {jaccard_index} ")

                aa_jaccard_index = calculate_aa_jaccard_index(fam1, fam2, fam_proteins, fam_metadata)
                if (aa_jaccard_index >= 0.5):
                    with open(log_file, 'a') as file:
                        file.write(f"aa_jaccard_index = {aa_jaccard_index}")
                    with open(output_file, 'a') as file:
                        file.write(f"{fam1},{fam2},{aa_jaccard_index}\n")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python calculate_fam_similarities.py <fam_proteins_file> <fam_metadata_file> <log_file> <output_file>")
        sys.exit(1)

    # Local: python ./bin/calculate_fam_similarities.py assets/test_data/test_output/tables/mgnifam_proteins.csv assets/test_data/test_output/tables/mgnifam.csv output/log.txt output/similarities.csv

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
