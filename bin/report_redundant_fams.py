#!/usr/bin/env python3

import sys
import os
import csv
import pandas as pd
import gc

def read_fam_to_size_map(fam_metadata_file):
    fam_to_size = {} # init dict

    # Open the CSV file and read it
    with open(fam_metadata_file, mode='r') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            # Store the id as the key and family_size as the value
            fam_to_size[int(row['id'])] = int(row['family_size'])
    
    return fam_to_size

def remove_redundant(df_redundant, df_similar, fam_to_size, redundant_fam_ids_file):
    # Open the file for writing redundant family IDs
    with open(redundant_fam_ids_file, 'w') as file:
        
        # Keep iterating while df_redundant has rows
        while not df_redundant.empty:
            print(f'Remaining rows: {len(df_redundant)}\n')
            # Get the last row of df_redundant
            last_row = df_redundant.iloc[-1]
            fam1 = last_row['fam1']
            fam2 = last_row['fam2']
            
            # Compare family sizes using the fam_to_size dictionary
            size1 = fam_to_size[fam1]
            size2 = fam_to_size[fam2]
            
            # Determine the redundant family (the one with the smaller size, or larger ID if sizes are equal)
            if size1 < size2:
                redundant_fam = fam1
            elif size1 > size2:
                redundant_fam = fam2
            else:  # sizes are the same
                redundant_fam = max(fam1, fam2)  # choose the one with the larger ID
            
            # Write the redundant family to the output file
            file.write(f"{int(redundant_fam)}\n")
            
            # Remove rows from df_redundant and df_similar where fam1 or fam2 matches redundant_fam
            df_redundant = df_redundant[(df_redundant['fam1'] != redundant_fam) & (df_redundant['fam2'] != redundant_fam)]
            df_similar = df_similar[(df_similar['fam1'] != redundant_fam) & (df_similar['fam2'] != redundant_fam)]
    
    return df_similar


def main(similarities_file, fam_metadata_file, redundant_fam_ids_file, non_redundant_similarities_file):
    fam_to_size = read_fam_to_size_map(fam_metadata_file)
    similarities_df = pd.read_csv(similarities_file, header=None, names=['fam1', 'fam2', 'similarity'])
    df_similar = similarities_df[similarities_df['similarity'] < 0.95]
    df_redundant = similarities_df[similarities_df['similarity'] >= 0.95]

    del similarities_df
    gc.collect()

    df_similar = remove_redundant(df_redundant, df_similar, fam_to_size, redundant_fam_ids_file)
    df_similar.to_csv(non_redundant_similarities_file, index=False)

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: ./bin/report_redundant_fams.py <similarities_file> <fam_metadata_file> <redundant_fam_ids_file> <non_redundant_similarities_file>")
        sys.exit(1)
    
    # Local: ./bin/report_redundant_fams.py assets/test_data/similarities.csv assets/test_data/mgnifam.csv output/redundancy_check/redundant_fam_ids.txt output/redundancy_check/non_redundant_similarities.csv 

    directory_path = os.path.dirname(sys.argv[3])
    os.makedirs(directory_path, exist_ok=True)
    directory_path = os.path.dirname(sys.argv[4])
    os.makedirs(directory_path, exist_ok=True)

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
