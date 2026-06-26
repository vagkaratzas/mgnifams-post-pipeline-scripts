#!/usr/bin/env python3

import argparse
import os
import pandas as pd

# Function to remove clusters with only one member (singletons)
def remove_singletons(clusters_df):
    # Group by cluster_rep and filter out clusters with size 1
    cluster_sizes = clusters_df.groupby('cluster_rep').size().reset_index(name='size')
    non_singletons = cluster_sizes[cluster_sizes['size'] > 1]['cluster_rep']
    
    # Filter the original clusters_df to keep only non-singleton clusters
    clusters_df = clusters_df[clusters_df['cluster_rep'].isin(non_singletons)]
    
    return clusters_df

# Function to filter mgnifam_proteins based on cluster representative, processing line by line
def filter_proteins_line_by_line(clusters_df, proteins_file):
    # Step 1: Split second column by "_" and keep the first part
    clusters_df['sequence_id'] = clusters_df['sequence_id'].str.split('_').str[0]

    # Step 2: Open proteins file and process line by line
    with open(proteins_file, 'r') as pf:
        # Skip header
        next(pf)
        for line in pf:
            # Read the protein data
            mgnifam_id, protein, region = line.strip().split(',')
            
            # Get the cluster rep for this protein if it exists
            seq_id = protein
            if seq_id in clusters_df['sequence_id'].values:
                # Get all cluster reps for this sequence ID
                cluster_reps = clusters_df.loc[clusters_df['sequence_id'] == seq_id, 'cluster_rep'].unique()
                # Filter out all rows in clusters_df with these cluster reps
                clusters_df = clusters_df[~clusters_df['cluster_rep'].isin(cluster_reps)]

    return clusters_df

# Function to report remaining cluster reps and their sizes
def report_clusters(clusters_df):
    # Group by cluster_rep and get sizes
    cluster_sizes = clusters_df.groupby('cluster_rep').size().reset_index(name='size')
    # Sorting by descending order of cluster size
    cluster_sizes = cluster_sizes.sort_values(by='size', ascending=False)
    return cluster_sizes

def main(clusters_file, proteins_file, output_file):
    # Step 1: Read input clusters file
    print("Reading input clusters tsv")
    clusters_df = pd.read_csv(clusters_file, sep='\t', header=None, names=['cluster_rep', 'sequence_id'])

    print("Removing singletons")
    # Step 2: Remove singletons
    clusters_df = remove_singletons(clusters_df)

    print("Filtering out clusters")
    # Step 3: Filter proteins line by line based on clusters
    clusters_df = filter_proteins_line_by_line(clusters_df, proteins_file)

    print("Ordering remaining clusters")
    # Step 4: Report remaining cluster representatives and their cluster sizes
    remaining_clusters = report_clusters(clusters_df)

    print("Writing out")
    # Step 5: Write result to the output CSV file
    remaining_clusters.to_csv(output_file, index=False, header=['rep_name', 'cluster_size'])

    # Confirmation message
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description='Process clustering and protein files and output remaining clusters.')
    parser.add_argument('clusters_file', type=str, help='Path to the clustering file (linclust_clusters.tsv)')
    parser.add_argument('proteins_file', type=str, help='Path to the proteins file (mgnifam_proteins.csv)')
    parser.add_argument('output_file', type=str, help='Path to save the output CSV file (rep_name, cluster_size)')
    
    args = parser.parse_args()
    # Examples
    # local: ./bin/identify_unmgnifamed_space.py assets/test_data/test_output/mmseqs/linclust_clusters.tsv assets/test_data/test_output/tables/mgnifam_proteins.csv output/unexplored_cluster_space.csv

    os.makedirs(os.path.dirname(args.output_file), exist_ok=True)

    # Run the main function
    main(args.clusters_file, args.proteins_file, args.output_file)
