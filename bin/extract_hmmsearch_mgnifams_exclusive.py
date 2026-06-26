#!/usr/bin/env python3

import argparse
import pandas as pd
import re
from collections import defaultdict

def parse_ranges(region_str):
    return [(int(start), int(end)) for start, end in re.findall(r'(\d+)-(\d+)', region_str)]

def merge_ranges(ranges):
    sorted_ranges = sorted(ranges)
    merged_ranges = []
    current_start, current_end = sorted_ranges[0]
    
    for start, end in sorted_ranges[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            merged_ranges.append((current_start, current_end))
            current_start, current_end = start, end
    merged_ranges.append((current_start, current_end))
    
    return merged_ranges

def get_total_coverage(merged_ranges):
    return sum(end - start + 1 for start, end in merged_ranges)

def range_difference(ranges1, ranges2):
    diff_ranges = []
    for start1, end1 in ranges1:
        current_start, current_end = start1, end1
        for start2, end2 in ranges2:
            if start2 <= current_end and end2 >= current_start:  # Overlap exists
                if current_start < start2:  # Add non-overlapping left portion
                    diff_ranges.append((current_start, start2 - 1))
                current_start = max(current_start, end2 + 1)
        if current_start <= current_end:
            diff_ranges.append((current_start, current_end))
    return diff_ranges

def main():
    parser = argparse.ArgumentParser(description="Compute exclusive annotations for File2.")
    parser.add_argument("file1", help="Path to the first annotated file")
    parser.add_argument("file2", help="Path to the second annotated file")
    parser.add_argument("output_file", help="Path to the output CSV file")

    args = parser.parse_args()

    # Read input files
    file1 = pd.read_csv(args.file1)
    file2 = pd.read_csv(args.file2)

    # Group regions by protein
    protein_ranges1 = defaultdict(list)
    protein_ranges2 = defaultdict(list)

    for _, row in file1.iterrows():
        protein_ranges1[row["protein"]].extend(parse_ranges(row["region"]))

    for _, row in file2.iterrows():
        protein_ranges2[row["protein"]].extend(parse_ranges(row["region"]))

    results = []

    for protein, ranges2 in protein_ranges2.items():
        merged_ranges2 = merge_ranges(ranges2)
        if protein in protein_ranges1:
            merged_ranges1 = merge_ranges(protein_ranges1[protein])
            diff_2_in_1 = range_difference(merged_ranges2, merged_ranges1)
        else:
            diff_2_in_1 = merged_ranges2  # All ranges are exclusive if protein not in File1

        for start, end in diff_2_in_1:
            results.append({
                "mgnifam_id": file2[file2["protein"] == protein]["mgnifam_id"].iloc[0],
                "protein": protein,
                "region": f"{start}-{end}",
                "annotation_length": end - start + 1
            })

    # Create DataFrame, sort by annotation_length in descending order, and write to CSV
    output_df = pd.DataFrame(results)
    output_df = output_df.sort_values(by="annotation_length", ascending=False)
    output_df.to_csv(args.output_file, index=False)

if __name__ == "__main__":
    main()

    # Local: ./bin/extract_hmmsearch_mgnifams_exclusive.py assets/test_data/annotations/pfam_annots_test.csv assets/test_data/annotations/mgnifam_annots.csv output/mgnifam_exclusive_annots.csv
