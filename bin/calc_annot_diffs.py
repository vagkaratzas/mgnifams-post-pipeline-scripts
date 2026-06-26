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
                # Move the current start to after the overlapping portion
                current_start = max(current_start, end2 + 1)
        # If there is any remaining range, add it
        if current_start <= current_end:
            diff_ranges.append((current_start, current_end))
    return diff_ranges

def range_overlap(ranges1, ranges2):
    overlap_ranges = []
    for start1, end1 in ranges1:
        for start2, end2 in ranges2:
            if start2 <= end1 and end2 >= start1:
                overlap_ranges.append((max(start1, start2), min(end1, end2)))
    return overlap_ranges

def compare_annotations(file1, file2):
    data1 = pd.read_csv(file1)
    data2 = pd.read_csv(file2)

    protein_ranges1 = defaultdict(list)
    protein_ranges2 = defaultdict(list)

    for _, row in data1.iterrows():
        protein_ranges1[row['protein']].extend(parse_ranges(row['region']))

    for _, row in data2.iterrows():
        protein_ranges2[row['protein']].extend(parse_ranges(row['region']))

    results = []

    for protein, ranges1 in protein_ranges1.items():
        merged_ranges1 = merge_ranges(ranges1)
        coverage1 = get_total_coverage(merged_ranges1)
        
        if protein in protein_ranges2:
            merged_ranges2 = merge_ranges(protein_ranges2[protein])
            coverage2 = get_total_coverage(merged_ranges2)
            
            # Unique and shared coverages
            diff_1_in_2 = range_difference(merged_ranges1, merged_ranges2)
            exclusive_coverage_1 = get_total_coverage(diff_1_in_2)
            
            diff_2_in_1 = range_difference(merged_ranges2, merged_ranges1)
            exclusive_coverage_2 = get_total_coverage(diff_2_in_1)
            
            overlap_ranges = range_overlap(merged_ranges1, merged_ranges2)
            overlap_coverage = get_total_coverage(overlap_ranges)

            results.append({
                "protein": protein,
                "total_coverage_file1": coverage1,
                "total_coverage_file2": coverage2,
                "exclusive_coverage_file1": exclusive_coverage_1,
                "exclusive_coverage_file2": exclusive_coverage_2,
                "overlap_coverage": overlap_coverage
            })
        else:
            results.append({
                "protein": protein,
                "total_coverage_file1": coverage1,
                "total_coverage_file2": 0,
                "exclusive_coverage_file1": coverage1,
                "exclusive_coverage_file2": 0,
                "overlap_coverage": 0
            })

    for protein, ranges2 in protein_ranges2.items():
        if protein not in protein_ranges1:
            merged_ranges2 = merge_ranges(ranges2)
            coverage2 = get_total_coverage(merged_ranges2)
            results.append({
                "protein": protein,
                "total_coverage_file1": 0,
                "total_coverage_file2": coverage2,
                "exclusive_coverage_file1": 0,
                "exclusive_coverage_file2": coverage2,
                "overlap_coverage": 0
            })

    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(description="Compare annotations between two files.")
    parser.add_argument("annot_file_1", help="Path to the first annotated file")
    parser.add_argument("annot_file_2", help="Path to the second annotated file")
    parser.add_argument("output_csv", help="Path to the output csv")
    parser.add_argument("output_stats_incr", help="Path to the output file with stats increase")

    args = parser.parse_args()

    # Example:
    # Local: ./bin/calc_annot_diffs.py assets/test_data/annotations/pfam_annots_test.csv assets/test_data/annotations/mgnifam_annots.csv output/annot_diff.csv output/stats_incr.txt

    print("Starting comparisons")
    comparison_df = compare_annotations(args.annot_file_1, args.annot_file_2)
    # print(comparison_df)
    print("Printing out CSV")
    comparison_df.to_csv(args.output_csv, index=False)

    print("Calculating the total amino acid increase and increase ratio")
    # Calculate sums of all relevant columns and increase ratio
    total_coverage_file1_sum = comparison_df["total_coverage_file1"].sum()
    total_coverage_file2_sum = comparison_df["total_coverage_file2"].sum()
    exclusive_coverage_file1_sum = comparison_df["exclusive_coverage_file1"].sum()
    exclusive_coverage_file2_sum = comparison_df["exclusive_coverage_file2"].sum()
    overlap_coverage_sum = comparison_df["overlap_coverage"].sum()
    # Increase ratio
    increase_ratio = ((exclusive_coverage_file2_sum + total_coverage_file1_sum) / total_coverage_file1_sum) * 100 if total_coverage_file1_sum else 0

    # Write the summed statistics to the specified output file
    with open(args.output_stats_incr, "w") as stats_file:
        stats_file.write("Summed Statistics:\n")
        stats_file.write("Total Coverage File 1: {}\n".format(total_coverage_file1_sum))
        stats_file.write("Total Coverage File 2: {}\n".format(total_coverage_file2_sum))
        stats_file.write("Exclusive Coverage File 1: {}\n".format(exclusive_coverage_file1_sum))
        stats_file.write("Exclusive Coverage File 2: {}\n".format(exclusive_coverage_file2_sum))
        stats_file.write("Overlap Coverage: {}\n".format(overlap_coverage_sum))
        stats_file.write("Amino Acid Level Increase: {} from initial {} total_coverage_file1\n".format(
            exclusive_coverage_file2_sum, total_coverage_file1_sum))
        stats_file.write("Increase Ratio: {:.2f}%\n".format(increase_ratio))

if __name__ == "__main__":
    main()
