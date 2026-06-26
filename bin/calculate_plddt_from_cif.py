#!/usr/bin/env python3

def extract_metrics_from_cif(cif_file):
    in_metric_local = False
    plddt_sum = 0
    count = 0

    # Open the CIF file and read line by line
    with open(cif_file, 'r') as f:
        for line in f:
            line = line.strip()

            # Detect when we are in the _ma_qa_metric_local section
            if line.startswith("loop_"):
                in_metric_local = False
            elif "_ma_qa_metric_local.ordinal_id" in line:
                in_metric_local = True
                continue

            # If we are inside the relevant block, process the lines
            if in_metric_local and len(line.split()) == 7:
                # Extract the pLDDT score (5th column)
                parts = line.split()
                plddt = float(parts[4])
                
                # Sum up the pLDDT values and increment the count
                plddt_sum += plddt
                count += 1

    # Calculate the average length (number of rows) and average pLDDT score
    if count > 0:
        average_length = count
        average_plddt = plddt_sum / count
    else:
        average_length = 0
        average_plddt = 0.0

    return average_length, average_plddt

# Example usage
cif_file = '/home/vangelis/Desktop/Projects/mgnifams/output/structures/cif/1.cif'
average_length, average_plddt = extract_metrics_from_cif(cif_file)
print(f"Average length: {average_length}")
print(f"Average pLDDT score: {average_plddt:.2f}")

# Call:
# ./bin/calculate_plddt_from_cif.py
