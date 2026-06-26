#!/bin/bash

# Navigate to the root directory containing the `books` folder
cd PANTHER19.0_hmmscoring/target/famlib/rel/PANTHER19.0_altVersion/hmmscoring/PANTHER19.0 || exit

# Output file to store all combined HMMs
output_file="panther_combined.hmm"

# Ensure the output file does not exist
rm -f "$output_file"

# Find and append all .hmm files
find books -type f -name "*.hmm" -exec cat {} + >> "$output_file"

echo "All .hmm files have been combined into $output_file"
