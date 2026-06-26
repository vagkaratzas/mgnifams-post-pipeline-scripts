#!/usr/bin/env python3

import os
import json
import argparse

# Function to process JSON files in the given folder
def process_json_files(folder_path):
    # Iterate over all files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith(".json"):
            file_path = os.path.join(folder_path, filename)
            basename = os.path.splitext(filename)[0]  # Get the basename without extension
            
            # Convert basename to an integer and format it to 10 digits
            basename_int = int(basename)
            formatted_basename = f"{basename_int:010d}"  # Ensure basename is 10 digits long
            
            # Open and read the JSON file
            with open(file_path, 'r') as f:
                data = json.load(f)

            # Iterate through the architecture_containers and domains
            for architecture in data.get("architecture_containers", []):
                for domain in architecture.get("domains", []):
                    if "MGnifam" in domain.get("name", ""):
                        # Modify the domain's id and name
                        domain['id'] = basename
                        domain['name'] = "MGnifam" + basename

                        # Modify the link to follow the new pattern
                        domain['link'] = f"http://mgnifams-demo.mgnify.org/details/?id=MGYF{formatted_basename}"

            # Save the modified data back to the file
            with open(file_path, 'w') as f:
                json.dump(data, f, indent=4)

            print(f"Processed {filename}")

# Main function to parse arguments
def main():
    parser = argparse.ArgumentParser(description="Process JSON files in a folder.")
    parser.add_argument("folder", help="Path to the folder containing JSON files")
    args = parser.parse_args()

    # Call the function to process the files
    process_json_files(args.folder)

if __name__ == "__main__":
    main()