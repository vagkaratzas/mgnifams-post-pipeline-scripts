#!/usr/bin/env python3
import os
import argparse
import shutil

def check_versions_file(subfolder_path):
    versions_file_path = os.path.join(subfolder_path, 'versions.yml')
    if os.path.isfile(versions_file_path):
        try:
            if "EBIMETAGENOMICS_MGNIFAMS:MGNIFAMS:GENERATE_NONREDUNDANT_FAMILIES:GENERATE_FAMILIES_PARALLEL:REFINE_FAMILIES_PARALLEL" in open(versions_file_path).read():
                return True
        except Exception as e:
            print(f"An error occurred while reading {versions_file_path}: {e}")
    return False

def check_ids_in_rf_folder(subfolder_path, ids):
    rf_folder_path = os.path.join(subfolder_path, 'rf')
    if os.path.isdir(rf_folder_path):
        rf_files = {os.path.splitext(f.name)[0] for f in os.scandir(rf_folder_path) if f.is_file()}
        return [id_ for id_ in ids if id_ in rf_files]
    return []

def copy_files(subfolder_path, hit_id, output_dir):
    file_paths = {
        "rf": f"rf/{hit_id}.txt",
        "hmm": f"hmm/{hit_id}.hmm",
        "msa_sto": f"msa_sto/{hit_id}.sto",
        "seed_msa_sto": f"seed_msa_sto/{hit_id}.sto",
        "domtblout": f"domtblout/{hit_id}.domtblout"
    }
    for folder, relative_path in file_paths.items():
        source_path = os.path.join(subfolder_path, relative_path)
        if os.path.isfile(source_path):
            dest_dir = os.path.join(output_dir, hit_id, folder)
            os.makedirs(dest_dir, exist_ok=True)
            dest_path = os.path.join(dest_dir, os.path.basename(source_path))
            shutil.copyfile(source_path, dest_path)

    # Process refined_families file
    refined_families_dir = os.path.join(subfolder_path, "refined_families")
    if os.path.isdir(refined_families_dir):
        refined_files = [f for f in os.listdir(refined_families_dir) if os.path.isfile(os.path.join(refined_families_dir, f))]
        if refined_files:
            refined_file_path = os.path.join(refined_families_dir, refined_files[0])
            refined_basename = os.path.splitext(refined_files[0])[0]
            
            refined_output_file = os.path.join(output_dir, hit_id, "refined_families.tsv")
            os.makedirs(os.path.dirname(refined_output_file), exist_ok=True)
            
            with open(refined_file_path, 'r') as infile, open(refined_output_file, 'w') as outfile:
                for line in infile:
                    parts = line.strip().split('\t')
                    if len(parts) == 2:
                        concatenated_name = f"{refined_basename}_{parts[0]}"
                        if concatenated_name == hit_id:
                            outfile.write(concatenated_name + "\t" + parts[1] + "\n")
    
    # Process family_metadata file
    family_metadata_dir = os.path.join(subfolder_path, "family_metadata")
    if os.path.isdir(family_metadata_dir):
        metadata_files = [f for f in os.listdir(family_metadata_dir) if os.path.isfile(os.path.join(family_metadata_dir, f))]
        if metadata_files:
            metadata_file_path = os.path.join(family_metadata_dir, metadata_files[0])
            metadata_basename = os.path.splitext(metadata_files[0])[0]
            
            metadata_output_file = os.path.join(output_dir, hit_id, "family_metadata.csv")
            os.makedirs(os.path.dirname(metadata_output_file), exist_ok=True)
            
            with open(metadata_file_path, 'r') as infile, open(metadata_output_file, 'w') as outfile:
                for line in infile:
                    parts = line.strip().split(',')
                    if len(parts) >= 4:
                        concatenated_name = f"{metadata_basename}_{parts[0]}"
                        if concatenated_name == hit_id:
                            outfile.write(concatenated_name + "," + parts[1] + "," + parts[2] + "," + parts[3] + "\n")

def list_valid_subfolders(folder_path, ids, output_dir):
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
                    hit_ids = check_ids_in_rf_folder(second_level_subfolder, ids)
                    if hit_ids:
                        valid_subfolders.append(second_level_subfolder)
                        for hit_id in hit_ids:
                            copy_files(second_level_subfolder, hit_id, output_dir)
    except Exception as e:
        print(f"An error occurred: {e}")
    return valid_subfolders

def main():
    parser = argparse.ArgumentParser(description="List valid second level subfolder paths of a given folder based on specific criteria.")
    parser.add_argument("folder_path", type=str, help="Path to the folder")
    parser.add_argument("ids_file", type=str, help="Path to the file containing line-separated list of IDs")
    parser.add_argument("output_dir", type=str, help="Path to the output directory where files will be copied")
    # Example call: ./bin/extract_tm_families.py /home/vangelis/Desktop/Projects/mgnifams/work /home/vangelis/Desktop/Projects/mgnifams/output/redundancy/tm/tm_ids.txt output
    
    args = parser.parse_args()
    folder_path = args.folder_path
    ids_file = args.ids_file
    output_dir = args.output_dir

    if not os.path.isfile(ids_file):
        print("The provided IDs file is not valid.")
        return

    with open(ids_file, 'r') as file:
        ids = [line.strip() for line in file.readlines()]

    if os.path.isdir(folder_path):
        valid_subfolders = list_valid_subfolders(folder_path, ids, output_dir)
        print("Valid second level subfolders containing the specific versions.yml and matching IDs in rf folder:")
        for subfolder in valid_subfolders:
            print(subfolder)
    else:
        print("The provided path is not a valid directory.")

if __name__ == "__main__":
    main()
