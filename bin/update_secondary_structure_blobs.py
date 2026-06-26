#!/usr/bin/env python3

import os
import argparse
import sqlite3

def update_mgnifam_table(folder_path, db_path):
    """
    Updates the `mgnifam` table with file names and blobs for JSON predictions.
    
    Args:
        folder_path (str): Path to the folder containing JSON files.
        db_path (str): Path to the SQLite database.
    """
    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Disable synchronous for faster writes (use with caution)
    cursor.execute("PRAGMA synchronous = OFF")
    # Begin a transaction
    cursor.execute("BEGIN TRANSACTION")

    # Iterate over all files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith('.json'):
            file_path = os.path.join(folder_path, filename)
            
            # Extract the ID from the filename (basename without extension)
            mgnifam_id = os.path.splitext(filename)[0]
            print(f"Processing {mgnifam_id}")

            # Check if the columns are already updated for this ID
            cursor.execute(
                """
                SELECT pred_secondary_structure_file, pred_secondary_structure_blob
                FROM mgnifam
                WHERE id = ?
                """, 
                (mgnifam_id,)
            )
            result = cursor.fetchone()
            
            if result and result[0] is not None and result[1] is not None:
                # Skip if already updated
                print(f"Skipping {mgnifam_id}, already updated.")
                continue

            # Read the JSON file as a blob
            with open(file_path, 'rb') as file:
                file_blob = file.read()

            # Update the database for the matching ID
            cursor.execute(
                """
                UPDATE mgnifam
                SET pred_secondary_structure_file = ?,
                    pred_secondary_structure_blob = ?
                WHERE id = ?
                """,
                (filename, file_blob, mgnifam_id)
            )

    # Commit changes and close the connection
    conn.commit()
    conn.close()

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Update mgnifam table with JSON predictions.")
    parser.add_argument("folder_path", type=str, help="Path to the folder containing JSON files.")
    parser.add_argument("db_path", type=str, help="Path to the SQLite database.")

    args = parser.parse_args()
    # Example: ./bin/update_secondary_structure_blobs.py /home/vangelis/Desktop/Projects/mgnifams-post-pipeline-scripts/output/s4pred/mgnifams-local /home/vangelis/Desktop/Projects/mgnifams-site/mgnifams_site/dbs/mgnifams.sqlite3
    
    # Run the update function
    update_mgnifam_table(args.folder_path, args.db_path)
