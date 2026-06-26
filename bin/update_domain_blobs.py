#!/usr/bin/env python3

import sqlite3
import sys
import os

def test_connection(conn): 
    try:
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM mgnifam")
        row_count = cursor.fetchone()[0]
        print("Connection successful! Number of rows in the table: ", row_count)
    except sqlite3.Error as e:
        print("Connection failed:", e)

def construct_file_path(base_dir, file_column):
    directory = "post-processing/domain_results"
    return os.path.join(base_dir, directory, file_column)

def read_file(file_path):
    try:
        with open(file_path, 'rb') as file:
            return file.read()
    except (OSError, IOError) as e:
        print(f"Error reading file {file_path}: {e}")
        return None

def update_blob_column(db_path, blob_data, row_id):
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        query = "UPDATE mgnifam SET domain_architecture_blob = ? WHERE id = ?"
        cursor.execute(query, (sqlite3.Binary(blob_data), row_id))
        conn.commit()
        conn.close()
        print(f"Updated domain_architecture_blob for row {row_id}")
    except sqlite3.Error as e:
        print(f"Failed to update domain_architecture_blob for row {row_id}: {e}")

def process_row(db_path, base_dir, row):
    row_id = row[0]
    file_column = row[1]

    if file_column is not None: 
        file_path = construct_file_path(base_dir, file_column)
        blob_data = read_file(file_path)
        
        if blob_data:
            update_blob_column(db_path, blob_data, row_id)
        else:
            print(f"Skipping update for row {row_id} due to read error")

def import_files(db_path, base_dir):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT id, domain_architecture_file FROM mgnifam")
    rows = cursor.fetchall()
    conn.close()

    for row in rows:
        process_row(db_path, base_dir, row)

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 append_domain_blobs_sqlite.py <db.sqlite3> <output_dir>")
        sys.exit(1)

    db_path = sys.argv[1]
    base_dir = sys.argv[2]

    # Import domain files
    import_files(db_path, base_dir)
    
if __name__ == "__main__":
    main()
