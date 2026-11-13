#!/usr/bin/env python3
import sqlite3
import pandas as pd
import argparse

def ensure_columns_exist(cursor):
    # Check if quality_rank and novelty_rank columns exist in the mgnifam table
    cursor.execute("PRAGMA table_info(mgnifam)")
    columns = {row[1] for row in cursor.fetchall()}

    # If quality_rank column doesn't exist, add it
    if 'quality_rank' not in columns:
        cursor.execute("ALTER TABLE mgnifam ADD COLUMN quality_rank INTEGER")

    # If novelty_rank column doesn't exist, add it
    if 'novelty_rank' not in columns:
        cursor.execute("ALTER TABLE mgnifam ADD COLUMN novelty_rank INTEGER")

def update_mgnifam_ranks(db_path, csv_path):
    # Load the CSV file into a DataFrame
    df = pd.read_csv(csv_path)

    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Ensure the columns exist before updating
    ensure_columns_exist(cursor)

    # Iterate over the DataFrame rows and update the database
    for index, row in df.iterrows():
        mgnifam_id = int(row['MGnify Family ID'])
        print(f"{mgnifam_id}-", end="")
        quality_stars = int(row['Quality Stars'])
        novelty_stars = int(row['Novelty Stars'])

        # Update the mgnifam table
        cursor.execute("""
            UPDATE mgnifam
            SET quality_rank = ?, novelty_rank = ?
            WHERE id = ?
        """, (quality_stars, novelty_stars, mgnifam_id))

    # Commit the changes and close the connection
    conn.commit()
    conn.close()

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Update mgnifam table with quality and novelty ranks from a CSV file.")
    parser.add_argument("db_path", help="Path to the SQLite database")
    parser.add_argument("csv_path", help="Path to the stars CSV file")

    # Parse command-line arguments
    args = parser.parse_args()

    # Update mgnifam ranks
    update_mgnifam_ranks(args.db_path, args.csv_path)

    # Test by querying the first line of the mgnifam table
    conn = sqlite3.connect(args.db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT id, novelty_rank, quality_rank FROM mgnifam LIMIT 1")
    result = cursor.fetchone()
    conn.close()

    # Print the result to the console
    print("First line of mgnifam table:", result)

if __name__ == "__main__":
    main()
