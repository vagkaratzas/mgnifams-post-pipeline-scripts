#!/usr/bin/env python3

import argparse
import os
import sqlite3

def create_subset_db_from_schema(schema_file, subset_db_path):
    with open(schema_file, 'r') as file:
        schema_sql = file.read()
    # Create a new subset database with the same schema
    with sqlite3.connect(subset_db_path) as conn:
        cursor = conn.cursor()
        cursor.executescript(schema_sql)
        conn.commit()

def insert_data(subset_conn, main_data, table_name):
    placeholders = ', '.join(['?'] * len(main_data[0]))
    query = f'INSERT INTO {table_name} VALUES ({placeholders})'
    subset_conn.executemany(query, main_data)

def subset_db(main_db_path, family_ids_file, subset_db_path, schema_file):
    # Read family IDs from the text file
    with open(family_ids_file, 'r') as f:
        family_ids = [line.strip() for line in f]

    # Create the new subset database by executing the schema SQL
    create_subset_db_from_schema(schema_file, subset_db_path)

    # Open connections to the main and subset databases
    with sqlite3.connect(main_db_path) as main_conn, sqlite3.connect(subset_db_path) as subset_conn:
        cursor = main_conn.cursor()

        # For each table, fetch and insert the corresponding rows for the family IDs
        tables = ['mgnifam', 'mgnifam_proteins', 'mgnifam_pfams', 'mgnifam_folds']
        for table in tables:
            if table == 'mgnifam':
                query = f"SELECT * FROM {table} WHERE id IN ({','.join(['?']*len(family_ids))})"
            else:
                query = f"SELECT * FROM {table} WHERE mgnifam_id IN ({','.join(['?']*len(family_ids))})"

            cursor.execute(query, family_ids)
            rows = cursor.fetchall()

            if rows:
                insert_data(subset_conn, rows, table)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subset an SQLite database by family IDs.")
    parser.add_argument("main_db", help="Path to the main SQLite database")
    parser.add_argument("family_ids_file", help="Path to the file containing family IDs")
    parser.add_argument("subset_db", help="Path to the new subset SQLite database")
    parser.add_argument("schema_file", help="Path to the schema file containing SQL")
    # Local: python ./bin/subset_mgnifams_sqlite.py assets/test_data/mgnifams.sqlite3 assets/test_data/family_ids.txt output/mgnifams_subset.sqlite3 assets/test_data/schema.sqlite

    args = parser.parse_args()
    if os.path.exists(args.subset_db):
        os.remove(args.subset_db)
    subset_db(args.main_db, args.family_ids_file, args.subset_db, args.schema_file)
