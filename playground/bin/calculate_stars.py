#!/usr/bin/env python3
import sqlite3
import argparse
import pandas as pd

def get_mgnifam_data(db_path):
    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Execute the query to select all IDs, plddt, family_size, and rep_length from the mgnifam table
    cursor.execute("SELECT id, plddt, family_size, rep_length FROM mgnifam")
    
    # Fetch all results
    mgnifam_data = cursor.fetchall()

    # Execute the query pfams table
    # tempalte_hmm_range format is x-y (z), function to parse that into coverage is (y-x+1) / z
    cursor.execute("""
        SELECT 
            mgnifam_id,
            MAX(
            (CAST(SUBSTR(template_hmm_range, INSTR(template_hmm_range, '-') + 1, 
                        INSTR(template_hmm_range, '(') - INSTR(template_hmm_range, '-') - 1) AS REAL) -
            CAST(SUBSTR(template_hmm_range, 1, INSTR(template_hmm_range, '-') - 1) AS REAL) + 1
            ) /
            CAST(SUBSTR(template_hmm_range, INSTR(template_hmm_range, '(') + 1, 
                        INSTR(template_hmm_range, ')') - INSTR(template_hmm_range, '(') - 1) AS REAL)
        ) AS max_pfam_coverage
        FROM 
            mgnifam_pfams
        WHERE 
            template_hmm_range LIKE '%-% (%' -- Ensure correct format
        GROUP BY 
            mgnifam_id
    """)
    mgnifam_pfams_res = cursor.fetchall()

    # Execute the query for folds table
    # Disregarding ESMAtlas 'MGYP' hits. We are MGnify !
    cursor.execute("""
        SELECT f.mgnifam_id, f.aligned_length, m.rep_length
        FROM mgnifam_folds f
        JOIN mgnifam m ON f.mgnifam_id = m.id
        WHERE f.target_structure NOT LIKE 'MGYP%' 
        GROUP BY f.mgnifam_id 
        HAVING MAX(f.aligned_length)
    """)
    mgnifam_folds_res = cursor.fetchall()

    # Close the connection
    conn.close()

    return mgnifam_data, mgnifam_pfams_res, mgnifam_folds_res

def calculate_quality_stars(plddt, family_size, rep_length):
    quality_stars = 0

    # Add stars based on plddt
    if (plddt is not None) and (plddt >= 70):
        quality_stars += 2
    elif (plddt is not None) and (plddt >= 50):
        quality_stars += 1

    # Add stars based on family_size
    if (family_size is not None) and (family_size >= 100):
        quality_stars += 1

    # Add stars based on rep_length
    if (rep_length is not None) and (rep_length >= 100):
        quality_stars += 1

    return quality_stars

def calculate_novelty_stars(mgnifam_id, mgnifam_pfams_res, mgnifam_folds_res):
    print(f"{mgnifam_id}-", end="")
    novelty_stars = 0

    # Calculate and add pfam stars:
    # Check if mgnifam_id is not in mgnifam_pfams_res
    pfam_coverage = next((row[1] for row in mgnifam_pfams_res if row[0] == mgnifam_id), None)
    
    if pfam_coverage is None:
        # Add 2 stars if no result was found for this mgnifam_id in mgnifam_pfams_res
        novelty_stars += 2
    else:
        # Add 1 star if max_pfam_coverage < 0.5
        if pfam_coverage < 0.5:
            novelty_stars += 1
    
    # Calculate and add fold stars:
    # Find the result for the current mgnifam_id
    result = next((row for row in mgnifam_folds_res if row[0] == mgnifam_id), None)

    if result is None:
        # Add 2 stars if no result was found for this mgnifam_id
        novelty_stars += 2
    else:
        aligned_length = result[1]
        rep_length = result[2]

        # Add 1 star if aligned_length / rep_length is less than 0.5
        if rep_length and (aligned_length / rep_length) < 0.5:
            novelty_stars += 1

    return novelty_stars

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Query mgnifam IDs and calculate quality and novelty stars from an SQLite database.")
    parser.add_argument("db_path", help="Path to the SQLite database")

    # Parse command-line arguments
    args = parser.parse_args()

    print("Executing queries...")
    # Get the mgnifam data and the sets of mgnifam_ids from the other tables
    mgnifam_data, mgnifam_pfams_res, mgnifam_folds_res = get_mgnifam_data(args.db_path)

    print("Calculating stars...")
    # Prepare the data for the DataFrame
    ids = [row[0] for row in mgnifam_data]
    plddts = [row[1] for row in mgnifam_data]
    family_sizes = [row[2] for row in mgnifam_data]
    rep_lengths = [row[3] for row in mgnifam_data]
    print("Quality stars...")
    quality_stars = [calculate_quality_stars(plddt, family_size, rep_length) for plddt, family_size, rep_length in zip(plddts, family_sizes, rep_lengths)]
    print("Novely stars...")
    novelty_stars = [calculate_novelty_stars(mgnifam_id, mgnifam_pfams_res, mgnifam_folds_res) for mgnifam_id in ids]

    # Create the DataFrame
    df = pd.DataFrame({
        'MGnify Family ID': ids,
        'Quality Stars': quality_stars,
        'Novelty Stars': novelty_stars
    })

    # Print the DataFrame
    print("\n")
    print(df)
    df.to_csv('output/stars.csv', index=False)

if __name__ == "__main__":
    main()
