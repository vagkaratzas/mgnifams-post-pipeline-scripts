#!/usr/bin/env python3

"""Load per-family domain architecture JSONs into the mgnifam SQLite blob column.

Input is the output directory of bin/parse_domain_architectures.py: one <family_id>.json per
family, including empty ones for families that got no annotated sequence.
"""

import argparse
import logging
import sqlite3
import time
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
log = logging.getLogger(__name__)


def _check_identifiers(connection, table, column):
    """Table and column are interpolated into SQL, so accept only names the schema really has."""
    tables = {row[0] for row in connection.execute(
        "SELECT name FROM sqlite_master WHERE type = 'table'")}
    if table not in tables:
        raise ValueError(f"No table named {table!r} in this database")

    columns = {row[1] for row in connection.execute(f"PRAGMA table_info({table})")}
    if column not in columns:
        raise ValueError(f"No column named {column!r} in table {table!r}")


def update_blobs(db, json_dir, table="mgnifam", column="domain_blob"):
    """Write <json_dir>/<id>.json into <column> for every row of <table>.

    Returns (rows updated, ids with no JSON file, JSON files with no row).
    """
    json_dir = Path(json_dir)
    connection = sqlite3.connect(db)

    try:
        _check_identifiers(connection, table, column)

        row_ids = [str(row[0]) for row in connection.execute(f"SELECT id FROM {table}")]
        updates = []
        missing = []

        for row_id in row_ids:
            path = json_dir / f"{row_id}.json"
            if path.is_file():
                updates.append((sqlite3.Binary(path.read_bytes()), row_id))
            else:
                missing.append(row_id)

        connection.executemany(f"UPDATE {table} SET {column} = ? WHERE id = ?", updates)
        connection.commit()
    finally:
        connection.close()

    on_disk = {path.stem for path in json_dir.glob("*.json")}
    orphans = sorted(on_disk - set(row_ids), key=lambda i: int(i) if i.isdigit() else 0)

    log.info(f"Updated {column} for {len(updates)} of {len(row_ids)} rows in {table}")
    if missing:
        log.warning(f"{len(missing)} rows had no JSON file, e.g. {missing[:5]}")
    if orphans:
        log.warning(f"{len(orphans)} JSON files had no matching row, e.g. {orphans[:5]}")

    return len(updates), missing, orphans


def main():
    parser = argparse.ArgumentParser(
        description="Load domain architecture JSONs into the mgnifam SQLite blob column.")
    parser.add_argument("--db", required=True, help="Path to the SQLite database")
    parser.add_argument("--json-dir", required=True,
                        help="Directory of <family_id>.json domain architectures")
    parser.add_argument("--table", default="mgnifam", help="Table to update")
    parser.add_argument("--column", default="domain_blob", help="Blob column to write")

    args = parser.parse_args()
    started = time.time()
    log.info("Starting update_domain_blobs")

    update_blobs(args.db, args.json_dir, args.table, args.column)

    log.info(f"update_domain_blobs complete in {time.time() - started:.1f}s")


if __name__ == "__main__":
    main()
