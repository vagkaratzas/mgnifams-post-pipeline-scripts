#!/usr/bin/env python3

import argparse
import gzip
import sqlite3
import sys
from pathlib import Path


def seed_msa_id(path):
    name = path.name
    if not name.endswith(".fas.gz"):
        raise ValueError(f"Seed MSA file must end with .fas.gz: {path}")

    family_id = name[: -len(".fas.gz")]
    if not family_id.isdigit():
        raise ValueError(f"Seed MSA filename must start with a numeric family id: {name}")
    return int(family_id)


def discover_seed_msa_files(seed_msa_dir):
    seed_msa_dir = Path(seed_msa_dir)
    if not seed_msa_dir.is_dir():
        raise NotADirectoryError(f"Seed MSA directory does not exist: {seed_msa_dir}")

    files_by_id = {}
    for path in sorted(seed_msa_dir.glob("*.fas.gz")):
        family_id = seed_msa_id(path)
        if family_id in files_by_id:
            raise ValueError(
                f"Duplicate seed MSA file for mgnifam id {family_id}: "
                f"{files_by_id[family_id]} and {path}"
            )
        files_by_id[family_id] = path

    if not files_by_id:
        raise ValueError(f"No .fas.gz seed MSA files found in {seed_msa_dir}")
    return files_by_id


def fetch_mgnifam_ids(conn):
    return {row[0] for row in conn.execute("SELECT id FROM mgnifam")}


def apply_fast_pragmas(conn):
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA locking_mode=EXCLUSIVE")
    conn.execute("PRAGMA cache_size=-2000000")


def restore_default_pragmas(conn):
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA locking_mode=NORMAL")
    conn.execute("PRAGMA journal_mode=DELETE")


def validate_id_matches(db_ids, files_by_id, allow_missing):
    file_ids = set(files_by_id)
    missing_file_ids = sorted(db_ids - file_ids)
    extra_file_ids = sorted(file_ids - db_ids)

    errors = []
    if missing_file_ids and not allow_missing:
        errors.append(
            "Missing seed MSA file(s) for mgnifam id(s): "
            + ", ".join(str(id_) for id_ in missing_file_ids)
        )
    if extra_file_ids:
        errors.append(
            "Seed MSA file(s) without mgnifam row: "
            + ", ".join(files_by_id[id_].name for id_ in extra_file_ids)
        )

    if errors:
        raise ValueError("\n".join(errors))


def read_decompressed_seed_msa(path):
    with gzip.open(path, "rb") as handle:
        return handle.read()


def update_seed_msa_blobs(db_path, seed_msa_dir, allow_missing=False, fast_pragmas=False):
    files_by_id = discover_seed_msa_files(seed_msa_dir)

    conn = sqlite3.connect(db_path)
    try:
        if fast_pragmas:
            apply_fast_pragmas(conn)

        db_ids = fetch_mgnifam_ids(conn)
        validate_id_matches(db_ids, files_by_id, allow_missing)

        update_ids = sorted(set(files_by_id) & db_ids)
        with conn:
            conn.executemany(
                "UPDATE mgnifam SET seed_msa_blob = ? WHERE id = ?",
                [
                    (sqlite3.Binary(read_decompressed_seed_msa(files_by_id[id_])), id_)
                    for id_ in update_ids
                ],
            )
    finally:
        if fast_pragmas:
            restore_default_pragmas(conn)
        conn.close()

    return len(update_ids)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Update mgnifam.seed_msa_blob from trimmed gzipped seed MSA files."
    )
    parser.add_argument("--db", required=True, help="SQLite database to update")
    parser.add_argument(
        "--seed-msa-dir",
        required=True,
        help="Directory containing trimmed seed MSA files named <id>.fas.gz",
    )
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="Update available files even when some mgnifam rows have no matching MSA file",
    )
    parser.add_argument(
        "--fast-pragmas",
        action="store_true",
        help="Apply faster, less durable SQLite PRAGMAs during this bulk update",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    try:
        updated = update_seed_msa_blobs(
            args.db, args.seed_msa_dir, args.allow_missing, args.fast_pragmas
        )
    except Exception as error:
        print(str(error), file=sys.stderr)
        return 1

    print(f"Updated {updated} seed_msa_blob value(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
