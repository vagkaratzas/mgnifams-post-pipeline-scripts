#!/usr/bin/env python3
"""Apply `bin/rebuild_full_msas.py`'s CSV to the `mgnifam` table of the website SQLite DB.

Reads the `mgnifam_updates.csv` that the rebuild wrote and updates the columns the
rebuild invalidated. Rows with `status` other than `rebuilt` are skipped: a family whose
alignment could not be regenerated must keep its released values rather than be blanked.

Previews by default. `--apply` is what writes, and it copies the database to `<db>.bak`
first unless that file already exists.

Columns
-------
`full_size`, `protein_rep`, `rep_region`, `rep_length`, `rep_sequence` -- everything the
re-recruitment can move. Not `consensus` or `converged`: both come from the HMM and the
seed loop, and the rebuild re-runs neither. Restrict further with `--columns`.

Usage
-----
    python bin/update_mgnifam_db_full_msa_bug.py --csv full_msa_fixed/mgnifam_updates.csv --sqlite mgnifams.sqlite3
    python bin/update_mgnifam_db_full_msa_bug.py --csv full_msa_fixed/mgnifam_updates.csv --sqlite mgnifams.sqlite3 --apply
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sqlite3
import sys
from pathlib import Path

TABLE = "mgnifam"
DB_COLUMNS = ("full_size", "protein_rep", "rep_region", "rep_length", "rep_sequence")
INTEGER_COLUMNS = {"full_size", "protein_rep", "rep_length"}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--csv", type=Path, help="mgnifam_updates.csv from bin/rebuild_full_msas.py")
    parser.add_argument("--sqlite", type=Path, help="website SQLite database")
    parser.add_argument(
        "--columns", nargs="+", default=list(DB_COLUMNS), choices=DB_COLUMNS, help="columns to write"
    )
    parser.add_argument("--apply", action="store_true", help="write; without it this only previews")
    parser.add_argument("--no_backup", action="store_true", help="skip the <db>.bak copy")
    parser.add_argument("--self_test", action="store_true", help="run the built-in check and exit")
    return parser.parse_args(argv)


def load_updates(path: Path, columns: list[str]) -> list[dict[str, object]]:
    """Read the rebuild CSV into one dict per family, coercing the integer columns.

    A missing column is an error rather than a silent skip: it means the CSV came from a
    rebuild that did not produce this value, and half-applying the update would leave the
    row internally inconsistent -- a `full_size` that no longer matches its `rep_*`.
    """
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        available = reader.fieldnames or []
        missing = [column for column in ["id", "status", *columns] if column not in available]
        if missing:
            raise SystemExit(f"{path}: missing column(s) {', '.join(missing)}")

        updates = []
        for line_number, row in enumerate(reader, 2):
            if row["status"] != "rebuilt":
                continue
            if not row["id"].isdigit():
                raise SystemExit(f"{path}:{line_number}: id {row['id']!r} is not a mgnifam id")
            values: dict[str, object] = {"id": int(row["id"])}
            for column in columns:
                value = row[column]
                if value == "":
                    raise SystemExit(f"{path}:{line_number}: {column} is empty for family {row['id']}")
                values[column] = int(value) if column in INTEGER_COLUMNS else value
            updates.append(values)
    return updates


def diff_against_db(
    connection: sqlite3.Connection, updates: list[dict[str, object]], columns: list[str]
) -> tuple[list[dict[str, object]], list[int], dict[str, int]]:
    """Split the updates into (present, missing ids) and count what each column changes."""
    selection = ", ".join(columns)
    current = {
        row[0]: dict(zip(columns, row[1:], strict=True))
        for row in connection.execute(f"SELECT id, {selection} FROM {TABLE}")
    }
    present, missing, changed = [], [], dict.fromkeys(columns, 0)
    for update in updates:
        existing = current.get(update["id"])
        if existing is None:
            missing.append(update["id"])
            continue
        present.append(update)
        for column in columns:
            if existing[column] != update[column]:
                changed[column] += 1
    return present, missing, changed


def main(argv: list[str] | None = None) -> int:
    options = parse_args(argv)
    if options.self_test:
        return self_test()
    for required in ("csv", "sqlite"):
        if getattr(options, required) is None:
            raise SystemExit(f"--{required} is required")

    updates = load_updates(options.csv, options.columns)
    with sqlite3.connect(options.sqlite) as connection:
        present, missing, changed = diff_against_db(connection, updates, options.columns)

    print(f"{len(updates)} rebuilt families in {options.csv}")
    if missing:
        print(f"  {len(missing)} not present in {TABLE}, skipping: {missing[:10]}")
    for column in options.columns:
        print(f"  {column}: {changed[column]} of {len(present)} rows change")

    if not options.apply:
        print("\npreview only; re-run with --apply to write")
        return 0
    if not present:
        print("\nnothing to write")
        return 0

    if not options.no_backup:
        backup = options.sqlite.with_suffix(options.sqlite.suffix + ".bak")
        if backup.exists():
            print(f"\nbackup {backup} already exists, keeping it")
        else:
            shutil.copy2(options.sqlite, backup)
            print(f"\nbacked up {options.sqlite} to {backup}")

    assignments = ", ".join(f"{column} = :{column}" for column in options.columns)
    with sqlite3.connect(options.sqlite) as connection:
        connection.executemany(f"UPDATE {TABLE} SET {assignments} WHERE id = :id", present)
    print(f"updated {len(present)} rows in {options.sqlite}")
    return 0


def self_test() -> int:
    """Apply a CSV to a throwaway database and assert what landed."""
    import tempfile

    with tempfile.TemporaryDirectory() as workdir:
        work = Path(workdir)
        database = work / "test.sqlite3"
        with sqlite3.connect(database) as connection:
            connection.execute(
                f"CREATE TABLE {TABLE} (id INTEGER PRIMARY KEY, full_size INTEGER, "
                "protein_rep INTEGER, rep_region TEXT, rep_length INTEGER, rep_sequence TEXT, "
                "consensus TEXT)"
            )
            connection.executemany(
                f"INSERT INTO {TABLE} VALUES (?, ?, ?, ?, ?, ?, ?)",
                [(1, 99, 111, "1-9", 9, "MKVLAAGIV", "keepme"), (2, 5, 222, "2-6", 5, "MKVLA", "keepme")],
            )

        updates = work / "mgnifam_updates.csv"
        with updates.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["id", "status", *DB_COLUMNS, "note"])
            writer.writeheader()
            writer.writerow({"id": "1", "status": "rebuilt", "full_size": "21", "protein_rep": "333",
                             "rep_region": "10-18", "rep_length": "9", "rep_sequence": "AAAAAAAAA", "note": ""})
            # A skipped family must keep everything it had.
            writer.writerow({"id": "2", "status": "skipped", "full_size": "", "protein_rep": "",
                             "rep_region": "", "rep_length": "", "rep_sequence": "", "note": "no hits"})
            # An id the database does not carry must be reported, not raise.
            writer.writerow({"id": "77", "status": "rebuilt", "full_size": "3", "protein_rep": "444",
                             "rep_region": "1-3", "rep_length": "3", "rep_sequence": "MKV", "note": ""})

        loaded = load_updates(updates, list(DB_COLUMNS))
        assert [row["id"] for row in loaded] == [1, 77], loaded
        assert loaded[0]["full_size"] == 21 and isinstance(loaded[0]["full_size"], int)

        assert main(["--csv", str(updates), "--sqlite", str(database)]) == 0
        with sqlite3.connect(database) as connection:
            assert connection.execute(f"SELECT full_size FROM {TABLE} WHERE id = 1").fetchone()[0] == 99, (
                "a preview must not write"
            )

        assert main(["--csv", str(updates), "--sqlite", str(database), "--apply"]) == 0
        with sqlite3.connect(database) as connection:
            row = connection.execute(
                f"SELECT full_size, protein_rep, rep_region, rep_length, rep_sequence, consensus "
                f"FROM {TABLE} WHERE id = 1"
            ).fetchone()
            assert row == (21, 333, "10-18", 9, "AAAAAAAAA", "keepme"), row
            untouched = connection.execute(
                f"SELECT full_size, rep_sequence FROM {TABLE} WHERE id = 2"
            ).fetchone()
            assert untouched == (5, "MKVLA"), untouched
        assert (work / "test.sqlite3.bak").exists(), "a write must leave a backup"

    print("self_test: ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
