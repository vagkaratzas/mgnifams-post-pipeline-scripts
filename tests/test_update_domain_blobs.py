import importlib.util
import json
import sqlite3
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "update_domain_blobs.py"


def load_module():
    spec = importlib.util.spec_from_file_location("update_domain_blobs", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mod = load_module()


def make_sqlite(path, ids):
    connection = sqlite3.connect(path)
    connection.execute("CREATE TABLE mgnifam (id INTEGER PRIMARY KEY, domain_blob BLOB)")
    connection.executemany("INSERT INTO mgnifam (id) VALUES (?)", [(i,) for i in ids])
    connection.commit()
    connection.close()


def write_json(json_dir, family_id, payload=None):
    json_dir.mkdir(parents=True, exist_ok=True)
    if payload is None:
        payload = {"architecture_containers": [{"architecture_text": str(family_id),
                                                "domains": []}]}
    (json_dir / f"{family_id}.json").write_text(json.dumps(payload))


def read_blob(db, family_id):
    connection = sqlite3.connect(db)
    blob = connection.execute("SELECT domain_blob FROM mgnifam WHERE id = ?",
                              (family_id,)).fetchone()[0]
    connection.close()
    return blob


def test_matching_json_is_written_into_the_blob_column(tmp_path):
    db = tmp_path / "db.sqlite3"
    json_dir = tmp_path / "json"
    make_sqlite(db, [1, 2])
    write_json(json_dir, 1)
    write_json(json_dir, 2)

    mod.update_blobs(db, json_dir)

    assert json.loads(read_blob(db, 1))["architecture_containers"][0]["architecture_text"] == "1"
    assert json.loads(read_blob(db, 2))["architecture_containers"][0]["architecture_text"] == "2"


def test_rows_without_a_json_are_left_alone_and_reported(tmp_path):
    db = tmp_path / "db.sqlite3"
    json_dir = tmp_path / "json"
    make_sqlite(db, [1, 2])
    write_json(json_dir, 1)

    updated, missing, orphans = mod.update_blobs(db, json_dir)

    assert (updated, missing, orphans) == (1, ["2"], [])
    assert read_blob(db, 2) is None


def test_json_files_without_a_row_are_reported(tmp_path):
    db = tmp_path / "db.sqlite3"
    json_dir = tmp_path / "json"
    make_sqlite(db, [1])
    write_json(json_dir, 1)
    write_json(json_dir, 99)

    updated, missing, orphans = mod.update_blobs(db, json_dir)

    assert (updated, missing, orphans) == (1, [], ["99"])


def test_an_unknown_column_is_rejected_rather_than_interpolated(tmp_path):
    db = tmp_path / "db.sqlite3"
    json_dir = tmp_path / "json"
    make_sqlite(db, [1])
    write_json(json_dir, 1)

    with pytest.raises(ValueError):
        mod.update_blobs(db, json_dir, column="domain_blob = 1; DROP TABLE mgnifam")
