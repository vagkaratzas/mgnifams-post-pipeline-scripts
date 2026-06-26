import gzip
import sqlite3
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "update_seed_msa_blobs.py"


def write_gz_text(path, text):
    with gzip.open(path, "wt") as handle:
        handle.write(text)


def create_db(path, family_ids):
    conn = sqlite3.connect(path)
    conn.execute(
        """
        CREATE TABLE mgnifam (
            id INTEGER PRIMARY KEY,
            seed_msa_blob BLOB DEFAULT NULL
        )
        """
    )
    conn.executemany("INSERT INTO mgnifam (id) VALUES (?)", [(id_,) for id_ in family_ids])
    conn.commit()
    conn.close()


def blob_for(db_path, family_id):
    conn = sqlite3.connect(db_path)
    value = conn.execute(
        "SELECT seed_msa_blob FROM mgnifam WHERE id = ?", (family_id,)
    ).fetchone()[0]
    conn.close()
    return value


def run_script(db_path, seed_msa_dir, *extra_args):
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--db",
            str(db_path),
            "--seed-msa-dir",
            str(seed_msa_dir),
            *extra_args,
        ],
        check=False,
        text=True,
        capture_output=True,
    )


def test_updates_seed_msa_blob_with_decompressed_fasta_bytes(tmp_path):
    db_path = tmp_path / "mgnifams.sqlite3"
    seed_msa_dir = tmp_path / "seed_msa_trimmed"
    seed_msa_dir.mkdir()
    create_db(db_path, [1, 2])

    family_1 = ">seq1/2-4\nBCD\n"
    family_2 = ">seq2/10-12\nXYZ\n"
    write_gz_text(seed_msa_dir / "1.fas.gz", family_1)
    write_gz_text(seed_msa_dir / "2.fas.gz", family_2)

    result = run_script(db_path, seed_msa_dir)

    assert result.returncode == 0, result.stderr
    assert blob_for(db_path, 1) == family_1.encode()
    assert blob_for(db_path, 2) == family_2.encode()
    assert "Updated 2 seed_msa_blob value(s)" in result.stdout


def test_fails_by_default_when_db_row_has_no_matching_trimmed_msa(tmp_path):
    db_path = tmp_path / "mgnifams.sqlite3"
    seed_msa_dir = tmp_path / "seed_msa_trimmed"
    seed_msa_dir.mkdir()
    create_db(db_path, [1, 2])
    write_gz_text(seed_msa_dir / "1.fas.gz", ">seq1\nAAA\n")

    result = run_script(db_path, seed_msa_dir)

    assert result.returncode != 0
    assert "Missing seed MSA file(s) for mgnifam id(s): 2" in result.stderr
    assert blob_for(db_path, 1) is None


def test_allow_missing_updates_available_files_only(tmp_path):
    db_path = tmp_path / "mgnifams.sqlite3"
    seed_msa_dir = tmp_path / "seed_msa_trimmed"
    seed_msa_dir.mkdir()
    create_db(db_path, [1, 2])
    write_gz_text(seed_msa_dir / "1.fas.gz", ">seq1\nAAA\n")

    result = run_script(db_path, seed_msa_dir, "--allow-missing")

    assert result.returncode == 0, result.stderr
    assert blob_for(db_path, 1) == b">seq1\nAAA\n"
    assert blob_for(db_path, 2) is None


def test_fast_pragmas_option_still_updates_seed_msa_blobs(tmp_path):
    db_path = tmp_path / "mgnifams.sqlite3"
    seed_msa_dir = tmp_path / "seed_msa_trimmed"
    seed_msa_dir.mkdir()
    create_db(db_path, [1])
    write_gz_text(seed_msa_dir / "1.fas.gz", ">seq1\nAAA\n")

    result = run_script(db_path, seed_msa_dir, "--fast-pragmas")

    assert result.returncode == 0, result.stderr
    assert blob_for(db_path, 1) == b">seq1\nAAA\n"


def test_fails_when_trimmed_msa_file_has_no_db_row(tmp_path):
    db_path = tmp_path / "mgnifams.sqlite3"
    seed_msa_dir = tmp_path / "seed_msa_trimmed"
    seed_msa_dir.mkdir()
    create_db(db_path, [1])
    write_gz_text(seed_msa_dir / "1.fas.gz", ">seq1\nAAA\n")
    write_gz_text(seed_msa_dir / "2.fas.gz", ">seq2\nBBB\n")

    result = run_script(db_path, seed_msa_dir)

    assert result.returncode != 0
    assert "Seed MSA file(s) without mgnifam row: 2.fas.gz" in result.stderr
    assert blob_for(db_path, 1) is None
