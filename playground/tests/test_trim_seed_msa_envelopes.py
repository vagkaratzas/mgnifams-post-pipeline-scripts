import gzip
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "trim_seed_msa_envelopes.py"


def write_gz_text(path, text):
    with gzip.open(path, "wt") as handle:
        handle.write(text)


def read_gz_text(path):
    with gzip.open(path, "rt") as handle:
        return handle.read()


def run_script(seed_msa_dir, rf_dir, output_dir):
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            str(seed_msa_dir),
            str(rf_dir),
            str(output_dir),
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr


def test_clips_terminal_rf_dots_and_recalculates_gap_aware_coordinates(tmp_path):
    seed_msa_dir = tmp_path / "seed_msa"
    rf_dir = tmp_path / "rf"
    output_dir = tmp_path / "trimmed"
    seed_msa_dir.mkdir()
    rf_dir.mkdir()

    write_gz_text(
        seed_msa_dir / "1.fas.gz",
        (
            ">seq1/10-16\n"
            "ABCD-EFG\n"
            ">seq2/20-25\n"
            "-AB-CDEF\n"
            ">full_length\n"
            "ABCDEFGH\n"
        ),
    )
    (rf_dir / "1.txt").write_text("..xxx.x.\n")

    run_script(seed_msa_dir, rf_dir, output_dir)

    assert read_gz_text(output_dir / "1.fas.gz") == (
        ">seq1/12-15\n"
        "CD-EF\n"
        ">seq2/21-24\n"
        "B-CDE\n"
        ">full_length/3-7\n"
        "CDEFG\n"
    )


def test_writes_unchanged_msa_when_rf_has_only_internal_dots(tmp_path):
    seed_msa_dir = tmp_path / "seed_msa"
    rf_dir = tmp_path / "rf"
    output_dir = tmp_path / "trimmed"
    seed_msa_dir.mkdir()
    rf_dir.mkdir()

    original_msa = ">seq1/5-7\nA-C\n>full_length\nABC\n"
    write_gz_text(seed_msa_dir / "2.fas.gz", original_msa)
    (rf_dir / "2.txt").write_text("x.x\n")

    run_script(seed_msa_dir, rf_dir, output_dir)

    assert sorted(path.name for path in output_dir.iterdir()) == ["2.fas.gz"]
    assert read_gz_text(output_dir / "2.fas.gz") == original_msa


def test_processes_all_seed_msas_by_matching_rf_basename(tmp_path):
    seed_msa_dir = tmp_path / "seed_msa"
    rf_dir = tmp_path / "rf"
    output_dir = tmp_path / "trimmed"
    seed_msa_dir.mkdir()
    rf_dir.mkdir()

    write_gz_text(seed_msa_dir / "family_a.fas.gz", ">a/1-4\nABCD\n")
    write_gz_text(seed_msa_dir / "family_b.fas.gz", ">b/10-13\nWXYZ\n")
    (rf_dir / "family_a.txt").write_text(".xx.\n")
    (rf_dir / "family_b.txt").write_text("xxxx\n")

    run_script(seed_msa_dir, rf_dir, output_dir)

    assert sorted(path.name for path in output_dir.iterdir()) == [
        "family_a.fas.gz",
        "family_b.fas.gz",
    ]
    assert read_gz_text(output_dir / "family_a.fas.gz") == ">a/2-3\nBC\n"
    assert read_gz_text(output_dir / "family_b.fas.gz") == ">b/10-13\nWXYZ\n"
