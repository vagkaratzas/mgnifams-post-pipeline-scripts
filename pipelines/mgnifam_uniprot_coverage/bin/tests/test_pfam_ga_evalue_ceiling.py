import gzip
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "pfam_ga_evalue_ceiling.py"


def row(target, query, acc, seq_e, dom_ie):
    """A --cut_ga domtbl row; only target, query, accession and E columns matter."""
    return " ".join(str(x) for x in [
        target, "-", 200, query, acc, 100, seq_e, 30.0, 0.0, 1, 1,
        dom_ie, dom_ie, 30.0, 0.0, 1, 100, 10, 60, 10, 60, 0.9, "desc",
    ])


def write_chunk(path, rows):
    with gzip.open(path, "wt") as fh:
        fh.write("# header\n" + "\n".join(rows) + "\n")


def run(argv):
    return subprocess.run([sys.executable, str(SCRIPT)] + argv,
                          capture_output=True, text=True)


def test_ceiling_rescales_each_chunk_by_its_own_size(tmp_path):
    # 2 full chunks of 10 records + a last chunk of 5: Z per chunk = 10, 10, 5
    d = tmp_path / "pfam"
    d.mkdir()
    write_chunk(d / "db_chunk_000001_pfam.domtbl.gz",
                [row("a", "ABC_tran", "PF00005.31", "1e-4", "1e-5")])
    write_chunk(d / "db_chunk_000002_pfam.domtbl.gz",
                [row("b", "ABC_tran", "PF00005.31", "2e-4", "1e-6")])
    write_chunk(d / "db_chunk_000003_pfam.domtbl.gz",
                [row("c", "Tiny", "PF99999.1", "1e-6", "1e-6")])
    out = tmp_path / "ceiling.tsv"
    proc = run(["--subset", str(d), "25", "--records-per-chunk", "10",
                "--target-z", "100", "--output", str(out)])
    assert proc.returncode == 0, proc.stderr
    rows = [line.split("\t") for line in out.read_text().splitlines()]
    assert rows[0] == ["family", "accession", "n_domains",
                       "max_seq_evalue", "max_dom_ievalue"]
    got = {r[0]: r for r in rows[1:]}
    # ABC_tran: chunk 2 gives the max full-sequence E, 2e-4 x 100/10 = 2e-3;
    # chunk 1 gives the max i-Evalue, 1e-5 x 100/10 = 1e-4
    assert got["ABC_tran"][1:] == ["PF00005.31", "2", "0.002", "0.0001"]
    # last chunk holds 25 - 2x10 = 5 records: 1e-6 x 100/5 = 2e-5
    assert got["Tiny"][1:] == ["PF99999.1", "1", "2e-05", "2e-05"]
    assert "families=2" in proc.stdout
    assert "max_seq_evalue > 0.001: 1 (50.0%)" in proc.stdout


def test_ceiling_refuses_a_partial_chunk_set(tmp_path):
    d = tmp_path / "pfam"
    d.mkdir()
    write_chunk(d / "db_chunk_000001_pfam.domtbl.gz",
                [row("a", "ABC_tran", "PF00005.31", "1e-4", "1e-5")])
    proc = run(["--subset", str(d), "25", "--records-per-chunk", "10",
                "--target-z", "100", "--output", str(tmp_path / "x.tsv")])
    assert proc.returncode != 0
    assert "expected 3 chunks" in proc.stderr
