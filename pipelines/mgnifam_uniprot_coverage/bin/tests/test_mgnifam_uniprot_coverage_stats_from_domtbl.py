import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest


SCRIPT = (Path(__file__).resolve().parents[1]
          / "mgnifam_uniprot_coverage_stats_from_domtbl.py")


def load_module():
    """Import the standalone script without requiring `bin` to be a package."""
    spec = importlib.util.spec_from_file_location("mgnifam_coverage", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_merge_intervals_counts_each_residue_once():
    merge = load_module().merge_intervals
    assert merge([]) == ([], 0)
    # overlapping, adjacent and disjoint spans, fed out of order
    assert merge([(10, 20), (1, 5)]) == ([(1, 5), (10, 20)], 16)
    assert merge([(1, 10), (5, 8)]) == ([(1, 10)], 10)
    assert merge([(1, 10), (11, 20)]) == ([(1, 20)], 20)


def test_subtract_intervals_leaves_only_exclusive_residues():
    sub = load_module().subtract_intervals
    # no mask: everything survives
    assert sub([(1, 10)], []) == ([(1, 10)], 10)
    # fully masked
    assert sub([(3, 6)], [(1, 10)]) == ([], 0)
    # masked in the middle, both flanks survive
    assert sub([(1, 20)], [(6, 10)]) == ([(1, 5), (11, 20)], 15)
    # mask overhangs each end
    assert sub([(5, 15)], [(1, 7), (12, 30)]) == ([(8, 11)], 4)
    # several spans against several mask spans, mask walked only forwards
    assert sub([(1, 5), (20, 30)], [(2, 3), (25, 40)]) == (
        [(1, 1), (4, 5), (20, 24)], 8)


def test_interval_is_masked_matches_subtraction():
    module = load_module()
    mask = ((10, 20), (30, 40))
    assert module.interval_is_masked(mask, 12, 18)
    assert module.interval_is_masked(mask, 10, 20)
    assert not module.interval_is_masked(mask, 9, 20)   # overhangs the left
    assert not module.interval_is_masked(mask, 10, 21)  # overhangs the right
    assert not module.interval_is_masked(mask, 21, 29)  # in the gap
    assert not module.interval_is_masked(mask, 1, 5)    # before everything
    # a span covered by two mask intervals it does not fit inside is not masked
    assert not module.interval_is_masked(mask, 15, 35)
    assert module.subtract_intervals([(15, 35)], list(mask))[1] > 0


def test_family_ids_round_trip_between_domtbl_and_lists():
    module = load_module()
    # domtbl query names for MGnifams are bare integers; the lists are padded
    assert module.bare_family("MGYF0000000444") == "444"
    assert module.bare_family("444") == "444"
    assert module.pretty_family("444") == "MGYF0000000444"
    # Pfam HMM names are not numeric and must pass through untouched
    assert module.pretty_family("1-cysPrx_C") == "1-cysPrx_C"


# --------------------------------------------------------------- CLI guards

def domtbl_row(target, tlen, query, ali_from, ali_to):
    """One hmmsearch --domtblout row; only the columns this script reads matter."""
    return " ".join(str(x) for x in [
        target, "-", tlen, query, "-", 100, "1e-30", 110.0, 0.0, 1, 1,
        "1e-30", "1e-30", 100.0, 0.0, 1, 100, ali_from, ali_to,
        ali_from, ali_to, 0.99, "a description with spaces",
    ])


@pytest.fixture
def chunk(tmp_path):
    """A one-chunk map output directory, plus the domtbl it came from."""
    domtbl = tmp_path / "uniprot_test_chunk_000001_mgnifams.domtbl"
    domtbl.write_text(
        "# a comment line\n"
        + domtbl_row("sp|P1|A", 200, "444", 10, 60) + "\n"
        + domtbl_row("sp|P2|B", 300, "7", 5, 105) + "\n"
    )
    outdir = tmp_path / "out"
    run(["map", "--domtbl", str(domtbl), "--outdir", str(outdir)])
    return domtbl, outdir


def run(argv, expect_ok=True):
    proc = subprocess.run([sys.executable, str(SCRIPT)] + argv,
                          capture_output=True, text=True)
    if expect_ok:
        assert proc.returncode == 0, proc.stderr
    return proc


def test_reduce_accepts_totals_as_plain_numbers(chunk):
    _domtbl, outdir = chunk
    run(["reduce", "--outdir", str(outdir),
         "--total-sequences", "1000", "--total-residues", "500000"])
    rows = (outdir / "reduced.tsv").read_text().splitlines()
    header = rows[0].split("\t")
    any_row = dict(zip(header, rows[1].split("\t")))
    assert any_row["category"] == "any"
    # 51 + 101 residues covered out of 500000
    assert any_row["n_residues"] == "152"
    assert float(any_row["pct_db_residues"]) == pytest.approx(100 * 152 / 500000)
    assert float(any_row["pct_db_sequences"]) == pytest.approx(100 * 2 / 1000)


def test_reduce_refuses_a_partial_chunk_set(chunk):
    _domtbl, outdir = chunk
    proc = run(["reduce", "--outdir", str(outdir), "--expect-chunks", "151"],
               expect_ok=False)
    assert proc.returncode != 0
    assert "expected 151 chunk summaries" in proc.stderr


def test_map_refuses_a_missing_mask(chunk, tmp_path):
    domtbl, _outdir = chunk
    proc = run(["map", "--domtbl", str(domtbl), "--outdir", str(tmp_path / "x"),
                "--mask", str(tmp_path / "does_not_exist.tsv.gz")],
               expect_ok=False)
    assert proc.returncode != 0
    # the whole point: it must not fall through to plain total coverage
    assert "refusing to report total coverage as exclusive" in proc.stderr
