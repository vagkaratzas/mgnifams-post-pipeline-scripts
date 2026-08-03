import importlib.util
from pathlib import Path


SCRIPT = (Path(__file__).resolve().parents[1] / "bin"
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


def test_family_ids_round_trip_between_domtbl_and_lists():
    module = load_module()
    # domtbl query names for MGnifams are bare integers; the lists are padded
    assert module.bare_family("MGYF0000000444") == "444"
    assert module.bare_family("444") == "444"
    assert module.pretty_family("444") == "MGYF0000000444"
    # Pfam HMM names are not numeric and must pass through untouched
    assert module.pretty_family("1-cysPrx_C") == "1-cysPrx_C"
