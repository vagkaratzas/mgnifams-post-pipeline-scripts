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
