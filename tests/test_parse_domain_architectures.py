import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "parse_domain_architectures.py"
FIXTURES = ROOT / "assets" / "test_data" / "domain_architecture"


def load_module():
    spec = importlib.util.spec_from_file_location("parse_domain_architectures", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mod = load_module()

# The clan lookups every chip test needs, matching clan_membership_dummy.csv.
FAMILY_TO_CLAN = {"100": "SF_1", "101": "SF_1", "102": "SF_1",
                  "200": "SF_2", "201": "SF_2", "202": "SF_2",
                  "300": "SF_3", "400": "SF_4"}
CLAN_TO_REP = {"SF_1": "100", "SF_2": "200", "SF_3": "300", "SF_4": "400"}
PFAM_MAPPING = {"PF00001": "Alpha domain", "PF00002": "Beta domain"}
BASE_URL = "http://mgnifams-demo.mgnify.org/details/"


def cluster(hits, fraction=0.5):
    return mod.cluster_hits(hits, FAMILY_TO_CLAN, fraction)


def chips(metadata, fraction=0.5):
    return mod.build_chips(metadata, FAMILY_TO_CLAN, CLAN_TO_REP, PFAM_MAPPING, BASE_URL, fraction)


# --- Phase 1: pure functions -------------------------------------------------


def test_load_clan_membership_maps_families_and_reps():
    family_to_clan, clan_to_rep = mod.load_clan_membership(FIXTURES / "clan_membership_dummy.csv")

    assert family_to_clan["102"] == "SF_1"
    assert family_to_clan["300"] == "SF_3"
    assert clan_to_rep["SF_2"] == "200"


def test_overlaps_is_true_just_above_half_the_shorter_hit():
    # 51 shared residues of a 100-long hit
    assert mod.overlaps((1, 100), (50, 149), 0.5) is True


def test_overlaps_is_false_at_exactly_half_the_shorter_hit():
    # 50 shared residues of a 100-long hit: the test is strictly greater
    assert mod.overlaps((1, 100), (51, 150), 0.5) is False


def test_cluster_hits_merges_same_clan_hits_on_the_same_region():
    groups = cluster([("100", 10, 100), ("101", 10, 100), ("102", 10, 100)])

    assert len(groups) == 1
    assert {hit[0] for hit in groups[0]} == {"100", "101", "102"}


def test_cluster_hits_keeps_same_clan_disjoint_hits_apart():
    groups = cluster([("100", 10, 100), ("101", 200, 290)])

    assert [[hit[0] for hit in group] for group in groups] == [["100"], ["101"]]


def test_cluster_hits_never_merges_different_clans():
    groups = cluster([("100", 10, 100), ("200", 10, 100)])

    assert len(groups) == 2


def test_cluster_hits_chains_through_a_shared_neighbour():
    # 200-201 and 201-202 each overlap by >50%, 200-202 does not: single linkage joins all three
    groups = cluster([("200", 1, 100), ("201", 40, 139), ("202", 80, 179)])

    assert len(groups) == 1
    assert {hit[0] for hit in groups[0]} == {"200", "201", "202"}


def test_unknown_family_never_merges_and_keeps_its_own_label():
    groups = cluster([("100", 10, 100), ("99999", 10, 100)])

    assert len(groups) == 2


def test_mgnifam_chip_labels_a_multi_family_group_with_its_clan():
    chip = mod.mgnifam_chip([("100", 10, 100), ("101", 10, 100)],
                            FAMILY_TO_CLAN, CLAN_TO_REP, BASE_URL)

    assert chip.id == "SF_1"
    assert chip.name == "MGnifam clan 1"
    assert chip.link == "http://mgnifams-demo.mgnify.org/details/MGYF0000000100"


def test_mgnifam_chip_labels_a_single_family_group_with_the_family():
    chip = mod.mgnifam_chip([("101", 10, 100)], FAMILY_TO_CLAN, CLAN_TO_REP, BASE_URL)

    assert chip.id == "101"
    assert chip.name == "MGnifam101"
    assert chip.link == "http://mgnifams-demo.mgnify.org/details/MGYF0000000101"


def test_mgnifam_chip_spans_the_whole_merged_group():
    chip = mod.mgnifam_chip([("100", 40, 139), ("101", 10, 100)],
                            FAMILY_TO_CLAN, CLAN_TO_REP, BASE_URL)

    assert (chip.start, chip.end) == (10, 139)


def test_pfam_chip_falls_back_to_the_accession_when_unmapped():
    chip = mod.pfam_chip(["PF09999", 1e-04, 45.0, 1, 30, 200, 230], PFAM_MAPPING)

    assert chip.id == "PF09999"
    assert chip.name == "PF09999"


def test_pfam_chip_uses_alignment_coordinates_not_hmm_coordinates():
    chip = mod.pfam_chip(["PF00001", 1e-05, 50.0, 1, 60, 200, 260], PFAM_MAPPING)

    assert (chip.start, chip.end) == (200, 260)


def test_string_to_hex_color_matches_the_previous_pipeline_output():
    # Asserted against assets/mgnifams_v2_results/output_db/domain_results/1.json
    name = "ATP synthase alpha/beta family, nucleotide-binding domain"

    assert mod.string_to_hex_color(name) == "#8b5115"


def test_decide_font_color_picks_white_on_a_dark_background():
    assert mod.decide_font_color("#310000") == "white"


def test_decide_font_color_picks_black_on_a_light_background():
    assert mod.decide_font_color("#ffffff") == "black"


# --- Phase 2: per-row assembly ----------------------------------------------


def test_pfam_sorts_after_the_mgnifam_when_its_alignment_starts_later():
    # Regression for the old parse_domains.py, which sorted Pfams by hmm_from (1) instead of
    # ali_from (200) and so rendered this Pfam first.
    key = mod.architecture_key(chips({
        "p": [["PF00001", 1e-05, 50.0, 1, 60, 200, 260]],
        "m": [["300", 1e-05, 50.0, 10, 100]],
    }))

    assert key == "300\tPF00001"


def test_equal_start_orders_the_shorter_chip_first():
    key = mod.architecture_key(chips({
        "m": [["100", 1e-05, 50.0, 10, 100], ["200", 1e-04, 45.0, 10, 50]],
    }))

    assert key == "200\t100"


def test_equal_start_and_length_orders_alphabetically_by_id():
    key = mod.architecture_key(chips({
        "p": [["PF00001", 1e-05, 50.0, 1, 60, 10, 100]],
        "m": [["300", 1e-05, 50.0, 10, 100]],
    }))

    assert key == "300\tPF00001"


def test_repeated_pfams_stay_as_separate_chips():
    key = mod.architecture_key(chips({
        "p": [["PF00002", 1e-05, 50.0, 1, 40, 10, 50],
              ["PF00002", 1e-06, 55.0, 1, 40, 100, 140]],
        "m": [["300", 1e-05, 50.0, 300, 380]],
    }))

    assert key == "PF00002\tPF00002\t300"


def test_fifty_and_fifty_one_percent_overlap_give_different_architectures():
    at_fifty = mod.architecture_key(chips({
        "m": [["100", 1e-05, 50.0, 1, 100], ["101", 1e-04, 45.0, 51, 150]],
    }))
    at_fifty_one = mod.architecture_key(chips({
        "m": [["100", 1e-05, 50.0, 1, 100], ["101", 1e-04, 45.0, 50, 149]],
    }))

    assert at_fifty == "100\t101"
    assert at_fifty_one == "SF_1"
    assert at_fifty != at_fifty_one


# --- Phase 3: streaming and counting -----------------------------------------


PROTEINS = FIXTURES / "proteins_dummy.csv.gz"
CLANS = FIXTURES / "clan_membership_dummy.csv"
MAPPING = FIXTURES / "pfam_mapping_dummy.tsv"


def count(metadata_rows, fraction=0.5):
    return mod.count_architectures(metadata_rows, FAMILY_TO_CLAN, CLAN_TO_REP, PFAM_MAPPING,
                                   BASE_URL, fraction, log_every=1000)


def test_rows_without_mgnifam_hits_are_skipped():
    rows = list(mod.iter_metadata(PROTEINS, use_prefilter=False))

    assert len(rows) == 11
    assert all(row.get("m") for row in rows)


def test_prefilter_and_plain_read_yield_the_same_rows():
    assert (list(mod.iter_metadata(PROTEINS, use_prefilter=True))
            == list(mod.iter_metadata(PROTEINS, use_prefilter=False)))


def test_family_is_credited_even_when_its_hit_merged_into_a_clan_chip():
    counters = count([{"m": [["100", 1e-05, 50.0, 10, 100],
                             ["101", 1e-04, 45.0, 10, 100],
                             ["102", 1e-03, 40.0, 10, 100]]}])

    assert counters["100"]["SF_1"] == 1
    assert counters["101"]["SF_1"] == 1
    assert counters["102"]["SF_1"] == 1


def test_family_hitting_twice_is_credited_once():
    counters = count([{"m": [["300", 1e-05, 50.0, 10, 100],
                             ["300", 1e-04, 45.0, 200, 290]]}])

    assert counters["300"] == {"300\t300": 1}


# --- Phase 4: output and validation ------------------------------------------


def test_architecture_json_matches_the_viewer_schema():
    payload = mod.architecture_json({"SF_1\tPF00001": 7}, 15, CLAN_TO_REP, PFAM_MAPPING, BASE_URL)
    container = payload["architecture_containers"][0]

    assert container["architecture_text"] == "7"
    assert container["domains"] == [
        {"id": "SF_1", "color": mod.string_to_hex_color("MGnifam clan 1"),
         "link": "http://mgnifams-demo.mgnify.org/details/MGYF0000000100",
         "name": "MGnifam clan 1", "font_color": mod.decide_font_color(
             mod.string_to_hex_color("MGnifam clan 1"))},
        {"id": "PF00001", "color": mod.string_to_hex_color("Alpha domain"),
         "link": "https://www.ebi.ac.uk/interpro/entry/pfam/PF00001",
         "name": "Alpha domain", "font_color": mod.decide_font_color(
             mod.string_to_hex_color("Alpha domain"))},
    ]


def test_architecture_json_keeps_the_top_n_in_descending_count_order():
    payload = mod.architecture_json({"100": 1, "101": 9, "102": 5}, 2,
                                    CLAN_TO_REP, PFAM_MAPPING, BASE_URL)

    assert [c["architecture_text"] for c in payload["architecture_containers"]] == ["9", "5"]


def test_families_without_hits_get_an_empty_file_and_are_reported(tmp_path, caplog):
    counters = {"100": {"SF_1": 1}}

    with caplog.at_level("WARNING"):
        missing = mod.write_outputs(counters, {"100", "400"}, tmp_path, 15,
                                    CLAN_TO_REP, PFAM_MAPPING, BASE_URL)

    assert missing == ["400"]
    assert json.loads((tmp_path / "400.json").read_text()) == {"architecture_containers": []}
    assert (tmp_path / "missing_families.txt").read_text() == "400\n"
    assert "1 expected famil" in caplog.text


def test_end_to_end_cli_run_over_the_dummy_input(tmp_path):
    import subprocess
    import sys

    subprocess.run([sys.executable, str(SCRIPT),
                    "--proteins", str(PROTEINS),
                    "--clan-membership", str(CLANS),
                    "--pfam-mapping", str(MAPPING),
                    "--output-dir", str(tmp_path)], check=True)

    assert sorted(p.stem for p in tmp_path.glob("*.json")) == [
        "100", "101", "102", "200", "201", "202", "300", "400"]

    top = json.loads((tmp_path / "300.json").read_text())["architecture_containers"][0]
    assert top["architecture_text"] == "2"
    assert [domain["id"] for domain in top["domains"]] == ["300", "PF00001"]

    assert json.loads((tmp_path / "400.json").read_text()) == {"architecture_containers": []}
    assert (tmp_path / "missing_families.txt").read_text() == "400\n"
