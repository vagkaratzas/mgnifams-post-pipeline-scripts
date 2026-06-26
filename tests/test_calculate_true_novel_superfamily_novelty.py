import importlib.util
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "bin"
    / "calculate_true_novel_superfamily_novelty.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location(
        "calculate_true_novel_superfamily_novelty", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_filters_true_novel_superfamilies_without_novel_filtered_hits(tmp_path):
    module = load_module()
    novel_filtered_csv = tmp_path / "novel_filtered.csv"
    stats_csv = tmp_path / "superfamily_statistics.csv"
    true_novel_txt = tmp_path / "true_novel_superfamilies.txt"

    novel_filtered_csv.write_text(
        "ID,Full size\n"
        "MGYF0000000001,100\n"
        "MGYF0000000003,200\n"
    )
    stats_csv.write_text(
        "Cluster Id,Family Rep Id,Family Size,Family Ids\n"
        "SF_with_hits,1,3,1;2;3\n"
        "SF_without_hits,4,2,4;5\n"
        "SF_not_true_novel,6,1,1\n"
    )
    true_novel_txt.write_text("SF_with_hits\nSF_without_hits\n")

    rows = module.true_novel_superfamily_novelty_scores(
        novel_filtered_csv,
        stats_csv,
        true_novel_txt,
    )

    assert [row["Cluster Id"] for row in rows] == ["SF_with_hits"]
    assert rows[0]["Novelty Score"] == "66.666667"
