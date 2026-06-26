import importlib.util
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "bin"
    / "extract_true_novel_superfamilies.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location(
        "extract_true_novel_superfamilies", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_selects_superfamilies_where_all_member_families_are_novel(tmp_path):
    module = load_module()
    novel_ids = tmp_path / "novel_ids.txt"
    stats_csv = tmp_path / "superfamily_statistics.csv"
    output_txt = tmp_path / "true_novel_superfamilies.txt"

    novel_ids.write_text("1\n2\n3\n8\n")
    stats_csv.write_text(
        "Cluster Id,Family Rep Id,Family Size,Family Ids\n"
        "SF_1,2,3,1;2;3\n"
        "SF_2,4,2,3;4\n"
        "Singleton_8,8,1,8\n"
    )

    novel_superfamilies = module.true_novel_superfamilies(stats_csv, novel_ids)
    module.write_superfamily_ids(novel_superfamilies, output_txt)

    assert novel_superfamilies == [
        module.Superfamily("SF_1", ["1", "2", "3"]),
        module.Superfamily("Singleton_8", ["8"]),
    ]
    assert output_txt.read_text() == "SF_1\nSingleton_8\n"

