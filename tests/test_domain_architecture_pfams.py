import csv
import importlib.util
import json
import sqlite3
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "annotate_novel_through_domain_architecture.py"


def load_module():
    spec = importlib.util.spec_from_file_location("annotate_novel_through_domain_architecture", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_sqlite(path: Path) -> None:
    connection = sqlite3.connect(path)
    connection.execute("CREATE TABLE mgnifam (id INTEGER PRIMARY KEY, domain_blob BLOB)")
    rows = [
        (
            1,
            {
                "architecture_containers": [
                    {
                        "domains": [
                            {"id": "1"},
                            {"id": "PF00002"},
                            {"id": "PF00001"},
                            {"id": "PF00002"},
                        ]
                    }
                ]
            },
        ),
        (2, {"architecture_containers": [{"domains": [{"id": "2"}]}]}),
        (3, {"architecture_containers": [{"domains": [{"id": "PF00003"}]}]}),
    ]
    connection.executemany(
        "INSERT INTO mgnifam (id, domain_blob) VALUES (?, ?)",
        [(family_id, json.dumps(payload).encode()) for family_id, payload in rows],
    )
    connection.commit()
    connection.close()


def test_extracts_unique_pfams_for_novel_mgnifams(tmp_path):
    module = load_module()
    sqlite_path = tmp_path / "mgnifams.sqlite3"
    novel_ids = tmp_path / "novel_ids.txt"
    make_sqlite(sqlite_path)
    novel_ids.write_text("1\n2\n99\n")

    family_pfams = module.novel_family_pfams(sqlite_path, novel_ids)

    assert family_pfams == {
        "1": ["PF00001", "PF00002"],
        "2": [],
    }


def test_marks_true_novel_clans_with_and_without_architecture_pfams(tmp_path):
    module = load_module()
    sqlite_path = tmp_path / "mgnifams.sqlite3"
    novel_ids = tmp_path / "novel_ids.txt"
    clan_membership = tmp_path / "clan_membership.csv"
    pfam_output = tmp_path / "family_pfams.csv"
    clan_output = tmp_path / "transient_clans.csv"
    true_novel_output = tmp_path / "true_novel_without_pfams.txt"
    make_sqlite(sqlite_path)
    novel_ids.write_text("1\n2\n3\n")
    clan_membership.write_text(
        "Cluster Id,Family Rep Id,Family Size,Family Ids\n"
        "SF_annotated,1,2,1;2\n"
        "Singleton_2,2,1,2\n"
        "SF_mixed_not_true_novel,4,2,1;4\n"
        "SF_other_annotated,3,1,3\n"
    )

    family_pfams = module.novel_family_pfams(sqlite_path, novel_ids)
    clan_rows = module.transiently_annotated_clans(clan_membership, novel_ids, family_pfams)
    module.write_family_pfams(family_pfams, pfam_output)
    module.write_clan_annotations(clan_rows, clan_output)
    module.write_true_novel_without_pfams(clan_rows, true_novel_output)

    assert clan_rows == [
        module.ClanAnnotation("SF_annotated", True, ["PF00001", "PF00002"]),
        module.ClanAnnotation("Singleton_2", False, []),
        module.ClanAnnotation("SF_other_annotated", True, ["PF00003"]),
    ]
    assert true_novel_output.read_text() == "Singleton_2\n"

    with open(pfam_output, newline="") as handle:
        assert list(csv.reader(handle)) == [
            ["MGnifam id", "Pfam annotations"],
            ["1", "PF00001;PF00002"],
            ["2", ""],
            ["3", "PF00003"],
        ]

    with open(clan_output, newline="") as handle:
        assert list(csv.reader(handle)) == [
            ["Clan Id", "Annotated", "Pfam annotations"],
            ["SF_annotated", "True", "PF00001;PF00002"],
            ["Singleton_2", "False", ""],
            ["SF_other_annotated", "True", "PF00003"],
        ]
