import csv
import gzip
import importlib.util
import json
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "calculate_annotation_stats.py"


def load_module():
    spec = importlib.util.spec_from_file_location("calculate_annotation_stats", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_csv_gz(path, rows):
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["mgyp", "sequence", "full_length", "cluster_size", "metadata"],
        )
        writer.writeheader()
        writer.writerows(rows)


def test_calculates_pfam_only_annotation_percentages_from_gz_csv(tmp_path):
    module = load_module()
    input_csv = tmp_path / "proteins.csv.gz"
    output_csv = tmp_path / "pfam_stats.csv"
    rows = [
        {
            "mgyp": "seq1",
            "sequence": "AAAAAA",
            "full_length": "false",
            "cluster_size": "1",
            "metadata": json.dumps(
                {
                    "p": [
                        ["PF00001", 1e-5, 20.0, 1, 3, 1, 3],
                        ["PF00002", 1e-6, 22.0, 1, 3, 3, 5],
                    ]
                }
            ),
        },
        {
            "mgyp": "seq2",
            "sequence": "AAAA",
            "full_length": "false",
            "cluster_size": "1",
            "metadata": "{}",
        },
    ]
    write_csv_gz(input_csv, rows)

    stats = module.calculate_stats(input_csv, ["p"], "pfam")
    module.write_stats(output_csv, stats)

    assert stats["label"] == "pfam"
    assert stats["annotation_keys"] == "p"
    assert stats["total_sequences"] == 2
    assert stats["annotated_sequences"] == 1
    assert stats["annotated_sequence_percentage"] == 50.0
    assert stats["total_amino_acids"] == 10
    assert stats["annotated_amino_acids"] == 5
    assert stats["annotated_amino_acid_percentage"] == 50.0

    with output_csv.open(newline="") as handle:
        written = next(csv.DictReader(handle))
    assert written["label"] == "pfam"
    assert written["annotated_amino_acid_percentage"] == "50.000000"


def test_calculates_combined_pfam_and_mgnifam_coverage_without_double_counting(tmp_path):
    module = load_module()
    input_csv = tmp_path / "proteins.csv"
    rows = [
        {
            "mgyp": "seq1",
            "sequence": "AAAAAA",
            "full_length": "false",
            "cluster_size": "1",
            "metadata": json.dumps(
                {
                    "p": [["PF00001", 1e-5, 20.0, 1, 3, 1, 3]],
                    "m": [["MGF00001", 1e-7, 35.0, 3, 6]],
                }
            ),
        }
    ]
    with input_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["mgyp", "sequence", "full_length", "cluster_size", "metadata"],
        )
        writer.writeheader()
        writer.writerows(rows)

    stats = module.calculate_stats(input_csv, ["p", "m"], "pfam_mgnifam")

    assert stats["total_sequences"] == 1
    assert stats["annotated_sequences"] == 1
    assert stats["total_amino_acids"] == 6
    assert stats["annotated_amino_acids"] == 6
    assert stats["annotated_amino_acid_percentage"] == 100.0
