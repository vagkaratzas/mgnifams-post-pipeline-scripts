import csv
import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "compare_annotation_stats.py"


def load_module():
    spec = importlib.util.spec_from_file_location("compare_annotation_stats", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_stats(path, label, sequence_pct, residue_pct, annotated_sequences, annotated_residues):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "label",
                "annotation_keys",
                "total_sequences",
                "annotated_sequences",
                "annotated_sequence_percentage",
                "total_amino_acids",
                "annotated_amino_acids",
                "annotated_amino_acid_percentage",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "label": label,
                "annotation_keys": "p",
                "total_sequences": 10,
                "annotated_sequences": annotated_sequences,
                "annotated_sequence_percentage": sequence_pct,
                "total_amino_acids": 100,
                "annotated_amino_acids": annotated_residues,
                "annotated_amino_acid_percentage": residue_pct,
            }
        )


def test_compares_annotation_stats_and_writes_absolute_and_relative_increases(tmp_path):
    module = load_module()
    before = tmp_path / "pfam.csv"
    after = tmp_path / "pfam_mgnifam.csv"
    output = tmp_path / "increase.csv"
    write_stats(before, "pfam", 30.0, 40.0, 3, 40)
    write_stats(after, "pfam_mgnifam", 50.0, 70.0, 5, 70)

    result = module.compare_stats(before, after)
    module.write_comparison(output, result)

    assert result["before_label"] == "pfam"
    assert result["after_label"] == "pfam_mgnifam"
    assert result["annotated_sequence_percentage_point_increase"] == 20.0
    assert result["annotated_sequence_relative_increase_percentage"] == 66.666667
    assert result["annotated_amino_acid_percentage_point_increase"] == 30.0
    assert result["annotated_amino_acid_relative_increase_percentage"] == 75.0

    with output.open(newline="") as handle:
        written = next(csv.DictReader(handle))
    assert written["annotated_amino_acid_percentage_point_increase"] == "30.000000"


def test_relative_increase_is_zero_when_before_percentage_is_zero(tmp_path):
    module = load_module()
    before = tmp_path / "pfam.csv"
    after = tmp_path / "pfam_mgnifam.csv"
    write_stats(before, "pfam", 0.0, 0.0, 0, 0)
    write_stats(after, "pfam_mgnifam", 20.0, 25.0, 2, 25)

    result = module.compare_stats(before, after)

    assert result["annotated_sequence_relative_increase_percentage"] == 0.0
    assert result["annotated_amino_acid_relative_increase_percentage"] == 0.0
