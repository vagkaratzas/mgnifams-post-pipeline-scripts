import csv
import gzip
import importlib.util
import json
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "append_mgnifams_annot.py"


def load_module():
    spec = importlib.util.spec_from_file_location("append_mgnifams_annot", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_parses_hmmsearch_domtblout_with_domain_filter_and_aligned_coordinates(tmp_path):
    module = load_module()
    domtbl = tmp_path / "hits.domtbl"
    domtbl.write_text(
        "# target name accession tlen query name accession qlen ...\n"
        "seq1 - 100 7 - 120 1e-20 50.0 0.0 1 1 2e-10 3e-10 45.0 0.0 4 60 8 64 2 70 0.95 -\n"
        "seq1 - 100 8 - 90 1e-20 40.0 0.0 1 1 0.2 0.3 10.0 0.0 1 20 10 30 5 35 0.50 -\n"
        "seq2 - 80 9 - 70 1e-5 25.0 0.0 1 1 1e-4 8e-4 20.0 0.0 2 40 11 49 9 55 0.80 -\n"
    )

    annotations = module.parse_domtblout(domtbl, domain_evalue_threshold=0.001)

    assert annotations == {
        "seq1": [["7", 3e-10, 45.0, 8, 64]],
        "seq2": [["9", 8e-4, 20.0, 11, 49]],
    }


def test_updates_gzipped_csv_and_replaces_existing_m_annotations(tmp_path):
    module = load_module()
    input_csv = tmp_path / "proteins.csv.gz"
    output_csv = tmp_path / "proteins_mgnifams.csv.gz"
    rows = [
        {
            "mgyp": "seq1",
            "sequence": "AAAA",
            "full_length": "false",
            "cluster_size": "1",
            "metadata": json.dumps({"p": [["PF1", 1e-6, 30.0, 1, 2, 1, 2]], "m": [["old", 1, 1, 1, 1]]}),
        },
        {
            "mgyp": "seq2",
            "sequence": "AAAA",
            "full_length": "false",
            "cluster_size": "1",
            "metadata": "{}",
        },
    ]
    with gzip.open(input_csv, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["mgyp", "sequence", "full_length", "cluster_size", "metadata"],
        )
        writer.writeheader()
        writer.writerows(rows)

    module.update_csv_with_annotations(
        input_csv,
        {"seq1": [["7", 3e-10, 45.0, 8, 64]]},
        output_csv,
    )

    with gzip.open(output_csv, "rt", newline="") as handle:
        updated = list(csv.DictReader(handle))

    assert json.loads(updated[0]["metadata"]) == {
        "p": [["PF1", 1e-06, 30.0, 1, 2, 1, 2]],
        "m": [["7", 3e-10, 45.0, 8, 64]],
    }
    assert json.loads(updated[1]["metadata"]) == {}
