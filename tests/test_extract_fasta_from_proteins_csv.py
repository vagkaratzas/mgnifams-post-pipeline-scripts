import csv
import gzip
import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "extract_fasta_from_proteins_csv.py"


def load_module():
    spec = importlib.util.spec_from_file_location("extract_fasta_from_proteins_csv", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_extracts_fasta_from_gzipped_protein_csv(tmp_path):
    module = load_module()
    input_csv = tmp_path / "proteins.csv.gz"
    output_fasta = tmp_path / "proteins.fasta"

    with gzip.open(input_csv, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["mgyp", "sequence", "full_length", "cluster_size", "metadata"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "mgyp": "seq1",
                "sequence": "MPEPTIDE",
                "full_length": "false",
                "cluster_size": "1",
                "metadata": "{}",
            }
        )

    module.extract_fasta(input_csv, output_fasta)

    assert output_fasta.read_text() == ">seq1\nMPEPTIDE\n"
