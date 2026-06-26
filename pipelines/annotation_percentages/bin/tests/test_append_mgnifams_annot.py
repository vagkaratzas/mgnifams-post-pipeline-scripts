import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


PIPELINE_DIR = Path(__file__).resolve().parents[1]
BIN_DIR = PIPELINE_DIR / "bin"


def load_script(script_name):
    script_path = BIN_DIR / script_name
    spec = importlib.util.spec_from_file_location(script_path.stem, script_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class AppendMgnifamsChunkTests(unittest.TestCase):
    def test_merges_annotations_from_multiple_domtbl_chunks(self):
        append_mgnifams = load_script("append_mgnifams_annot.py")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_csv = tmpdir / "proteins.csv"
            output_csv = tmpdir / "proteins_mgnifams.csv"
            chunk_a = tmpdir / "chunk_a.domtbl"
            chunk_b = tmpdir / "chunk_b.domtbl"

            with input_csv.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["mgyp", "metadata"])
                writer.writeheader()
                writer.writerow(
                    {
                        "mgyp": "MGYP0001",
                        "metadata": json.dumps({"p": [["PF00001"]], "m": [["old"]]}),
                    }
                )
                writer.writerow(
                    {
                        "mgyp": "MGYP0002",
                        "metadata": json.dumps({"m": [["old"]]}),
                    }
                )

            chunk_a.write_text(
                "MGYP0001 - - MGYFAM0001 - - - - - - - - 1e-10 50.0 - - - 3 24 - - -\n"
            )
            chunk_b.write_text(
                "MGYP0001 - - MGYFAM0002 - - - - - - - - 1e-20 60.0 - - - 30 52 - - -\n"
            )

            annotation_data = append_mgnifams.parse_domtblouts([chunk_a, chunk_b])
            append_mgnifams.update_csv_with_annotations(
                input_csv, annotation_data, output_csv
            )

            with output_csv.open(newline="") as handle:
                rows = {row["mgyp"]: json.loads(row["metadata"]) for row in csv.DictReader(handle)}

            self.assertEqual(
                rows["MGYP0001"]["m"],
                [
                    ["MGYFAM0001", 1e-10, 50.0, 3, 24],
                    ["MGYFAM0002", 1e-20, 60.0, 30, 52],
                ],
            )
            self.assertEqual(rows["MGYP0001"]["p"], [["PF00001"]])
            self.assertNotIn("m", rows["MGYP0002"])

if __name__ == "__main__":
    unittest.main()
