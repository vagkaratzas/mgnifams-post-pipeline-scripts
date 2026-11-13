#!/usr/bin/env python3
"""Tests for contig_lookup.py. Run: python paper/tests/test_contig_lookup.py

Fixture: family 7 has members 1, 2, 3; family 9 has member 2.
  protein 1 on contig 100 (5 kb, ERZ1, v5.0) and contig 101 (9 kb, ERZ2, v4.1) -> 101, the longest
  protein 2 on contig 200 (3 kb, ERZ1) -> in both families
  protein 3 on no contig -> not placed
Genes: contig 101 holds protein 1 and two neighbours; contig 200 holds protein 2; contig 999 is not a member contig.
"""

import gzip
import io
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import contig_lookup as c  # noqa: E402


def sto(ids):
    return "# STOCKHOLM 1.0\n" + "".join(f"{i}/1-10   ACDEFGHIKL\n" for i in ids) + "#=GC RF  xxxxxxxxxx\n//\n"


class ContigLookup(unittest.TestCase):
    def test_run(self):
        import pyarrow as pa
        import pyarrow.parquet as pq
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            for fam, ids in ((7, (1, 2, 3)), (9, (2,))):
                with gzip.open(d / f"{fam}.sto.gz", "wt") as f:
                    f.write(sto(ids))
            occ = [  # protein, contig, length, start, end, strand, assembly, name, version, truncation, study, biome
                (1, 100, 5000, 10, 300, 1, "ERZ1", "c100", 5.0, "11", "MGYS1", "root:A"),
                (1, 101, 9000, 400, 700, -1, "ERZ2", "c101", 4.1, "01", "MGYS2", "root:B"),
                (2, 200, 3000, 1, 600, 1, "ERZ1", "c200", 5.0, "00", "MGYS1", "root:A"),
            ]
            OCC = c.OCCURRENCE + c.EXTRA
            pq.write_table(pa.table({k: [r[i] for r in occ] for i, k in enumerate(OCC)}), d / "clusters.parquet",
                           row_group_size=1)
            genes = [(101, 1, 300, 1, 50), (101, 400, 700, -1, 1), (101, 800, 900, -1, 51),
                     (200, 1, 600, 1, 2), (999, 1, 90, 1, 60)]
            pq.write_table(pa.table({k: [g[i] for g in genes] for i, k in enumerate(c.GENE_COLUMNS)}),
                           d / "meta.parquet", row_group_size=2)
            out = io.StringIO()
            with redirect_stdout(out):
                c.main([str(d / "7.sto.gz"), str(d / "9.sto.gz"), "--clusters-metadata", str(d / "clusters.parquet"),
                        "--out", str(d / "m.tsv"), "--metadata", str(d / "meta.parquet"),
                        "--genes", str(d / "g.tsv.gz")])
            rows = [l.split("\t") for l in (d / "m.tsv").read_text().splitlines()]
            self.assertEqual(rows[0], list(c.MEMBER_FIELDS))
            self.assertEqual([r[:6] for r in rows[1:]], [
                ["7", "MGYP000000000001", "101", "ERZ2", "c101", "4.1"],
                ["7", "MGYP000000000002", "200", "ERZ1", "c200", "5.0"],
                ["9", "MGYP000000000002", "200", "ERZ1", "c200", "5.0"],
            ])
            self.assertEqual(rows[1][6:], ["9000", "400", "700", "-1", "01", "2", "MGYS2", "root:B"])
            g = gzip.open(d / "g.tsv.gz", "rt").read().splitlines()
            self.assertEqual(len(g), 5)  # header + 4 genes; contig 999 left out
            self.assertNotIn("999", "".join(g))
            text = out.getvalue()
            self.assertIn("7\t3\t2\t2\t2\t4.1:1 5.0:1\t6000.0", text)
            self.assertIn("2 proteins placed, 2 contigs, 2 assemblies", text)


if __name__ == "__main__":
    unittest.main()
