#!/usr/bin/env python3
"""Tests for architectures.py. Run: python paper/tests/test_architectures.py"""

import gzip
import io
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import architectures as a  # noqa: E402

# model of 10 match columns; protein 1: one copy split in two rows (model 1-5, 6-10), full; protein 2: one partial row
# (model 3-8); protein 3: two full copies (tandem); protein 4: one full row
STO = """# STOCKHOLM 1.0
1/11-15     ACDEF.....
1/31-35     .....GHIKL
2/5-10      ..DEFGHI..
3/1-10      ACDEFGHIKL
3/21-30     ACDEFGHIKL
4/101-110   ACDEFGHIKL
#=GC RF     xxxxxxxxxx
//
"""


def write_msa(folder, family="7"):
    path = Path(folder) / f"{family}.sto.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(STO)
    return path


class Copies(unittest.TestCase):
    def test_split_copy_counts_once_and_ends(self):
        with tempfile.TemporaryDirectory() as d:
            family, copies = a.mgnifam_copies(write_msa(d))
        self.assertEqual(family, "7")
        self.assertEqual(copies[1], [(11, 35, True, True, True)])
        self.assertEqual(copies[2], [(5, 10, False, False, False)])
        self.assertEqual(len(copies[3]), 2)


class Resolve(unittest.TestCase):
    clans = {"PF02518": "CL0025", "PF14501": "CL0025", "PF00512": "CL0001"}

    def test_same_clan_overlap_keeps_best(self):
        hits = [("PF14501", 1e-5, 30.0, 300, 420), ("PF02518", 1e-20, 90.0, 310, 430)]
        self.assertEqual([h[0] for h in a.resolve(hits, self.clans)], ["PF02518"])

    def test_other_clan_or_no_overlap_kept(self):
        hits = [("PF02518", 1e-20, 90.0, 310, 430), ("PF00512", 1e-9, 40.0, 400, 460),
                ("PF14501", 1e-5, 30.0, 431, 520)]
        self.assertEqual(sorted(h[0] for h in a.resolve(hits, self.clans)), ["PF00512", "PF02518", "PF14501"])

    def test_family_without_clan_is_its_own(self):
        hits = [("PF03382", 1e-3, 20.0, 100, 130), ("PF03382", 1e-6, 25.0, 120, 150), ("PF03382", 1e-6, 25.0, 151, 180)]
        self.assertEqual([(h[3], h[4]) for h in a.resolve(hits, {})], [(120, 150), (151, 180)])


class Architectures(unittest.TestCase):
    def test_order_and_tokens(self):
        copies = [(11, 35, True, True, True)]
        pfams = [("PF00528", 1e-9, 50.0, 60, 200), ("PF16296", 1e-9, 50.0, 1, 9)]
        key, mgnifam = a.architecture("7", copies, pfams)
        self.assertEqual(key, ("PF16296", "7", "PF00528"))
        self.assertEqual(mgnifam, [(True, True, True)])

    def test_aggregate_rank_and_cover(self):
        with tempfile.TemporaryDirectory() as d:
            family, copies = a.mgnifam_copies(write_msa(d))
        pfams = {4: [("PF00528", 1e-9, 50.0, 200, 300)]}
        rows = a.aggregate(family, copies, pfams, {})
        self.assertEqual([(r["architecture"], r["proteins"]) for r in rows],
                         [(("7",), 2), (("7", "7"), 1), (("7", "PF00528"), 1)])
        self.assertEqual(rows[0]["ends"], [[1, 1, 1]])  # of 2 proteins, 1 full, 1 curved N, 1 curved C
        self.assertEqual(a.rank_at([r["proteins"] for r in rows], 0.95), 3)
        self.assertEqual(a.rank_at([r["proteins"] for r in rows], 0.5), 1)


class Parquet(unittest.TestCase):
    def test_pfam_hits_and_main(self):
        import pyarrow as pa
        import pyarrow.parquet as pq
        with tempfile.TemporaryDirectory() as d:
            msa = write_msa(d)
            table = pa.table({"protein_id": [4, 4, 99], "pfam_accession": [528, 14501, 1],
                              "i_evalue": [1e-9, 1e-3, 1e-9], "score": [50.0, 20.0, 50.0], "hmm_from": [1, 1, 1],
                              "hmm_to": [9, 9, 9], "env_from": [200, 210, 1], "env_to": [300, 290, 9]})
            parquet = Path(d) / "pfam.parquet"
            pq.write_table(table, parquet, row_group_size=2)
            expected = {4: [("PF00528", 1e-9, 50.0, 200, 300), ("PF14501", 1e-3, 20.0, 210, 290)]}
            self.assertEqual(a.pfam_hits(parquet, [4], cpus=2), expected)
            self.assertEqual(a.pfam_hits(parquet, [4]), expected)
            self.assertEqual(a._scan_row_group(1), [])  # row group of protein 99 only: skipped by its min / max
            clans = Path(d) / "clans.tsv.gz"
            with gzip.open(clans, "wt") as handle:
                handle.write("PF00528\tCL0404\tABC-membrane\tBPD_transp_1\tBinding-protein-dependent\n"
                             "PF14501\tCL0404\tABC-membrane\tOther\tx\n")
            out = Path(d) / "arch.tsv.gz"
            stdout = io.StringIO()
            sys.stdout, saved = stdout, sys.stdout
            try:
                a.main([str(msa), "--pfam", str(parquet), "--clans", str(clans), "--out", str(out)])
            finally:
                sys.stdout = saved
            lines = gzip.open(out, "rt").read().splitlines()
        self.assertEqual(lines[0].split("\t"), list(a.FIELDS))
        last = dict(zip(a.FIELDS, lines[-1].split("\t")))
        self.assertEqual((last["architecture"], last["names"], last["proteins"]), ("7 PF00528", "7 BPD_transp_1", "1"))
        self.assertIn("7\t4\t3\t3\t", stdout.getvalue())  # family, proteins, architectures, rank at 95%


if __name__ == "__main__":
    unittest.main()
