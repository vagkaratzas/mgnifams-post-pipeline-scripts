#!/usr/bin/env python3
"""Tests for neighbourhood.py. Run: python paper/tests/test_neighbourhood.py

Fixture (family 25, members 1-4; protein 4 is on no contig):

contig 100 (20 kb), member 1 on +, longest of its two contigs (101 is 5 kb):
    <-60   50 PF00005->  51 PF00496->  52 PF00528->  [1]->  70 (rep 3)->  <-71 PF00005   72 PF00528->
     -4         -3            -2            -1        0          +1           +2             +3
    directon: -3..+1 (stopped by 60 and 71 on the other strand) -> complete
    census: 4 genes and < 10 kb upstream -> incomplete
contig 200 (3 kb), member 2 on -, read right to left:
    [2]<-  80 PF00005<-  81->        directon: 80 at -1, then the contig edge downstream -> incomplete
contig 300 (12 kb), member 3 on + with 5 same-strand genes on each side (gene 302 has rep 1) -> both complete
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
import neighbourhood as n  # noqa: E402

STO = """# STOCKHOLM 1.0
1/1-10   ACDEFGHIKL
2/1-10   ACDEFGHIKL
3/1-10   ACDEFGHIKL
4/1-10   ACDEFGHIKL
#=GC RF  xxxxxxxxxx
//
"""
CLANS = "PF00496\tCL0177\tPBP\tSBP_bac_5\t.\nPF00005\tCL0023\tP-loop_NTPase\tABC_tran\t.\nPF00528\t\t\tBPD_transp_1\t.\n"
# (contig, contig_length, start, end, strand, protein, cluster_rep)
GENES = [
    (100, 20000, 500, 900, -1, 60, 60), (100, 20000, 1000, 1800, 1, 50, 50), (100, 20000, 2000, 3000, 1, 51, 51),
    (100, 20000, 3100, 4900, 1, 52, 52), (100, 20000, 5000, 5900, 1, 1, 1), (100, 20000, 6000, 6900, 1, 70, 3),
    (100, 20000, 7000, 7500, -1, 71, 71), (100, 20000, 8000, 8500, 1, 72, 72),
    (101, 5000, 100, 1000, 1, 1, 1), (101, 5000, 1100, 2000, 1, 90, 90),
    (200, 3000, 100, 1000, -1, 2, 2), (200, 3000, 1100, 2000, -1, 80, 80), (200, 3000, 2100, 2900, 1, 81, 81),
] + [(300, 12000, k * 1000 + 100, k * 1000 + 900, 1, 300 + k if k != 5 else 3, 300 + k if k != 2 else 1)
     for k in range(11)] + [(900, 4000, 100, 900, 1, 9, 9)]
PFAM = [(50, 5), (51, 496), (52, 528), (71, 5), (72, 528), (80, 5)]


def write_fixture(d):
    import pyarrow as pa
    import pyarrow.parquet as pq
    d = Path(d)
    with gzip.open(d / "25.sto.gz", "wt") as h:
        h.write(STO)
    with gzip.open(d / "clans.tsv.gz", "wt") as h:
        h.write(CLANS)
    cols = ["contig_id", "contig_length", "start_position", "end_position", "strand", "protein_id", "cluster_rep"]
    rows = [dict(zip(cols, g)) for g in GENES]
    genes = pa.Table.from_pylist(rows[::-1])  # unsorted, as the real file
    pq.write_table(genes, d / "metadata.parquet", row_group_size=3)
    reps = sorted((r for r in rows if r["protein_id"] == r["cluster_rep"] or r["protein_id"] in (1, 2, 3)),
                  key=lambda r: r["protein_id"])
    pq.write_table(pa.Table.from_pylist([{k: r[k] for k in n.OCCURRENCE} for r in reps]), d / "clusters.parquet",
                   row_group_size=2)
    pq.write_table(pa.Table.from_pylist([{"protein_id": p, "pfam_accession": a} for p, a in PFAM]),
                   d / "pfam.parquet", row_group_size=2)
    return d


def contig(c):
    return sorted(g[2:5] + g[5:7] for g in GENES if g[0] == c)


class Windows(unittest.TestCase):
    def test_plus_strand_directon_stops_at_strand_switch(self):
        genes = contig(100)
        directon, census, d_complete, c_complete = n.windows(genes, 4, 20000)
        self.assertEqual([(o, genes[j][3]) for o, j in sorted(directon)], [(-3, 50), (-2, 51), (-1, 52), (1, 70)])
        self.assertTrue(d_complete)
        self.assertEqual(sorted(genes[j][3] for _, j in census), [50, 51, 52, 60, 70, 71, 72])
        self.assertFalse(c_complete)  # 4 genes and 4,999 bp upstream

    def test_minus_strand_reads_right_to_left_and_edge_is_incomplete(self):
        genes = contig(200)
        directon, _, d_complete, _ = n.windows(genes, 0, 3000)
        self.assertEqual([(o, genes[j][3]) for o, j in directon], [(-1, 80)])
        self.assertFalse(d_complete)

    def test_five_same_strand_genes_complete_both_windows(self):
        genes = contig(300)
        directon, census, d_complete, c_complete = n.windows(genes, 5, 12000)
        self.assertEqual(sorted(o for o, _ in directon), [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
        self.assertTrue(d_complete and c_complete)
        self.assertEqual(len(census), 10)

    def test_census_respects_bp_cap(self):
        genes = [(1, 100, 1, 1, 1), (20000, 20100, 1, 2, 2), (40000, 40100, 1, 3, 3)]
        _, census, _, _ = n.windows(genes, 1, 50000)
        self.assertEqual(census, [])


class Names(unittest.TestCase):
    def test_pfam_member_and_dark(self):
        self.assertEqual(n.gene_names((0, 0, 1, 51, 51), {51: ("PF00496",)}, {1, 3}, "25"), ["PF00496"])
        self.assertEqual(n.gene_names((0, 0, 1, 70, 3), {}, {1, 3}, "25"), ["MGnifam_25"])
        self.assertEqual(n.gene_names((0, 0, 1, 60, 60), {}, {1, 3}, "25"), ["cluster:60"])

    def test_expand_clan_and_family(self):
        self.assertEqual(n.expand("CL0177", {"CL0177": {"PF00496", "PF01547"}}, "25"), {"PF00496", "PF01547"})
        self.assertEqual(n.expand("PF00528,family", {}, "25"), {"PF00528", "MGnifam_25"})
        with self.assertRaises(ValueError):
            n.expand("CL9999", {}, "25")


class EndToEnd(unittest.TestCase):
    def run_main(self, cpus):
        with tempfile.TemporaryDirectory() as d:
            d = write_fixture(d)
            out = io.StringIO()
            with redirect_stdout(out):
                n.main([str(d / "25.sto.gz"), "--clusters-metadata", str(d / "clusters.parquet"),
                        "--metadata", str(d / "metadata.parquet"), "--pfam", str(d / "pfam.parquet"),
                        "--clans", str(d / "clans.tsv.gz"), "--out", str(d / "p.tsv.gz"),
                        "--proteins", str(d / "r.tsv.gz"), "--cpus", str(cpus)])
            with gzip.open(d / "p.tsv.gz", "rt") as h:
                partners = {r.split("\t")[0]: r.rstrip("\n").split("\t") for r in h}
            with gzip.open(d / "r.tsv.gz", "rt") as h:
                proteins = {r.split("\t")[0]: r.rstrip("\n").split("\t") for r in h}
            return out.getvalue(), partners, proteins

    def test_counts_targets_and_files(self):
        for cpus in (1, 2):
            stdout, partners, proteins = self.run_main(cpus)
            self.assertIn("# 3 member proteins on a contig (1 not found), 3 contigs", stdout)
            self.assertIn("directon window complete 2, census window complete 1", stdout)
            self.assertIn("CL0177\t1\t50\t1\n", stdout)
            self.assertIn("PF00005\t1\t50\t2\n", stdout)  # protein 2 has it, but its window is cut by the edge
            self.assertIn("PF00528,family\t2\t100\t2\n", stdout)  # protein 1 (PF00528), protein 3 (member rep)
            self.assertIn("all 3 targets\t1\t50\t\n", stdout)
            p = partners["PF00528"]
            self.assertEqual(p[2:5] + [p[8], p[9]], ["1", "1", "50", "-1", "100"])
            self.assertEqual(partners["PF00005"][2:4], ["2", "1"])
            # MGnifam_25: 2 of 15 window genes (4 + 1 + 10) vs 2 of 19 non-member genes on contigs 100, 200, 300
            self.assertAlmostEqual(float(partners["MGnifam_25"][11]), (2 / 15) / (2 / 19))
            self.assertEqual(proteins["1"][1], "100")  # longest contig, not 101
            self.assertEqual(proteins["1"][7], "PF00005@-3 PF00496@-2 PF00528@-1 MGnifam_25@+1")
            self.assertEqual(proteins["2"][3:7], ["-1", "3", "0", "0"])
            self.assertNotIn("cluster:90", partners)


if __name__ == "__main__":
    unittest.main(verbosity=1)
