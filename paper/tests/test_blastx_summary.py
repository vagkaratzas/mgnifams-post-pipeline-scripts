#!/usr/bin/env python3
"""Tests for blastx_summary.py. Run: python paper/tests/test_blastx_summary.py

Fixture contig 1 (3,000 bp), member M on 1001-1600, family 7 with hit proteins P (Malassezia) and Q (Candidatus X):
    genes:   [g0 1-900]   [M 1001-1600]   [g2 1701-2400]
    P HSPs:   f+1 100-800 (aa 1-200)   f+2 1050-1500 (aa 190-340, 98%)   f+3 1750-2300 (aa 330-500)
    Q HSP:    f+2 1100-1400 (60%, lower bit score)
    R HSP:    f+2 1000-1600 (99%, best bit score, but R is a hit of another family: ignored)
-> best on M = P at 98%; P spans 3 frames and 3 genes; P residues covered 1-500 of 500.
Member N on contig 2 has no HSP -> member_hsp "no".
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import blastx_summary as b  # noqa: E402


def hsp(acc, pident, qs, qe, ss, se, frame, bits, slen=500):
    return dict(acc=acc, pident=pident, lo=min(qs, qe), hi=max(qs, qe), sstart=ss, send=se, qframe=frame,
                bitscore=bits, evalue=1e-10, slen=slen)


class Summary(unittest.TestCase):
    def test_helpers(self):
        self.assertEqual(b.covered([(1, 10), (5, 20), (30, 30)]), 21)
        self.assertEqual(b.genus("Malassezia globosa (strain X)"), "Malassezia")
        self.assertEqual(b.genus("Candidatus Nanosynbacter lyticus"), "Candidatus Nanosynbacter")
        self.assertEqual(b.acc("tr|A8PR02|A8PR02_MALGO"), "A8PR02")

    def test_summarise(self):
        hsps = {"1": [hsp("P", 90, 100, 800, 1, 200, 1, 300), hsp("P", 98, 1050, 1500, 190, 340, 2, 280),
                      hsp("P", 92, 1750, 2300, 330, 500, 3, 320), hsp("Q", 60, 1100, 1400, 10, 110, 2, 90),
                      hsp("R", 99, 1000, 1600, 1, 200, 2, 900)]}
        genes = {"1": [(1, 900, 1, "g0"), (1001, 1600, 1, "M"), (1701, 2400, 1, "g2")]}
        org = {"P": "Malassezia globosa", "Q": "Candidatus X y"}
        m = dict(family="7", protein="M", contig_id="1", contig_length="3000", start="1001", end="1600")
        r = b.summarise(m, hsps, genes, {"P", "Q"}, org)
        self.assertEqual((r["member_hsp"], r["hit"], r["pident"], r["genus"]), ("yes", "P", "98.0", "Malassezia"))
        self.assertEqual((r["frames"], r["genes_under_hit"], r["hit_coverage"]), (3, 3, "1.000"))
        n = dict(family="7", protein="N", contig_id="2", contig_length="900", start="1", end="300")
        self.assertEqual(b.summarise(n, hsps, genes, {"P"}, org)["member_hsp"], "no")
        pic = b.draw("1", m, hsps, genes, "P", width=30)
        self.assertIn("#", pic.splitlines()[1])
        self.assertEqual(len(pic.splitlines()), 5)  # title, genes, the 3 HSPs of P only
        rev = {"1": [(1, 900, -1, "g0"), (1001, 1600, 1, "M"), (1701, 2400, 1, "g2")]}  # g0 on the other strand
        self.assertEqual(b.summarise(m, hsps, rev, {"P", "Q"}, org)["genes_under_hit"], 2)
        touch = {"1": genes["1"] + [(1490, 1700, 1, "t")]}  # same strand, 11 bp under the f+2 HSP: not counted
        self.assertEqual(b.summarise(m, hsps, touch, {"P", "Q"}, org)["genes_under_hit"], 3)


if __name__ == "__main__":
    unittest.main()
