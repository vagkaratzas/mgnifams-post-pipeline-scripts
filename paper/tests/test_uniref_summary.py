#!/usr/bin/env python3
"""Tests for uniref_summary.py. Run: python paper/tests/test_uniref_summary.py

Fixture contig 1 (3,000 bp), genes g0 1-900 (+), M 1001-1600 (+, the member), g2 1701-2400 (-):
    f+1 100-800   Malassezia restricta 80%      -> g0 Malassezia
    f+2 1050-1500 Malassezia globosa 97%        -> best on M (higher bits than the Basidiomycota HSP)
    f+3 1020-1580 Basidiomycota 99%, fewer bits
    f+1 1750-2300 Escherichia coli              -> on the + strand: not g2's (-) hit
    f-1 1800-2300 Bathycoccus prasinos          -> g2 Bathycoccus
-> member Malassezia 97%, same genus; contig 3 genes, 3 hit, majority Malassezia (2 genes).
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import uniref_summary as u  # noqa: E402


def hsp(qs, qe, frame, bits, pident, taxon, taxid="1"):
    return dict(lo=qs, hi=qe, qframe=frame, bitscore=bits, pident=pident, evalue=1e-20, sseqid=f"UniRef90_{bits}",
                stitle=f"UniRef90_{bits} some protein n=2 Tax={taxon} TaxID={taxid} RepID=X")


HSPS = [hsp(100, 800, 1, 300, 80, "Malassezia restricta"), hsp(1050, 1500, 2, 280, 97, "Malassezia globosa", "76773"),
        hsp(1020, 1580, 3, 200, 99, "Basidiomycota"), hsp(1750, 2300, 1, 150, 60, "Escherichia coli"),
        hsp(1800, 2300, -1, 120, 70, "Bathycoccus prasinos")]
GENES = [(1, 900, 1), (1001, 1600, 1), (1701, 2400, -1)]


class UnirefSummary(unittest.TestCase):
    def test_tax(self):
        self.assertEqual(u.tax("UniRef90_A n=1 Tax=Malassezia sp. CBS 1 TaxID=99 RepID=A"), ("Malassezia sp. CBS 1", "99"))
        self.assertEqual(u.tax("no taxon"), ("", ""))

    def test_member_and_contig(self):
        m = dict(family="7", protein="M", contig_id="1", contig_length="3000", start="1600", end="1001", strand="1")
        r = u.member_row(m, "Skin", "Malassezia", HSPS, GENES)
        self.assertEqual((r["member_hit"], r["pident"], r["genus"], r["same_genus"], r["taxid"]),
                         ("yes", "97.0", "Malassezia", "yes", "76773"))
        self.assertEqual((r["genes"], r["genes_hit"], r["contig_genus"], r["contig_genus_genes"],
                          r["contig_expected_genes"]), (3, 3, "Malassezia", 2, 2))
        none = u.member_row(dict(m, protein="N", contig_id="2"), "Skin", "Malassezia", [], GENES)  # contig 2: no HSP
        self.assertEqual((none["member_hit"], none["contig_genus"], none["genes_hit"]), ("no", "", 0))
        line = u.group_line("7", [r, none]).split("\t")
        self.assertEqual(line[1:4], ["2", "1", "1"])
        self.assertEqual(line[4:6], ["1/1", "0/0"])  # 97% in the >=95 bin, none in >=90
        self.assertEqual(line[-2], "1/2")  # contig 1 majority Malassezia, contig 2 none


if __name__ == "__main__":
    unittest.main()
