#!/usr/bin/env python3
"""Regenerate the domain-architecture test fixtures.

Each protein row below targets one branch of bin/parse_domain_architectures.py; the comment on
the row names the branch and the architecture it must produce. Run from anywhere:

    python assets/test_data/domain_architecture/make_fixtures.py
"""

import csv
import gzip
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

# Clans. SF_1 and SF_2 are multi-family, SF_3 is a singleton, SF_4 never gets a hit so it
# exercises the missing-family validation path.
CLANS = [
    ("SF_1", "100", "3", "100;101;102"),
    ("SF_2", "200", "3", "200;201;202"),
    ("SF_3", "300", "1", "300"),
    ("SF_4", "400", "1", "400"),
    ("SF_5", "500", "1", "500"),
]

# Three columns: accession, name, clan. PF00003 and PF00004 share CL0003 so they can merge;
# PF00005 has a blank clan so it never can; PF09999 is absent entirely, exercising the fallback.
PFAMS = [
    ("PF00001", "Alpha domain", "CL0001"),
    ("PF00002", "Beta domain", "CL0002"),
    ("PF00003", "Gamma domain", "CL0003"),
    ("PF00004", "Delta domain", "CL0003"),
    ("PF00005", "Epsilon domain", ""),
]

SEQ = "MKVLATTVGSKQWDEHRPLY"

# (mgyp, metadata, expected architecture) - the third element is documentation, not written out.
ROWS = [
    # no "m" at all: must be dropped by the prefilter
    (1, {"p": [["PF00001", 1e-05, 50.0, 1, 60, 10, 70]]}, None),
    # three same-clan hits on the same region -> one clan chip
    (2, {"m": [["100", 1e-05, 50.0, 10, 100],
               ["101", 1e-04, 45.0, 10, 100],
               ["102", 1e-03, 40.0, 10, 100]]}, "SF_1"),
    # two same-clan hits on disjoint regions -> two family chips
    (3, {"m": [["100", 1e-05, 50.0, 10, 100],
               ["101", 1e-04, 45.0, 200, 290]]}, "100\t101"),
    # two different-clan hits that overlap -> both shown, never merged
    (4, {"m": [["100", 1e-05, 50.0, 10, 100],
               ["200", 1e-04, 45.0, 20, 110]]}, "100\t200"),
    # bug 1 regression: hmm_from is 1 but ali_from is 200, so the Pfam must sort AFTER the MGnifam
    (5, {"p": [["PF00001", 1e-05, 50.0, 1, 60, 200, 260]],
         "m": [["300", 1e-05, 50.0, 10, 100]]}, "300\tPF00001"),
    # repeated Pfam, plus an accession missing from the mapping
    (6, {"p": [["PF00002", 1e-05, 50.0, 1, 40, 10, 50],
               ["PF00002", 1e-06, 55.0, 1, 40, 100, 140],
               ["PF09999", 1e-04, 45.0, 1, 30, 200, 230]],
         "m": [["300", 1e-05, 50.0, 300, 380]]}, "PF00002\tPF00002\tPF09999\t300"),
    # exactly 50% overlap (50 of 100) -> NOT merged, because the test is strictly greater
    (7, {"m": [["100", 1e-05, 50.0, 1, 100],
               ["101", 1e-04, 45.0, 51, 150]]}, "100\t101"),
    # 51% overlap (51 of 100) -> merged
    (8, {"m": [["100", 1e-05, 50.0, 1, 100],
               ["101", 1e-04, 45.0, 50, 149]]}, "SF_1"),
    # equal start, different length, different clans -> shorter first
    (9, {"m": [["100", 1e-05, 50.0, 10, 100],
               ["200", 1e-04, 45.0, 10, 50]]}, "200\t100"),
    # equal start, equal length -> alphabetical id ("300" < "PF00001")
    (10, {"p": [["PF00001", 1e-05, 50.0, 1, 60, 10, 100]],
          "m": [["300", 1e-05, 50.0, 10, 100]]}, "300\tPF00001"),
    # single linkage: 200-201 and 201-202 overlap >50%, 200-202 does not -> one group
    (11, {"m": [["200", 1e-05, 50.0, 1, 100],
                ["201", 1e-04, 45.0, 40, 139],
                ["202", 1e-03, 40.0, 80, 179]]}, "SF_2"),
    # same family twice on disjoint regions -> two chips, but credited once
    (12, {"m": [["300", 1e-05, 50.0, 10, 100],
                ["300", 1e-04, 45.0, 200, 290]]}, "300\t300"),

    # Pfam clan collapsing. All of these sit on family 500 alone, so every other family's JSON
    # stays byte-identical to the pre-clan output.
    # two same-clan Pfams overlapping 48/52 -> one clan chip
    (13, {"p": [["PF00003", 1e-05, 50.0, 1, 60, 1, 52],
                ["PF00004", 1e-04, 45.0, 1, 60, 5, 82]],
          "m": [["500", 1e-05, 50.0, 300, 380]]}, "CL0003\t500"),
    # two different-clan Pfams on the same region -> never merged
    (14, {"p": [["PF00001", 1e-05, 50.0, 1, 60, 1, 52],
                ["PF00003", 1e-04, 45.0, 1, 60, 1, 52]],
          "m": [["500", 1e-05, 50.0, 300, 380]]}, "PF00001\tPF00003\t500"),
    # two same-clan Pfams on disjoint regions -> two chips, repeat preserved
    (15, {"p": [["PF00003", 1e-05, 50.0, 1, 60, 1, 52],
                ["PF00004", 1e-04, 45.0, 1, 60, 200, 260]],
          "m": [["500", 1e-05, 50.0, 300, 380]]}, "PF00003\tPF00004\t500"),
    # a blank-clan Pfam overlapping a clanned one -> never merged, and no warning
    (16, {"p": [["PF00003", 1e-05, 50.0, 1, 60, 1, 52],
                ["PF00005", 1e-04, 45.0, 1, 60, 1, 52]],
          "m": [["500", 1e-05, 50.0, 300, 380]]}, "PF00003\tPF00005\t500"),
    # exactly 50% Pfam overlap -> NOT merged
    (17, {"p": [["PF00003", 1e-05, 50.0, 1, 60, 1, 100],
                ["PF00004", 1e-04, 45.0, 1, 60, 51, 150]],
          "m": [["500", 1e-05, 50.0, 300, 380]]}, "PF00003\tPF00004\t500"),
    # 51% Pfam overlap -> merged
    (18, {"p": [["PF00003", 1e-05, 50.0, 1, 60, 1, 100],
                ["PF00004", 1e-04, 45.0, 1, 60, 50, 149]],
          "m": [["500", 1e-05, 50.0, 300, 380]]}, "CL0003\t500"),
]


def main():
    with (HERE / "clan_membership_dummy.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Cluster Id", "Family Rep Id", "Family Size", "Family Ids"])
        writer.writerows(CLANS)

    with (HERE / "pfam_mapping_dummy.tsv").open("w") as handle:
        handle.write("pfam_id\tname\tclan_id\n")
        for accession, name, clan in PFAMS:
            handle.write(f"{accession}\t{name}\t{clan}\n")

    with gzip.open(HERE / "proteins_dummy.csv.gz", "wt", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["mgyp", "sequence", "full_length", "cluster_size", "metadata"])
        for mgyp, metadata, _ in ROWS:
            writer.writerow([mgyp, SEQ, "false", 1, json.dumps(metadata, separators=(",", ":"))])

    print(f"wrote {len(CLANS)} clans, {len(PFAMS)} pfams, {len(ROWS)} proteins to {HERE}")


if __name__ == "__main__":
    main()
