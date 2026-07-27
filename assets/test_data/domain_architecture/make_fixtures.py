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
]

# PF09999 is deliberately absent so the mapping fallback is exercised.
PFAMS = [
    ("PF00001", "Alpha domain"),
    ("PF00002", "Beta domain"),
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
]


def main():
    with (HERE / "clan_membership_dummy.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Cluster Id", "Family Rep Id", "Family Size", "Family Ids"])
        writer.writerows(CLANS)

    with (HERE / "pfam_mapping_dummy.tsv").open("w") as handle:
        for accession, name in PFAMS:
            handle.write(f"{accession}\t{name}\n")

    with gzip.open(HERE / "proteins_dummy.csv.gz", "wt", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["mgyp", "sequence", "full_length", "cluster_size", "metadata"])
        for mgyp, metadata, _ in ROWS:
            writer.writerow([mgyp, SEQ, "false", 1, json.dumps(metadata, separators=(",", ":"))])

    print(f"wrote {len(CLANS)} clans, {len(PFAMS)} pfams, {len(ROWS)} proteins to {HERE}")


if __name__ == "__main__":
    main()
