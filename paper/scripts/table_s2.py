"""Table S2: the biome-exclusive MGnifams with their best UniProtKB hit (stdlib only)."""
import argparse
import csv
import gzip
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
HEADER = [
    "MGnifam", "Biome (end of its only lineage)", "Member proteins with a biome",
    "Full-alignment regions", "Representative length (aa)", "pLDDT", "pTM",
    "Pfam, FunFam, PDB or AlphaFoldDB match", "UniProtKB proteins hit",
    "Best UniProtKB hit", "Best hit description", "Best hit organism", "Best hit E-value",
]


def best_hits(path):
    """Family -> (proteins, E, accession, description, organism); lowest E across SwissProt / TrEMBL rows."""
    best = {}
    for fam, _subset, proteins, e, target, desc in csv.reader(open(path), delimiter="\t"):
        e = float(e)
        if fam in best and best[fam][1] <= e:
            continue
        organism = re.search(r" OS=(.*?) OX=", desc)
        best[fam] = (int(proteins), e, target.split("|")[1], desc.split(" OS=")[0],
                     organism.group(1) if organism else "")
    return best


def build(exclusive, hits, catalogue):
    cat = {r["id"]: r for r in csv.DictReader(gzip.open(catalogue, "rt"))}
    best = best_hits(hits)
    rows = []
    for fam, biome, members in csv.reader(open(exclusive), delimiter="\t"):
        c = cat[fam]  # KeyError = the exclusive set and the catalogue are from different releases
        h = best.get(fam)
        rows.append([fam, biome, members, c["full_size"], c["rep_length"], c["plddt"], c["ptm"],
                     "yes" if c["any_match"] == "1" else "no",
                     *([h[0], h[2], h[3], h[4], "0" if h[1] == 0 else f"{h[1]:.1e}"] if h else [0, "", "", "", ""])])
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exclusive", default=HERE / "data/biome_exclusive_1.1.tsv")
    p.add_argument("--hits", default=HERE / "data/uniprot_hits_1.1.tsv")
    p.add_argument("--catalogue", default=HERE.parent / "data/catalogue_1.1.csv.gz")
    p.add_argument("--out", default="supplementary_tables/table-s2_biome-exclusive-families.tsv")
    a = p.parse_args()
    rows = build(a.exclusive, a.hits, a.catalogue)
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(HEADER)
        w.writerows(rows)
    print(f"{len(rows)} families, {sum(1 for r in rows if r[8])} with a UniProtKB hit, "
          f"{sum(1 for r in rows if r[7] == 'yes')} with a match -> {a.out}")


if __name__ == "__main__":
    main()
