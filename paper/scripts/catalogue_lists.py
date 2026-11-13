"""Table 2 annotation-category lists from the Figure 3 catalogue (paper/data/catalogue_1.1.csv.gz).

Usage: python -I catalogue_lists.py CATALOGUE.csv.gz OUTDIR
Writes OUTDIR/novel.txt (no match to any resource: any_match == 0) and
OUTDIR/novel_structure_any.txt (no Pfam sequence / profile or FunFam match: family_match == 0),
one MGYF%010d per line, the names and form of the UniProtKB coverage pipeline's --lists_dir.
On the 1.0 catalogue they must give 7,266 and 17,883 (Table 2) and equal the published run's lists.
"""
import csv
import gzip
import os
import sys


def lists(rows):
    novel = sorted(f"MGYF{int(r['id']):010d}" for r in rows if r["any_match"] == "0")
    no_family = sorted(f"MGYF{int(r['id']):010d}" for r in rows if r["family_match"] == "0")
    return {"novel": novel, "novel_structure_any": no_family}


def main(path, outdir):
    with gzip.open(path, "rt") as fh:
        out = lists(list(csv.DictReader(fh)))
    os.makedirs(outdir, exist_ok=True)
    for name, ids in out.items():
        with open(os.path.join(outdir, f"{name}.txt"), "w") as fh:
            fh.write("".join(f"{i}\n" for i in ids))
        print(f"{outdir}/{name}.txt: {len(ids)}")


if __name__ == "__main__":
    if sys.argv[1:] == ["--selftest"]:
        r = lists([{"id": "7", "any_match": "0", "family_match": "0"}, {"id": "8", "any_match": "1", "family_match": "0"}])
        assert r == {"novel": ["MGYF0000000007"], "novel_structure_any": ["MGYF0000000007", "MGYF0000000008"]}
        print("selftest ok")
    else:
        main(*sys.argv[1:])
