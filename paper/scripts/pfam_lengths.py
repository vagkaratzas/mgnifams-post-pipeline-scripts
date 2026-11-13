"""Model lengths of the Pfam families in the size baseline (Results catalogue paragraph 5).

Usage: python -I pfam_lengths.py <pfam_2026_07_family_sizes.tsv> <Pfam-A.hmm.dat.gz> [floor]

Joins the per-accession domain counts from pfam_recount.sh with the Pfam `#=GF ML`
model lengths, keeps families with at least `floor` domains (default 6, the smallest
1.1 MGnifam), and reports the length range, the median and IQR of those with 75-2,000
match states, and how many each bound excludes.
"""
import gzip
import re
import statistics
import sys


def accession(raw):
    """'PF00005', 'PF00005.30' or '5' -> 'PF00005'."""
    raw = raw.strip().split(".")[0]
    return raw if raw.startswith("PF") else f"PF{int(raw):05d}"


def read_lengths(path):
    lengths, acc = {}, None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#=GF AC"):
                acc = accession(line.split()[2])
            elif line.startswith("#=GF ML"):
                lengths[acc] = int(line.split()[2])
    return lengths


def read_sizes(path):
    sizes = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if not re.fullmatch(r"\d+", fields[1]):
                continue  # header
            sizes[accession(fields[0])] = int(fields[1])
    return sizes


def summarise(sizes, lengths, floor):
    kept = [a for a, n in sizes.items() if n >= floor]
    missing = [a for a in kept if a not in lengths]
    lens = [lengths[a] for a in kept if a in lengths]
    inside = [x for x in lens if 75 <= x <= 2000]
    q1, med, q3 = statistics.quantiles(inside, n=4, method="inclusive")
    return {
        "families": len(kept),
        "missing_length": len(missing),
        "length_range": (min(lens), max(lens)),
        "below_75": sum(x < 75 for x in lens),
        "above_2000": sum(x > 2000 for x in lens),
        "in_75_2000": len(inside),
        "median": med,
        "iqr": (q1, q3),
    }


def selftest():
    assert accession("5") == "PF00005" and accession("PF00528.27") == "PF00528"
    sizes = {"PF00001": 10, "PF00002": 6, "PF00003": 5, "PF00004": 100}
    lengths = {"PF00001": 80, "PF00002": 200, "PF00003": 300, "PF00004": 50}
    r = summarise(sizes, lengths, 6)
    assert r["families"] == 3 and r["below_75"] == 1 and r["in_75_2000"] == 2, r


if __name__ == "__main__":
    selftest()
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    floor = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    for key, value in summarise(read_sizes(sys.argv[1]), read_lengths(sys.argv[2]), floor).items():
        print(f"{key}\t{value}")
