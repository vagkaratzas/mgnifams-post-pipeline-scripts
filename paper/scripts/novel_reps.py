"""Novel families of a release and their NMPFamsDB hits.

  set      families with no match (catalogue any_match 0), in no clan, and no Pfam domain on any member
           protein (pfam_per_family with_pfam 0; only updated families have that column, so the others
           are left out and counted). Writes family, rep_id, rep_region, rep_length.
  reps     FASTA of each representative region, sequence from the MGnify Proteins API; its length must
           equal the catalogue rep_length, or the run stops.
  summary  best hit per family from an hmmsearch --domtblout (NMPFamsDB HMMs as queries, reps as targets);
           a hit counts if one of its domains has i-Evalue < 0.001.
"""
import argparse
import csv
import gzip
import json
import sys
import time
import urllib.error
import urllib.request

API = "https://www.ebi.ac.uk/metagenomics/proteins/api/v1/protein/"


def rows(path, delim):
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as f:
        yield from csv.DictReader(f, delimiter=delim)


def novel(catalogue, clans, pfam, families):
    """(selected rows, no-match no-clan families without Pfam data)."""
    in_clan = {r["family_id"] for r in rows(clans, ",")}
    with_pfam = {r["family"]: int(r["with_pfam"]) for r in rows(pfam, "\t")}
    length = {r["id"]: r["rep_length"] for r in rows(catalogue, ",") if r["any_match"] == "0"}
    left = {f for f in length if f not in in_clan}
    no_data = sorted(left - with_pfam.keys(), key=int)
    keep = {f for f in left if with_pfam.get(f) == 0}
    out = [dict(family=r["family_id"], rep_id=r["rep_id"], rep_region=r["rep_region"], rep_length=length[r["family_id"]])
           for r in rows(families, "\t") if r["family_id"] in keep]
    assert len(out) == len(keep), "families.tsv misses a selected family"
    return sorted(out, key=lambda r: int(r["family"])), no_data


def cut(seq, region):
    """1-based inclusive region; '-' = the whole protein."""
    if region == "-":
        return seq
    a, b = map(int, region.split("-"))
    return seq[a - 1:b]


def fetch(mgyp, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(API + mgyp, timeout=60) as r:
                return json.load(r)["sequence"]
        except (urllib.error.URLError, TimeoutError) as e:
            err = e
            time.sleep(5 * (i + 1))
    raise RuntimeError(f"{mgyp}: {err}")


def best_hits(domtbl, dom_e=0.001):
    """{family: (nmpfam, seq E, seq bits, hmm from, hmm to, qlen)}: the top full-sequence score among hits
    with at least one domain of i-Evalue < dom_e; HMM span of that hit's best-scoring domain."""
    doms = {}
    for line in open(domtbl):
        if not line.startswith("#"):
            c = line.split()
            doms.setdefault((c[0], c[3]), []).append(c)
    best = {}
    for (fam, nmp), rows in doms.items():
        if not any(float(c[12]) < dom_e for c in rows):
            continue
        d = max(rows, key=lambda c: float(c[13]))
        hit = (nmp, float(rows[0][6]), float(rows[0][7]), int(d[15]), int(d[16]), int(d[5]))
        if fam not in best or hit[2] > best[fam][2]:
            best[fam] = hit
    return best


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("set")
    s.add_argument("--catalogue", default="paper/data/catalogue_1.1.csv.gz")
    s.add_argument("--clans", default="tmp/clans_1.1/1.1_clans.csv")
    s.add_argument("--pfam", default="tmp/pfam_per_family_out/pfam_per_family_25590863.tsv.gz")
    s.add_argument("--families", default="tmp/families_1.1.tsv.gz")
    s.add_argument("--out", default="paper/data/novel_1.1.tsv")
    r = sub.add_parser("reps")
    r.add_argument("--set", default="paper/data/novel_1.1.tsv")
    r.add_argument("--out", default="paper/data/novel_reps_1.1.faa")
    m = sub.add_parser("summary")
    m.add_argument("domtbl")
    m.add_argument("--set", default="paper/data/novel_1.1.tsv")
    t = sub.add_parser("table")
    t.add_argument("--set", default="paper/data/novel_1.1.tsv")
    t.add_argument("--domtbl", default="paper/data/novel_1.1_nmpfamsdb2.domtbl")
    t.add_argument("--catalogue", default="paper/data/catalogue_1.1.csv.gz")
    t.add_argument("--out", default="supplementary_tables/table-s3_novel-families.tsv")
    a = p.parse_args()

    if a.cmd == "table":
        cat = {r["id"]: r for r in rows(a.catalogue, ",")}
        best = best_hits(a.domtbl)
        with open(a.out, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["MGnifam", "Representative", "Representative region", "Representative length (aa)",
                        "Full-alignment regions", "pLDDT", "Best NMPFamsDB2 hit", "E-value", "Bit score",
                        "NMPFamsDB2 HMM positions matched (best domain)"])
            for row in rows(a.set, "\t"):
                c, h = cat[row["family"]], best.get(row["family"])
                hit = [h[0], f"{h[1]:.2g}", h[2], f"{h[3]}-{h[4]} of {h[5]}"] if h else ["", "", "", ""]
                w.writerow([row["family"], row["rep_id"], row["rep_region"], row["rep_length"], c["full_size"], c["plddt"], *hit])
        print(f"table -> {a.out}; {len(best)} with a hit")
        return

    if a.cmd == "set":
        out, no_data = novel(a.catalogue, a.clans, a.pfam, a.families)
        with open(a.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t")
            w.writeheader()
            w.writerows(out)
        print(f"{len(out)} novel families -> {a.out}; {len(no_data)} with no match and no clan lack Pfam data (not updated)")
    elif a.cmd == "reps":
        with open(a.out, "w") as f:
            for row in rows(a.set, "\t"):
                seq = cut(fetch(row["rep_id"]), row["rep_region"])
                if len(seq) != int(row["rep_length"]):
                    sys.exit(f"{row['family']}: region {row['rep_region']} gives {len(seq)} aa, catalogue says {row['rep_length']}")
                f.write(f">{row['family']} {row['rep_id']}/{row['rep_region']}\n{seq}\n")
        print(f"reps -> {a.out}")
    else:
        fams = [row["family"] for row in rows(a.set, "\t")]
        best = best_hits(a.domtbl)
        bits = sorted(b[2] for b in best.values())
        print(f"{len(best)} of {len(fams)} families with a hit; best-hit bits: "
              f"<20 {sum(b < 20 for b in bits)}, 20-30 {sum(20 <= b < 30 for b in bits)}, >=30 {sum(b >= 30 for b in bits)}")
        print("family\tnmpfam\tseq_E\tbits\thmm_from\thmm_to\thmm_len")
        for fam, h in sorted(best.items(), key=lambda kv: -kv[1][2]):
            print(fam, *h, sep="\t")


def selftest():
    assert cut("ABCDEF", "2-4") == "BCD" and cut("ABC", "-") == "ABC"
    import os, tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".domtbl", delete=False) as f:
        f.write("# c\n7 - 110 F1 - 139 0.0008 12.4 0 1 1 6e-6 0.001 12.1 0 76 121 28 72 23 89 0.84 -\n"
                "7 - 110 F2 - 50 1e-9 40.0 0 1 1 1e-9 1e-9 40 0 1 50 1 50 1 50 0.9 -\n"
                "8 - 110 F3 - 99 1e-20 70.0 0 1 2 1e-3 1.0 7.0 0 1 10 1 10 1 12 0.8 -\n"
                "8 - 110 F3 - 99 1e-20 70.0 0 2 2 1e-25 1e-20 60.0 0 20 90 30 100 25 105 0.9 -\n"
                "8 - 110 F4 - 99 1e-22 80.0 0 1 2 1e-5 0.01 30.0 0 1 40 1 40 1 40 0.8 -\n"
                "8 - 110 F4 - 99 1e-22 80.0 0 2 2 1e-4 0.5 25.0 0 50 90 50 90 50 90 0.8 -\n")
    # F1 / F4: no domain with i-Evalue < 0.001; F3: span of its best domain, not the first row
    assert best_hits(f.name) == {"7": ("F2", 1e-9, 40.0, 1, 50, 50), "8": ("F3", 1e-20, 70.0, 20, 90, 99)}
    os.unlink(f.name)
    print("selftest ok")


if __name__ == "__main__":
    selftest() if sys.argv[1:] == ["--selftest"] else main()
