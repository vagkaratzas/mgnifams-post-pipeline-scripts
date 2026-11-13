#!/usr/bin/env python3
"""Organism of the members' genes and contigs from the blastx of the contigs against all of UniRef90 2026_02.

Input: blastx_uniref.sh output (OUTFMT + stitle, gzipped or not); contig_lookup.py members.tsv and genes.tsv.gz;
uniprot_hit_taxa.tsv (the family's expected genus = the genus of its best UniProtKB hit, lowest E, as in Table S2); the biome of each
family (paper/data/biome_exclusive_1.1.tsv: family, biome, n).

A gene's taxon = the UniRef90 cluster's common taxon (stitle Tax=) of the best HSP (bit score) that covers it by
>= 30 bp on its strand (blastx_summary.on_gene). Genus = first word of Tax (two for Candidatus): a Tax above genus
rank (e.g. Basidiomycota) never equals the expected genus.
Per member: the member gene's taxon and identity. Per contig: the most common genus over its genes with a hit.
Per family and per biome: members in the expected genus at each identity bin, and contigs whose majority is it.

Usage: uniref_summary.py blastx_uniref.tsv.gz members.tsv genes.tsv.gz uniprot_hit_taxa.tsv biomes.tsv --out x.tsv
"""

import argparse
import csv
import gzip
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python -I drops the script dir
from blastx_summary import BINS, genus, on_gene, read_hsps  # noqa: E402

TAX = re.compile(r"Tax=(.+?) TaxID=(\d+)")
FIELDS = ("family", "biome", "protein", "contig_id", "contig_length", "expected", "member_hit", "pident", "evalue",
          "uniref", "tax", "taxid", "genus", "same_genus", "genes", "genes_hit", "contig_genus", "contig_genus_genes",
          "contig_expected_genes")


def tax(stitle):
    m = TAX.search(stitle or "")
    return (m.group(1), m.group(2)) if m else ("", "")


def best_on(hsps, g):
    on = [h for h in hsps if on_gene(h, g)]
    return max(on, key=lambda h: h["bitscore"]) if on else None


def contig_view(hsps, genes, expected):
    """Genes on the contig, genes with a hit, majority genus and its gene count, genes in the expected genus."""
    gen = Counter()
    for g in genes:
        b = best_on(hsps, g)
        if b:
            gen[genus(tax(b.get("stitle"))[0])] += 1
    top, n = gen.most_common(1)[0] if gen else ("", 0)
    return len(genes), sum(gen.values()), top, n, gen[expected]


def member_row(m, biome, expected, hsps, genes):
    g = (int(m["start"]), int(m["end"]), int(m["strand"]))
    g = (min(g[:2]), max(g[:2]), g[2])
    row = dict(family=m["family"], biome=biome, protein=m["protein"], contig_id=m["contig_id"],
               contig_length=m["contig_length"], expected=expected, member_hit="no", pident="", evalue="", uniref="",
               tax="", taxid="", genus="", same_genus="")
    b = best_on(hsps, g)
    if b:
        t, tid = tax(b.get("stitle"))
        row.update(member_hit="yes", pident=f"{b['pident']:.1f}", evalue=f"{b['evalue']:.2g}", uniref=b["sseqid"],
                   tax=t, taxid=tid, genus=genus(t), same_genus="yes" if genus(t) == expected else "no")
    v = contig_view(hsps, genes, expected)
    row.update(zip(("genes", "genes_hit", "contig_genus", "contig_genus_genes", "contig_expected_genes"), v))
    return row


def group_line(label, rows):
    """members, with hit, in expected genus; per identity bin in/total; contigs with the expected majority; others."""
    hit = [r for r in rows if r["member_hit"] == "yes"]
    same = [r for r in hit if r["same_genus"] == "yes"]
    bins = defaultdict(lambda: [0, 0])
    for r in hit:
        b = next(b for b in BINS if float(r["pident"]) >= b)
        bins[b][0] += r["same_genus"] == "yes"
        bins[b][1] += 1
    contigs = {r["contig_id"]: r for r in rows}.values()
    major = sum(r["contig_genus"] == r["expected"] for r in contigs)
    other = Counter(r["genus"] for r in hit if r["same_genus"] == "no").most_common(3)
    return "\t".join([label, str(len(rows)), str(len(hit)), str(len(same))]
                     + [f"{bins[b][0]}/{bins[b][1]}" for b in BINS]
                     + [f"{major}/{len(contigs)}", ", ".join(f"{g or '-'} {n}" for g, n in other)])


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for a in ("blastx", "members", "genes", "taxa", "biomes"):
        p.add_argument(a)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)
    hsps = read_hsps(a.blastx)
    genes = defaultdict(list)
    for r in csv.DictReader(gzip.open(a.genes, "rt"), delimiter="\t"):
        s, e = sorted((int(r["start"]), int(r["end"])))
        genes[r["contig_id"]].append((s, e, int(r["strand"])))
    best = {}  # family -> (E, genus); not the most common genus: 12231 has 10 Fusarium hits but its best is Bathycoccus
    for fam, _subset, _target, e, org, _taxid in csv.reader(open(a.taxa), delimiter="\t"):
        best[fam] = min(best.get(fam, (float("inf"), "")), (float(e), genus(org)))
    expected = {f: g for f, (_e, g) in best.items()}
    biome = {r[0]: r[1] for r in csv.reader(open(a.biomes), delimiter="\t")}
    rows = [member_row(m, biome.get(m["family"], ""), expected.get(m["family"], ""), hsps.get(m["contig_id"], []),
                       genes.get(m["contig_id"], []))
            for m in csv.DictReader(open(a.members), delimiter="\t")]
    with open(a.out, "w") as f:
        f.write("\t".join(FIELDS) + "\n")
        f.writelines("\t".join(str(r[k]) for k in FIELDS) + "\n" for r in rows)

    head = ("members\twith_hit\tin_expected\t" + "\t".join(f">={b}" if b else "<50" for b in BINS)
            + "\tcontig_majority\tother_genera")
    print("family\tbiome\texpected\t" + head)
    for fam in sorted({r["family"] for r in rows}, key=int):
        fr = [r for r in rows if r["family"] == fam]
        print(group_line(f"{fam}\t{fr[0]['biome']}\t{fr[0]['expected']}", fr))
    print("\nbiome\t" + head)
    for b in sorted({r["biome"] for r in rows}):
        print(group_line(b, [r for r in rows if r["biome"] == b]))


if __name__ == "__main__":
    main()
