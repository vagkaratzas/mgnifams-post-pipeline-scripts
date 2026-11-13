#!/usr/bin/env python3
"""Per member and per family: how the family's UniProtKB hit proteins align to the members' source contigs.

Input: DIAMOND blastx of the contigs (query) against the hit proteins (subject), --outfmt 6 with the columns of
OUTFMT; contig_lookup.py members.tsv and genes.tsv.gz; uniprot_hit_taxa.tsv (family, subset, target, E, organism,
taxid). Only HSPs to the member's own family's hit proteins count.

Per member: the best HSP (bit score) that overlaps the member gene: identity, hit protein, organism. Then the HSPs
of that same hit protein anywhere on the contig: frames, gene calls on the same strand they cover by >= 30 bp, and the hit protein's residues
covered. Several gene calls in different frames under one hit protein = a gene split by the gene caller (exons).
Per family: identity distribution of the member-gene HSPs.

Usage: blastx_summary.py blastx.tsv members.tsv genes.tsv.gz uniprot_hit_taxa.tsv --out per_member.tsv
       [--draw CONTIG_ID ...]   (ASCII picture of a contig: gene calls and HSPs)
"""

import argparse
import csv
import gzip
import statistics
from collections import Counter, defaultdict

OUTFMT = "qseqid sseqid pident length qstart qend sstart send evalue bitscore qframe slen qlen".split()
BINS = (95, 90, 80, 70, 50, 0)
MIN_GENE_BP = 30  # a gene is under the hit only if an HSP covers 30 bp of it (adjacent bacterial genes overlap by 1-4 bp)
FIELDS = ("family", "protein", "contig_id", "contig_length", "member_hsp", "pident", "evalue", "hit", "organism",
          "genus", "frames", "genes_under_hit", "hit_coverage")


def acc(sseqid):
    return sseqid.split("|")[1] if "|" in sseqid else sseqid


def genus(organism):
    w = organism.split()
    return " ".join(w[:2]) if w and w[0] == "Candidatus" else (w[0] if w else "")


def read_hsps(path):
    hsps = defaultdict(list)  # contig_id -> HSPs
    for row in csv.reader((gzip.open if path.endswith(".gz") else open)(path, "rt"), delimiter="\t"):
        h = dict(zip(OUTFMT + ["stitle"], row))  # stitle: the UniRef90 run only
        for k in ("length", "qstart", "qend", "sstart", "send", "qframe", "slen", "qlen"):
            h[k] = int(h[k])
        h["pident"], h["evalue"], h["bitscore"] = float(h["pident"]), float(h["evalue"]), float(h["bitscore"])
        h["acc"] = acc(h["sseqid"])
        h["lo"], h["hi"] = sorted((h["qstart"], h["qend"]))
        hsps[h["qseqid"]].append(h)
    return hsps


def overlaps(a_lo, a_hi, b_lo, b_hi):
    return a_lo <= b_hi and b_lo <= a_hi


def on_gene(h, g):
    """HSP h covers gene g = (start, end, strand, ...) by >= MIN_GENE_BP on the same strand."""
    return min(h["hi"], g[1]) - max(h["lo"], g[0]) + 1 >= MIN_GENE_BP and (h["qframe"] > 0) == (g[2] == 1)


def covered(intervals):
    total, end = 0, 0
    for lo, hi in sorted(intervals):
        if hi > end:
            total += hi - max(lo, end + 1) + 1
            end = hi
    return total


def summarise(member, hsps, genes, fam_hits, organism):
    """One output row for a member (dict from members.tsv)."""
    lo, hi = sorted((int(member["start"]), int(member["end"])))
    mine = [h for h in hsps.get(member["contig_id"], []) if h["acc"] in fam_hits]
    on_member = [h for h in mine if overlaps(h["lo"], h["hi"], lo, hi)]
    row = dict(family=member["family"], protein=member["protein"], contig_id=member["contig_id"],
               contig_length=member["contig_length"], member_hsp="no", pident="", evalue="", hit="", organism="",
               genus="", frames="", genes_under_hit="", hit_coverage="")
    if not on_member:
        return row
    best = max(on_member, key=lambda h: h["bitscore"])
    same = [h for h in mine if h["acc"] == best["acc"]]
    under = {(g[0], g[1]) for g in genes.get(member["contig_id"], [])
             if any(on_gene(h, g) for h in same)}
    cov = covered([sorted((h["sstart"], h["send"])) for h in same]) / best["slen"]
    org = organism.get(best["acc"], "")
    row.update(member_hsp="yes", pident=f"{best['pident']:.1f}", evalue=f"{best['evalue']:.2g}", hit=best["acc"],
               organism=org, genus=genus(org), frames=len({h["qframe"] for h in same}), genes_under_hit=len(under),
               hit_coverage=f"{cov:.3f}")
    return row


def draw(contig, member, hsps, genes, hit, width=100):
    """ASCII: gene calls (member = #) and the HSPs of the member's best hit protein, scaled to the contig."""
    n = int(member["contig_length"])
    scale = lambda x: min(width - 1, (x - 1) * width // n)  # noqa: E731
    lines = [f"contig {contig} ({n:,} bp), member {member['protein']} (family {member['family']})"]
    gl = [" "] * width
    for s, e, strand, prot in genes.get(contig, []):
        ch = "#" if prot == member["protein"] else (">" if strand == 1 else "<")
        for x in range(scale(s), scale(e) + 1):
            gl[x] = ch
    lines.append("genes  |" + "".join(gl) + "|")
    for h in sorted((h for h in hsps.get(contig, []) if h["acc"] == hit), key=lambda h: h["lo"]):
        hl = [" "] * width
        for x in range(scale(h["lo"]), scale(h["hi"]) + 1):
            hl[x] = "="
        lines.append(f"f{h['qframe']:+d}    |" + "".join(hl) + f"| {h['acc']} {h['pident']:.0f}% "
                     f"aa {min(h['sstart'], h['send'])}-{max(h['sstart'], h['send'])}/{h['slen']} E {h['evalue']:.1g}")
    return "\n".join(lines)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("blastx")
    p.add_argument("members")
    p.add_argument("genes")
    p.add_argument("taxa")
    p.add_argument("--out", required=True)
    p.add_argument("--draw", nargs="*", default=[], help="contig ids to draw")
    a = p.parse_args(argv)
    hsps = read_hsps(a.blastx)
    genes = defaultdict(list)
    for r in csv.DictReader(gzip.open(a.genes, "rt"), delimiter="\t"):
        genes[r["contig_id"]].append((int(r["start"]), int(r["end"]), int(r["strand"]), r["protein"]))
    fam_hits, organism = defaultdict(set), {}
    for fam, _subset, target, _e, org, _taxid in csv.reader(open(a.taxa), delimiter="\t"):
        fam_hits[fam].add(acc(target))
        organism[acc(target)] = org
    members = list(csv.DictReader(open(a.members), delimiter="\t"))
    rows = [summarise(m, hsps, genes, fam_hits[m["family"]], organism) for m in members]
    with open(a.out, "w") as f:
        f.write("\t".join(FIELDS) + "\n")
        f.writelines("\t".join(str(r[k]) for k in FIELDS) + "\n" for r in rows)

    print("family\tmembers\twith_hsp\tmedian_pident\t" + "\t".join(f">={b}" if b else "<50" for b in BINS)
          + "\tgenera\tmulti_gene\tmulti_frame")
    for fam in sorted({r["family"] for r in rows}, key=int):
        fr = [r for r in rows if r["family"] == fam]
        hit = [r for r in fr if r["member_hsp"] == "yes"]
        ids = [float(r["pident"]) for r in hit]
        bins = Counter(next(b for b in BINS if x >= b) for x in ids)
        genera = Counter(r["genus"] for r in hit).most_common(3)
        print(f"{fam}\t{len(fr)}\t{len(hit)}\t{statistics.median(ids) if ids else 'NA'}\t"
              + "\t".join(str(bins[b]) for b in BINS)
              + f"\t{', '.join(f'{g} {n}' for g, n in genera)}\t{sum(r['genes_under_hit'] > 1 for r in hit)}"
              f"\t{sum(r['frames'] > 1 for r in hit)}")
    by_contig = {r["contig_id"]: (m, r) for m, r in zip(members, rows)}
    for c in a.draw:
        m, r = by_contig[c]
        print(f"\n{r['family']} {r['protein']}: best {r['hit']} {r['organism']} {r['pident']}%")
        print(draw(c, m, hsps, genes, r["hit"]))


if __name__ == "__main__":
    main()
