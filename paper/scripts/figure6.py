#!/usr/bin/env python3
"""Figure 6: the genes next to the MGnifam 243 representative, with their domains.

Genes are drawn to scale on the contig and point in the coding direction. Domains are drawn inside their gene at
their protein positions (envelope), with Pfam's ends: curved when the model is matched to within 5% of its end,
jagged otherwise (classify_multi_hits.ends). Box shapes come from Figure 4.

Inputs (data/): the contig (MGnify API), gene calls (contig_lookup.py), Pfam 38.2 domtbl (hmmsearch
--cut_ga of the translated genes), the MGnifam region (full MSA), Pfam names. Each gene is translated from the
contig; its length must equal the domtbl's, or the run stops.
Usage (from the repository root):
  python paper/scripts/figure6.py --out figures/figure-6_mgnifam-243-gene-neighbourhood
  --translate FAA writes the translated genes (the hmmsearch input); --selftest checks the coordinate maths.
"""
import argparse
import csv
import itertools
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE))
from classify_multi_hits import ends  # noqa: E402

DATA = HERE.parent / "data"
WIDTH_MM, HEIGHT_MM, MARGIN_MM = 180, 40, 14   # margin: room for the end genes' accessions
GENE_H_MM, HEAD_MM, DOMAIN_INSET_MM = 5.0, 2.2, 0.5
Y_GENE = 31                      # gene centre line
LEGEND_Y, LEGEND_ROW_MM, LEGEND_COL1_MM = 15, 3.6, 62
COLORS = {"MGYF0000000243": "#D55E00", "PF13472": "#009E73", "PF25309": "#0072B2", "PF01183": "#CC79A7"}
GENE_FILL, LINE, INK, MUTED = "#E6E6E6", "#888888", "#000000", "#555555"

_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {"".join(c): _AA[i] for i, c in enumerate(itertools.product(_BASES, repeat=3))}


def tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def translate(nt):
    """Protein without its stop codon; an internal stop or a missing one stops the run."""
    protein = "".join(CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))
    if not protein.endswith("*") or "*" in protein[:-1]:
        sys.exit(f"bad ORF: {protein}")
    return protein[:-1]


def gene_protein(contig, start, end, strand):
    nt = contig[start - 1:end]
    return translate(nt[::-1].translate(str.maketrans("ACGT", "TGCA")) if strand < 0 else nt)


def nt_span(start, end, strand, aa_from, aa_to):
    """Contig positions of protein residues aa_from..aa_to (1-based, inclusive)."""
    if strand > 0:
        return start + 3 * (aa_from - 1), start + 3 * aa_to - 1
    return end - 3 * aa_to + 1, end - 3 * (aa_from - 1)


def read_domtbl(path):
    """{protein: [(accession, env_from, env_to, hmm_from, hmm_to, model_length)]} and {protein: length}."""
    domains, lengths = {}, {}
    for line in open(path):
        if line.startswith("#"):
            continue
        c = line.split()
        lengths[c[0]] = int(c[2])
        domains.setdefault(c[0], []).append((c[4].split(".")[0], int(c[19]), int(c[20]), int(c[15]), int(c[16]), int(c[5])))
    return domains, lengths


def load():
    contig_lines = open(DATA / "contig_243.fasta").read().split("\n")
    header, contig = contig_lines[0], "".join(contig_lines[1:])
    genes = [dict(g, start=int(g["start"]), end=int(g["end"]), strand=int(g["strand"])) for g in tsv(DATA / "genes_243.tsv")]
    domains, lengths = read_domtbl(DATA / "pfam38.2_243.domtbl")
    for m in tsv(DATA / "mgnifam_243.tsv"):
        domains.setdefault(m["protein"], []).append(
            (m["accession"], int(m["seq_from"]), int(m["seq_to"]), int(m["hmm_from"]), int(m["hmm_to"]), int(m["model_length"])))
    for g in genes:
        g["protein_seq"] = gene_protein(contig, g["start"], g["end"], g["strand"])
        g["length"] = len(g["protein_seq"])
        if g["protein"] in lengths and lengths[g["protein"]] != g["length"]:
            sys.exit(f"{g['protein']}: translated {g['length']} aa, domtbl says {lengths[g['protein']]}")
        g["domains"] = sorted(domains.get(g["protein"], []), key=lambda d: d[1])
    names = {r["accession"]: r for r in tsv(DATA / "pfam38.2_names.tsv")}
    return header, len(contig), genes, names


def gene_path(x0, x1, y, strand):
    from matplotlib.path import Path as MPath
    h, head = GENE_H_MM / 2, min(HEAD_MM, 0.45 * (x1 - x0))
    if strand > 0:
        pts = [(x0, y - h), (x1 - head, y - h), (x1, y), (x1 - head, y + h), (x0, y + h)]
    else:
        pts = [(x1, y - h), (x0 + head, y - h), (x0, y), (x0 + head, y + h), (x1, y + h)]
    return MPath(pts + [pts[0]], closed=True)


def draw(genes, names, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import PathPatch, Rectangle
    from figure_architectures import FONT_STACK, PREVIEW_DPI, TEXT_SIZE_PT, outline

    plt.rcParams.update({"font.family": "serif", "font.serif": FONT_STACK, "font.size": TEXT_SIZE_PT,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    lo, hi = min(g["start"] for g in genes) - 60, max(g["end"] for g in genes) + 60
    scale = (WIDTH_MM - 2 * MARGIN_MM) / (hi - lo)
    x = lambda bp: MARGIN_MM + (bp - lo) * scale  # noqa: E731

    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
    ax = fig.add_axes([0, 0, 1, 1], xlim=(0, WIDTH_MM), ylim=(0, HEIGHT_MM))
    ax.set_axis_off()
    ax.plot([x(lo), x(hi)], [Y_GENE, Y_GENE], color=LINE, lw=0.6, zorder=0)
    seen = []
    for g in genes:
        x0, x1 = x(g["start"]), x(g["end"] + 1)
        arrow = PathPatch(gene_path(x0, x1, Y_GENE, g["strand"]), facecolor=GENE_FILL, edgecolor=LINE, lw=0.5, zorder=1)
        ax.add_patch(arrow)
        for acc, a, b, hf, ht, ml in g["domains"]:
            n_end, c_end = ends(hf, ht, ml)
            left, right = (n_end, c_end) if g["strand"] > 0 else (c_end, n_end)
            s, e = nt_span(g["start"], g["end"], g["strand"], a, b)
            h = GENE_H_MM / 2 - DOMAIN_INSET_MM
            box = PathPatch(outline(x(s), x(e + 1), Y_GENE - h, Y_GENE + h, left, right),
                            facecolor=COLORS[acc], edgecolor="none", zorder=2)
            ax.add_patch(box)
            box.set_clip_path(arrow)
            label = names[acc]["name"] if acc in names else acc
            ax.text((x(s) + x(e + 1)) / 2, Y_GENE + GENE_H_MM / 2 + 1.2, label, ha="center", va="bottom")
            if acc not in seen:
                seen.append(acc)
        mid = (x0 + x1) / 2
        ax.text(mid, Y_GENE - GENE_H_MM / 2 - 1.2, f"{g['length']} aa", ha="center", va="top", color=MUTED)
        ax.text(mid, Y_GENE - GENE_H_MM / 2 - 4.6, g["protein"], ha="center", va="top", color=MUTED, fontsize=TEXT_SIZE_PT - 1)

    # 1 kb scale bar, bottom right
    bar_y = LEGEND_Y - 2 * LEGEND_ROW_MM
    ax.plot([x(hi) - 1000 * scale, x(hi)], [bar_y, bar_y], color=INK, lw=0.8)
    ax.text(x(hi) - 500 * scale, bar_y - 0.9, "1 kb", ha="center", va="top")

    # legend: MGnifam first, then Pfam in gene order, then grey
    entries = [(COLORS[a], a if a not in names else f"{names[a]['name']} ({a}), {names[a]['description']}") for a in seen]
    entries.append((GENE_FILL, "No Pfam domain"))
    for i, (color, text) in enumerate(entries):
        cx, cy = 4 + (i % 2) * LEGEND_COL1_MM, LEGEND_Y - (i // 2) * LEGEND_ROW_MM
        ax.add_patch(Rectangle((cx, cy - 1.1), 3.2, 2.2, facecolor=color, edgecolor=LINE if color == GENE_FILL else "none", lw=0.5))
        ax.text(cx + 4.4, cy, text, va="center")
    for ext in ("pdf", "png"):
        fig.savefig(f"{out}.{ext}", dpi=PREVIEW_DPI)
    print(f"figure -> {out}.pdf, .png")


def selftest():
    assert nt_span(100, 108, 1, 1, 3) == (100, 108) and nt_span(100, 108, 1, 2, 2) == (103, 105)
    assert nt_span(100, 108, -1, 1, 1) == (106, 108) and nt_span(100, 108, -1, 3, 3) == (100, 102)
    assert gene_protein("ATGAAATAA", 1, 9, 1) == "MK" and gene_protein("TTATTTCAT", 1, 9, -1) == "MK"
    print("selftest ok")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="figures/figure-6_mgnifam-243-gene-neighbourhood")
    p.add_argument("--translate", help="write the translated genes to this FASTA and stop")
    a = p.parse_args()
    header, contig_length, genes, names = load()
    if a.translate:
        with open(a.translate, "w") as f:
            f.writelines(f">{g['protein']}\n{g['protein_seq']}\n" for g in genes)
        print(f"translated -> {a.translate}")
        return
    print(header, f"{contig_length} bp")
    for g in genes:
        doms = ", ".join(f"{d[0]} {d[1]}-{d[2]} (model {d[3]}-{d[4]} of {d[5]}, {'/'.join(ends(*d[3:]))})" for d in g["domains"])
        print(g["protein"], g["start"], g["end"], g["strand"], f"{g['length']} aa", doms or "-", sep="\t")
    draw(genes, names, a.out)


if __name__ == "__main__":
    selftest() if sys.argv[1:] == ["--selftest"] else main()
