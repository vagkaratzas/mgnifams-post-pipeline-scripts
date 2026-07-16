#!/usr/bin/env python3
"""
Publication-ready gene-neighbourhood (contig) figures.

Reads the ``*_maps.txt`` written by ``synteny_census.py`` and emits one 600-dpi
PNG per contig::

    # contig 88320196
    [79-486(-) 914919837 | cluster:969294708 | 136aa] -- [5032-5583(+) ... <<ANCHOR]

Usage::

    python bin/make_contig_figure.py --input input/partners_maps.txt --output figures/
"""

import argparse
import io
import re
from pathlib import Path

import cairosvg
from PIL import Image, ImageChops

PAD = 2                           # px of white kept around the cropped content
DPI = 600

GENE_RE = re.compile(
    r"\[(\d+)-(\d+)\(([+-?])\)\s+(\S+)\s+\|\s+([^|]+?)\s+\|\s+(\d+)aa(\s+<<\w+)?\]"
)

# Okabe-Ito derived, colour-vision-deficiency safe. First entry is reserved for
# the anchor, the rest are cycled over the Pfam labels of a contig.
FILL = ["#E69F00", "#009E73", "#0072B2", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442"]
EDGE = ["#9A6A00", "#00654A", "#004C77", "#8C4E72", "#2E7BA6", "#8C3C00", "#A99C1E"]
UNASSIGNED_FILL, UNASSIGNED_EDGE = "#C9CDD2", "#8A9099"

FONT = "Helvetica Neue, Helvetica, Arial, DejaVu Sans, sans-serif"
INK, MUTED = "#1A1D21", "#6B7280"

W = 1180
X0, X1 = 74, 1116                 # plot area in px
Y_ARROW, AH = 150, 32             # arrow vertical centre / height
Y_ID, Y_AA, Y_COORD = 108, 121, 182
Y_AXIS, Y_LEG = 244, 300


def parse_maps(text):
    """Yield (contig_id, [gene, ...]) for every '# contig <id>' block."""
    contig = None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("# contig"):
            contig = line.split()[-1]
        elif line.startswith("[") and contig:
            genes = []
            for start, end, strand, mgyp, annot, aa, mark in GENE_RE.findall(line):
                annot = annot.strip()
                genes.append(dict(
                    start=int(start), end=int(end),
                    strand={"+": 1, "-": -1}.get(strand, 1),
                    mgyp=mgyp, aa=int(aa),
                    pfams=None if annot.startswith("cluster:") else annot.replace(",", " · "),
                    anchor=mark.strip() == "<<ANCHOR",
                ))
            yield contig, genes
            contig = None


def colour_genes(genes):
    """Assign (fill, edge) per gene and build the legend for one contig."""
    legend, palette = [], {}
    for g in genes:
        if g["anchor"]:
            key, name = "anchor", g["pfams"] or "MGnifam anchor"
        elif g["pfams"]:
            key, name = g["pfams"], g["pfams"]
        else:
            key, name = "none", "unassigned"
        if key not in palette:
            if key == "anchor":
                palette[key] = (FILL[0], EDGE[0])
            elif key == "none":
                palette[key] = (UNASSIGNED_FILL, UNASSIGNED_EDGE)
            else:
                i = 1 + len([k for k in palette if k not in ("anchor", "none")]) % (len(FILL) - 1)
                palette[key] = (FILL[i], EDGE[i])
            legend.append((palette[key], name))
        g["fill"], g["edge"] = palette[key]
    return legend


def place(rows, xm, width):
    """Pick a stagger row whose last label ends before this one starts.

    Returns the row index, or None when every row is taken — dense contigs have
    more genes than legible labels, so the crowded ones are dropped.
    """
    for i, right in enumerate(rows):
        if xm - width / 2 >= right:
            rows[i] = xm + width / 2 + 4
            return i
    return None


def arrow(xa, xb, strand, yc=Y_ARROW, h=AH):
    """Pentagon gene arrow; xa<xb are the genomic bounds in px."""
    t, b = yc - h / 2, yc + h / 2
    head = min(15.0, (xb - xa) * 0.42)
    if strand > 0:
        p = [(xa, t), (xb - head, t), (xb, yc), (xb - head, b), (xa, b)]
    else:
        p = [(xb, t), (xa + head, t), (xa, yc), (xa + head, b), (xb, b)]
    return " ".join(f"{px:.2f},{py:.2f}" for px, py in p)


def tick_step(span):
    """A round tick step giving roughly 6 major ticks over `span` bp."""
    for step in (100, 250, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000):
        if span / step <= 8:
            return step
    return 200000


def render(contig, genes):
    lo, hi = min(g["start"] for g in genes), max(g["end"] for g in genes)
    pad = max((hi - lo) * 0.03, 60)
    bp_lo, bp_hi = lo - pad, hi + pad
    sc = (X1 - X0) / (bp_hi - bp_lo)
    x = lambda bp: X0 + (bp - bp_lo) * sc

    legend = colour_genes(genes)
    H = Y_LEG + ((len(legend) + 1) // 2) * 21 + 46
    rows = [0.0, 0.0]                 # right edge of the last label in each stagger row
    # the anchor claims a row before the greedy left-to-right pass, so a crowded
    # neighbourhood can never be the reason its label goes missing
    for g in genes:
        if g["anchor"]:
            g["row"] = place(rows, (x(g["start"]) + x(g["end"])) / 2, 118)

    s = []
    a = s.append
    a(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
      f'viewBox="0 0 {W} {H}" font-family="{FONT}">')
    a(f'<rect width="{W}" height="{H}" fill="#FFFFFF"/>')

    # --- title ---------------------------------------------------------------
    anchor = next((g for g in genes if g["anchor"]), None)
    a(f'<text x="{X0}" y="34" font-size="16" font-weight="700" fill="{INK}">'
      f'Contig {contig}</text>')
    who = (f'<tspan font-weight="600" fill="{INK}">MGYP{int(anchor["mgyp"]):012d}</tspan>'
           if anchor else "—")
    a(f'<text x="{X0}" y="54" font-size="11.5" fill="{MUTED}">'
      f'Gene neighbourhood of the MGnifams anchor {who}'
      f' — {len(genes)} predicted CDS, {lo:,}–{hi:,} bp</text>')
    a(f'<line x1="{X0}" y1="68" x2="{X1}" y2="68" stroke="#E3E6EA" stroke-width="1"/>')

    # --- backbone ------------------------------------------------------------
    a(f'<line x1="{X0}" y1="{Y_ARROW}" x2="{X1}" y2="{Y_ARROW}" '
      f'stroke="#9AA1A9" stroke-width="1.6"/>')

    # --- genes ---------------------------------------------------------------
    for g in genes:
        xa, xb = x(g["start"]), x(g["end"])
        if g["anchor"]:  # soft highlight band behind the anchor gene
            a(f'<rect x="{xa-7:.1f}" y="76" width="{xb-xa+14:.1f}" height="147" rx="5" '
              f'fill="#E69F00" fill-opacity="0.10"/>')

        a(f'<polygon points="{arrow(xa, xb, g["strand"])}" fill="{g["fill"]}" '
          f'stroke="{g["edge"]}" stroke-width="{1.6 if g["anchor"] else 1.0}" '
          f'stroke-linejoin="round"/>')

        if g["pfams"] and (xb - xa) > 92:
            a(f'<text x="{(xa+xb)/2:.1f}" y="{Y_ARROW+3.6:.1f}" font-size="10.5" '
              f'font-weight="700" fill="#FFFFFF" text-anchor="middle" '
              f'letter-spacing="0.2">{g["pfams"]}</text>')

        xm = (xa + xb) / 2
        # 118px ≈ the widest MGYP label at 9.8px
        row = g["row"] if "row" in g else place(rows, xm, 118)
        if row is not None:
            dy = row * 26
            a(f'<text x="{xm:.1f}" y="{Y_ID-dy}" font-size="9.8" font-weight="600" '
              f'fill="{INK}" text-anchor="middle">MGYP{int(g["mgyp"]):012d}</text>')
            a(f'<text x="{xm:.1f}" y="{Y_AA-dy}" font-size="9" fill="{MUTED}" '
              f'text-anchor="middle">{g["aa"]} aa · '
              f'{"+" if g["strand"] > 0 else "−"} strand</text>')
            a(f'<text x="{xm:.1f}" y="{Y_COORD+row*13}" font-size="8.6" fill="{MUTED}" '
              f'text-anchor="middle">{g["start"]:,}–{g["end"]:,}</text>')

        if g["anchor"]:
            a(f'<text x="{xm:.1f}" y="218" font-size="9.6" '
              f'font-weight="700" fill="#9A6A00" text-anchor="middle" '
              f'letter-spacing="0.6">▲ ANCHOR</text>')

    # --- axis ----------------------------------------------------------------
    a(f'<line x1="{X0}" y1="{Y_AXIS}" x2="{X1}" y2="{Y_AXIS}" '
      f'stroke="{MUTED}" stroke-width="1"/>')
    step = tick_step(bp_hi - bp_lo)
    first = int(bp_lo // step + 1) * step
    for bp in range(first, int(bp_hi) + 1, step):
        a(f'<line x1="{x(bp):.1f}" y1="{Y_AXIS}" x2="{x(bp):.1f}" y2="{Y_AXIS+5}" '
          f'stroke="{MUTED}" stroke-width="1"/>')
        a(f'<text x="{x(bp):.1f}" y="{Y_AXIS+17}" font-size="9" fill="{MUTED}" '
          f'text-anchor="middle">{bp/1000:g} kb</text>')
    for bp in range(first - step, int(bp_hi) + 1, step // 4):
        if bp % step and bp >= bp_lo:
            a(f'<line x1="{x(bp):.1f}" y1="{Y_AXIS}" x2="{x(bp):.1f}" '
              f'y2="{Y_AXIS+2.5}" stroke="#B9BEC5" stroke-width="0.7"/>')
    # scale bar, one tick step wide
    sb0, sb1 = X0, X0 + step * sc
    a(f'<line x1="{sb0:.1f}" y1="{Y_AXIS+34}" x2="{sb1:.1f}" y2="{Y_AXIS+34}" '
      f'stroke="{INK}" stroke-width="2.2"/>')
    a(f'<text x="{(sb0+sb1)/2:.1f}" y="{Y_AXIS+29}" font-size="8.8" fill="{INK}" '
      f'text-anchor="middle" font-weight="600">{step/1000:g} kb</text>')

    # --- legend --------------------------------------------------------------
    cols = [X0 + 160, X0 + 480]
    for i, ((fill, edge), name) in enumerate(legend):
        cx, cy = cols[i % 2], Y_LEG + (i // 2) * 21
        a(f'<polygon points="{arrow(cx, cx+26, 1, cy-3.5, 11)}" fill="{fill}" '
          f'stroke="{edge}" stroke-width="0.9" stroke-linejoin="round"/>')
        a(f'<text x="{cx+33}" y="{cy}" font-size="9.4" fill="{INK}" '
          f'font-weight="600">{name}</text>')

    a('</svg>')
    return "\n".join(s)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, type=Path,
                   help="*_maps.txt written by synteny_census.py")
    p.add_argument("--output", required=True, type=Path,
                   help="output folder for the per-contig PNGs")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=True)
    n = 0
    for contig, genes in parse_maps(args.input.read_text()):
        if not genes:
            continue
        png = cairosvg.svg2png(bytestring=render(contig, genes).encode(),
                               scale=DPI / 96)
        img = Image.open(io.BytesIO(png)).convert("RGB")
        # the canvas is sized for the widest possible contig, so trim it back to
        # what was actually drawn
        x0, t, r, b = ImageChops.invert(img).getbbox()
        img.crop((max(x0 - PAD, 0), max(t - PAD, 0),
                  min(r + PAD, img.width), min(b + PAD, img.height))
                 ).save(args.output / f"contig_{contig}.png", dpi=(DPI, DPI))
        n += 1
    print(f"wrote {n} PNGs to {args.output}")


if __name__ == "__main__":
    main()
