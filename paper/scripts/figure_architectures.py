#!/usr/bin/env python3
"""Domain-architecture figure of one MGnifam (Figure 4 for MGYF0000000025).

A: the most common architectures, schematic N- to C-terminal (not to scale), until 95% of the member proteins or
15 rows, then one row for all the others; per row the proteins, % of the family and % with every MGnifam copy
full length. MGnifam boxes have curved or jagged ends (Pfam's convention, within 5% of the model end) by the
majority of the row's copies; Pfam boxes are plain. B: cumulative % of proteins against architecture rank.

Input: architectures_<job>.tsv.gz from architectures.py. Style as figure3.py: Times 8 pt, 180 mm, Wong palette,
PDF (deliverable) + 600 dpi PNG.
Usage (from the repository root):
  python paper/scripts/figure_architectures.py \
      --rows paper/data/architectures_1.1.tsv.gz --family 25 \
      --out figures/figure-4_mgnifam-25-architectures
  --selftest checks the row cut and the run compression only.
"""
import argparse
import csv
import gzip
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import PathPatch  # noqa: E402
from matplotlib.ticker import LogLocator  # noqa: E402
from matplotlib.path import Path as MPath  # noqa: E402

TEXT_SIZE_PT = 8
FONT_STACK = ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"]
WIDTH_MM = 180
ARCH_MM = 128               # panel A width, of which the numbers take NUMBERS_MM
NUMBERS_MM = 36
ROW_MM = 6.0                # row pitch
BOX_H_MM = 3.6
BOX_W_MM = 13.0            # widest box; narrower when the longest drawn row needs it
BOX_GAP_MM = 1.2
TOOTH_MM = 0.9              # depth of a jagged end
COVER, MAX_ROWS, MAX_TOKENS = 0.95, 15, 7
PREVIEW_DPI = 600
WONG = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#F0E442", "#000000"]
MGNIFAM_COLOR = "#D55E00"   # vermillion: the MGnifam, as the dark red of the 1.0 figure
OTHER_PFAM = "#BBBBBB"      # Pfams past the seventh share one grey (labelled by accession in the box)
INK, MUTED, LINE = "#000000", "#555555", "#888888"


def read_rows(path, family):
    with gzip.open(path, "rt") as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t") if r["family"] == family]
    if not rows:
        sys.exit(f"family {family} not in {path}")
    return rows


def cut(rows):
    """Rows drawn: until COVER of the proteins or MAX_ROWS."""
    shown = []
    for r in rows:
        shown.append(r)
        if float(r["cum_pct"]) >= 100 * COVER or len(shown) == MAX_ROWS:
            break
    return shown


def compress(tokens, ends):
    """[(token, run length, ends of the run's first and last copy)]; a run is merged only when the row is long."""
    items = []
    mgnifam = iter(ends)
    for t in tokens:
        e = next(mgnifam) if not t.startswith("PF") else None
        if len(tokens) > MAX_TOKENS and items and items[-1][0] == t:
            items[-1] = (t, items[-1][1] + 1, (items[-1][2][0], e) if e else None)
        else:
            items.append((t, 1, (e, e) if e else None))
    return items


def n_boxes(tokens):
    """Boxes compress() draws for a row."""
    return len(tokens) if len(tokens) <= MAX_TOKENS else 1 + sum(a != b for a, b in zip(tokens, tokens[1:]))


def outline(x0, x1, y0, y1, left, right):
    """Closed box path within [x0, x1]; left / right in {'curved', 'jagged', 'flat'}. Ends are drawn inside the
    box, so neighbouring boxes never overlap and the gap between them stays visible."""
    h = y1 - y0
    inset = {"curved": h / 2, "jagged": TOOTH_MM, "flat": 0}
    a, b = x0 + inset[left], x1 - inset[right]
    pts = [(a, y0), (b, y0)]
    if right == "curved":
        pts += [(b + h / 2 * (1 - (2 * i / 8 - 1) ** 2) ** 0.5, y0 + h * i / 8) for i in range(1, 8)]
    elif right == "jagged":
        pts += [(b + (TOOTH_MM if i % 2 else 0), y0 + h * i / 4) for i in range(1, 4)]
    pts += [(b, y1), (a, y1)]
    if left == "curved":
        pts += [(a - h / 2 * (1 - (2 * i / 8 - 1) ** 2) ** 0.5, y1 - h * i / 8) for i in range(1, 8)]
    elif left == "jagged":
        pts += [(a - (TOOTH_MM if i % 2 else 0), y1 - h * i / 4) for i in range(1, 4)]
    return MPath(pts + [(a, y0)], closed=True)


def end_style(counts, n):
    """('curved' | 'jagged', ...) for the N and C end of one MGnifam token, by majority of its n copies."""
    full, n_curved, c_curved = map(int, counts.split("/"))
    return ("curved" if n_curved * 2 > n else "jagged", "curved" if c_curved * 2 > n else "jagged")


def pct(value):
    return "<0.1" if 0 < value < 0.05 else f"{value:.1f}"


def draw(rows, family, out):
    plt.rcParams.update({"font.family": "serif", "font.serif": FONT_STACK, "font.size": TEXT_SIZE_PT,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    shown = cut(rows)
    total = sum(int(r["proteins"]) for r in rows)
    rest = rows[len(shown):]
    colors = {}
    for r in shown:
        for t in r["architecture"].split():
            if t.startswith("PF") and t not in colors:
                colors[t] = WONG[len(colors)] if len(colors) < len(WONG) else OTHER_PFAM
    names = {t: n for r in shown for t, n in zip(r["architecture"].split(), r["names"].split()) if t.startswith("PF")}
    n_rows = len(shown) + bool(rest)
    legend_rows = (len(names) + 2) // 3
    height_mm = max(60, 9 + ROW_MM * n_rows + 4 + 4.2 * legend_rows)  # panel B needs ~40 mm
    mm = 1 / 25.4
    fig = plt.figure(figsize=(WIDTH_MM * mm, height_mm * mm))
    ax = fig.add_axes([0, 0, ARCH_MM / WIDTH_MM, 1])
    ax.set_xlim(0, ARCH_MM)
    ax.set_ylim(height_mm, 0)
    ax.axis("off")
    ax.text(0, 2.5, "A", fontweight="bold", va="center")
    x_num = ARCH_MM - NUMBERS_MM
    longest = max(n_boxes(r["architecture"].split()) for r in shown)
    box_w = min(BOX_W_MM, (x_num - 6) / longest - BOX_GAP_MM)
    columns = (x_num + 11, x_num + 23, x_num + 35)
    for x, label in zip(columns, ("Proteins", "% family", "% full")):
        ax.text(x, 6, label, ha="right", va="center", color=MUTED)
    for i, r in enumerate(shown):
        y = 9 + ROW_MM * (i + 0.5)
        n = int(r["proteins"])
        tokens = r["architecture"].split()
        items = compress(tokens, [end_style(e, n) for e in r["mgnifam_ends"].split(",") if e])
        x = 3
        ax.plot([1.5, x + len(items) * (box_w + BOX_GAP_MM)], [y, y], color=LINE, lw=0.6, zorder=1)
        for token, run, e in items:
            pfam = token.startswith("PF")
            left, right = ("flat", "flat") if pfam else (e[0][0], e[1][1])
            ax.add_patch(PathPatch(outline(x, x + box_w, y - BOX_H_MM / 2, y + BOX_H_MM / 2, left, right),
                                   facecolor=colors.get(token, MGNIFAM_COLOR), edgecolor="white", lw=0.4, zorder=2))
            fill = colors.get(token, MGNIFAM_COLOR)
            ink = INK if fill in ("#E69F00", "#56B4E9", "#F0E442", OTHER_PFAM, "#CC79A7", "#009E73") else "white"
            ax.text(x + box_w / 2, y, token + (f" ×{run}" if run > 1 else ""), ha="center", va="center",
                    color=ink, zorder=3)
            x += box_w + BOX_GAP_MM
        for xr, value in zip(columns, (f"{n:,}", pct(float(r["pct"])), pct(float(r["all_full_pct"])))):
            ax.text(xr, y, value, ha="right", va="center")
    if rest:
        y = 9 + ROW_MM * (len(shown) + 0.5)
        n = sum(int(r["proteins"]) for r in rest)
        ax.text(3, y, f"{len(rest):,} other architectures", va="center", color=MUTED)
        ax.text(columns[0], y, f"{n:,}", ha="right", va="center")
        ax.text(columns[1], y, pct(100 * n / total), ha="right", va="center")
    y = 9 + ROW_MM * n_rows + 4
    # the MGnifam first: its boxes carry the short number, the legend its MGYF accession
    legend = {family: f"MGYF{int(family):010d}", **names}
    for i, (acc, name) in enumerate(legend.items()):
        cx, cy = 3 + (i % 3) * 38, y + 4.2 * (i // 3)
        ax.add_patch(plt.Rectangle((cx, cy - 1.2), 4, 2.4, facecolor=colors.get(acc, MGNIFAM_COLOR), edgecolor="none"))
        ax.text(cx + 5.5, cy, f"{acc} {name}", va="center")

    bx = fig.add_axes([(ARCH_MM + 12) / WIDTH_MM, 13 / height_mm,
                       (WIDTH_MM - ARCH_MM - 15) / WIDTH_MM, 1 - 22 / height_mm])
    ranks = list(range(1, len(rows) + 1))
    bx.plot(ranks, [float(r["cum_pct"]) for r in rows], color=MGNIFAM_COLOR, lw=1.0)
    at = next(i for i, r in enumerate(rows, 1) if float(r["cum_pct"]) >= 100 * COVER)
    bx.axhline(100 * COVER, color=LINE, lw=0.5, ls=(0, (3, 2)))
    bx.plot([at], [float(rows[at - 1]["cum_pct"])], "o", ms=3, color=MGNIFAM_COLOR)
    # above the plot area, where it never meets the curve
    bx.text(1, 1.02, f"95% at rank {at}", transform=bx.transAxes, ha="right", va="bottom")
    bx.set_xscale("log")
    bx.xaxis.set_major_locator(LogLocator(base=10, numticks=12))
    bx.set_xlim(0.8, max(len(rows), 10) * 1.3)
    bx.set_ylim(min(50, float(rows[0]["cum_pct"]) - 5), 101)
    bx.set_xlabel("Architecture rank")
    bx.set_ylabel("Cumulative % of proteins")
    for side in ("top", "right"):
        bx.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        bx.spines[side].set_linewidth(0.5)
    bx.tick_params(width=0.5, length=2)
    bx.grid(axis="y", lw=0.3, color="#DDDDDD")
    fig.text((ARCH_MM + 2) / WIDTH_MM, 1 - 2.5 / height_mm, "B", fontweight="bold", va="center")
    for ext in ("pdf", "png"):
        fig.savefig(f"{out}.{ext}", dpi=PREVIEW_DPI)
    print(f"{family}: {total:,} proteins, {len(rows):,} architectures, drawn {len(shown)} rows "
          f"({shown[-1]['cum_pct'][:6]}%), others {len(rest):,}, 95% at rank {at}")
    for r in shown:
        print(f"  {r['rank']}\t{r['proteins']}\t{r['pct'][:6]}\t{r['all_full_pct'][:6]}\t{r['names']}\t{r['mgnifam_ends']}")


def selftest():
    rows = [dict(cum_pct=str(c)) for c in (50, 80, 95.5, 99, 100)]
    assert len(cut(rows)) == 3
    assert len(cut([dict(cum_pct="1")] * 20)) == MAX_ROWS
    ends = [("jagged", "curved")] * 9
    items = compress(["PF1"] + ["7"] * 9, ends)
    assert items == [("PF1", 1, None), ("7", 9, (("jagged", "curved"), ("jagged", "curved")))], items
    assert compress(["7", "7"], [("curved", "curved"), ("jagged", "jagged")])[1] == ("7", 1, (("jagged", "jagged"),) * 2)
    assert n_boxes(["PF1"] + ["7"] * 9) == 2 and n_boxes(["7", "7"]) == 2
    for left in ("curved", "jagged", "flat"):
        for right in ("curved", "jagged", "flat"):
            xs = [x for x, _ in outline(10, 20, 0, 3.6, left, right).vertices]
            assert 10 - 1e-9 <= min(xs) and max(xs) <= 20 + 1e-9, (left, right)  # ends stay inside the box
    assert end_style("1/3/1", 4) == ("curved", "jagged")
    print("selftest ok")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--rows")
    parser.add_argument("--family")
    parser.add_argument("--out")
    parser.add_argument("--selftest", action="store_true")
    a = parser.parse_args()
    if a.selftest:
        return selftest()
    draw(read_rows(a.rows, a.family), a.family, a.out)


if __name__ == "__main__":
    main()
