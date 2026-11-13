"""Figure 5 (SF_1556 structure): the ESMFold model of the 1.1_clan_3863 representative, MGYP009160527382.

a) model coloured by region, b) the same view by pLDDT, c) the regions to scale on the 635-aa protein.
Regions: the representative's rows in the 1.1 full MSAs of MGYF0000004542 and MGYF0000026748 (classify_multi_hits.py rows) and
DUF285 from Pfam 38.2 (GA). PyMOL renders a and b (render_structure.py, system Python); matplotlib draws c.

    python paper/scripts/figure5.py [--out PREFIX]
    python paper/scripts/figure5.py --selftest
"""
import argparse
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE))
from classify_multi_hits import ends  # noqa: E402

DATA = HERE.parent / "data"
PDB = DATA / "esmfold_MGYP009160527382.pdb"
FAMILIES = {"4542": "MGYF0000004542", "26748": "MGYF0000026748"}
MODEL_LENGTH = {"4542": 106, "26748": 113}   # LENG in mgnifams_hmm.lib.gz (model release 1.0, the HMMs of 1.1)
MGNIFAM_COLOR, PFAM_COLOR, PARTNER_COLOR = "#D55E00", "#0072B2", "#009E73"   # Wong, as Figures 4 and B1
PLDDT_LEGEND = [("#0053D6", "Very high (pLDDT > 90)"), ("#65CBF3", "Confident (90 > pLDDT > 70)"),
                ("#FFDB13", "Low (70 > pLDDT > 50)"), ("#FF7D45", "Very low (pLDDT < 50)")]
WIDTH_MM, TOP_MM, BOTTOM_MM = 180, 82, 42   # a and b with their legends; c
RENDER_PX = (1400, 1100)
# v_x0 of the 2026-11-12 view search: cmd.orient, no turn. Replace with structure_stats.json "view" to freeze another.
VIEW = None
SYSTEM_PYTHON = ["env", "PATH=/usr/bin:/bin", "/usr/bin/python3", "-I"]


def read_regions():
    """{family: [(from, to, n_end, c_end, copy, hmm_from, hmm_to)]} from the classify_multi_hits.py rows of the representative."""
    out = {}
    with open(DATA / "mgnifam_regions.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            out.setdefault(r["family"], []).append((int(r["seq_from"]), int(r["seq_to"]), r["n_end"], r["c_end"], int(r["copy"]),
                                                    int(r["hmm_from"]), int(r["hmm_to"])))
    return out


def read_domtbl(path):
    """[(name, acc, env_from, env_to, hmm_from, hmm_to, model_length)] in protein order."""
    rows = []
    for line in open(path):
        if line.startswith("#"):
            continue
        c = line.split()
        rows.append((c[3], c[4].split(".")[0], int(c[19]), int(c[20]), int(c[15]), int(c[16]), int(c[5])))
    return sorted(rows, key=lambda r: r[2])


def copy_unions(regions):
    """Per copy number, the span covered by that copy in any family: [(from, to)]."""
    spans = {}
    for rows in regions.values():
        for f, t, _, _, k, *_ in rows:
            a, b = spans.get(k, (f, t))
            spans[k] = (min(a, f), max(b, t))
    return [spans[k] for k in sorted(spans)]


def lanes(boxes):
    """Lane index per box so that overlapping boxes (sorted by start) never share a lane."""
    last, out = [], []
    for f, t in boxes:
        for i, e in enumerate(last):
            if f > e:
                last[i] = t
                out.append(i)
                break
        else:
            last.append(t)
            out.append(len(last) - 1)
    return out


def crop(img):
    """The render without its white margin (any channel below 250 counts as ink)."""
    ink = (img[..., :3] < 250 / 255).any(axis=2)
    rows, cols = ink.any(axis=1).nonzero()[0], ink.any(axis=0).nonzero()[0]
    return img[rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1]


def render(segments, contact_of, tmp):
    spec = {"pdb": str(PDB), "outdir": tmp, "width": RENDER_PX[0], "height": RENDER_PX[1], "view": VIEW,
            "segments": segments, "contact_of": contact_of, "contact_cutoff": 8.0, "partner_color": PARTNER_COLOR}
    Path(tmp, "spec.json").write_text(json.dumps(spec))
    subprocess.run(SYSTEM_PYTHON + [str(HERE / "render_structure.py"), str(Path(tmp, "spec.json"))],
                   check=True, stdout=subprocess.DEVNULL)
    return json.loads(Path(tmp, "structure_stats.json").read_text())


def draw(out, regions, pfam, stats, tmp, length):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.image as mimg
    import matplotlib.pyplot as plt
    from matplotlib.patches import PathPatch, Rectangle
    from figure_architectures import FONT_STACK, INK, LINE, MUTED, PREVIEW_DPI, TEXT_SIZE_PT, outline

    plt.rcParams.update({"font.family": "serif", "font.serif": FONT_STACK, "font.size": TEXT_SIZE_PT,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    H = TOP_MM + BOTTOM_MM
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, H / 25.4))
    ax = fig.add_axes([0, 0, 1, 1], xlim=(0, WIDTH_MM), ylim=(0, H))
    ax.set_axis_off()

    # a, b: the two renders side by side (same crop, so the views stay identical), legends under each
    legend_top = BOTTOM_MM + 16
    imgs = [mimg.imread(Path(tmp, f"{n}.png")) for n in ("panel_a", "panel_b")]
    img = crop(imgs[0])
    box_w, box_h = WIDTH_MM / 2 - 6, H - legend_top - 6
    pw = min(box_w, box_h * img.shape[1] / img.shape[0])
    ph = pw * img.shape[0] / img.shape[1]
    ink = (imgs[0][..., :3] < 250 / 255).any(axis=2)
    r0, c0 = ink.any(axis=1).nonzero()[0][0], ink.any(axis=0).nonzero()[0][0]
    for i, full in enumerate(imgs):
        x0 = 4 + i * (WIDTH_MM / 2) + (box_w - pw) / 2
        part = full[r0:r0 + img.shape[0], c0:c0 + img.shape[1]]
        ax.imshow(part, extent=(x0, x0 + pw, legend_top + 2, legend_top + 2 + ph), zorder=1)
        ax.text(2 + i * (WIDTH_MM / 2), H - 2, "ab"[i], fontweight="bold", fontsize=TEXT_SIZE_PT + 2, va="top")
    partner = ", ".join(f"{f}–{t}" for f, t in stats["partner"])
    spans = copy_unions(regions)
    leg_a = [(MGNIFAM_COLOR, "Family region (" + ", ".join(f"{f}–{t}" for f, t in spans) + ")"),
             (PFAM_COLOR, "DUF285 (PF03382)"), (PARTNER_COLOR, f"Packs against the first family region ({partner})"),
             ("#CCCCCC", "No family region or Pfam domain")]
    for col, entries in ((0, leg_a), (1, PLDDT_LEGEND)):
        for j, (color, text) in enumerate(entries):
            cx, cy = 6 + col * (WIDTH_MM / 2), legend_top - j * 3.4
            ax.add_patch(Rectangle((cx, cy - 1.0), 3.0, 2.0, facecolor=color, edgecolor="none"))
            ax.text(cx + 4.2, cy, text, va="center")

    # c: regions to scale
    ax.text(2, BOTTOM_MM - 1, "c", fontweight="bold", fontsize=TEXT_SIZE_PT + 2, va="top")
    left, right = 34, WIDTH_MM - 6
    ax.text(left, BOTTOM_MM - 1.6, f"MGYP009160527382 ({length} aa)", va="top", fontweight="bold")
    x = lambda r: left + (r - 1) * (right - left) / length  # noqa: E731
    # boxes carry their model span inside; row labels give the model length
    rows = [(FAMILIES[k], f"model 1–{MODEL_LENGTH[k]}", [(f, t, n, c, hf, ht) for f, t, n, c, _, hf, ht in v], MGNIFAM_COLOR)
            for k, v in sorted(regions.items(), key=lambda kv: int(kv[0]))]
    pf = [(f, t, *ends(hf, ht, ml), hf, ht) for _, _, f, t, hf, ht, ml in pfam]
    rows.append((f"{pfam[0][0]} (Pfam 38.2)", f"model 1–{pfam[0][6]}", pf, PFAM_COLOR))
    box_h, pitch, y = 3.6, 4.4, BOTTOM_MM - 9
    for label, sub, boxes, color in rows:
        ln = lanes([(f, t) for f, t, *_ in boxes])
        depth = max(ln) + 1
        ly = y - (depth - 1) * pitch / 2
        ax.text(left - 2, ly + 0.9, label, ha="right", va="center")
        ax.text(left - 2, ly - 2.0, sub, ha="right", va="center", color=MUTED, fontsize=TEXT_SIZE_PT - 1)
        ax.plot([x(1), x(length + 1)], [y, y], color=LINE, lw=0.6, zorder=0)
        for (f, t, n, c, hf, ht), lane in zip(boxes, ln):
            yy = y - lane * pitch
            ax.add_patch(PathPatch(outline(x(f), x(t + 1), yy - box_h / 2, yy + box_h / 2, n, c), facecolor=color, edgecolor="none"))
            ax.text((x(f) + x(t + 1)) / 2, yy, f"{hf}–{ht}", ha="center", va="center", color="white",
                    fontsize=TEXT_SIZE_PT - 1)
        y -= pitch * depth + 2.4
    for tick in [1] + list(range(100, length, 100)) + [length]:
        ax.plot([x(tick), x(tick)], [y + 1.6, y + 2.4], color=INK, lw=0.5)
        ax.text(x(tick), y + 1.2, str(tick), ha="center", va="top", color=MUTED)
    ax.plot([x(1), x(length)], [y + 2.4, y + 2.4], color=INK, lw=0.5)
    for ext in ("pdf", "png"):
        fig.savefig(f"{out}.{ext}", dpi=PREVIEW_DPI)
    print(f"figure -> {out}.pdf, .png")


def selftest():
    assert lanes([(1, 10), (5, 20), (21, 30)]) == [0, 1, 0]
    assert lanes([(349, 469), (415, 509), (576, 635)]) == [0, 1, 0]
    assert copy_unions({"a": [(2, 136, "", "", 1), (277, 341, "", "", 2)], "b": [(1, 152, "", "", 1), (274, 347, "", "", 2)]}) == [(1, 152), (274, 347)]
    subprocess.run(SYSTEM_PYTHON + [str(HERE / "render_structure.py"), "--selftest"], check=True)
    print("selftest ok")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="figures/figure-5_sf1556-structure")
    a = p.parse_args()
    seq = "".join(l.strip() for l in open(DATA / "MGYP009160527382.fasta") if not l.startswith(">"))
    regions = read_regions()
    pfam = read_domtbl(DATA / "pfam38.2_sf1556.domtbl")
    spans = copy_unions(regions)
    segments = [[f"family region {i + 1}", f, t, MGNIFAM_COLOR] for i, (f, t) in enumerate(spans)]
    segments += [[f"{n} {f}-{t}", f, t, PFAM_COLOR] for n, _, f, t, _, _, _ in pfam]
    with tempfile.TemporaryDirectory() as tmp:
        stats = render(segments, list(spans[0]), tmp)
        assert stats["length"] == len(seq), (stats["length"], len(seq))
        print(f"MGYP009160527382 {len(seq)} aa, mean pLDDT {stats['mean_plddt']}")
        for k, v in sorted(regions.items(), key=lambda kv: int(kv[0])):
            print(FAMILIES[k], f"(model {MODEL_LENGTH[k]})", "; ".join(f"{f}-{t} model {hf}-{ht} ({n}/{c})" for f, t, n, c, _, hf, ht in v))
        for n, acc, f, t, hf, ht, ml in pfam:
            print(n, acc, f"{f}-{t}", f"model {hf}-{ht} of {ml}", "/".join(ends(hf, ht, ml)))
        print("mean pLDDT per segment:", stats["segments"])
        print("packs against family region 1:", stats["partner"])
        for n, r in stats["contacts"].items():
            print(f"CA within 8 A of {n}:", ", ".join(f"{f}-{t}" for f, t in r))
        draw(a.out, regions, pfam, stats, tmp, len(seq))


if __name__ == "__main__":
    selftest() if sys.argv[1:] == ["--selftest"] else main()
