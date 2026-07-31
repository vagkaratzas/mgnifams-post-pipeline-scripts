import importlib.util
import re
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin" / "make_contig_figure.py"
spec = importlib.util.spec_from_file_location("make_contig_figure", BIN)
mcf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mcf)

MAPS = """# contig 19137674
[1-456(-) 874089832 | cluster:874089832 | 152aa <<ANCHOR]

# contig 88320196
[79-486(-) 914919837 | cluster:969294708 | 136aa] -- [5032-5583(+) 914271416 | \
cluster:914271416 | 184aa <<ANCHOR] -- [5634-7595(+) 949579749 | PF13472 | 654aa] -- \
[7626-9374(+) 919921044 | PF01183,PF25309 | 583aa]
"""


def test_parse_maps():
    blocks = list(mcf.parse_maps(MAPS))
    assert [c for c, _ in blocks] == ["19137674", "88320196"]

    genes = blocks[1][1]
    assert len(genes) == 4
    assert genes[0] == dict(start=79, end=486, strand=-1, mgyp="914919837", aa=136,
                            pfams=None, anchor=False)
    assert genes[1]["anchor"] and genes[1]["strand"] == 1
    assert genes[3]["pfams"] == "PF01183 · PF25309"
    assert not genes[3]["anchor"]


def test_colour_genes():
    genes = list(mcf.parse_maps(MAPS))[1][1]
    legend = mcf.colour_genes(genes)

    assert genes[1]["fill"] == mcf.FILL[0]                    # anchor
    assert genes[0]["fill"] == mcf.UNASSIGNED_FILL            # cluster: only
    assert genes[2]["fill"] not in (mcf.FILL[0], mcf.UNASSIGNED_FILL)
    assert genes[2]["fill"] != genes[3]["fill"]               # distinct Pfam sets
    assert [name for _, name in legend] == [
        "unassigned", "MGnifam anchor", "PF13472", "PF01183 · PF25309"]


def test_place_staggers_then_drops():
    rows = [0.0, 0.0]
    assert mcf.place(rows, 50, 100) == 0
    assert mcf.place(rows, 60, 100) == 1     # overlaps row 0
    assert mcf.place(rows, 70, 100) is None  # both rows taken
    assert mcf.place(rows, 400, 100) == 0    # clear of everything again


def test_tick_step():
    assert mcf.tick_step(500) == 100
    assert mcf.tick_step(14000) == 2000
    assert mcf.tick_step(1e9) == 200000


def long_contig(n, anchor_at):
    """n genes, each with its own Pfam, as parse_maps would have produced them."""
    return [dict(start=i * 1000 + 1, end=i * 1000 + 301, strand=1, aa=100,
                 mgyp=str(900000000 + i), pfams="PF%05d" % i, anchor=(i == anchor_at))
            for i in range(n)]


def test_crop_to_anchor():
    genes = long_contig(800, 400)
    kept = mcf.crop_to_anchor(genes, 10)
    assert len(kept) == 21
    assert [g["mgyp"] for g in kept] == [str(900000000 + i) for i in range(390, 411)]

    assert mcf.crop_to_anchor(genes, 0) is genes             # 0 = whole contig
    short = long_contig(5, 2)
    assert mcf.crop_to_anchor(short, 10) is short            # already inside the window
    assert len(mcf.crop_to_anchor(long_contig(800, 3), 10)) == 14   # anchor near the edge


def test_canvas_stays_inside_cairo_size_limit():
    """cairo refuses a surface over 32767 px; the legend used to blow past it."""
    for flank in (10, 0):
        svg = mcf.render("X", mcf.crop_to_anchor(long_contig(800, 400), flank))
        height = int(re.search(r'height="(\d+)"', svg).group(1))
        assert height * mcf.DPI / 96 < 32767
    assert "further annotations" in svg      # flank=0 kept every gene, legend spilled


def test_render_is_svg_covering_every_gene():
    contig, genes = list(mcf.parse_maps(MAPS))[1]
    svg = mcf.render(contig, genes)

    assert svg.startswith("<svg") and svg.endswith("</svg>")
    assert "Contig 88320196" in svg
    assert "▲ ANCHOR" in svg
    for g in genes:
        assert f'MGYP{int(g["mgyp"]):012d}' in svg
