"""PyMOL half of the SF_1556 structure figure: renders panels a (regions) and b (pLDDT) and measures contacts.

Runs on the system Python that carries PyMOL, not conda (see SKILL.md); figure5.py calls it:
    env PATH=/usr/bin:/bin /usr/bin/python3 -I render_structure.py SPEC.json
SPEC: {"pdb", "outdir", "width", "height", "view" (18 floats or null), "segments": [[name, from, to, hex]],
       "contact_of": [from, to], "contact_cutoff": 8.0}
Writes outdir/panel_a.png, outdir/panel_b.png and outdir/structure_stats.json.
"""
import json
import sys

from pymol import cmd, stored

PLDDT_BANDS = [(90, 101, "0x0053D6"), (70, 90, "0x65CBF3"), (50, 70, "0xFFDB13"), (0, 50, "0xFF7D45")]  # AlphaFold DB


def ranges(residues, join=3):
    """Sorted residue numbers -> [[from, to]], merging gaps of up to `join` residues."""
    out = []
    for r in sorted(residues):
        if out and r - out[-1][1] <= join:
            out[-1][1] = r
        else:
            out.append([r, r])
    return out


def main(spec_path):
    spec = json.load(open(spec_path))
    cmd.load(spec["pdb"], "m")
    cmd.hide("everything")
    cmd.show("cartoon")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)

    stored.ca = []
    cmd.iterate("m and name CA", "stored.ca.append((int(resi), b))")
    plddt = dict(stored.ca)

    def near(a, b):
        """CA within the cutoff of a-b, outside it and its 8 flanking residues (chain neighbours are trivially close)."""
        stored.r = set()
        cmd.iterate(f"(m and name CA and not resi {max(1, a - 8)}-{b + 8}) within {spec['contact_cutoff']} of (m and name CA and resi {a}-{b})",
                    "stored.r.add(int(resi))")
        return set(stored.r)

    contacts = {n: ranges(near(f, t)) for n, f, t, _ in spec["segments"] if n.startswith("family region")}
    coloured = {i for _, f, t, _ in spec["segments"] for i in range(f, t + 1)}
    partner = [r for r in ranges(near(*spec["contact_of"]) - coloured) if r[1] - r[0] >= 4]

    cmd.color("0xCCCCCC", "m")
    for _, f, t, hexcol in spec["segments"]:
        cmd.color(hexcol.replace("#", "0x"), f"m and resi {f}-{t}")
    for f, t in partner:
        cmd.color(spec["partner_color"].replace("#", "0x"), f"m and resi {f}-{t}")
    cmd.orient("m")
    if spec.get("view"):
        cmd.set_view(spec["view"])
    w, h = spec["width"], spec["height"]
    cmd.png(f"{spec['outdir']}/panel_a.png", width=w, height=h, ray=1)
    for lo, hi, hexcol in PLDDT_BANDS:
        cmd.color(hexcol, f"m and b > {lo - 0.001} and b < {hi}")
    cmd.png(f"{spec['outdir']}/panel_b.png", width=w, height=h, ray=1)

    stats = {"length": len(plddt), "mean_plddt": round(sum(plddt.values()) / len(plddt), 2),
             "segments": {n: round(sum(plddt[i] for i in range(f, t + 1)) / (t - f + 1), 1) for n, f, t, _ in spec["segments"]},
             "partner": partner, "contacts": contacts, "view": list(cmd.get_view())}
    json.dump(stats, open(f"{spec['outdir']}/structure_stats.json", "w"), indent=1)


if __name__ == "__main__":
    if sys.argv[1:] == ["--selftest"]:
        assert ranges([5, 1, 2, 3, 9, 10]) == [[1, 5], [9, 10]] and ranges([1, 6]) == [[1, 1], [6, 6]]
        print("selftest ok")
    else:
        main(sys.argv[1])
