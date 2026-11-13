#!/usr/bin/env python3
"""Build the MGnifams family-generation figure as an .excalidraw.svg
(clean SVG + embedded Excalidraw payload, so it stays editable in the
VS Code Excalidraw editor)."""
import base64, json, random, zlib, sys

FS = 19          # box font size
FS_S = 16        # edge-label font size
LH = 1.25
random.seed(7)

els = []
def nonce(): return random.randint(1, 2**31)

def base(t, x, y, w, h, **kw):
    d = dict(type=t, version=1, versionNonce=nonce(), isDeleted=False,
             id=kw.pop("id"), fillStyle="solid", strokeWidth=2,
             strokeStyle="solid", roughness=0, opacity=100, angle=0,
             x=x, y=y, strokeColor="#1e1e1e", backgroundColor="transparent",
             width=w, height=h, seed=nonce(), groupIds=[], frameId=None,
             roundness=None, boundElements=[], updated=1740000000000,
             link=None, locked=False)
    d.update(kw)
    return d

def txt(id, x, y, text, size=FS, family=2, color="#1e1e1e",
        container=None, align="center", bold=False):
    lines = text.split("\n")
    w = max(len(l) for l in lines) * size * 0.52
    h = len(lines) * size * LH
    e = base("text", x, y, w, h, id=id, strokeColor=color,
             fontSize=size, fontFamily=family, text=text, textAlign=align,
             verticalAlign="middle", containerId=container,
             originalText=text, lineHeight=LH, autoResize=True)
    if bold: e["_bold"] = True
    els.append(e)
    return e

def node(id, shape, x, y, w, h, text, stroke="#1e1e1e", bg="#ffffff",
         fill="solid", size=FS, family=2, tcolor=None, radius=True):
    e = base(shape, x, y, w, h, id=id, strokeColor=stroke,
             backgroundColor=bg, fillStyle=fill,
             roundness={"type": 3} if (radius and shape == "rectangle") else None)
    els.append(e)
    t = txt(id + "_t", x + 10, y + h / 2 - (len(text.split("\n")) * size * LH) / 2,
            text, size=size, family=family, color=tcolor or stroke, container=id)
    e["boundElements"].append({"type": "text", "id": t["id"]})
    return e

def arrow(id, pts, label=None, start=None, end=None, dashed=False, seg=None):
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    x0, y0 = pts[0]
    e = base("arrow", x0, y0, max(xs) - min(xs), max(ys) - min(ys), id=id,
             points=[[p[0] - x0, p[1] - y0] for p in pts],
             lastCommittedPoint=None,
             startBinding={"elementId": start, "focus": 0, "gap": 4} if start else None,
             endBinding={"elementId": end, "focus": 0, "gap": 4} if end else None,
             startArrowhead=None, endArrowhead="arrow", elbowed=False,
             strokeStyle="dashed" if dashed else "solid")
    els.append(e)
    for eid in (start, end):
        if eid:
            for n in els:
                if n["id"] == eid:
                    n["boundElements"].append({"id": id, "type": "arrow"})
    if label:
        # midpoint of the longest segment
        if seg is None:
            segs = [(abs(pts[i+1][0]-pts[i][0]) + abs(pts[i+1][1]-pts[i][1]), i)
                    for i in range(len(pts)-1)]
            _, seg = max(segs)
        i = seg
        mx = (pts[i][0] + pts[i+1][0]) / 2
        my = (pts[i][1] + pts[i+1][1]) / 2
        t = txt(id + "_l", 0, 0, label, size=FS_S, color="#343a40")
        t["x"] = mx - t["width"] / 2
        t["y"] = my - t["height"] / 2
        t["_arrowlabel"] = True
    return e

# ---------------------------------------------------------------- counts
# counts from the release 1.0 run outputs
N = dict(seed="63,780", lost="2", built="58,062", discarded="5,716", nohit="244", member="207", short="5,091", long="174",
         converged="3,686", capped="31,773", redundant="22,603", final="35,459")  # lost: in unique_reps, in neither pooled file

# ---------------------------------------------------------------- layout
LX, LW = 280, 430          # left column
RX, RW = 940, 480          # right column
LC, RC = LX + LW/2, RX + RW/2
TX, TW = 735, 120          # lane discard terminals
INK, GREY = "#1e1e1e", "#868e96"
RED, GREEN, BLUE = "#e03131", "#2f9e44", "#1971c2"
BG_STEP, BG_IN, BG_OUT, BG_DEC = "#f8f9fa", "#fff5f5", "#ebfbee", "#fff9db"

for hid, cx, title in (("h0", LC, "Refinement (rounds 1–3)"), ("h1", RC, "Final pass, then across families")):
    h = txt(hid, 0, 18, title, size=21, bold=True)
    h["x"] = cx - h["width"] / 2

node("L1", "rectangle", LX, 60, LW, 100,
     f"seed cluster  (≥ 25 sequences)\nMMseqs2 linclust --min-seq-id 0.5 -c 0.9\nn = {N['seed']}  ({N['lost']} gave no output)",
     stroke=RED, bg=BG_IN, tcolor=RED, size=16)
node("L2", "rectangle", LX, 200, LW, 84,
     "align seed-cluster sequences\nFAMSA (pyfamsa)", bg=BG_STEP, size=17)
node("L3", "rectangle", LX, 330, LW, 84,
     "build profile HMM\nhmmbuild --fast", bg=BG_STEP, size=17)
node("L4", "rectangle", LX, 455, LW, 110,
     "search MGnify-UR\nhmmsearch -E 0.001  (full-sequence E-value)\nkeep domains with envelope ≥ 0.9 × HMM length,\ncut to the envelope",
     bg=BG_STEP, size=15)
node("L5", "diamond", LX, 605, LW, 116,
     "new sequence\nrecruited?", bg=BG_DEC, size=17)
node("L6", "rectangle", LX, 765, LW, 84,
     "align recruited regions to the HMM\nhmmalign", bg=BG_STEP, size=16)
node("L6b", "diamond", LX, 890, LW, 116,
     "75 ≤ representative\nlength ≤ 2,000?", bg=BG_DEC, size=16)
node("L7", "rectangle", LX, 1050, LW, 120,
     "reduce redundancy, trim ends\ntrimAl representative filter: identity ≤ 0.8,\nat most 2,000 sequences; cut leading and\ntrailing columns with ≤ 50% occupancy",
     bg=BG_STEP, size=14)
node("T1", "rectangle", TX, 482, TW, 56, "discard\nno hits", stroke=RED, bg=BG_IN, tcolor=RED, size=13)
node("T2", "rectangle", TX, 920, TW, 56, "discard\nlength", stroke=RED, bg=BG_IN, tcolor=RED, size=13)

node("R10", "rectangle", RX, 200, RW, 150,
     "final pass\n1. family HMM: hmmbuild --hand on the current seed MSA\n2. hmmsearch -E 0.001 with the last --fast HMM:\nevery domain, no envelope-length filter",
     bg=BG_STEP, size=14)
node("R11", "diamond", RX, 395, RW, 116,
     "≥ 90% of seed-cluster\nsequences recruited?", bg=BG_DEC, size=16)
node("R12", "rectangle", RX, 555, RW, 104,
     "align to the last --fast HMM\nhmmalign  →  full MSA\nrepresentative = first region of top-scoring sequence",
     bg=BG_STEP, size=15)
node("R13", "diamond", RX, 700, RW, 116,
     "75 ≤ representative\nlength ≤ 2,000?", bg=BG_DEC, size=16)
node("R14", "rectangle", RX, 855, RW, 94,
     f"family (before redundancy removal)\nfamily HMM + seed MSA + full MSA\nn = {N['built']}", stroke=BLUE, bg="#e7f5ff", tcolor=BLUE, size=16)
node("D", "rectangle", 1500, 380, 190, 240,
     f"discard seed cluster\nn = {N['discarded']}\n(both columns)\n\nno hits: {N['nohit']}\n< 90% recruited: {N['member']}\nrepresentative\n  < 75: {N['short']}\n  > 2,000: {N['long']}",
     stroke=RED, bg=BG_IN, tcolor=RED, size=15)

node("B1", "rectangle", RX, 990, RW, 84,
     "after all seed clusters: identical representative\nregion → keep the family with more sequences",
     bg=BG_STEP, size=14)
node("B2", "rectangle", RX, 1115, RW, 104,
     "library search\nhmmsearch -E 0.001: all family HMMs vs all representatives;\nenvelope ≥ 0.9 × HMM length; no self-matches",
     bg=BG_STEP, size=14)
node("B3", "diamond", RX, 1260, RW, 130,
     "Jaccard index of member\nproteins ≥ 0.95?", bg=BG_DEC, size=16)
node("B5", "rectangle", 1500, 1275, 190, 100,
     f"remove smaller family\n(ties: identifier)\nn = {N['redundant']}\n(incl. identical repr.)", stroke=RED, bg=BG_IN, tcolor=RED, size=13)
node("B4", "rectangle", RX, 1435, RW, 80,
     f"MGnifams release 1.0\nn = {N['final']}",  # no clans: they start with release 1.1
     stroke=GREEN, bg=BG_OUT, tcolor=GREEN, size=15)

# ---------------------------------------------------------------- arrows
arrow("A1", [(LC, 160), (LC, 200)], "sequences", "L1", "L2")
arrow("A2", [(LC, 284), (LC, 330)], "seed MSA (round 1)", "L2", "L3")
arrow("A3", [(LC, 414), (LC, 455)], "profile HMM", "L3", "L4")
arrow("A4", [(LC, 565), (LC, 605)], "recruited regions", "L4", "L5")
arrow("A5", [(LC, 721), (LC, 765)], "yes", "L5", "L6")
arrow("A6", [(LC, 849), (LC, 890)], "alignment", "L6", "L6b")
arrow("A6b", [(LC, 1006), (LC, 1050)], "yes", "L6b", "L7")
arrow("A7", [(LX, 1110), (100, 1110), (100, 372), (LX, 372)],
      "next seed MSA\n(round < 3;\nnot realigned)", "L7", "L3", seg=1)
arrow("A8", [(LX + LW, 663), (875, 663), (875, 245), (RX, 245)],
      f"no\n(converged:\n{N['converged']})", "L5", "R10", seg=0)
arrow("A9", [(LX + LW, 1110), (905, 1110), (905, 305), (RX, 305)],
      f"after\nround 3:\n{N['capped']}", "L7", "R10", seg=1)
arrow("AT1", [(LX + LW, 510), (TX, 510)], None, "L4", "T1")
arrow("AT2", [(LX + LW, 948), (TX, 948)], None, "L6b", "T2")
arrow("A11", [(RC, 350), (RC, 395)], "recruited regions", "R10", "R11")
arrow("A12", [(RC, 511), (RC, 555)], "yes", "R11", "R12")
arrow("A13", [(RC, 659), (RC, 700)], "full MSA", "R12", "R13")
arrow("A14", [(RC, 816), (RC, 860)], "yes", "R13", "R14")
arrow("A15", [(RX + RW, 275), (1465, 275), (1465, 420), (1500, 420)], "no hits", "R10", "D", seg=1)
arrow("A16", [(RX + RW, 453), (1500, 453)], "no", "R11", "D")
arrow("A17", [(RX + RW, 758), (1465, 758), (1465, 580), (1500, 580)], "no", "R13", "D", seg=1)
arrow("A18", [(RC, 944), (RC, 990)], "all families", "R14", "B1")
arrow("A19", [(RC, 1074), (RC, 1115)], None, "B1", "B2")
arrow("A20", [(RC, 1219), (RC, 1260)], "matching pairs", "B2", "B3")
arrow("A21", [(RX + RW, 1325), (1500, 1325)], "yes", "B3", "B5")
arrow("A22", [(RC, 1390), (RC, 1435)], "no", "B3", "B4")

# ---------------------------------------------------------------- write
scene = {"type": "excalidraw", "version": 2,
         "source": "https://marketplace.visualstudio.com/items?itemName=pomdtr.excalidraw-editor",
         "elements": [{k: v for k, v in e.items() if not k.startswith("_")} for e in els],
         "appState": {"gridSize": None, "viewBackgroundColor": "#ffffff"},
         "files": {}}

out = sys.argv[1]
open(out.replace(".excalidraw.svg", ".excalidraw"), "w").write(json.dumps(scene, indent=2))

# ---- minimal SVG render (roughness 0 => plain shapes) -------------------
W, H = 1720, 1550
def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

body = []
for e in els:
    t = e["type"]
    if t == "rectangle":
        body.append(f'<rect x="{e["x"]}" y="{e["y"]}" width="{e["width"]}" height="{e["height"]}" '
                    f'rx="{12 if e["roundness"] else 0}" fill="{e["backgroundColor"]}" '
                    f'stroke="{e["strokeColor"]}" stroke-width="2"/>')
    elif t == "diamond":
        x, y, w, h = e["x"], e["y"], e["width"], e["height"]
        body.append(f'<polygon points="{x+w/2},{y} {x+w},{y+h/2} {x+w/2},{y+h} {x},{y+h/2}" '
                    f'fill="{e["backgroundColor"]}" stroke="{e["strokeColor"]}" stroke-width="2"/>')
    elif t == "arrow":
        pts = " ".join(f'{e["x"]+p[0]},{e["y"]+p[1]}' for p in e["points"])
        dash = ' stroke-dasharray="8 6"' if e["strokeStyle"] == "dashed" else ""
        body.append(f'<polyline points="{pts}" fill="none" stroke="{e["strokeColor"]}" '
                    f'stroke-width="2" marker-end="url(#ah)"{dash}/>')

for e in els:
    if e["type"] != "text":
        continue
    lines = e["text"].split("\n")
    fsz = e["fontSize"]
    if e.get("containerId"):
        c = next(n for n in els if n["id"] == e["containerId"])
        cx = c["x"] + c["width"] / 2
        cy = c["y"] + c["height"] / 2 - (len(lines) - 1) * fsz * LH / 2 + fsz * 0.35
        anchor = "middle"
    elif e.get("_arrowlabel"):
        cx = e["x"] + e["width"] / 2
        cy = e["y"] + fsz * 0.9
        anchor = "middle"
        body.append(f'<rect x="{e["x"]-5}" y="{e["y"]-2}" width="{e["width"]+10}" '
                    f'height="{e["height"]+4}" fill="#ffffff"/>')
    else:
        cx, cy, anchor = e["x"], e["y"] + fsz * 0.9, e["textAlign"]
        if anchor == "center":
            cx = e["x"] + e["width"] / 2
        anchor = {"left": "start", "center": "middle", "right": "end"}[anchor]
    for i, ln in enumerate(lines):
        weight = "700" if e.get("_bold") else "600" if (e.get("containerId") and c["type"] != "diamond" and i == 0) else "400"
        body.append(f'<text x="{cx}" y="{cy + i*fsz*LH}" font-family="Nunito, Helvetica, Arial, sans-serif" '
                    f'font-size="{fsz}" font-weight="{weight}" fill="{e["strokeColor"]}" '
                    f'text-anchor="{anchor}">{esc(ln)}</text>')

raw = json.dumps(scene, separators=(",", ":")).encode()
comp = zlib.compress(raw, 9)
payload = base64.b64encode(json.dumps({
    "version": "1", "encoding": "bstring", "compressed": True,
    "encoded": comp.decode("latin-1")}, ensure_ascii=False).encode("latin-1")).decode()

svg = (f'<svg version="1.1" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
       f'width="{W}" height="{H}">\n  <!-- svg-source:excalidraw -->\n'
       f'  <!-- payload-type:application/vnd.excalidraw+json --><!-- payload-version:2 -->'
       f'<!-- payload-start -->{payload}<!-- payload-end -->\n'
       f'  <defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
       f'markerHeight="7" orient="auto-start-reverse">'
       f'<path d="M 0 0 L 10 5 L 0 10 z" fill="#1e1e1e"/></marker></defs>\n'
       f'  <rect width="{W}" height="{H}" fill="#ffffff"/>\n  '
       + "\n  ".join(body) + "\n</svg>\n")
open(out, "w").write(svg)
print("wrote", out)
