#!/usr/bin/env python3
"""
synteny_census.py
=================
Gene-neighbourhood / synteny analysis for an MGnifams family, against the MGnify proteins
parquet store (~6 billion proteins) queried with DuckDB.

Input: a text file of family member protein_ids (the MSA members). Output: what the family
recurrently sits next to, and whether the arrangement is conserved.

  1. protein_id -> contig_id                (metadata parquet)  [skipped with --contigs]
  2. every gene on those contigs            (metadata parquet)
  3. Pfam annotations for all those genes   (pfam parquet)
  4. per contig: sort genes, locate the family gene(s) = ANCHORS, take the window around
     each (within --window-genes ranks AND --window-bp bp), and rank the partners

PARTNERS
--------
A partner is anything recurring in the anchors' windows:
  * a Pfam accession, or
  * a cluster_rep, for genes with NO Pfam -- 'dark' neighbours are still nameable, and a
    recurrent unannotated ORF next to a novel family is often the interesting result.

For each partner, over INDEPENDENT contigs (dereplicated by the anchor's cluster_rep -- 300
copies of one over-sampled genome are one observation, not 300):

  offset      signed, strand-aware gene offset. -1 = immediately upstream of the anchor,
              +1 = immediately downstream, reading in the anchor's direction. Unsigned
              distance would make 'always upstream' indistinguishable from 'either side',
              which is the difference between synteny and mere co-occurrence.
  conserved   fraction of independent contigs showing the MODAL arrangement (same offset,
              same relative strand). This is the synteny statistic: co-occurrence says the
              gene is near, conservation says it is in the same place every time. A lineage
              with several anchors votes its own modal arrangement, so a multi-anchor contig
              is neither ignored nor able to outvote a single-anchor lineage.
  frac        positives / (positives + confident negatives), where a confident negative is
              an anchor whose FULL window fits on the contig and lacks the partner. Anchors
              truncated by a contig edge cannot testify to absence and are excluded.
  enrichment  partner frequency in the windows vs on the same contigs at large. ~1 means
              'merely abundant here' (transposases, ribosomal proteins), not 'adjacent'.

--targets is OPTIONAL and only needed to test a named hypothesis (e.g. "is it the GH25
endolysin cassette?"): it adds a POSITIVE/NEGATIVE verdict treating the listed Pfams as ONE
hypothesis (a hit on any of them counts). Everything above is reported either way.

Notes
-----
* IN-lists are passed as registered relations (hash semi-joins), not giant SQL literals:
  far cheaper on a 6B-row scan and injection-safe.
* Column types are read from the parquet schema and input ids CAST to match, so ids may be
  integers or strings.
* Pfam accessions are canonicalised, so 'PF01183', 'pf01183.7' and the bare integer 1183
  (how the pfam parquet stores them) compare equal.
* Contig edges come from metadata.contig_length when present, else from the span of the
  called genes.
* The metadata is occurrence-level: one protein_id can sit on hundreds of contigs.

Usage
-----
    python synteny_census.py ids.txt \
        --metadata /path/to/mgy_proteins_metadata.parquet \
        --pfam     /path/to/mgy_proteins_pfam.parquet \
        [--contigs contig_ids.txt] [--targets PF01183,PF25309] \
        [--window-genes 5] [--window-bp 10000] [--id-strip-prefix MGYP] \
        [--outdir output/synteny] [--top 25] [--threads 16] [--memory-limit 64GB]

Tests: tests/test_synteny_census.py (pytest).
"""

import argparse
import logging
import os
import re
import statistics
import sys
from collections import Counter, defaultdict

import duckdb
import pandas as pd


LOG = logging.getLogger("synteny")


# ----------------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------------
def setup_logging(log_path):
    """Console (stderr) + log.txt, same messages in both."""
    LOG.setLevel(logging.INFO)
    LOG.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%Y-%m-%d %H:%M:%S")
    for h in (logging.StreamHandler(sys.stderr), logging.FileHandler(log_path, mode="w")):
        h.setFormatter(fmt)
        LOG.addHandler(h)


def read_ids(path, strip_prefix=None):
    """One id per line; drops '/start-end' MSA suffixes, blanks, comments, dups."""
    ids, seen = [], set()
    with open(path) as fh:
        for line in fh:
            tok = line.strip()
            if not tok or tok.startswith("#"):
                continue
            tok = tok.split()[0].split("/", 1)[0]      # 'MGYP..../1-176' -> 'MGYP....'
            if strip_prefix and tok.startswith(strip_prefix):
                tok = tok[len(strip_prefix):]
            if tok not in seen:
                seen.add(tok)
                ids.append(tok)
    return ids


def parquet_schema(con, path):
    df = con.execute("DESCRIBE SELECT * FROM read_parquet(?) LIMIT 0", [path]).df()
    return dict(zip(df["column_name"], df["column_type"]))


def find_col(schema, candidates):
    norm = {c.lower().replace(" ", "").replace("_", ""): c for c in schema}
    for cand in candidates:
        k = cand.lower().replace(" ", "").replace("_", "")
        if k in norm:
            return norm[k]
    return None


def sanitise_ids_for_type(ids, dtype, what="id"):
    """If the column is integer, keep only numeric ids and normalise them (warn on drops).

    Normalising matters: stripping 'MGYP' off 'MGYP000108638007' leaves '000108638007',
    which SQL casts fine but never equals the '108638007' the parquet round-trips back,
    so the later string comparisons all miss.
    """
    if "INT" in dtype.upper():
        numeric = [i for i in ids if re.fullmatch(r"-?\d+", i)]
        if len(numeric) < len(ids):
            LOG.warning("%d %s(s) are non-numeric but the column is %s; dropped. "
                        "(Did you mean --id-strip-prefix MGYP?)", len(ids) - len(numeric),
                        what, dtype)
        # int() drops the leading zeros; dict.fromkeys re-dedups ids that only differed by them
        return list(dict.fromkeys(str(int(i)) for i in numeric))
    return ids


def strand_to_int(s):
    m = {"+": 1, "-": -1, "1": 1, "-1": -1, "plus": 1, "minus": -1,
         "forward": 1, "reverse": -1}
    return m.get(str(s).strip().lower(), 0) if s is not None else 0


def pfam_key(acc):
    """Canonical PFxxxxx. Accepts 'PF01183', 'PF01183.2', 'pf1183', 1183, '1183'."""
    s = str(acc).strip().upper().split(".", 1)[0]
    if s.startswith("PF"):
        s = s[2:]
    return f"PF{int(s):05d}" if s.isdigit() else s


def bp_gap(a_start, a_end, b_start, b_end):
    """Intergenic gap; 0 if the genes overlap."""
    if a_end < b_start:
        return b_start - a_end
    if b_end < a_start:
        return a_start - b_end
    return 0


def aa_len(start, end):
    return max(0, (int(end) - int(start) + 1) // 3)


def partner_keys(gene):
    """What this neighbour is called: its Pfams, or its cluster_rep if it has none."""
    if gene["pfams"]:
        return sorted(gene["pfams"])
    return [f"cluster:{gene['cluster_rep']}"]


# ----------------------------------------------------------------------------------
# DuckDB queries
# ----------------------------------------------------------------------------------
def q1_protein_to_contig(con, meta_path, ids, id_type):
    con.register("ids_rel", pd.DataFrame({"pid": ids}))
    return con.execute(f"""
        SELECT DISTINCT m.protein_id, m.contig_id
        FROM read_parquet(?) m
        JOIN ids_rel i ON m.protein_id = CAST(i.pid AS {id_type})
    """, [meta_path]).df()


def q2_contig_neighbourhoods(con, meta_path, contig_ids, contig_type, has_len):
    con.register("contig_rel", pd.DataFrame({"cid": [str(c) for c in contig_ids]}))
    length_col = "m.contig_length" if has_len else "NULL AS contig_length"
    return con.execute(f"""
        SELECT m.protein_id, m.contig_id, m.contig_name, m.cluster_rep,
               m.start_position, m.end_position, m.strand, {length_col}
        FROM read_parquet(?) m
        JOIN contig_rel c ON m.contig_id = CAST(c.cid AS {contig_type})
    """, [meta_path]).df()


def q3_pfam(con, pfam_path, protein_ids, pid_type):
    con.register("pid_rel", pd.DataFrame({"pid": [str(p) for p in protein_ids]}))
    return con.execute(f"""
        SELECT p.*
        FROM read_parquet(?) p
        JOIN pid_rel x ON p.protein_id = CAST(x.pid AS {pid_type})
    """, [pfam_path]).df()


# ----------------------------------------------------------------------------------
# per-anchor window
# ----------------------------------------------------------------------------------
def anchor_window(anchor, genes, window_genes, window_bp):
    """Neighbours within --window-genes ranks AND --window-bp bp, with signed offsets.

    AND, not OR: with OR the bp cap could never exclude anything, so a gene 2 positions
    away but 76 kb down a sparse contig would count as 'adjacent'.

    The offset is signed in the ANCHOR's reading direction, so -1 is always 'immediately
    upstream' whether the anchor is on the + or the - strand. Without that, a partner
    always upstream and a partner on either side look identical.
    """
    r = anchor["rank"]
    ori = anchor["strand_i"] or 1                   # unknown strand: read left-to-right
    clen = anchor.get("contig_length")
    if clen and clen > 0:
        left, right = anchor["start"] - 1, clen - anchor["end"]
    else:                                           # no contig_length: use called-gene span
        left = anchor["start"] - min(g["start"] for g in genes)
        right = max(g["end"] for g in genes) - anchor["end"]
    # in the anchor's frame, 'upstream' is the left edge on +, the right edge on -
    up_span, down_span = (left, right) if ori == 1 else (right, left)

    window = []
    for g in genes:
        if g["rank"] == r:
            continue
        gap = bp_gap(anchor["start"], anchor["end"], g["start"], g["end"])
        if abs(g["rank"] - r) <= window_genes and gap <= window_bp:
            window.append(dict(g,
                               offset=(g["rank"] - r) * ori,
                               gap_bp=gap,
                               same_strand=g["strand_i"] == anchor["strand_i"]))

    if not any(g["rank"] != r for g in genes):
        status = "ISOLATED"                 # lone gene on a fragment: says nothing
    elif up_span >= window_bp and down_span >= window_bp:
        status = "FULL"                     # whole window on-contig: a real negative
    elif up_span >= window_bp or down_span >= window_bp:
        status = "PARTIAL"                  # one side truncated by a contig edge
    else:
        status = "EDGE"                     # both sides truncated
    return window, status, int(up_span), int(down_span)


# ----------------------------------------------------------------------------------
# partner ranking (the main result)
# ----------------------------------------------------------------------------------
def rank_partners(anchors, neigh, targets):
    """anchors: list of dicts with 'rep', 'status', 'window'. Returns the partner table."""
    assessable_reps = {a["rep"] for a in anchors if a["status"] == "FULL"}

    reps, contigs, n_anch = defaultdict(set), defaultdict(set), Counter()
    gaps, aas = defaultdict(list), defaultdict(list)
    arrangements = defaultdict(lambda: defaultdict(Counter))
    win_hits, win_total = Counter(), 0

    for a in anchors:
        win_total += len(a["window"])
        seen = set()
        for g in a["window"]:
            for key in partner_keys(g):
                win_hits[key] += 1
                gaps[key].append(g["gap_bp"])
                aas[key].append(aa_len(g["start"], g["end"]))
                if key in seen:
                    continue                       # count each anchor once per partner
                seen.add(key)
                n_anch[key] += 1
                reps[key].add(a["rep"])
                contigs[key].add(a["contig_id"])
                # one ballot per anchor; the lineage's vote is the mode of its ballots, so a
                # contig carrying several anchors is neither dropped nor able to outvote a
                # single-anchor lineage
                arrangements[key][a["rep"]][(g["offset"], bool(g["same_strand"]))] += 1

    # background: partner frequency across ALL genes on the same contigs
    bg_total = len(neigh)
    bg_hits = Counter()
    for _, g in neigh.iterrows():
        for key in partner_keys(g):
            bg_hits[key] += 1

    rows = []
    for key in n_anch:
        pos = reps[key]
        neg = assessable_reps - pos                # FULL window, partner absent
        ballots = [c.most_common(1)[0][0] for c in arrangements[key].values()]
        modal, n_modal = Counter(ballots).most_common(1)[0]
        obs = win_hits[key] / win_total if win_total else 0.0
        bg = bg_hits[key] / bg_total if bg_total else 0.0
        rows.append({
            "partner": key,
            "type": "pfam" if key.startswith("PF") else "unannotated_cluster",
            "n_indep": len(pos),
            "n_contigs": len(contigs[key]),
            "n_anchors": n_anch[key],
            "frac_indep": round(len(pos) / (len(pos) + len(neg)), 3) if (pos or neg) else "",
            "offset": modal[0],
            "same_strand": modal[1],
            "frac_conserved": round(n_modal / len(ballots), 3),
            "median_gap_bp": int(statistics.median(gaps[key])),
            "median_aa_len": int(statistics.median(aas[key])),
            "enrichment": round(obs / bg, 2) if bg else float("inf"),
            "is_target": key in targets,
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["n_indep", "frac_conserved", "enrichment"],
                            ascending=False).reset_index(drop=True)
    return df, assessable_reps


def ascii_map(genes, anchor_pids, targets):
    parts = []
    for g in genes:
        pf = ",".join(sorted(g["pfams"])) if g["pfams"] else f"cluster:{g['cluster_rep']}"
        if g["protein_id"] in anchor_pids:
            mark = " <<ANCHOR"
        elif g["pfams"] & targets:
            mark = " <<TARGET"
        else:
            mark = ""
        strand = {1: "+", -1: "-"}.get(g["strand_i"], "?")
        parts.append(f"[{g['start']}-{g['end']}({strand}) {g['protein_id']} | "
                     f"{pf} | {aa_len(g['start'], g['end'])}aa{mark}]")
    return " -- ".join(parts)


# ----------------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ids", nargs="?", help="text file of family protein_ids (one per line)")
    ap.add_argument("--metadata", help="mgy_proteins_metadata.parquet")
    ap.add_argument("--pfam", help="mgy_proteins_pfam.parquet")
    ap.add_argument("--contigs", default=None,
                    help="text file of contig_ids already known to carry the family; skips "
                         "the protein_id -> contig_id scan (step 1)")
    ap.add_argument("--targets", default=None,
                    help="OPTIONAL. Pfams to test as one hypothesis, comma-sep (e.g. "
                         "PF01183,PF25309). Adds a POSITIVE/NEGATIVE verdict. Omit it and "
                         "the partner ranking still tells you what the family sits next to.")
    ap.add_argument("--window-genes", type=int, default=5,
                    help="neighbours within this many genes of the anchor (default 5)")
    ap.add_argument("--window-bp", type=int, default=10000,
                    help="AND within this many bp of the anchor (default 10000)")
    ap.add_argument("--id-strip-prefix", default=None,
                    help="prefix to strip from input ids (e.g. MGYP) if the parquet stores "
                         "bare integers")
    ap.add_argument("--outdir", default="output/synteny",
                    help="output directory (created if missing); log.txt is written here")
    ap.add_argument("--out-prefix", default="synteny", help="basename prefix for results")
    ap.add_argument("--top", type=int, default=25, help="partners to print (default 25)")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--memory-limit", default=None)
    args = ap.parse_args(argv)
    if not (args.ids and args.metadata and args.pfam):
        ap.error("ids, --metadata and --pfam are required")
    return args


def run(args):
    targets = {pfam_key(t) for t in (args.targets or "").split(",") if t.strip()}

    for p in [args.metadata, args.pfam, args.ids] + ([args.contigs] if args.contigs else []):
        if not os.path.exists(p):
            sys.exit(f"ERROR: file not found: {p}")

    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(os.path.join(args.outdir, "log.txt"))
    prefix = os.path.join(args.outdir, args.out_prefix)
    LOG.info("START | window=%d genes AND %d bp | targets=%s | out=%s/",
             args.window_genes, args.window_bp, ",".join(sorted(targets)) or "none",
             args.outdir)

    con = duckdb.connect()
    if args.threads:
        con.execute(f"PRAGMA threads={args.threads}")
    if args.memory_limit:
        con.execute(f"PRAGMA memory_limit='{args.memory_limit}'")

    # ---- schema ------------------------------------------------------------------
    meta_schema = parquet_schema(con, args.metadata)
    pfam_schema = parquet_schema(con, args.pfam)
    for req in ("protein_id", "contig_id", "cluster_rep", "start_position",
                "end_position", "strand"):
        if req not in meta_schema:
            sys.exit(f"ERROR: metadata parquet lacks required column '{req}'. "
                     f"Found: {list(meta_schema)}")
    pid_type, contig_type = meta_schema["protein_id"], meta_schema["contig_id"]
    has_len = "contig_length" in meta_schema
    pfam_pid_col = find_col(pfam_schema, ["protein_id"])
    pfam_acc_col = find_col(pfam_schema, ["pfam_accession", "pfam_acc", "pfam_id", "pfam",
                                          "accession"])
    if not pfam_pid_col or not pfam_acc_col:
        sys.exit(f"ERROR: no protein_id/pfam-accession columns in pfam parquet. "
                 f"Found: {list(pfam_schema)}")
    LOG.info("schema OK | protein_id=%s contig_id=%s contig_length=%s", pid_type, contig_type,
             "yes" if has_len else "NO (edges from gene span)")

    ids = sanitise_ids_for_type(read_ids(args.ids, args.id_strip_prefix), pid_type,
                                "protein_id")
    if not ids:
        sys.exit("ERROR: no usable protein_ids in input.")

    # ---- step 1: protein_id -> contig_id -----------------------------------------
    if args.contigs:
        contig_ids = sanitise_ids_for_type(read_ids(args.contigs), contig_type, "contig_id")
        if not contig_ids:
            sys.exit("ERROR: no usable contig_ids in --contigs.")
        LOG.info("[1/3] SKIPPED (--contigs): %d contigs given, %d input protein_ids",
                 len(contig_ids), len(ids))
    else:
        LOG.info("[1/3] mapping %d protein_ids -> contigs ...", len(ids))
        hit = q1_protein_to_contig(con, args.metadata, ids, pid_type)
        if hit.empty:
            sys.exit("ERROR: none of the input ids are in the metadata parquet.")
        missing = set(ids) - set(hit["protein_id"].astype(str))
        if missing:
            LOG.warning("%d id(s) not in metadata (e.g. %s)", len(missing),
                        sorted(missing)[:5])
        contig_ids = sorted(hit["contig_id"].astype(str).unique())
        # step 1 is the slow scan; cache it now so a later failure doesn't cost it again
        with open(f"{prefix}_contigs.txt", "w") as fh:
            fh.write("\n".join(contig_ids) + "\n")
        LOG.info("[1/3] done: %d protein_ids on %d contigs (one protein can sit on many); "
                 "reusable via --contigs %s_contigs.txt",
                 hit["protein_id"].nunique(), len(contig_ids), prefix)

    # ---- step 2: neighbourhoods --------------------------------------------------
    LOG.info("[2/3] fetching every gene on %d contigs ...", len(contig_ids))
    neigh = q2_contig_neighbourhoods(con, args.metadata, contig_ids, contig_type, has_len)
    neigh = neigh.drop_duplicates(subset=["contig_id", "protein_id"])
    if neigh.empty:
        sys.exit("ERROR: no genes found on the requested contigs.")
    neigh["pid_str"] = neigh["protein_id"].astype(str)
    anchor_pids = set(neigh.loc[neigh["pid_str"].isin(set(ids)), "pid_str"])
    if not anchor_pids:
        sys.exit("ERROR: none of the input protein_ids are on the given contigs.")
    LOG.info("[2/3] done: %d genes on %d contigs; %d/%d family proteins found",
             len(neigh), neigh["contig_id"].nunique(), len(anchor_pids), len(ids))

    # ---- step 3: Pfams -----------------------------------------------------------
    LOG.info("[3/3] fetching Pfams for %d neighbourhood genes ...", neigh["pid_str"].nunique())
    pf = q3_pfam(con, args.pfam, sorted(neigh["pid_str"].unique()), pfam_schema[pfam_pid_col])
    pf_map = defaultdict(set)
    # iterate columns, not rows: iterrows upcasts a mixed int/float row to float64, turning
    # protein_id 917027859 into the string "917027859.0" which never matches neigh.pid_str.
    for pid, acc in zip(pf[pfam_pid_col], pf[pfam_acc_col]):
        pf_map[str(pid)].add(pfam_key(acc))
    LOG.info("[3/3] done: %d Pfam rows for %d annotated genes", len(pf), len(pf_map))

    # ---- windows -----------------------------------------------------------------
    LOG.info("building neighbourhoods ...")
    neigh["strand_i"] = neigh["strand"].map(strand_to_int)
    neigh["start"] = neigh["start_position"].astype("int64")
    neigh["end"] = neigh["end_position"].astype("int64")
    neigh["pfams"] = neigh["pid_str"].map(lambda p: pf_map.get(p, set()))

    anchors, maps = [], []
    for cid, grp in neigh.groupby("contig_id"):
        grp = grp.sort_values("start").reset_index(drop=True)
        clen = int(grp["contig_length"].iloc[0]) if has_len and pd.notna(
            grp["contig_length"].iloc[0]) else None
        genes = [{
            "rank": i, "protein_id": g["pid_str"], "contig_id": cid,
            "contig_name": g.get("contig_name"), "cluster_rep": g["cluster_rep"],
            "start": int(g["start"]), "end": int(g["end"]), "contig_length": clen,
            "strand_i": int(g["strand_i"]), "pfams": set(g["pfams"]),
        } for i, (_, g) in enumerate(grp.iterrows())]

        on_contig = [g for g in genes if g["protein_id"] in anchor_pids]
        for a in on_contig:
            window, status, up, down = anchor_window(a, genes, args.window_genes,
                                                     args.window_bp)
            anchors.append({
                "protein_id": a["protein_id"], "contig_id": cid,
                "contig_name": a["contig_name"], "rep": str(a["cluster_rep"]),
                "strand": a["strand_i"], "contig_length": clen,
                "n_genes_on_contig": len(genes), "n_genes_in_window": len(window),
                "upstream_bp": up, "downstream_bp": down, "status": status,
                "window": window,
            })
        if on_contig:
            maps.append((cid, ascii_map(genes, anchor_pids, targets)))

    # ---- partners ----------------------------------------------------------------
    pdf, assessable_reps = rank_partners(anchors, neigh, targets)
    n_indep_total = len({a["rep"] for a in anchors})
    counts = Counter(a["status"] for a in anchors)
    LOG.info("%d anchors on %d contigs (%d independent lineages) | FULL=%d PARTIAL=%d "
             "EDGE=%d ISOLATED=%d", len(anchors), len({a['contig_id'] for a in anchors}),
             n_indep_total, counts["FULL"], counts["PARTIAL"], counts["EDGE"],
             counts["ISOLATED"])
    LOG.info("ranked %d partners (%d Pfam, %d unannotated clusters)", len(pdf),
             int((pdf["type"] == "pfam").sum()) if not pdf.empty else 0,
             int((pdf["type"] == "unannotated_cluster").sum()) if not pdf.empty else 0)

    # ---- write -------------------------------------------------------------------
    adf = pd.DataFrame([{
        k: v for k, v in a.items() if k != "window"
    } | {
        "window_partners": ";".join(f"{k}@{g['offset']:+d}"
                                    for g in sorted(a["window"], key=lambda g: g["offset"])
                                    for k in partner_keys(g)),
        "target_hit": ";".join(sorted({k for g in a["window"] for k in partner_keys(g)
                                       if k in targets})),
    } for a in anchors])
    adf.to_csv(f"{prefix}_anchors.tsv", sep="\t", index=False)
    pdf.to_csv(f"{prefix}_partners.tsv", sep="\t", index=False)
    with open(f"{prefix}_maps.txt", "w") as fh:
        for cid, m in maps:
            fh.write(f"# contig {cid}\n{m}\n\n")

    # ---- report ------------------------------------------------------------------
    lines = [
        "=" * 78,
        "SYNTENY CENSUS  --  what does this family sit next to?",
        "=" * 78,
        f"family proteins (input) ........ {len(ids)}",
        f"anchors found .................. {len(anchors)} occurrences of "
        f"{len(anchor_pids)} proteins",
        f"contigs ........................ {len({a['contig_id'] for a in anchors})}",
        f"independent lineages ........... {n_indep_total}  (distinct anchor cluster_rep)",
        "-" * 78,
        f"FULL     window fully on-contig, absence is informative .... {counts['FULL']}",
        f"PARTIAL  one side truncated by a contig edge ............... {counts['PARTIAL']}",
        f"EDGE     both sides truncated ............................... {counts['EDGE']}",
        f"ISOLATED lone gene on the fragment, says nothing ........... {counts['ISOLATED']}",
        f"-> {len(assessable_reps)} independent lineages are ASSESSABLE (the denominator)",
        "-" * 78,
    ]

    if pdf.empty:
        lines.append("No neighbours in any window. Every anchor is alone on its contig -- "
                     "the assemblies are too fragmented to say anything about synteny.")
    else:
        top = pdf.head(args.top)
        lines += [
            f"TOP {len(top)} PARTNERS of {len(pdf)}  (ranked by independent lineages, then "
            "conservation)",
            "",
            f"{'partner':<18}{'indep':>6}{'frac':>7}{'offset':>7}{'str':>5}{'consv':>7}"
            f"{'gap_bp':>8}{'aa':>6}{'enrich':>8}",
        ]
        for _, r in top.iterrows():
            star = "*" if r["is_target"] else " "
            lines.append(f"{star}{r['partner']:<17}{r['n_indep']:>6}{str(r['frac_indep']):>7}"
                         f"{r['offset']:>+7d}{'same' if r['same_strand'] else 'opp':>5}"
                         f"{r['frac_conserved']:>7}{r['median_gap_bp']:>8}"
                         f"{r['median_aa_len']:>6}{r['enrichment']:>8.1f}")
        lines += [
            "",
            "indep  independent lineages (distinct anchor cluster_rep) with this partner in",
            "       the window. THE column: 300 copies of one genome are one observation.",
            "frac   indep / (indep + confident absences), counting only FULL-window anchors",
            "       as absences -- an edge-truncated window cannot testify to absence.",
            "offset modal signed gene offset, read in the ANCHOR's direction: -1 = directly",
            "       upstream, +1 = directly downstream. str = same/opposite strand.",
            "consv  fraction of independent lineages showing that exact arrangement. THIS is",
            "       synteny: high indep + high consv = same gene, same place, every time.",
            "       High indep + low consv = it is merely nearby (co-occurrence).",
            "enrich vs the partner's frequency on the same contigs at large. ~1 = merely",
            "       abundant here (transposases, ribosomal), not specifically adjacent.",
            "aa     median length; unannotated ~50-150 aa neighbours are holin/spanin-sized.",
        ]
        if targets:
            lines += ["", "-" * 78]
            hit_reps = {a["rep"] for a in anchors
                        if any(k in targets for g in a["window"] for k in partner_keys(g))}
            neg_reps = assessable_reps - hit_reps
            denom = len(hit_reps) + len(neg_reps)
            frac = f"{len(hit_reps)/denom:.0%}" if denom else "n/a"
            lines += [
                f"HYPOTHESIS {','.join(sorted(targets))} (any of them counts as a hit)",
                f"  independent lineages with the target adjacent = {len(hit_reps)}",
                f"  confident absences (FULL window, no target)   = {len(neg_reps)}",
                f"  VERDICT: {len(hit_reps)}/{denom} = {frac} of assessable lineages",
            ]
        lines.append("-" * 78)
        lines.append("NOTE: partners were ranked on the same data that suggested them. Treat "
                     "the top hit as a")
        lines.append("      hypothesis to confirm on held-out members, not as a test that "
                     "has already passed.")

    lines.append("=" * 78)
    report = "\n".join(lines)
    print("\n" + report)
    with open(f"{prefix}_report.txt", "w") as fh:
        fh.write(report + "\n")
    LOG.info("wrote %s_{partners,anchors}.tsv, _maps.txt, _report.txt", prefix)
    LOG.info("DONE")
    return pdf, adf


def main():
    run(parse_args())


if __name__ == "__main__":  # pragma: no cover
    main()
