#!/usr/bin/env python3
"""
synteny_census.py
=================
Recurrent gene-neighbourhood (synteny) analysis for an MGnifams family against the
MGnify proteins DuckDB/parquet store (~6 billion proteins).

Given a plain-text list of family member protein_ids (e.g. exported from a seed or full
MSA), it:

  1. maps each protein_id  ->  contig_id                  (metadata parquet, query 1)
     -- SKIPPED when --contigs is supplied (the contig_ids are already known)
  2. pulls EVERY gene on those contigs (the neighbourhoods) (metadata parquet, query 2)
  3. attaches Pfam annotations to all neighbourhood genes   (pfam parquet,     query 3)
  4. orders genes per contig, locates the family 'anchor' gene(s), and extracts the
     surrounding window (within --window-genes ranks AND --window-bp bp), then runs in
     one of two modes:

TWO MODES
---------
DISCOVERY (no --targets) -- "I don't know the architecture I'm looking for yet."
    Ranks every Pfam appearing in the anchors' windows by how many INDEPENDENT contigs
    (distinct cluster_rep) it neighbours, plus an enrichment ratio against that Pfam's
    background frequency on the same contigs. Enrichment demotes Pfams that are merely
    abundant in these communities (transposases, ribosomal proteins) rather than
    specifically adjacent to the family. Read the table, pick candidates, re-run in TEST
    mode. This is the entry point for a new family, and it costs nothing extra: every
    gene on the contigs is already in memory.

TEST (--targets PF...) -- "I have a hypothesis, hold it to a number."
    Classifies each contig POSITIVE / NEGATIVE(confident) / NEGATIVE(edge-ambiguous) /
    AMBIGUOUS_EDGE / ISOLATED, dereplicates positives by cluster_rep (independence),
    records strand co-directionality, and prints a numerator/denominator VERDICT. The
    neighbour-Pfam table is still written, so you always see what else is around.

Query 2 returns the *whole* contig (all genes), which is what a neighbourhood requires --
restricting it to the family protein_ids would leave no neighbours to analyse.

Notes / assumptions
-------------------
* IN-lists are passed as *registered relations* (hash semi-joins) rather than giant SQL
  literals: far cheaper on a 6B-row scan and injection-safe.
* protein_id / contig_id column types are read from the parquet schema and the input ids
  are CAST to match, so the script works whether ids are integers or strings.
* Pfam accessions are canonicalised to PFxxxxx, so 'PF02566', 'pf2566' and the bare
  integer 2566 (how the pfam parquet stores them) all compare equal.
* Contig edges come from metadata.contig_length when present, else from the span of the
  called genes. An anchor closer than --window-bp to a contig edge has a truncated
  window: its 'no-hit' side is ambiguous, never a confident negative.
* Performance depends on the parquet layout. If these scans are slow, partitioning or
  sorting the metadata parquet by contig_id (and/or protein_id) enables row-group
  pruning; pass --threads / --memory-limit to tune DuckDB.

Usage
-----
    # DISCOVERY: what does this family sit next to?
    python synteny_census.py ids.txt \
        --metadata /path/to/mgy_proteins_metadata.parquet \
        --pfam     /path/to/mgy_proteins_pfam.parquet

    # TEST: is it the GH25 endolysin cassette, and in how many independent contigs?
    python synteny_census.py ids.txt \
        --metadata /path/to/mgy_proteins_metadata.parquet \
        --pfam     /path/to/mgy_proteins_pfam.parquet \
        --targets PF01183,PF25309 [--context-pfams PF13472] \
        [--contigs contig_ids.txt] \
        [--window-genes 5] [--window-bp 10000] \
        [--small-orf-aa 120] [--id-strip-prefix MGYP] [--top 25] \
        [--outdir output/synteny] [--out-prefix synteny] \
        [--threads 8] [--memory-limit 32GB]

    python synteny_census.py --self-test      # synthetic end-to-end check
"""

import argparse
import logging
import os
import re
import sys
from collections import defaultdict

import duckdb
import pandas as pd


LOG = logging.getLogger("synteny")


# ----------------------------------------------------------------------------------
# small helpers
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
            tok = tok.split()[0].split("/", 1)[0]          # 'MGYP..../1-176' -> 'MGYP....'
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


def is_int_type(dtype):
    return "INT" in dtype.upper() and "POINT" not in dtype.upper()


def sanitise_ids_for_type(ids, dtype, what="id"):
    """If the column is integer, keep only numeric ids (warn on drops)."""
    if is_int_type(dtype):
        good = [i for i in ids if re.fullmatch(r"-?\d+", i)]
        dropped = len(ids) - len(good)
        if dropped:
            LOG.warning("%d %s(s) are non-numeric but the column is %s; dropped.",
                        dropped, what, dtype)
        return good
    return ids


def strand_to_int(s):
    m = {"+": 1, "-": -1, "1": 1, "-1": -1, "plus": 1, "minus": -1, "forward": 1, "reverse": -1}
    if s is None:
        return 0
    return m.get(str(s).strip().lower(), 0)


def pfam_key(acc):
    """Canonical PFxxxxx. Accepts 'PF01183', 'PF01183.2', 'pf1183', 1183, '1183'."""
    s = str(acc).strip().upper().split(".", 1)[0]
    if s.startswith("PF"):
        s = s[2:]
    return f"PF{int(s):05d}" if s.isdigit() else s


def bp_gap(a_start, a_end, b_start, b_end):
    if a_end < b_start:
        return b_start - a_end
    if b_end < a_start:
        return a_start - b_end
    return 0  # overlapping


def aa_len(start, end):
    return max(0, (int(end) - int(start) + 1) // 3)


# ----------------------------------------------------------------------------------
# DuckDB queries
# ----------------------------------------------------------------------------------
def q1_protein_to_contig(con, meta_path, ids, id_type):
    """protein_id -> contig_id  (skipped entirely when --contigs is given)."""
    con.register("ids_rel", pd.DataFrame({"pid": ids}))
    sql = f"""
        SELECT DISTINCT m.protein_id AS protein_id, m.contig_id AS contig_id
        FROM read_parquet(?) m
        JOIN ids_rel i ON m.protein_id = CAST(i.pid AS {id_type})
    """
    return con.execute(sql, [meta_path]).df()


def q2_contig_neighbourhoods(con, meta_path, contig_ids, contig_type, has_contig_length):
    """contig_id -> every gene on the contig (the neighbourhood)."""
    con.register("contig_rel", pd.DataFrame({"cid": [str(c) for c in contig_ids]}))
    length_col = "m.contig_length" if has_contig_length else "NULL AS contig_length"
    sql = f"""
        SELECT m.protein_id, m.contig_id, m.contig_name, m.cluster_rep,
               m.start_position, m.end_position, m.strand, {length_col}
        FROM read_parquet(?) m
        JOIN contig_rel c ON m.contig_id = CAST(c.cid AS {contig_type})
    """
    return con.execute(sql, [meta_path]).df()


def q3_pfam(con, pfam_path, protein_ids, pid_type):
    """protein_id -> Pfam rows for every neighbourhood gene."""
    con.register("pid_rel", pd.DataFrame({"pid": [str(p) for p in protein_ids]}))
    sql = f"""
        SELECT p.*
        FROM read_parquet(?) p
        JOIN pid_rel x ON p.protein_id = CAST(x.pid AS {pid_type})
    """
    return con.execute(sql, [pfam_path]).df()


# ----------------------------------------------------------------------------------
# neighbourhood assembly + per-anchor classification
# ----------------------------------------------------------------------------------
def classify_anchor(anchor, genes_sorted, targets, ctx, window_genes, window_bp, small_orf_aa):
    """genes_sorted: list of gene dicts on the contig, ascending by start, with 'rank'.

    Returns (row, window) where window is the in-window neighbour genes, each annotated
    with dist_genes/dist_bp/same_strand -- the discovery pass aggregates over it.
    """
    r = anchor["rank"]
    clen = anchor.get("contig_length")
    if clen and clen > 0:
        left_span = anchor["start"] - 1
        right_span = clen - anchor["end"]
    else:                                    # no contig_length: fall back to called-gene span
        left_span = anchor["start"] - min(g["start"] for g in genes_sorted)
        right_span = max(g["end"] for g in genes_sorted) - anchor["end"]

    neighbours = [g for g in genes_sorted if g["rank"] != r]

    # The window: a neighbour must be within --window-genes ranks AND --window-bp bp.
    # AND, not OR: with OR the bp cap can never exclude anything, so a gene 2 positions
    # away but 76 kb down a sparse contig would count as 'adjacent' -- not a cassette.
    window = []
    for g in neighbours:
        gap = bp_gap(anchor["start"], anchor["end"], g["start"], g["end"])
        gdist = abs(g["rank"] - r)
        if gdist <= window_genes and gap <= window_bp:
            window.append(dict(g, dist_bp=gap, dist_genes=gdist,
                               same_strand=g["strand_i"] == anchor["strand_i"]))

    hits = sorted((g for g in window if g["pfams"] & targets), key=lambda g: g["dist_bp"])

    # cassette context within the window (context pfams + small no-pfam ORFs)
    ctx_present, small_orfs = set(), 0
    for g in window:
        ctx_present |= (g["pfams"] & ctx)
        if not g["pfams"] and aa_len(g["start"], g["end"]) <= small_orf_aa:
            small_orfs += 1

    left_observed = left_span >= window_bp
    right_observed = right_span >= window_bp

    if hits:
        status = "POSITIVE"
    elif not neighbours:
        status = "ISOLATED"
    elif left_observed and right_observed:
        status = "NEGATIVE"                # confident: full window is on-contig
    elif left_observed or right_observed:
        status = "NEGATIVE_PARTIAL"        # one side truncated by a contig edge
    else:
        status = "AMBIGUOUS_EDGE"          # both sides truncated

    nearest = hits[0] if hits else None
    row = {
        "protein_id": anchor["protein_id"],
        "contig_id": anchor["contig_id"],
        "contig_name": anchor.get("contig_name"),
        "cluster_rep": anchor["cluster_rep"],
        "anchor_strand": anchor["strand_i"],
        "n_genes_on_contig": len(genes_sorted),
        "n_genes_in_window": len(window),
        "left_span_bp": int(left_span),
        "right_span_bp": int(right_span),
        "status": status,
        "target_pfam": ";".join(sorted(nearest["pfams"] & targets)) if nearest else "",
        "target_protein_id": nearest["protein_id"] if nearest else "",
        "dist_genes": nearest["dist_genes"] if nearest else "",
        "dist_bp": nearest["dist_bp"] if nearest else "",
        "target_same_strand": nearest["same_strand"] if nearest else "",
        "context_pfams": ";".join(sorted(ctx_present)),
        "small_orf_candidates": small_orfs,
    }
    return row, window


def discover_neighbour_pfams(windows, neigh):
    """Rank the Pfams that recur in the anchors' windows -- the 'I don't know what I'm
    looking for yet' pass. Use it to choose --targets, then re-run for a hard verdict.

    windows: list of (anchor_gene, in_window_neighbour_genes).
    neigh:   every gene on every anchor-bearing contig -- the local background.

    n_indep_contigs (distinct anchor cluster_rep) is the column that matters: 500 hits from
    one over-sampled genome is one observation. enrichment contrasts the Pfam's frequency
    among window genes with its frequency across all genes on the same contigs, which
    demotes Pfams that are merely abundant in these communities (transposases, ribosomal
    proteins) rather than specifically adjacent to the family.
    """
    import statistics

    n_anchors_total = len(windows)
    contigs, reps, n_anchors = defaultdict(set), defaultdict(set), defaultdict(int)
    dist_genes, dist_bp, same_strand = defaultdict(list), defaultdict(list), defaultdict(list)
    win_genes_with_pfam, win_genes_total = defaultdict(int), 0

    for anchor, window in windows:
        win_genes_total += len(window)
        seen_here = set()
        for g in window:
            for pf in g["pfams"]:
                win_genes_with_pfam[pf] += 1
                dist_genes[pf].append(g["dist_genes"])
                dist_bp[pf].append(g["dist_bp"])
                same_strand[pf].append(bool(g["same_strand"]))
                if pf not in seen_here:            # count each anchor once per Pfam
                    seen_here.add(pf)
                    n_anchors[pf] += 1
                    contigs[pf].add(anchor["contig_id"])
                    reps[pf].add(str(anchor["cluster_rep"]))

    # local background: Pfam frequency across ALL genes on the same contigs
    bg_total = len(neigh)
    bg_with_pfam = defaultdict(int)
    for pfams in neigh["pfams"]:
        for pf in pfams:
            bg_with_pfam[pf] += 1

    rows = []
    for pf in win_genes_with_pfam:
        obs = win_genes_with_pfam[pf] / win_genes_total if win_genes_total else 0.0
        bg = bg_with_pfam[pf] / bg_total if bg_total else 0.0
        rows.append({
            "pfam": pf,
            "n_indep_contigs": len(reps[pf]),
            "n_contigs": len(contigs[pf]),
            "n_anchors": n_anchors[pf],
            "frac_anchors": round(n_anchors[pf] / n_anchors_total, 4) if n_anchors_total else 0,
            "median_dist_genes": statistics.median(dist_genes[pf]),
            "median_dist_bp": int(statistics.median(dist_bp[pf])),
            "frac_same_strand": round(sum(same_strand[pf]) / len(same_strand[pf]), 3),
            "window_gene_freq": round(obs, 5),
            "contig_bg_freq": round(bg, 5),
            "enrichment": round(obs / bg, 2) if bg else float("inf"),
        })
    ddf = pd.DataFrame(rows)
    if not ddf.empty:
        ddf = ddf.sort_values(["n_indep_contigs", "enrichment", "n_anchors"],
                              ascending=False).reset_index(drop=True)
    return ddf


def ascii_map(genes_sorted, anchor_ids, targets, ctx):
    parts = []
    for g in genes_sorted:
        tag = g["protein_id"]
        pf = ",".join(sorted(g["pfams"])) if g["pfams"] else "no Pfam"
        mark = ""
        if g["protein_id"] in anchor_ids:
            mark = " <<ANCHOR"
        elif g["pfams"] & targets:
            mark = " <<TARGET"
        elif g["pfams"] & ctx:
            mark = " <context>"
        strand = "+" if g["strand_i"] == 1 else "-" if g["strand_i"] == -1 else "?"
        parts.append(f"[{g['start']}-{g['end']}({strand}) {tag} | {pf}{mark}]")
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
                    help="text file of contig_ids (one per line) already known to carry the "
                         "family. Skips the protein_id -> contig_id lookup (query 1).")
    ap.add_argument("--targets", default=None,
                    help="decisive Pfam accessions to test adjacency against (comma-sep, "
                         "e.g. PF01183,PF25309). OMIT for DISCOVERY MODE: no hypothesis is "
                         "tested, instead every Pfam recurring in the anchors' windows is "
                         "ranked so you can see what the family actually sits next to.")
    ap.add_argument("--context-pfams", default="",
                    help="reported-but-not-decisive Pfams (e.g. GDSL PF13472). Comma-sep.")
    ap.add_argument("--top", type=int, default=25,
                    help="how many neighbour Pfams to print in discovery mode")
    ap.add_argument("--window-genes", type=int, default=5)
    ap.add_argument("--window-bp", type=int, default=10000)
    ap.add_argument("--small-orf-aa", type=int, default=120,
                    help="max length (aa) of a no-Pfam ORF to flag as holin/spanin candidate")
    ap.add_argument("--id-strip-prefix", default=None,
                    help="prefix to strip from ids (e.g. MGYP) if parquet stores bare numbers")
    ap.add_argument("--outdir", default="output/synteny",
                    help="output directory (created if missing); log.txt is written here")
    ap.add_argument("--out-prefix", default="synteny", help="basename prefix for result files")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--memory-limit", default=None)
    ap.add_argument("--self-test", action="store_true",
                    help="run a synthetic end-to-end check and exit")
    args = ap.parse_args(argv)
    if not args.self_test and not (args.ids and args.metadata and args.pfam):
        ap.error("ids, --metadata and --pfam are required (unless --self-test)")
    return args


def run(args):
    targets = {pfam_key(t) for t in (args.targets or "").split(",") if t.strip()}
    ctx = {pfam_key(t) for t in (args.context_pfams or "").split(",") if t.strip()}
    discovery_only = not targets

    inputs = [args.metadata, args.pfam, args.ids] + ([args.contigs] if args.contigs else [])
    for p in inputs:
        if not os.path.exists(p):
            sys.exit(f"ERROR: file not found: {p}")

    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(os.path.join(args.outdir, "log.txt"))
    prefix = os.path.join(args.outdir, args.out_prefix)

    LOG.info("START synteny census | mode=%s targets=%s context=%s window=%dgenes/%dbp",
             "DISCOVERY (no --targets: ranking neighbour Pfams)" if discovery_only else "TEST",
             ",".join(sorted(targets)) or "-", ",".join(sorted(ctx)) or "-",
             args.window_genes, args.window_bp)
    LOG.info("outputs -> %s/", args.outdir)

    con = duckdb.connect()
    if args.threads:
        con.execute(f"PRAGMA threads={args.threads}")
    if args.memory_limit:
        con.execute(f"PRAGMA memory_limit='{args.memory_limit}'")

    # ---- schema discovery -------------------------------------------------------
    meta_schema = parquet_schema(con, args.metadata)
    pfam_schema = parquet_schema(con, args.pfam)
    for req in ("protein_id", "contig_id", "cluster_rep", "start_position",
                "end_position", "strand"):
        if req not in meta_schema:
            sys.exit(f"ERROR: metadata parquet is missing required column '{req}'. "
                     f"Found: {list(meta_schema)}")
    pid_type_meta = meta_schema["protein_id"]
    contig_type = meta_schema["contig_id"]
    has_contig_length = "contig_length" in meta_schema
    pfam_pid_col = find_col(pfam_schema, ["protein_id"])
    pfam_acc_col = find_col(pfam_schema, ["pfam_id", "pfam_acc", "pfam", "pfam_accession",
                                          "accession"])
    pfam_name_col = find_col(pfam_schema, ["pfam_name", "name", "description", "pfam_desc"])
    if pfam_pid_col is None or pfam_acc_col is None:
        sys.exit("ERROR: could not identify protein_id/pfam-accession columns in pfam parquet. "
                 f"Found: {list(pfam_schema)}")
    pid_type_pfam = pfam_schema[pfam_pid_col]
    LOG.info("schemas OK | protein_id=%s contig_id=%s contig_length=%s pfam_acc=%s",
             pid_type_meta, contig_type, "yes" if has_contig_length else "NO (gene-span edges)",
             pfam_acc_col)

    # ---- input ids --------------------------------------------------------------
    ids = read_ids(args.ids, args.id_strip_prefix)
    ids = sanitise_ids_for_type(ids, pid_type_meta, "protein_id")
    if not ids:
        sys.exit("ERROR: no usable protein_ids in input.")

    # ---- STEP 1: protein_id -> contig_id (skipped when --contigs is given) -------
    if args.contigs:
        contig_ids = sanitise_ids_for_type(read_ids(args.contigs), contig_type, "contig_id")
        if not contig_ids:
            sys.exit("ERROR: no usable contig_ids in --contigs file.")
        LOG.info("[1/3] SKIPPED protein_id->contig_id lookup: %d contig_ids supplied via "
                 "--contigs (%d input protein_ids)", len(contig_ids), len(ids))
    else:
        LOG.info("[1/3] mapping %d input protein_ids -> contig_ids ...", len(ids))
        anchors = q1_protein_to_contig(con, args.metadata, ids, pid_type_meta)
        if anchors.empty:
            sys.exit("ERROR: none of the input ids were found in the metadata parquet.")
        found = set(anchors["protein_id"].astype(str))
        missing = [i for i in ids if i not in found]
        if missing:
            LOG.warning("%d id(s) not found in metadata (e.g. %s)", len(missing), missing[:5])
        contig_ids = sorted(anchors["contig_id"].astype(str).unique())
        LOG.info("[1/3] done: %d protein_ids on %d distinct contigs", len(found), len(contig_ids))

    # ---- STEP 2: neighbourhoods -------------------------------------------------
    LOG.info("[2/3] fetching all genes on %d contigs ...", len(contig_ids))
    neigh = q2_contig_neighbourhoods(con, args.metadata, contig_ids, contig_type,
                                     has_contig_length)
    neigh = neigh.drop_duplicates(subset=["contig_id", "protein_id"])
    if neigh.empty:
        sys.exit("ERROR: no genes found on the requested contigs.")
    id_set = set(ids)
    neigh["pid_str"] = neigh["protein_id"].astype(str)
    anchor_ids = set(neigh.loc[neigh["pid_str"].isin(id_set), "pid_str"])
    if not anchor_ids:
        sys.exit("ERROR: none of the input protein_ids are present on the given contigs.")
    LOG.info("[2/3] done: %d genes on %d contigs; %d anchors matched",
             len(neigh), neigh["contig_id"].nunique(), len(anchor_ids))
    if args.contigs:
        absent = len(id_set) - len(anchor_ids)
        if absent:
            LOG.warning("%d input protein_id(s) are not on any supplied contig", absent)

    # ---- STEP 3: Pfam annotations -----------------------------------------------
    LOG.info("[3/3] fetching Pfam annotations for %d neighbourhood genes ...",
             neigh["pid_str"].nunique())
    all_pids = sorted(neigh["pid_str"].unique())
    pf = q3_pfam(con, args.pfam, all_pids, pid_type_pfam)
    pf_map = defaultdict(set)
    name_map = {}
    for _, row in pf.iterrows():
        acc = pfam_key(row[pfam_acc_col])
        pf_map[str(row[pfam_pid_col])].add(acc)
        if pfam_name_col:
            name_map[acc] = row[pfam_name_col]
    LOG.info("[3/3] done: %d Pfam rows for %d annotated genes", len(pf), len(pf_map))

    # ---- build neighbourhoods ---------------------------------------------------
    LOG.info("building neighbourhoods and classifying anchors ...")
    neigh = neigh.copy()
    neigh["strand_i"] = neigh["strand"].map(strand_to_int)
    neigh["start"] = neigh["start_position"].astype("int64")
    neigh["end"] = neigh["end_position"].astype("int64")
    neigh["pfams"] = neigh["pid_str"].map(lambda p: pf_map.get(p, set()))

    per_anchor, maps, windows = [], [], []
    for cid, grp in neigh.groupby("contig_id"):
        grp = grp.sort_values("start").reset_index(drop=True)
        clen = grp["contig_length"].iloc[0] if has_contig_length else None
        clen = int(clen) if pd.notna(clen) else None
        genes = []
        for rank, (_, g) in enumerate(grp.iterrows()):
            genes.append({
                "rank": rank, "protein_id": g["pid_str"], "contig_id": cid,
                "contig_name": g.get("contig_name"), "cluster_rep": g["cluster_rep"],
                "start": int(g["start"]), "end": int(g["end"]), "contig_length": clen,
                "strand_i": int(g["strand_i"]), "pfams": set(g["pfams"]),
            })
        contig_anchors = [x for x in genes if x["protein_id"] in anchor_ids]
        for a in contig_anchors:
            row, window = classify_anchor(a, genes, targets, ctx, args.window_genes,
                                          args.window_bp, args.small_orf_aa)
            per_anchor.append(row)
            windows.append((a, window))
        if contig_anchors:
            maps.append((cid, ascii_map(genes, anchor_ids, targets, ctx)))

    adf = pd.DataFrame(per_anchor)
    LOG.info("classified %d anchors on %d contigs", len(adf), adf["contig_id"].nunique())

    # ---- discovery: which Pfams actually recur next to the family? ---------------
    # Always run: it costs nothing (every gene is already in memory) and it is the only
    # way to pick --targets when the architecture is not known up front.
    ddf = discover_neighbour_pfams(windows, neigh)
    ddf.to_csv(f"{prefix}_neighbour_pfams.tsv", sep="\t", index=False)
    LOG.info("discovery: %d distinct Pfams in the anchor windows -> %s_neighbour_pfams.tsv",
             len(ddf), prefix)

    # ---- aggregate to per-contig (a contig is POSITIVE if any anchor is) --------
    order = ["POSITIVE", "NEGATIVE", "NEGATIVE_PARTIAL", "AMBIGUOUS_EDGE", "ISOLATED"]

    contig_rows = []
    for cid, sub in adf.groupby("contig_id"):
        s = set(sub["status"])
        st = next(o for o in order if o in s)
        pos = sub[sub["status"] == "POSITIVE"]
        rep_row = pos.iloc[0] if not pos.empty else sub.iloc[0]
        contig_rows.append({
            "contig_id": cid, "status": st,
            "cluster_rep": rep_row["cluster_rep"],
            "target_pfam": rep_row["target_pfam"],
            "dist_genes": rep_row["dist_genes"], "dist_bp": rep_row["dist_bp"],
            "target_same_strand": rep_row["target_same_strand"],
            "context_pfams": rep_row["context_pfams"],
            "small_orf_candidates": rep_row["small_orf_candidates"],
        })
    cdf = pd.DataFrame(contig_rows)

    # ---- verdict numbers --------------------------------------------------------
    n_input = len(ids)
    n_contigs = len(cdf)
    counts = cdf["status"].value_counts().to_dict()
    n_pos, n_neg, n_negp, n_amb, n_iso = (counts.get(o, 0) for o in order)

    assess_strict = n_pos + n_neg
    assess_lenient = n_pos + n_neg + n_negp
    frac_strict = n_pos / assess_strict if assess_strict else float("nan")
    frac_lenient = n_pos / assess_lenient if assess_lenient else float("nan")

    pos_contigs = cdf[cdf["status"] == "POSITIVE"]
    indep_pos = pos_contigs["cluster_rep"].astype(str).nunique()
    indep_assess = cdf[cdf["status"].isin(["POSITIVE", "NEGATIVE"])]["cluster_rep"].astype(str).nunique()
    same_strand_pos = (pos_contigs["target_same_strand"] == True).sum()  # noqa: E712

    # ---- write outputs ----------------------------------------------------------
    adf.to_csv(f"{prefix}_anchors.tsv", sep="\t", index=False)
    cdf.to_csv(f"{prefix}_contigs.tsv", sep="\t", index=False)
    with open(f"{prefix}_maps.txt", "w") as fh:
        for cid, m in maps:
            st = cdf.loc[cdf["contig_id"] == cid, "status"].iloc[0]
            fh.write(f"# contig {cid}  [{st}]\n{m}\n\n")
    LOG.info("wrote %s_anchors.tsv, %s_contigs.tsv, %s_maps.txt", prefix, prefix, prefix)

    # ---- verdict ----------------------------------------------------------------
    def pct(x):
        return "n/a" if x != x else f"{x*100:.0f}%"

    header = ("SYNTENY DISCOVERY  (no --targets: ranking what the family sits next to)"
              if discovery_only else
              "SYNTENY CENSUS VERDICT  (target Pfams: " + ", ".join(sorted(targets)) + ")")
    lines = [
        "=" * 72,
        header,
        "=" * 72,
        f"input protein_ids ............ {n_input}",
        f"mapped to a contig ........... {len(anchor_ids)}",
        f"distinct contigs ............. {n_contigs}",
        "-" * 72,
    ]

    if discovery_only:
        n_assessable = n_contigs - n_iso
        lines += [
            f"contigs with >=1 neighbour (usable) .. {n_assessable}",
            f"ISOLATED (lone gene on fragment) ..... {n_iso}",
            "-" * 72,
        ]
        if ddf.empty:
            lines.append("No Pfam-annotated neighbours in any window -- nothing to rank. The "
                         "family's neighbourhoods are unannotated (or all anchors are "
                         "ISOLATED); synteny cannot be assessed from Pfam alone here.")
        else:
            top = ddf.head(args.top)
            lines.append(f"TOP {len(top)} NEIGHBOUR PFAMS  (of {len(ddf)}; ranked by independent "
                         "contigs, then enrichment)")
            lines.append("")
            lines.append(f"{'pfam':<10}{'indep':>6}{'contigs':>8}{'anchors':>8}{'%anch':>7}"
                         f"{'d_genes':>8}{'d_bp':>8}{'strand':>7}{'enrich':>8}")
            for _, r in top.iterrows():
                lines.append(f"{r['pfam']:<10}{r['n_indep_contigs']:>6}{r['n_contigs']:>8}"
                             f"{r['n_anchors']:>8}{r['frac_anchors']*100:>6.0f}%"
                             f"{r['median_dist_genes']:>8.0f}{r['median_dist_bp']:>8}"
                             f"{r['frac_same_strand']:>7.2f}{r['enrichment']:>8.1f}")
            lines += [
                "",
                "indep   = distinct anchor cluster_reps -- THE column to read. Many hits from "
                "one",
                "          over-sampled genome are one observation, not many.",
                "enrich  = how much more often the Pfam sits in the window than on the contigs",
                "          at large. ~1 means 'merely abundant here', not 'adjacent to you'.",
                "strand  = fraction co-directional with the anchor (operon-like cassettes are",
                "          usually co-directional).",
                "-" * 72,
                "NEXT: pick candidates from the table, then re-run with a hypothesis to get a",
                "      hard numerator/denominator verdict, e.g.:",
                f"      --targets {','.join(ddf.head(2)['pfam'])}",
            ]
        lines.append("=" * 72)
        report = "\n".join(lines)
        print("\n" + report)
        with open(f"{prefix}_verdict.txt", "w") as fh:
            fh.write(report + "\n")
        LOG.info("DONE (discovery) | %d neighbour Pfams ranked; %d/%d contigs usable",
                 len(ddf), n_contigs - n_iso, n_contigs)
        return cdf

    lines += [
        f"POSITIVE  (target adjacent) .......... {n_pos}",
        f"NEGATIVE  (confident, full window) ... {n_neg}",
        f"NEGATIVE  (edge-ambiguous) ........... {n_negp}",
        f"AMBIGUOUS (both sides truncated) ..... {n_amb}",
        f"ISOLATED  (lone gene on fragment) .... {n_iso}",
        "-" * 72,
        f"strict  positive fraction  = {n_pos}/{assess_strict} = {pct(frac_strict)}"
        "   (confident calls only)",
        f"lenient positive fraction  = {n_pos}/{assess_lenient} = {pct(frac_lenient)}"
        "   (edge-ambiguous counted as negatives)",
        f"independent positives (distinct cluster_rep) = {indep_pos}"
        f"   of {indep_assess} independent assessable",
    ]
    if n_pos:
        lines.append(f"co-directional with anchor (same strand)     = {same_strand_pos}/{n_pos}")
    lines.append("-" * 72)

    # heuristic call -- raw numbers above are what matter; this is a summary aid.
    if assess_strict == 0:
        call = "INCONCLUSIVE: no confidently assessable contigs (mostly contig edges)."
    elif frac_strict >= 0.5 and indep_pos >= 3:
        call = ("STRONG: target adjacency is recurrent AND independent "
                f"({indep_pos} distinct clusters) -- the synteny claim is well supported.")
    elif n_pos >= 3 and indep_pos >= 2:
        call = ("MODERATE: recurrent across >=2 independent clusters, but either the "
                "fraction or the independence count is limited -- report exact numbers.")
    elif n_pos >= 1:
        call = ("WEAK: adjacency seen but in few and/or non-independent contigs "
                "-- not enough to claim conserved synteny.")
    else:
        call = "NEGATIVE: no target adjacency among assessable contigs."
    lines.append("VERDICT: " + call)
    if indep_pos < n_pos:
        lines.append("  note: some positives share a cluster_rep (near-identical) -- cite the "
                     "independent count, not the raw positive count.")
    lines.append("=" * 72)

    report = "\n".join(lines)
    print("\n" + report)
    with open(f"{prefix}_verdict.txt", "w") as fh:
        fh.write(report + "\n")
    LOG.info("DONE | POSITIVE=%d NEGATIVE=%d NEGATIVE_PARTIAL=%d AMBIGUOUS_EDGE=%d ISOLATED=%d",
             n_pos, n_neg, n_negp, n_amb, n_iso)
    return cdf


# ----------------------------------------------------------------------------------
# self-test: the real test parquets have one gene per contig, so they cannot exercise
# POSITIVE / NEGATIVE. Build a synthetic store that does.
# ----------------------------------------------------------------------------------
def self_test():
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        # contig 1: anchor 100 next to target gene 101 (PF01183)       -> POSITIVE
        # contig 2: anchor 200, neighbours but no target, long contig  -> NEGATIVE
        # contig 3: anchor 300 alone on the contig                     -> ISOLATED
        # contig 4: anchor 400 next to PF01183 again (other cluster)   -> POSITIVE
        # PF09999 sits on contig 2 in-window and on contig 4 far away (out of window),
        # so discovery must rank PF01183 (2 independent contigs) above it (1).
        meta = pd.DataFrame([
            (1, 100, 1_000, 1_500, 1, 100, 50_000),
            (1, 101, 2_000, 2_600, 1, 101, 50_000),
            (2, 200, 20_000, 20_900, -1, 200, 50_000),
            (2, 201, 30_500, 31_000, 1, 201, 50_000),
            (3, 300, 500, 900, 1, 300, 1_200),
            (4, 400, 1_000, 1_500, 1, 400, 90_000),
            (4, 401, 3_000, 3_600, 1, 401, 90_000),
            (4, 402, 80_000, 80_600, 1, 402, 90_000),
        ], columns=["contig_id", "protein_id", "start_position", "end_position",
                    "strand", "cluster_rep", "contig_length"])
        meta["contig_name"] = "contig_" + meta["contig_id"].astype(str)
        pfam = pd.DataFrame([(101, 1183), (201, 9999), (401, 1183), (402, 9999)],
                            columns=["protein_id", "pfam_accession"])

        meta_p = os.path.join(tmp, "meta.parquet")
        pfam_p = os.path.join(tmp, "pfam.parquet")
        ids_p = os.path.join(tmp, "ids.txt")
        contigs_p = os.path.join(tmp, "contigs.txt")
        meta.to_parquet(meta_p)
        pfam.to_parquet(pfam_p)
        with open(ids_p, "w") as fh:
            fh.write("MGYP100\nMGYP200\nMGYP300\nMGYP400\n")   # exercises --id-strip-prefix
        with open(contigs_p, "w") as fh:
            fh.write("1\n2\n3\n4\n")

        outdir = os.path.join(tmp, "out")
        base = ["--metadata", meta_p, "--pfam", pfam_p, "--id-strip-prefix", "MGYP",
                "--window-bp", "10000", "--outdir", outdir]
        tgt = base + ["--targets", "PF01183"]

        cdf = run(parse_args([ids_p] + tgt))
        got = dict(cdf.set_index("contig_id")["status"])
        assert got == {1: "POSITIVE", 2: "NEGATIVE", 3: "ISOLATED", 4: "POSITIVE"}, got
        assert cdf.loc[cdf.contig_id == 1, "target_pfam"].iloc[0] == "PF01183"

        # --contigs must skip query 1 and give an identical answer
        cdf2 = run(parse_args([ids_p] + tgt + ["--contigs", contigs_p]))
        assert dict(cdf2.set_index("contig_id")["status"]) == got

        # discovery mode: no --targets, the recurring neighbour must surface on its own
        run(parse_args([ids_p] + base))
        ddf = pd.read_csv(os.path.join(outdir, "synteny_neighbour_pfams.tsv"), sep="\t")
        assert ddf.iloc[0]["pfam"] == "PF01183", ddf
        assert ddf.iloc[0]["n_indep_contigs"] == 2, ddf   # contigs 1 and 4, distinct reps
        # PF09999's far-away copy on contig 4 is out of window -> only 1 independent contig
        assert ddf.set_index("pfam").loc["PF09999", "n_indep_contigs"] == 1, ddf
        # in-window enrichment must beat the contig background
        assert ddf.iloc[0]["enrichment"] > 1.0, ddf

        # integer pfam accession 1183 must match the 'PF01183' target string
        assert pfam_key(1183) == pfam_key("PF01183.7") == "PF01183"

    print("\nself-test OK: POSITIVE/NEGATIVE/ISOLATED, --contigs skip, discovery ranking, "
          "pfam id normalisation")


def main():
    args = parse_args()
    if args.self_test:
        logging.basicConfig(level=logging.INFO)
        self_test()
        return
    run(args)


if __name__ == "__main__":
    main()
