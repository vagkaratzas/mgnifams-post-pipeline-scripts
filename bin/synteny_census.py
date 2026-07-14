#!/usr/bin/env python3
"""
synteny_census.py
=================
Recurrent gene-neighbourhood (synteny) analysis for an MGnifams family against the
MGnify proteins DuckDB/parquet store (~6 billion proteins).

Given a plain-text list of family member protein_ids (e.g. exported from a seed or full
MSA), it:

  1. maps each protein_id  ->  contig_id                  (metadata parquet, query 1)
  2. pulls EVERY gene on those contigs (the neighbourhoods) (metadata parquet, query 2)
  3. attaches Pfam annotations to all neighbourhood genes   (pfam parquet,     query 3)
  4. orders genes per contig, locates the family 'anchor' gene(s), extracts the
     surrounding window, and tests adjacency to a configurable TARGET Pfam set
     (default: the GH25 endolysin cassette PF01183 / PF25309)
  5. classifies each contig as POSITIVE / NEGATIVE(confident) / NEGATIVE(edge-ambiguous)
     / ISOLATED, dereplicates positives by cluster_rep (independence), records strand
     co-directionality, and prints a numerator/denominator VERDICT.

The three SQL steps are the ones supplied by the user, adapted so that query 2 returns
the *whole* contig (all genes), which is what a neighbourhood requires -- the literal
`AND protein_id IN (...)` in the original query 2 would return only the family genes and
leave no neighbours to analyse.

Notes / assumptions
-------------------
* IN-lists are passed as *registered relations* (hash semi-joins) rather than giant SQL
  literals: far cheaper on a 6B-row scan and injection-safe.
* protein_id / contig_id column types are read from the parquet schema and the input ids
  are CAST to match, so the script works whether ids are integers or strings.
* Contig boundaries are approximated by the span of called genes (we have gene coords,
  not contig length); an anchor closer than --window-bp to the nearest called-gene edge
  is treated as sitting at a possible assembly/annotation break -> its 'no-hit' side is
  ambiguous, never a confident negative.
* Performance depends on the parquet layout. If these scans are slow, partitioning or
  sorting the metadata parquet by contig_id (and/or protein_id) enables row-group
  pruning; pass --threads / --memory-limit to tune DuckDB.

Usage
-----
    python synteny_census.py ids.txt \
        --metadata /path/to/mgy_proteins_metadata.parquet \
        --pfam     /path/to/mgy_proteins_pfam.parquet \
        [--targets PF01183,PF25309] [--context-pfams PF13472] \
        [--window-genes 5] [--window-bp 10000] \
        [--small-orf-aa 120] [--id-strip-prefix MGYP] \
        [--out-prefix synteny] [--threads 8] [--memory-limit 32GB]
"""

import argparse
import os
import re
import sys
from collections import defaultdict

import duckdb
import pandas as pd


# ----------------------------------------------------------------------------------
# small helpers
# ----------------------------------------------------------------------------------
def log(msg):
    print(msg, file=sys.stderr, flush=True)


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


def sanitise_ids_for_type(ids, dtype):
    """If the column is integer, keep only numeric ids (warn on drops)."""
    if is_int_type(dtype):
        good = [i for i in ids if re.fullmatch(r"-?\d+", i)]
        dropped = len(ids) - len(good)
        if dropped:
            log(f"  ! {dropped} id(s) are non-numeric but the column is {dtype}; dropped.")
        return good
    return ids


def strand_to_int(s):
    m = {"+": 1, "-": -1, "1": 1, "-1": -1, "plus": 1, "minus": -1, "forward": 1, "reverse": -1}
    if s is None:
        return 0
    return m.get(str(s).strip().lower(), 0)


def strip_pfam_version(acc):
    return str(acc).split(".", 1)[0].upper().strip()


def bp_gap(a_start, a_end, b_start, b_end):
    if a_end < b_start:
        return b_start - a_end
    if b_end < a_start:
        return a_start - b_end
    return 0  # overlapping


def aa_len(start, end):
    return max(0, (int(end) - int(start) + 1) // 3)


# ----------------------------------------------------------------------------------
# DuckDB queries (the three user-supplied steps, combined)
# ----------------------------------------------------------------------------------
def q1_protein_to_contig(con, meta_path, ids, id_type):
    """protein_id -> contig_id  (returns the anchor rows that were found)."""
    con.register("ids_rel", pd.DataFrame({"pid": ids}))
    sql = f"""
        SELECT DISTINCT m.protein_id AS protein_id, m.contig_id AS contig_id
        FROM read_parquet(?) m
        JOIN ids_rel i ON m.protein_id = CAST(i.pid AS {id_type})
    """
    return con.execute(sql, [meta_path]).df()


def q2_contig_neighbourhoods(con, meta_path, contig_ids, contig_type):
    """contig_id -> every gene on the contig (the neighbourhood)."""
    con.register("contig_rel", pd.DataFrame({"cid": [str(c) for c in contig_ids]}))
    sql = f"""
        SELECT m.protein_id, m.contig_id, m.contig_name, m.cluster_rep,
               m.start_position, m.end_position, m.strand
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
    """genes_sorted: list of gene dicts on the contig, ascending by start, with 'rank'."""
    r = anchor["rank"]
    contig_left = min(g["start"] for g in genes_sorted)
    contig_right = max(g["end"] for g in genes_sorted)
    left_span = anchor["start"] - contig_left
    right_span = contig_right - anchor["end"]

    neighbours = [g for g in genes_sorted if g["rank"] != r]

    # find target genes within the window (gene-rank OR bp)
    hits = []
    for g in neighbours:
        if not (g["pfams"] & targets):
            continue
        gap = bp_gap(anchor["start"], anchor["end"], g["start"], g["end"])
        gdist = abs(g["rank"] - r)
        if gdist <= window_genes or gap <= window_bp:
            hits.append((gap, gdist, g))
    hits.sort(key=lambda t: t[0])

    # cassette context within the window (context pfams + small no-pfam ORFs)
    ctx_present, small_orfs = set(), 0
    for g in neighbours:
        gap = bp_gap(anchor["start"], anchor["end"], g["start"], g["end"])
        if abs(g["rank"] - r) <= window_genes or gap <= window_bp:
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

    nearest = hits[0][2] if hits else None
    return {
        "protein_id": anchor["protein_id"],
        "contig_id": anchor["contig_id"],
        "contig_name": anchor.get("contig_name"),
        "cluster_rep": anchor["cluster_rep"],
        "anchor_strand": anchor["strand_i"],
        "n_genes_on_contig": len(genes_sorted),
        "left_span_bp": int(left_span),
        "right_span_bp": int(right_span),
        "status": status,
        "target_pfam": ";".join(sorted(nearest["pfams"] & targets)) if nearest else "",
        "target_protein_id": nearest["protein_id"] if nearest else "",
        "dist_genes": hits[0][1] if hits else "",
        "dist_bp": hits[0][0] if hits else "",
        "target_same_strand": (nearest["strand_i"] == anchor["strand_i"]) if nearest else "",
        "context_pfams": ";".join(sorted(ctx_present)),
        "small_orf_candidates": small_orfs,
    }


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
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ids", help="text file of family protein_ids (one per line)")
    ap.add_argument("--metadata", required=True, help="mgy_proteins_metadata.parquet")
    ap.add_argument("--pfam", required=True, help="mgy_proteins_pfam.parquet")
    ap.add_argument("--targets", default="PF01183,PF25309",
                    help="decisive Pfam accessions (comma-sep). Default: GH25 endolysin.")
    ap.add_argument("--context-pfams", default="PF13472",
                    help="reported-but-not-decisive Pfams (e.g. GDSL). Comma-sep.")
    ap.add_argument("--window-genes", type=int, default=5)
    ap.add_argument("--window-bp", type=int, default=10000)
    ap.add_argument("--small-orf-aa", type=int, default=120,
                    help="max length (aa) of a no-Pfam ORF to flag as holin/spanin candidate")
    ap.add_argument("--id-strip-prefix", default=None,
                    help="prefix to strip from ids (e.g. MGYP) if parquet stores bare numbers")
    ap.add_argument("--out-prefix", default="synteny")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--memory-limit", default=None)
    args = ap.parse_args()

    targets = {strip_pfam_version(t) for t in args.targets.split(",") if t.strip()}
    ctx = {strip_pfam_version(t) for t in args.context_pfams.split(",") if t.strip()}

    for p in (args.metadata, args.pfam, args.ids):
        if not os.path.exists(p):
            sys.exit(f"ERROR: file not found: {p}")

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
    pfam_pid_col = find_col(pfam_schema, ["protein_id"])
    pfam_acc_col = find_col(pfam_schema, ["pfam_id", "pfam_acc", "pfam", "pfam_accession", "accession"])
    pfam_name_col = find_col(pfam_schema, ["pfam_name", "name", "description", "pfam_desc"])
    if pfam_pid_col is None or pfam_acc_col is None:
        sys.exit(f"ERROR: could not identify protein_id/pfam-accession columns in pfam parquet. "
                 f"Found: {list(pfam_schema)}")
    pid_type_pfam = pfam_schema[pfam_pid_col]

    # ---- input ids --------------------------------------------------------------
    ids = read_ids(args.ids, args.id_strip_prefix)
    ids = sanitise_ids_for_type(ids, pid_type_meta)
    if not ids:
        sys.exit("ERROR: no usable protein_ids in input.")
    log(f"[1/3] {len(ids)} input protein_ids -> querying contigs ...")

    # ---- QUERY 1 ----------------------------------------------------------------
    anchors = q1_protein_to_contig(con, args.metadata, ids, pid_type_meta)
    found = set(anchors["protein_id"].astype(str))
    missing = [i for i in ids if i not in found]
    if missing:
        log(f"  ! {len(missing)} id(s) not found in metadata (e.g. {missing[:5]})")
    if anchors.empty:
        sys.exit("ERROR: none of the input ids were found in the metadata parquet.")
    anchor_ids = set(anchors["protein_id"].astype(str))
    contig_ids = sorted(anchors["contig_id"].astype(str).unique())
    log(f"      {len(anchor_ids)} anchors on {len(contig_ids)} distinct contigs.")

    # ---- QUERY 2 ----------------------------------------------------------------
    log("[2/3] fetching all genes on those contigs ...")
    neigh = q2_contig_neighbourhoods(con, args.metadata, contig_ids, contig_type)
    log(f"      {len(neigh)} genes retrieved across {neigh['contig_id'].nunique()} contigs.")

    # ---- QUERY 3 ----------------------------------------------------------------
    log("[3/3] fetching Pfam annotations for neighbourhood genes ...")
    all_pids = sorted(neigh["protein_id"].astype(str).unique())
    pf = q3_pfam(con, args.pfam, all_pids, pid_type_pfam)
    pf_map = defaultdict(set)
    name_map = {}
    for _, row in pf.iterrows():
        acc = strip_pfam_version(row[pfam_acc_col])
        pf_map[str(row[pfam_pid_col])].add(acc)
        if pfam_name_col:
            name_map[acc] = row[pfam_name_col]
    log(f"      {len(pf)} Pfam rows for {len(pf_map)} annotated genes.")

    # ---- build neighbourhoods ---------------------------------------------------
    neigh = neigh.copy()
    neigh["strand_i"] = neigh["strand"].map(strand_to_int)
    neigh["start"] = neigh["start_position"].astype("int64")
    neigh["end"] = neigh["end_position"].astype("int64")
    neigh["pid_str"] = neigh["protein_id"].astype(str)
    neigh["pfams"] = neigh["pid_str"].map(lambda p: pf_map.get(p, set()))

    per_anchor = []
    maps = []
    for cid, grp in neigh.groupby("contig_id"):
        grp = grp.sort_values("start").reset_index(drop=True)
        genes = []
        for rank, (_, g) in enumerate(grp.iterrows()):
            genes.append({
                "rank": rank, "protein_id": g["pid_str"], "contig_id": cid,
                "contig_name": g.get("contig_name"), "cluster_rep": g["cluster_rep"],
                "start": int(g["start"]), "end": int(g["end"]),
                "strand_i": int(g["strand_i"]), "pfams": set(g["pfams"]),
            })
        contig_anchor_ids = [x["protein_id"] for x in genes if x["protein_id"] in anchor_ids]
        for a in genes:
            if a["protein_id"] not in anchor_ids:
                continue
            res = classify_anchor(a, genes, targets, ctx,
                                  args.window_genes, args.window_bp, args.small_orf_aa)
            per_anchor.append(res)
        if contig_anchor_ids:
            maps.append((cid, ascii_map(genes, anchor_ids, targets, ctx)))

    adf = pd.DataFrame(per_anchor)

    # ---- aggregate to per-contig (a contig is POSITIVE if any anchor is) --------
    def contig_status(sub):
        s = set(sub["status"])
        if "POSITIVE" in s:
            return "POSITIVE"
        if "NEGATIVE" in s:
            return "NEGATIVE"
        if "NEGATIVE_PARTIAL" in s:
            return "NEGATIVE_PARTIAL"
        if "AMBIGUOUS_EDGE" in s:
            return "AMBIGUOUS_EDGE"
        return "ISOLATED"

    contig_rows = []
    for cid, sub in adf.groupby("contig_id"):
        st = contig_status(sub)
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
    n_pos = counts.get("POSITIVE", 0)
    n_neg = counts.get("NEGATIVE", 0)
    n_negp = counts.get("NEGATIVE_PARTIAL", 0)
    n_amb = counts.get("AMBIGUOUS_EDGE", 0)
    n_iso = counts.get("ISOLATED", 0)

    assess_strict = n_pos + n_neg
    assess_lenient = n_pos + n_neg + n_negp
    frac_strict = n_pos / assess_strict if assess_strict else float("nan")
    frac_lenient = n_pos / assess_lenient if assess_lenient else float("nan")

    pos_contigs = cdf[cdf["status"] == "POSITIVE"]
    indep_pos = pos_contigs["cluster_rep"].astype(str).nunique()
    indep_assess = cdf[cdf["status"].isin(["POSITIVE", "NEGATIVE"])]["cluster_rep"].astype(str).nunique()
    same_strand_pos = (pos_contigs["target_same_strand"] == True).sum()  # noqa: E712

    # ---- write outputs ----------------------------------------------------------
    adf.to_csv(f"{args.out_prefix}_anchors.tsv", sep="\t", index=False)
    cdf.to_csv(f"{args.out_prefix}_contigs.tsv", sep="\t", index=False)
    with open(f"{args.out_prefix}_maps.txt", "w") as fh:
        for cid, m in maps:
            st = cdf.loc[cdf["contig_id"] == cid, "status"].iloc[0]
            fh.write(f"# contig {cid}  [{st}]\n{m}\n\n")

    # ---- verdict ----------------------------------------------------------------
    def pct(x):
        return "n/a" if x != x else f"{x*100:.0f}%"

    print("\n" + "=" * 72)
    print("SYNTENY CENSUS VERDICT  (target Pfams: " + ", ".join(sorted(targets)) + ")")
    print("=" * 72)
    print(f"input protein_ids ............ {n_input}")
    print(f"mapped to a contig ........... {len(anchor_ids)}")
    print(f"distinct contigs ............. {n_contigs}")
    print("-" * 72)
    print(f"POSITIVE  (target adjacent) .......... {n_pos}")
    print(f"NEGATIVE  (confident, full window) ... {n_neg}")
    print(f"NEGATIVE  (edge-ambiguous) ........... {n_negp}")
    print(f"AMBIGUOUS (both sides truncated) ..... {n_amb}")
    print(f"ISOLATED  (lone gene on fragment) .... {n_iso}")
    print("-" * 72)
    print(f"strict  positive fraction  = {n_pos}/{assess_strict} = {pct(frac_strict)}"
          f"   (confident calls only)")
    print(f"lenient positive fraction  = {n_pos}/{assess_lenient} = {pct(frac_lenient)}"
          f"   (edge-ambiguous counted as negatives)")
    print(f"independent positives (distinct cluster_rep) = {indep_pos}"
          f"   of {indep_assess} independent assessable")
    if n_pos:
        print(f"co-directional with anchor (same strand)     = {same_strand_pos}/{n_pos}")
    print("-" * 72)

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
    print("VERDICT:", call)
    if indep_pos < n_pos:
        print("  note: some positives share a cluster_rep (near-identical) -- cite the "
              "independent count, not the raw positive count.")
    print("=" * 72)
    print(f"\nWrote: {args.out_prefix}_anchors.tsv, {args.out_prefix}_contigs.tsv, "
          f"{args.out_prefix}_maps.txt")


if __name__ == "__main__":
    main()