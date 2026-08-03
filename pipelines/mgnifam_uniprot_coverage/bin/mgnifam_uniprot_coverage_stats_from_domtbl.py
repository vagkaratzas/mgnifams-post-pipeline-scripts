#!/usr/bin/env python3
"""
MGnifams x UniProt coverage statistics, at sequence and residue level,
partitioned by family category lists.

Designed for TrEMBL scale: chunks are processed independently (one job per
chunk), then reduced. Because hmmsearch chunks partition the TARGET database,
every target sequence appears in exactly one chunk -- so per-target interval
merging is fully local and needs no cross-chunk reduction.

Stage 1 (map, parallel over chunks):
    python mgnifam_uniprot_coverage_stats_from_domtbl.py map \
        --domtbl uniprot_trembl_chunk_000001_mgnifams.domtbl.gz \
        --lists /path/to/lists \
        --outdir /path/to/out \
        --min-score 25

    Emits, per chunk:
      <stem>.pertarget.tsv.gz  target, tlen, category, n_res, intervals
      <stem>.summary.tsv       category, n_targets, n_residues, n_families,
                               n_targets_unannotated
      <stem>.families.tsv.gz   category, family  (for correct cross-chunk union)
      logs/chunk_NNNNNN.log    progress log

    --lists may be omitted (e.g. for the Pfam passes), in which case only the
    "any" pseudo-category is computed.

Stage 1b (map, MGnifam-EXCLUSIVE):
    Run the Pfam pass for a chunk first, then feed its per-target file back as
    a mask for the MGnifams pass over the SAME chunk:

    python mgnifam_uniprot_coverage_stats_from_domtbl.py map \
        --domtbl uniprot_trembl_chunk_000001_mgnifams.domtbl.gz \
        --mask   uniprot_trembl_chunk_000001_pfam.pertarget.tsv.gz \
        --lists /path/to/lists --outdir /path/to/out_exclusive

    Every count then excludes residues Pfam already explains, and
    n_targets_unannotated counts the sequences Pfam missed entirely.

Stage 2 (reduce):
    python mgnifam_uniprot_coverage_stats_from_domtbl.py reduce \
        --outdir /path/to/out --lists /path/to/lists \
        --totals-csv annotation_percentage_increase.csv

    --totals-csv supplies the whole-database denominators -- sequences and
    residues including everything with no hit at all, which this script never
    sees -- turning coverage into a quotable percentage-point gain.
    --total-sequences/--total-residues give the same two numbers directly.

    Those denominators are valid only when EVERY chunk of the search is present
    in --outdir, so pair them with --expect-chunks N, which fails the run rather
    than quietly understating every percentage.

Notes
-----
* Uses envelope coordinates (env from/to) by default; --no-env switches to
  alignment coordinates for a stricter sensitivity check.
* MGnifam HMM names in the domtbl are bare integers ("444"), while the
  category lists carry the padded form ("MGYF0000000444"). Both sides are
  normalised to the bare form, so lookups match.
* Bit score is preferred over E-value for thresholding: per-chunk E-values
  are only comparable across chunks if hmmsearch was given a common -Z/--domZ
  (the MGnifams passes were, `-Z 149234636`; the Pfam passes were not, they
  used --cut_ga instead). Bit score is safe either way.
* Category lists are NOT assumed disjoint. A family may appear in several;
  reduce reports the pairwise overlaps so this is explicit rather than
  assumed. The pseudo-category "any" is the union of all MGnifams hits and
  is the denominator for MGnifam-exclusive residue calculations.
* Per-category residues may sum to more than "any", because families from
  different lists can cover overlapping spans of the same target. Reduce
  reports this excess.
"""

import argparse
import bisect
import csv
import glob
import gzip
import logging
import os
import re
import sys
import time
from collections import defaultdict


# ---------------------------------------------------------------- utilities

def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def chunk_stem(path):
    """uniprot_trembl_chunk_000001_mgnifams.domtbl.gz -> that name, sans ext."""
    stem = os.path.basename(path)
    for suffix in (".domtbl.gz", ".domtbl", ".tsv.gz", ".tsv"):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return os.path.splitext(stem)[0]


def chunk_label(path):
    """Extract 'chunk_000001' for the log filename; fall back to the stem."""
    m = re.search(r"chunk[_-](\d+)", os.path.basename(path))
    return "chunk_%s" % m.group(1) if m else chunk_stem(path)


def setup_logging(outdir, label, verbose=False):
    logdir = os.path.join(outdir, "logs")
    os.makedirs(logdir, exist_ok=True)
    logpath = os.path.join(logdir, "%s.log" % label)

    logger = logging.getLogger("mgnifam_coverage")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    fh = logging.FileHandler(logpath, mode="w")
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    logger.info("log file: %s", logpath)
    return logger


def merge_intervals(ivs):
    """Merge closed 1-based intervals. Returns (merged_list, total_residues)."""
    if not ivs:
        return [], 0
    ivs.sort()
    out = []
    cur_s, cur_e = ivs[0]
    for s, e in ivs[1:]:
        if s <= cur_e + 1:
            if e > cur_e:
                cur_e = e
        else:
            out.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    out.append((cur_s, cur_e))
    return out, sum(e - s + 1 for s, e in out)


def bare_family(fam):
    """MGYF0000000444 -> 444, so list IDs match the raw domtbl query names.

    Normalising the lists (30k ids, once) rather than the rows (10^8+, per
    row) keeps the hot loop free of string work.
    """
    if fam.startswith("MGYF"):
        fam = fam[4:]
    return fam.lstrip("0") or "0"


def pretty_family(fam):
    """Inverse of bare_family for report output; non-numeric names pass through
    unchanged, so the Pfam passes keep their real HMM names."""
    return "MGYF%010d" % int(fam) if fam.isdigit() else fam


def subtract_intervals(ivs, mask):
    """ivs minus mask, both merged/sorted. Returns (list, total_residues).

    This is what makes a coverage figure *MGnifam-exclusive*: `mask` is the
    Pfam-covered span of the same target, so what survives is residues only
    MGnifams explains.
    """
    if not mask:
        return ivs, sum(e - s + 1 for s, e in ivs)
    out = []
    total = 0
    j = 0
    for s, e in ivs:
        while j < len(mask) and mask[j][1] < s:
            j += 1
        cur = s
        k = j
        while k < len(mask) and mask[k][0] <= e:
            ms, me = mask[k]
            if ms > cur:
                out.append((cur, ms - 1))
                total += ms - cur
            cur = max(cur, me + 1)
            k += 1
        if cur <= e:
            out.append((cur, e))
            total += e - cur + 1
    return out, total


def interval_is_masked(mask, s, e):
    """True if [s, e] lies wholly inside one mask interval.

    mask is merged, so a span that survives nowhere must sit inside a single
    interval -- one bisect answers it. Used to decide whether a family
    contributes any exclusive residue at all, without accumulating per-family
    intervals for every target.
    """
    i = bisect.bisect_left(mask, (s + 1,)) - 1
    return i >= 0 and mask[i][1] >= e


def load_mask(path, logger=None):
    """Read a .pertarget.tsv.gz from another pass into {target: merged spans}.

    Only the 'any' rows are used -- that is the full covered span of that
    pass (e.g. all of Pfam) for the target.
    """
    mask = {}
    with open_maybe_gz(path) as fh:
        fh.readline()  # header
        for ln in fh:
            target, _tlen, cat, _nres, spans = ln.rstrip("\n").split("\t")
            if cat != "any" or not spans:
                continue
            mask[target] = tuple(
                tuple(int(x) for x in span.split("-")) for span in spans.split(",")
            )
    if logger:
        logger.info("mask: %d targets from %s", len(mask), path)
    return mask


def load_totals(csv_path, logger=None):
    """(total_sequences, total_residues) of the searched database.

    Reads the `*_before` totals from an annotation_percentage_increase.csv
    produced by compare_annotation_stats.py -- the whole database, including
    sequences with no hit at all, which this script never sees.
    """
    with open(csv_path) as fh:
        row = next(csv.DictReader(fh))
    seqs = int(row["total_sequences_before"])
    res = int(row["total_amino_acids_before"])
    if logger:
        logger.info("db totals from %s: %d sequences, %d residues",
                    csv_path, seqs, res)
    return seqs, res


def load_lists(listdir, logger=None):
    """Return {category: set(bare_family_ids)}. Category name = filename stem."""
    cats = {}
    for path in sorted(glob.glob(os.path.join(listdir, "*.txt"))):
        name = os.path.splitext(os.path.basename(path))[0]
        with open(path) as fh:
            cats[name] = {bare_family(ln.strip()) for ln in fh if ln.strip()}
        if logger:
            logger.info("list %-22s %7d families", name, len(cats[name]))
    if not cats:
        sys.exit("no *.txt lists found in %s" % listdir)
    return cats


def build_lookup(cats):
    """Invert to {mgyf_id: tuple(categories)} for O(1) row-time lookup."""
    lut = defaultdict(list)
    for cat, ids in cats.items():
        for i in ids:
            lut[i].append(cat)
    return {k: tuple(v) for k, v in lut.items()}


# --------------------------------------------------------------------- map

def run_map(args):
    label = chunk_label(args.domtbl)
    logger = setup_logging(args.outdir, label, args.verbose)
    t0 = time.time()
    logger.info("START map  domtbl=%s", args.domtbl)
    logger.info("coords=%s  min_score=%s  max_ievalue=%s",
                "env" if args.env else "ali", args.min_score, args.max_ievalue)

    if args.lists:
        cats = load_lists(args.lists, logger)
        lut = build_lookup(cats)
        cat_names = sorted(cats) + ["any"]
    else:
        logger.info("no --lists given; computing 'any' only")
        cats, lut, cat_names = {}, {}, ["any"]

    mask = {}
    if args.mask:
        # A mask path that does not resolve would silently downgrade this run to
        # plain total coverage while every output still says "exclusive". That is
        # the one failure mode that quietly produces a wrong published number.
        if not os.path.isfile(args.mask) or os.path.getsize(args.mask) == 0:
            sys.exit("--mask %s is missing or empty; refusing to report total "
                     "coverage as exclusive" % args.mask)
        mask = load_mask(args.mask, logger)
        if mask:
            logger.info("EXCLUSIVE mode: reporting only residues absent from the mask")
        else:
            # legitimate for a chunk the mask pass never hit, but loud either way
            logger.error("mask %s parsed to 0 targets: this chunk's output is "
                         "identical to total coverage", args.mask)

    cov = defaultdict(lambda: defaultdict(list))
    tlens = {}
    fams_seen = defaultdict(set)

    lo_idx, hi_idx = (19, 20) if args.env else (17, 18)
    n_rows = n_kept = n_short = 0
    every = args.log_every

    logger.info("parsing rows...")
    with open_maybe_gz(args.domtbl) as fh:
        for line in fh:
            if line.startswith("#"):  # blank lines fall through to the field check
                continue
            n_rows += 1
            if n_rows % every == 0:
                logger.info("  rows=%d kept=%d targets=%d elapsed=%.1fs",
                            n_rows, n_kept, len(cov), time.time() - t0)

            f = line.split(None, 21)
            if len(f) < 21:
                n_short += 1
                continue

            if args.min_score is not None and float(f[13]) < args.min_score:
                continue
            if args.max_ievalue is not None and float(f[12]) > args.max_ievalue:
                continue
            n_kept += 1

            target = f[0]
            query = f[3]
            iv = (int(f[lo_idx]), int(f[hi_idx]))

            if target not in tlens:
                tlens[target] = int(f[2])

            in_cats = lut.get(query, ())
            slot = cov[target]
            slot["any"].append(iv)
            for c in in_cats:
                slot[c].append(iv)

            # a family counts as "seen" only if this hit leaves residues the
            # mask does not already explain
            tmask = mask.get(target)
            if tmask is not None and interval_is_masked(tmask, iv[0], iv[1]):
                continue
            fams_seen["any"].add(query)
            for c in in_cats:
                fams_seen[c].add(query)

    logger.info("parsed rows=%d kept=%d malformed=%d targets=%d elapsed=%.1fs",
                n_rows, n_kept, n_short, len(cov), time.time() - t0)
    if n_short:
        logger.warning("%d rows had fewer than 21 fields and were skipped",
                       n_short)

    stem = chunk_stem(args.domtbl)
    os.makedirs(args.outdir, exist_ok=True)

    n_targets = defaultdict(int)
    n_residues = defaultdict(int)
    n_unannotated = defaultdict(int)  # only meaningful under --mask
    n_overlong = 0

    logger.info("merging intervals and writing per-target detail...")
    detail = os.path.join(args.outdir, stem + ".pertarget.tsv.gz")
    with gzip.open(detail, "wt") as out:
        out.write("target\ttlen\tcategory\tn_res\tintervals\n")
        for target, bycat in cov.items():
            tlen = tlens[target]
            tmask = mask.get(target, ())
            for cat, ivs in bycat.items():
                merged, nres = merge_intervals(ivs)
                if nres > tlen:
                    n_overlong += 1
                if mask:
                    merged, nres = subtract_intervals(merged, tmask)
                    if not nres:  # nothing exclusive here; do not count target
                        continue
                    if not tmask:
                        # the mask pass never touched this sequence at all, so
                        # MGnifams takes it from unannotated to annotated
                        n_unannotated[cat] += 1
                n_targets[cat] += 1
                n_residues[cat] += nres
                spans = ",".join("%d-%d" % (s, e) for s, e in merged)
                out.write("%s\t%d\t%s\t%d\t%s\n" % (target, tlen, cat, nres, spans))

    if n_overlong:
        logger.error("%d (target, category) pairs report more covered residues "
                     "than the target length -- interval merging is wrong",
                     n_overlong)
    else:
        logger.info("sanity check passed: no coverage exceeds target length")

    summary = os.path.join(args.outdir, stem + ".summary.tsv")
    with open(summary, "w") as out:
        out.write("category\tn_targets\tn_residues\tn_families"
                  "\tn_targets_unannotated\n")
        for cat in cat_names:
            out.write("%s\t%d\t%d\t%d\t%d\n" % (
                cat, n_targets[cat], n_residues[cat], len(fams_seen[cat]),
                n_unannotated[cat]))
            logger.info("  %-22s targets=%-10d residues=%-12d families=%-7d "
                        "unannotated_before=%d",
                        cat, n_targets[cat], n_residues[cat],
                        len(fams_seen[cat]), n_unannotated[cat])

    fampath = os.path.join(args.outdir, stem + ".families.tsv.gz")
    with gzip.open(fampath, "wt") as out:
        out.write("category\tfamily\n")
        for cat in cat_names:
            for fam in sorted(pretty_family(f) for f in fams_seen[cat]):
                out.write("%s\t%s\n" % (cat, fam))

    logger.info("wrote %s", detail)
    logger.info("wrote %s", summary)
    logger.info("wrote %s", fampath)
    logger.info("DONE map in %.1fs", time.time() - t0)


# ------------------------------------------------------------------ reduce

def run_reduce(args):
    logger = setup_logging(args.outdir, "reduce", args.verbose)
    t0 = time.time()

    summaries = sorted(glob.glob(os.path.join(args.outdir, "*.summary.tsv")))
    if not summaries:
        sys.exit("no *.summary.tsv found in %s" % args.outdir)
    if args.expect_chunks is not None and len(summaries) != args.expect_chunks:
        sys.exit("expected %d chunk summaries in %s but found %d: an incomplete "
                 "chunk set still reduces cleanly, it just understates every "
                 "whole-database percentage" %
                 (args.expect_chunks, args.outdir, len(summaries)))
    logger.info("reducing %d chunk summaries", len(summaries))

    tot_targets = defaultdict(int)
    tot_residues = defaultdict(int)
    tot_unannotated = defaultdict(int)
    for path in summaries:
        with open(path) as fh:
            fh.readline()  # header; readline (not next) tolerates an empty file
            for ln in fh:
                cat, nt, nr, _nf, nu = ln.rstrip("\n").split("\t")
                tot_targets[cat] += int(nt)
                tot_residues[cat] += int(nr)
                tot_unannotated[cat] += int(nu)

    # family counts must be UNIONED across chunks, never summed: the same
    # family can hit targets in many chunks.
    fam_union = defaultdict(set)
    fampaths = sorted(glob.glob(os.path.join(args.outdir, "*.families.tsv.gz")))
    logger.info("unioning families across %d chunks", len(fampaths))
    for path in fampaths:
        with gzip.open(path, "rt") as fh:
            fh.readline()  # header
            for ln in fh:
                cat, fam = ln.rstrip("\n").split("\t")
                fam_union[cat].add(fam)

    cats = load_lists(args.lists, logger) if args.lists else {}
    any_res = tot_residues.get("any", 0)
    any_tgt = tot_targets.get("any", 0)
    if args.totals_csv:
        db_tgt, db_res = load_totals(args.totals_csv, logger)
    elif args.total_sequences and args.total_residues:
        db_tgt, db_res = args.total_sequences, args.total_residues
        logger.info("db totals given directly: %d sequences, %d residues",
                    db_tgt, db_res)
    else:
        db_tgt = db_res = 0
    if db_res:
        logger.warning("pct_db_* columns treat these %d chunk(s) as the WHOLE "
                       "database -- they are only quotable if every chunk of "
                       "the search is present in %s", len(summaries), args.outdir)

    out_path = os.path.join(args.outdir, "reduced.tsv")
    with open(out_path, "w") as out:
        hdr = ("category\tfamilies_in_list\tfamilies_with_hits\tpct_families_hit"
               "\tn_targets\tpct_targets\tn_residues\tpct_residues")
        if db_res:
            hdr += ("\tpct_db_sequences\tpct_db_residues"
                    "\tn_seqs_unannotated_before\tpp_gain_sequences")
        out.write(hdr + "\n")
        print(hdr)
        for cat in sorted(tot_residues, key=lambda c: (c == "any", c)):
            in_list = str(len(cats[cat])) if cat in cats else ""
            hit = len(fam_union.get(cat, ()))
            if cat in cats and cats[cat]:
                pct_fam = "%.1f" % (100.0 * hit / len(cats[cat]))
            else:
                pct_fam = ""
            pct_t = 100.0 * tot_targets[cat] / any_tgt if any_tgt else 0.0
            pct_r = 100.0 * tot_residues[cat] / any_res if any_res else 0.0
            row = "%s\t%s\t%d\t%s\t%d\t%.2f\t%d\t%.2f" % (
                cat, in_list, hit, pct_fam, tot_targets[cat], pct_t,
                tot_residues[cat], pct_r)
            if db_res:
                # the figures that can be quoted as percentage-point gains:
                # shares of the WHOLE database, unhit sequences included.
                # pct_db_residues under --mask IS the residue-level pp gain;
                # the sequence-level one counts only sequences the mask pass
                # missed entirely, which is the stricter claim.
                row += "\t%.4f\t%.4f\t%d\t%.4f" % (
                    100.0 * tot_targets[cat] / db_tgt,
                    100.0 * tot_residues[cat] / db_res,
                    tot_unannotated[cat],
                    100.0 * tot_unannotated[cat] / db_tgt)
            out.write(row + "\n")
            print(row)

    # cross-category residue overlap: how far the parts exceed the whole
    parts = sum(v for k, v in tot_residues.items() if k != "any")
    if any_res and parts:
        logger.info("per-category residues sum to %d vs 'any' %d (ratio %.3f)",
                    parts, any_res, parts / float(any_res))
        logger.info("categories are not disjoint; report overlaps explicitly")

    # pairwise list overlaps, so novel vs novel_structure_any is resolved
    # empirically rather than assumed
    if cats:
        ov_path = os.path.join(args.outdir, "list_overlaps.tsv")
        names = sorted(cats)
        with open(ov_path, "w") as out:
            out.write("category_a\tcategory_b\tn_a\tn_b\tn_shared"
                      "\tpct_of_a\tpct_of_b\trelation\n")
            for i, a in enumerate(names):
                for b in names[i + 1:]:
                    sa, sb = cats[a], cats[b]
                    shared = len(sa & sb)
                    pa = 100.0 * shared / len(sa) if sa else 0.0
                    pb = 100.0 * shared / len(sb) if sb else 0.0
                    if shared == 0:
                        rel = "disjoint"
                    elif sa <= sb:
                        rel = "%s subset of %s" % (a, b)
                    elif sb <= sa:
                        rel = "%s subset of %s" % (b, a)
                    else:
                        rel = "partial"
                    out.write("%s\t%s\t%d\t%d\t%d\t%.1f\t%.1f\t%s\n" % (
                        a, b, len(sa), len(sb), shared, pa, pb, rel))
                    logger.info("overlap %-20s %-20s shared=%-7d %s",
                                a, b, shared, rel)
        logger.info("wrote %s", ov_path)

    logger.info("wrote %s", out_path)
    logger.info("DONE reduce in %.1fs", time.time() - t0)


# -------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    sub.required = True

    m = sub.add_parser("map", help="process one domtbl chunk")
    m.add_argument("--domtbl", required=True)
    m.add_argument("--lists", default=None,
                   help="directory of *.txt ID lists; omit for Pfam passes")
    m.add_argument("--outdir", required=True)
    m.add_argument("--min-score", type=float, default=None,
                   help="minimum domain bit score (preferred over E-value)")
    m.add_argument("--max-ievalue", type=float, default=None,
                   help="maximum independent E-value; comparable across chunks "
                        "only if hmmsearch was given a common -Z/--domZ")
    m.add_argument("--env", action=argparse.BooleanOptionalAction, default=True,
                   help="envelope (default) rather than alignment coordinates; "
                        "--no-env for ali coords")
    m.add_argument("--mask", default=None,
                   help="a .pertarget.tsv.gz from another pass (e.g. the Pfam "
                        "one for the SAME chunk); every count becomes exclusive "
                        "of the residues it covers")
    m.add_argument("--log-every", type=int, default=5000000,
                   help="log progress every N rows read")
    m.add_argument("--verbose", action="store_true")
    m.set_defaults(func=run_map)

    r = sub.add_parser("reduce", help="aggregate per-chunk outputs")
    r.add_argument("--outdir", required=True,
                   help="directory containing the map-stage outputs")
    r.add_argument("--lists", default=None)
    r.add_argument("--totals-csv", default=None,
                   help="annotation_percentage_increase.csv; its *_before "
                        "totals give the whole-database denominators, so "
                        "coverage can be quoted as a percentage-point gain")
    r.add_argument("--total-sequences", type=int, default=None,
                   help="whole-database sequence count; alternative to "
                        "--totals-csv, must be paired with --total-residues")
    r.add_argument("--total-residues", type=int, default=None,
                   help="whole-database residue count")
    r.add_argument("--expect-chunks", type=int, default=None,
                   help="fail unless exactly N chunk summaries are present; "
                        "the guard against quoting a partial chunk set")
    r.add_argument("--verbose", action="store_true")
    r.set_defaults(func=run_reduce)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
