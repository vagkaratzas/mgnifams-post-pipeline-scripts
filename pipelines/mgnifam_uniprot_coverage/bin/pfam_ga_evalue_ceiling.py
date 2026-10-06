#!/usr/bin/env python3
"""
How permissive are Pfam gathering thresholds, in E-value terms, on UniProtKB?

The Pfam passes ran `hmmsearch --cut_ga` without -Z, so each chunk's E-values use
Z = the number of sequences in that chunk. This rescales them onto one search
space (--target-z, e.g. the whole of UniProtKB) and reports, per Pfam family, the
largest full-sequence E-value and domain i-Evalue among its GA-passing hits. Those
maxima approximate the E-value each family's GA corresponds to (a lower bound:
only realised hits are seen), and can be set against the single MGnifams cutoff.

    pfam_ga_evalue_ceiling.py \\
        --subset /path/swissprot/hmmsearch_pfams 575503 \\
        --subset /path/trembl/hmmsearch_pfams 149234636 \\
        --target-z 149810139 --threads 8 --output pfam_ga_evalue_ceiling.tsv

Chunk sizes come from the chunking: every chunk holds --records-per-chunk
sequences except the last, which holds the remainder of the subset total.
"""

import argparse
import glob
import gzip
import math
import os
import re
import statistics
import sys
from concurrent.futures import ProcessPoolExecutor


def scan_chunk(job):
    """Per family in one chunk: [accession, n_domains, max seq E, max dom i-E], rescaled."""
    path, scale = job
    out = {}
    with gzip.open(path, "rt") if path.endswith(".gz") else open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split(None, 13)
            if len(f) < 13:
                continue
            seq_e, dom_ie = float(f[6]) * scale, float(f[12]) * scale
            hit = out.get(f[3])
            if hit is None:
                out[f[3]] = [f[4], 1, seq_e, dom_ie]
            else:
                hit[1] += 1
                hit[2] = max(hit[2], seq_e)
                hit[3] = max(hit[3], dom_ie)
    return out


def chunk_jobs(directory, total, per_chunk, target_z):
    files = sorted(glob.glob(os.path.join(directory, "*_pfam.domtbl*")))
    n = math.ceil(total / per_chunk)
    index = {int(re.search(r"chunk_(\d+)", os.path.basename(p)).group(1)): p for p in files}
    if sorted(index) != list(range(1, n + 1)):
        sys.exit("%s: expected %d chunks numbered 1..%d for %d sequences, found %d"
                 % (directory, n, n, total, len(files)))
    return [(index[i], target_z / (per_chunk if i < n else total - (n - 1) * per_chunk))
            for i in range(1, n + 1)]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subset", nargs=2, action="append", required=True,
                    metavar=("PFAM_DOMTBL_DIR", "TOTAL_SEQUENCES"))
    ap.add_argument("--records-per-chunk", type=int, default=1000000,
                    help="fasta_records_per_chunk of the search run")
    ap.add_argument("--target-z", type=float, required=True,
                    help="search space to express E-values in")
    ap.add_argument("--cutoff", type=float, default=0.001,
                    help="the MGnifams E-value cutoff to compare against")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    jobs = [job for d, total in args.subset
            for job in chunk_jobs(d, int(total), args.records_per_chunk, args.target_z)]
    fams = {}
    with ProcessPoolExecutor(args.threads) as pool:
        for part in pool.map(scan_chunk, jobs):
            for fam, (acc, n, seq_e, dom_ie) in part.items():
                hit = fams.setdefault(fam, [acc, 0, 0.0, 0.0])
                hit[1] += n
                hit[2] = max(hit[2], seq_e)
                hit[3] = max(hit[3], dom_ie)

    with open(args.output, "w") as out:
        out.write("family\taccession\tn_domains\tmax_seq_evalue\tmax_dom_ievalue\n")
        for fam in sorted(fams):
            acc, n, seq_e, dom_ie = fams[fam]
            out.write("%s\t%s\t%d\t%.3g\t%.3g\n" % (fam, acc, n, seq_e, dom_ie))

    print("chunks=%d families=%d target_z=%g" % (len(jobs), len(fams), args.target_z))
    for col, name in ((2, "max_seq_evalue"), (3, "max_dom_ievalue")):
        vals = [h[col] for h in fams.values()]
        if not vals:
            continue
        print("%s median=%.3g" % (name, statistics.median(vals)))
        for cut in (args.cutoff, args.cutoff * 10):
            k = sum(v > cut for v in vals)
            print("%s > %g: %d (%.1f%%)" % (name, cut, k, 100 * k / len(vals)))


if __name__ == "__main__":
    main()
