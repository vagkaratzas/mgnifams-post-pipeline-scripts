#!/usr/bin/env python3
"""
Pfam-on-MGnify90 statistics for the MGnifams manuscript.

Single pass, self-parallelising. Shards the input by byte offset and runs
one worker per shard, then reduces in-process. Only the "p" array of the
metadata JSON is parsed; "s" and "b" are never touched.

  python3 pfam_stats.py annotations.csv -j 32
  python3 pfam_stats.py shards/*.csv.gz -j 32        # pre-split / gzipped
  python3 pfam_stats.py annotations.csv -j 32 --out-prefix pfam_v38

Byte-offset sharding assumes no embedded newlines inside quoted fields
(true for sequence and metadata columns). Gzipped inputs cannot be
offset-sharded, so each .gz file is handled by one worker -- pass several.
"""
import sys
import os
import csv
import json
import gzip
import argparse
import multiprocessing as mp
from collections import defaultdict, Counter

csv.field_size_limit(1 << 31)

MGNIFAM_MIN = 29          # smallest MGnifam, from the >=25-member seed rule
MGNIFAM_MAX = 1_515_677   # largest MGnifam observed (reported, never applied)
COV_BINS = 20

COLS = {'mgyp': 0, 'sequence': 1, 'full_length': 2, 'cluster_size': 3, 'metadata': 4}


def extract_p(meta):
    """Bracket-match the "p" array out of the metadata JSON.

    Avoids parsing "s", which dominates record size for large clusters.
    """
    i = meta.find('"p":')
    if i == -1:
        return []
    j = meta.find('[', i)
    if j == -1:
        return []
    depth = 0
    for k in range(j, len(meta)):
        c = meta[k]
        if c == '[':
            depth += 1
        elif c == ']':
            depth -= 1
            if depth == 0:
                return json.loads(meta[j:k + 1])
    return []


def merge_len(intervals):
    """Total length of the union of inclusive 1-based intervals."""
    if not intervals:
        return 0
    intervals.sort()
    total = 0
    cs, ce = intervals[0]
    for s, e in intervals[1:]:
        if s <= ce + 1:
            if e > ce:
                ce = e
        else:
            total += ce - cs + 1
            cs, ce = s, e
    return total + ce - cs + 1


def blank():
    return {
        'n_seqs': 0, 'n_residues': 0, 'n_full_length': 0, 'n_with_pfam': 0,
        'n_residues_pfam': 0, 'n_cluster_members': 0,
        'seq_len_hist': Counter(), 'ndom_hist': Counter(), 'nfam_hist': Counter(),
        'cov_hist': Counter(), 'fam_seqs': Counter(), 'fam_doms': Counter(),
        'fam_res': Counter(), 'fam_hmm_max': defaultdict(int),
    }


def consume(lines, acc):
    """Accumulate statistics from an iterable of CSV lines."""
    i_seq, i_meta = COLS['sequence'], COLS['metadata']
    i_full, i_csize = COLS['full_length'], COLS['cluster_size']
    for row in csv.reader(lines):
        if not row or row[0] == 'mgyp':      # header, wherever it lands
            continue
        seq = row[i_seq]
        L = len(seq)
        acc['n_seqs'] += 1
        acc['n_residues'] += L
        acc['seq_len_hist'][L] += 1
        if row[i_full].strip().lower() == 'true':
            acc['n_full_length'] += 1
        try:
            acc['n_cluster_members'] += int(row[i_csize])
        except (ValueError, IndexError):
            pass

        hits = extract_p(row[i_meta])
        if not hits:
            continue

        acc['n_with_pfam'] += 1
        acc['ndom_hist'][len(hits)] += 1

        per_fam = defaultdict(list)
        all_iv = []
        for h in hits:
            acc_id = h[0].split('.')[0]
            hmm_to, s_from, s_to = int(h[4]), int(h[5]), int(h[6])
            if s_to < s_from:
                s_from, s_to = s_to, s_from
            per_fam[acc_id].append((s_from, s_to))
            all_iv.append((s_from, s_to))
            acc['fam_doms'][acc_id] += 1
            if hmm_to > acc['fam_hmm_max'][acc_id]:
                acc['fam_hmm_max'][acc_id] = hmm_to

        acc['nfam_hist'][len(per_fam)] += 1
        for acc_id, ivs in per_fam.items():
            acc['fam_seqs'][acc_id] += 1
            acc['fam_res'][acc_id] += merge_len(ivs)

        covered = merge_len(all_iv)
        acc['n_residues_pfam'] += covered
        if L:
            acc['cov_hist'][min(COV_BINS - 1, int(covered / L * COV_BINS))] += 1
    return acc


def worker(task):
    """Process one shard: (path, start_byte, end_byte) or (path, None, None)."""
    path, start, end = task
    acc = blank()
    if path.endswith('.gz'):
        with gzip.open(path, 'rt', newline='') as fh:
            consume(fh, acc)
    else:
        with open(path, 'rb') as fh:
            if start:
                fh.seek(start)
                fh.readline()            # discard partial line
            buf = []
            pos = fh.tell()
            for raw in fh:
                buf.append(raw.decode('utf-8', 'replace'))
                pos += len(raw)
                if len(buf) >= 50000:
                    consume(buf, acc)
                    buf = []
                if end is not None and pos >= end:
                    break
            if buf:
                consume(buf, acc)
    acc['fam_hmm_max'] = dict(acc['fam_hmm_max'])
    return acc


def combine(parts):
    tot = blank()
    tot['fam_hmm_max'] = {}
    for p in parts:
        for k in ('n_seqs', 'n_residues', 'n_full_length', 'n_with_pfam',
                  'n_residues_pfam', 'n_cluster_members'):
            tot[k] += p[k]
        for k in ('seq_len_hist', 'ndom_hist', 'nfam_hist', 'cov_hist',
                  'fam_seqs', 'fam_doms', 'fam_res'):
            tot[k].update(p[k])
        for a, v in p['fam_hmm_max'].items():
            if v > tot['fam_hmm_max'].get(a, 0):
                tot['fam_hmm_max'][a] = v
    return tot


def quant_values(vals, qs=(0.25, 0.5, 0.75)):
    if not vals:
        return {q: float('nan') for q in qs}
    v = sorted(vals)
    n = len(v)
    out = {}
    for q in qs:
        pos = q * (n - 1)
        lo = int(pos)
        hi = min(lo + 1, n - 1)
        out[q] = v[lo] + (v[hi] - v[lo]) * (pos - lo)
    return out


def quant_hist(hist, qs=(0.25, 0.5, 0.75)):
    keys = sorted(hist)
    total = sum(hist.values())
    if not total:
        return {q: float('nan') for q in qs}
    out = {}
    for q in qs:
        target, seen = q * (total - 1), 0
        for k in keys:
            if seen + hist[k] > target:
                out[q] = float(k)
                break
            seen += hist[k]
        else:
            out[q] = float(keys[-1])
    return out


def report(t, prefix):
    ns, nr = t['n_seqs'], t['n_residues']
    ann = t['n_with_pfam']
    P = print
    P('=== MGnify90 baseline (self-consistency checks) ===')
    P(f'  sequences                    {ns:>18,}   [expect 717,738,164]')
    P(f'  residues                     {nr:>18,}   [expect 121,126,699,551]')
    P(f'  sequences with >=1 Pfam hit  {ann:>18,}   [expect 302,709,789]  '
      f'({100*ann/ns:.2f}%)')
    P(f'  Pfam-covered residues        {t["n_residues_pfam"]:>18,}   '
      f'[expect 48,750,452,731]  ({100*t["n_residues_pfam"]/nr:.2f}%)')
    P(f'  full-length sequences        {t["n_full_length"]:>18,}   '
      f'({100*t["n_full_length"]/ns:.1f}%)')
    if t['n_cluster_members']:
        P(f'  underlying cluster members   {t["n_cluster_members"]:>18,}')

    q = quant_hist(t['seq_len_hist'])
    P(f'\n  MGnify90 sequence length     median {q[0.5]:.0f} aa '
      f'(IQR {q[0.25]:.0f}-{q[0.75]:.0f})   [MGnifam rep median 212 aa]')

    P('\n=== Pfam families on MGnify90 ===')
    m = [v for v in t['fam_seqs'].values() if v]
    q = quant_values(m)
    P(f'  families with >=1 match      {len(m):,}')
    P(f'    PFAM_MED_SIZE              {q[0.5]:,.0f} sequences')
    P(f'    PFAM_IQR_SIZE              {q[0.25]:,.0f}-{q[0.75]:,.0f}')
    P(f'    range                      {min(m):,} - {max(m):,}')

    f = [v for v in m if v >= MGNIFAM_MIN]
    q = quant_values(f)
    P(f'\n  families with >={MGNIFAM_MIN} matches      {len(f):,}   '
      f'<-- like-for-like with MGnifams')
    P(f'    PFAM_MED_SIZE  (floored)   {q[0.5]:,.0f} sequences')
    P(f'    PFAM_IQR_SIZE  (floored)   {q[0.25]:,.0f}-{q[0.75]:,.0f}')

    over = sum(1 for v in m if v > MGNIFAM_MAX)
    P(f'\n  larger than largest MGnifam  {over:,} families exceed {MGNIFAM_MAX:,}'
      f'  ({100*over/len(m):.1f}% of matched Pfam families)')
    P('    [descriptive only -- no upper cap applied; the MGnifam maximum is an')
    P('     observed value, not a construction constraint, so censoring Pfam')
    P('     there would condition on an outcome and understate Pfam sizes]')

    P('\n=== Domain architecture context ===')
    q = quant_hist(t['ndom_hist'])
    P(f'  Pfam domains per annotated seq  median {q[0.5]:.0f} '
      f'(IQR {q[0.25]:.0f}-{q[0.75]:.0f})')
    P(f'    single domain                 {t["ndom_hist"][1]:,} '
      f'({100*t["ndom_hist"][1]/ann:.1f}% of annotated)')
    P(f'    single distinct family        {t["nfam_hist"][1]:,} '
      f'({100*t["nfam_hist"][1]/ann:.1f}% of annotated)')

    P('\n=== Per-sequence Pfam coverage fraction (annotated sequences) ===')
    for b in range(COV_BINS):
        c = t['cov_hist'][b]
        if c:
            P(f'    {b/COV_BINS:.2f}-{(b+1)/COV_BINS:.2f}  {c:>18,} '
              f'({100*c/ann:5.2f}%)')
    partial = sum(t['cov_hist'][b] for b in range(COV_BINS - 4))
    P(f'  annotated but <80% covered     {partial:,} '
      f'({100*partial/ann:.1f}% of annotated)')
    P('  [supports the residue-level gain exceeding the sequence-level gain]')

    fn = f'{prefix}_family_sizes.tsv'
    with open(fn, 'w') as out:
        out.write('pfam_acc\tn_sequences\tn_domains\tn_residues\tmax_hmm_to\n')
        for a in sorted(t['fam_seqs'], key=lambda x: -t['fam_seqs'][x]):
            out.write(f'{a}\t{t["fam_seqs"][a]}\t{t["fam_doms"][a]}\t'
                      f'{t["fam_res"][a]}\t{t["fam_hmm_max"].get(a, 0)}\n')
    P(f'\nper-family table -> {fn}')
    P(f'raw counters     -> {prefix}_counters.json  (re-report without re-scanning)')
    with open(f'{prefix}_counters.json', 'w') as out:
        json.dump({k: (dict(v) if isinstance(v, (Counter, dict)) else v)
                   for k, v in t.items()}, out)


def read_pfam_leng(path):
    """Model length per Pfam accession from Pfam-A.hmm (LENG) or .hmm.dat (#=GF ML)."""
    leng, acc = {}, None
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('ACC '):
                acc = line.split()[1].split('.')[0]
            elif line.startswith('#=GF AC'):
                acc = line.split()[2].split('.')[0]
            elif line.startswith('LENG ') and acc:
                leng[acc] = int(line.split()[1])
            elif line.startswith('#=GF ML') and acc:
                leng[acc] = int(line.split()[2])
    return leng


def report_leng(t, leng):
    P = print
    P('\n=== Pfam model lengths (from Pfam-A.hmm) ===')
    allv = list(leng.values())
    q = quant_values(allv)
    P(f'  all Pfam-A models            {len(allv):,}')
    P(f'    PFAM_MED_LENG              {q[0.5]:,.0f} match states '
      f'(IQR {q[0.25]:,.0f}-{q[0.75]:,.0f})')

    w = [v for v in allv if 75 <= v <= 2000]
    q = quant_values(w)
    P(f'\n  within 75-2,000 match states {len(w):,} '
      f'({100*len(w)/len(allv):.1f}% of Pfam-A)')
    P(f'    PFAM_MED_LENG_75_2000      {q[0.5]:,.0f} match states '
      f'(IQR {q[0.25]:,.0f}-{q[0.75]:,.0f})')
    P(f'    excluded below 75          {sum(1 for v in allv if v < 75):,}')
    P(f'    excluded above 2,000       {sum(1 for v in allv if v > 2000):,}')

    matched = {a: leng[a] for a in t['fam_seqs'] if a in leng}
    miss = [a for a in t['fam_seqs'] if a not in leng]
    if matched:
        q = quant_values(list(matched.values()))
        P(f'\n  matching MGnify90            {len(matched):,}  '
          f'<-- most comparable denominator')
        P(f'    median                     {q[0.5]:,.0f} match states '
          f'(IQR {q[0.25]:,.0f}-{q[0.75]:,.0f})')
        wm = [v for v in matched.values() if 75 <= v <= 2000]
        if wm:
            q = quant_values(wm)
            P(f'    within 75-2,000 ({len(wm):,})  {q[0.5]:,.0f} match states '
              f'(IQR {q[0.25]:,.0f}-{q[0.75]:,.0f})')
    if miss:
        P(f'  !! {len(miss):,} matched accessions absent from the HMM file '
          f'(release mismatch?): {", ".join(sorted(miss)[:5])}')

    sat = [a for a in matched if t['fam_hmm_max'].get(a, 0) >= matched[a]]
    if matched:
        P(f'\n  saturation check: {len(sat):,}/{len(matched):,} '
          f'({100*len(sat)/len(matched):.1f}%) of matched families have an '
          f'observed match reaching the final match state')
        P('    [confirms max_hmm_to is a tight lower bound, but LENG above is')
        P('     read from the HMM file and is the value to cite]')


def build_tasks(paths, jobs):
    tasks = []
    for p in paths:
        if p.endswith('.gz'):
            tasks.append((p, None, None))
            continue
        size = os.path.getsize(p)
        n = max(1, min(jobs, size // (32 << 20) or 1))
        step = size // n
        for i in range(n):
            tasks.append((p, i * step, size if i == n - 1 else (i + 1) * step))
    return tasks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('inputs', nargs='+')
    ap.add_argument('-j', '--jobs', type=int, default=os.cpu_count())
    ap.add_argument('--out-prefix', default='pfam')
    ap.add_argument('--pfam-hmm', help='Pfam-A.hmm or Pfam-A.hmm.dat, for model lengths')
    a = ap.parse_args()

    tasks = build_tasks(a.inputs, a.jobs)
    sys.stderr.write(f'{len(tasks)} shards across {a.jobs} workers\n')
    if a.jobs == 1:
        parts = [worker(t) for t in tasks]
    else:
        with mp.Pool(a.jobs) as pool:
            parts = pool.map(worker, tasks, chunksize=1)
    tot = combine(parts)
    report(tot, a.out_prefix)
    if a.pfam_hmm:
        report_leng(tot, read_pfam_leng(a.pfam_hmm))


if __name__ == '__main__':
    main()