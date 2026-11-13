#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_pfam_recount
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=mgnifams_pfam_recount_%j.log
# Pfam baseline recount on 2026_07: Pfam family size = n_domains (Pfam 38.2 hits) on the 2026_07
# cluster reps, the base the 1.1 search sliced. Same statistics as the 1.0 baseline on MGnify90
# (pfam_stats.py: floor 29, quantiles method inclusive, count larger than the largest MGnifam).
# Read-only on the release; writes pfam_2026_07_family_sizes.tsv in the submit dir. sbatch pfam_recount.sh
: "${DB:?set DB, see paper/README.md}"
: "${P11:?set P11, see paper/README.md}"
: "${CACHE:?set CACHE, see paper/README.md}"
IMG=${IMG:-$(ls "$CACHE"/*e415e86f361373d5bf62a587c79daada9fa797ce076c764f1ccadb32ba2c3861* 2>/dev/null | head -1)}
RUN=${RUN:-singularity exec ${BIND:+-B $BIND} $IMG}   # RUN=env for a local test
echo "== image: ${IMG:-NOT FOUND (set IMG=)}"
echo "== README on mgy_clusters"; grep -n -i -A4 "### mgy_clusters" "$DB/README.md" | head -20
echo "== 1.1 full_size: min, max"
MINMAX=$(zcat "$P11" | awk -F'\t' 'NR > 1 { if (n++ == 0 || $9 < mn) mn = $9; if ($9 > mx) mx = $9 } END { print mn, mx }'); echo "$MINMAX"
$RUN python3 - "$DB" ${MINMAX#* } ${SLURM_CPUS_PER_TASK:-8} <<'PY'
# Pure pyarrow (the image has no numpy), same join as the pipeline's extract_unannotated_parquet_slices.py:
# per Pfam row group, the cluster reps in its protein_id range, then is_in.
import sys, statistics, time
from collections import Counter
from multiprocessing import Pool
import pyarrow.compute as pc, pyarrow.dataset as ds, pyarrow.parquet as pq
db, largest, cpus = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
MIN = 29  # floor of the 1.0 baseline (smallest MGnifam by the seed rule)
PFAM = f"{db}/mgy_proteins_pfam.parquet"

def count(rg):
    t = pq.ParquetFile(PFAM).read_row_group(rg, columns=["protein_id", "pfam_accession"])
    if t.num_rows == 0:
        return Counter(), 0, 0
    ids = t.column("protein_id"); b = pc.min_max(ids); rep = ds.field("cluster_rep")
    # shortcut: min/max prunes cluster row groups only if mgy_clusters is sorted by cluster_rep (as the 1.1 run relied on)
    reps = ds.dataset(f"{db}/mgy_clusters.parquet").to_table(
        columns=["cluster_rep"], filter=(rep >= b["min"]) & (rep <= b["max"])).column(0).combine_chunks()
    acc = t.column("pfam_accession").filter(pc.is_in(ids, value_set=reps))
    vc = pc.value_counts(acc)
    return Counter(dict(zip(vc.field("values").to_pylist(), vc.field("counts").to_pylist()))), t.num_rows, len(acc)

if __name__ == "__main__":
    t0 = time.time(); n = pq.ParquetFile(PFAM).num_row_groups
    print(f"pfam row groups {n:,}, {cpus} workers", flush=True)
    cnt, rows, kept = Counter(), 0, 0
    with Pool(cpus) as pool:
        for i, (c, r, k) in enumerate(pool.imap_unordered(count, range(n)), 1):
            cnt.update(c); rows += r; kept += k
            if i % 200 == 0 or i == n:
                print(f"  {i:,}/{n:,} row groups, {time.time() - t0:.0f}s", flush=True)
    print(f"pfam rows {rows:,}, on cluster reps {kept:,}")
    with open("pfam_2026_07_family_sizes.tsv", "w") as out:
        out.write("pfam\tn_domains\n")
        for a in sorted(cnt): out.write(f"PF{a:05d}\t{cnt[a]}\n")
    s = sorted(x for x in cnt.values() if x >= MIN)
    q = statistics.quantiles(s, n=4, method="inclusive")
    print(f"Pfam families hitting reps {len(cnt):,}; with >= {MIN} domains {len(s):,}, "
          f"median {q[1]:,} (IQR {q[0]:,}-{q[2]:,}); larger than the largest MGnifam ({largest:,}): {sum(x > largest for x in s)}; "
          f"larger than 1,515,677 (1.0 largest): {sum(x > 1515677 for x in s)}")
PY
