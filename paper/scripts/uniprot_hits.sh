#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_uniprot_hits
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=mgnifams_uniprot_hits_%j.log
# UniProtKB hits per MGnifam (are the biome-exclusive families found in UniProtKB?). Read-only.
# Same library and rule as the paper's UniProtKB coverage: E-values rescaled to Z = 149,810,139,
# a domain row counts if rescaled sequence E and domain c-E are both <= 0.001.
# The HMMs are unchanged in 1.1 (update keeps seeds), so the 1.0 domtbls hold for 1.1 families.
# Check: families with a UniProtKB hit must be 34,144 (paper, Table 1 text).
# Output: uniprot_hits_per_family.tsv (family, subset, proteins, best E, best target, its description).
set -euo pipefail
: "${U:?set U, see paper/README.md}"
OUT=${OUT:-${SLURM_SUBMIT_DIR:-$PWD}/uniprot_hits_per_family.tsv}
TMPD=$(mktemp -d -p "${SLURM_SUBMIT_DIR:-$PWD}" uniprot_hits_parts.XXXX)
export TMPD

one() {  # $1 = domtbl.gz; search Z from the original -Z/--domZ (E-value rescale)
  case $1 in */swissprot/*) ss=swissprot z=575503 ;; *) ss=trembl z=149234636 ;; esac
  zcat "$1" | awk -v s="$(awk -v z=$z 'BEGIN{printf "%.17g", 149810139/z}')" -v ss=$ss 'BEGIN{OFS="\t"; OFMT="%.10g"}
    /^#/ {next}
    { e = $7 * s; if (e > 1e-3 || $12 * s > 1e-3) next
      k = $4 SUBSEP $1; if (!(k in seen)) { seen[k] = 1; n[$4]++ }
      if (!($4 in best) || e < best[$4]) { best[$4] = e; bt[$4] = $1; d = $23; for (i = 24; i <= NF; i++) d = d " " $i; bd[$4] = d } }
    END { for (q in n) print q, ss, n[q], best[q], bt[q], bd[q] }' > "$TMPD/$(basename "$1").tsv"
}
export -f one

ls "$U"/swissprot/hmmsearch_mgnifams/*.domtbl.gz "$U"/trembl/hmmsearch_mgnifams/*.domtbl.gz | tee "$TMPD/files.txt" | wc -l
xargs -P "${SLURM_CPUS_PER_TASK:-4}" -n 1 bash -c 'one "$0"' < "$TMPD/files.txt"
ls "$TMPD"/*.domtbl.gz.tsv | wc -l

# merge chunks: sum proteins, keep the best E per (family, subset)
cat "$TMPD"/*.domtbl.gz.tsv | awk -F'\t' 'BEGIN{OFS="\t"}
  { k = $1 OFS $2; n[k] += $3; if (!(k in b) || $4 + 0 < b[k] + 0) { b[k] = $4; t[k] = $5; d[k] = $6 } }
  END { for (k in n) print k, n[k], b[k], t[k], d[k] }' | sort -k1,1n -k2,2 > "$OUT"
rm -r "$TMPD"

echo "rows: $(wc -l < "$OUT")  -> $OUT"
echo "families with a UniProtKB hit (paper 34,144): $(cut -f1 "$OUT" | sort -u | wc -l)"
echo "families with a SwissProt hit: $(awk -F'\t' '$2=="swissprot"' "$OUT" | wc -l)"
echo "families with a TrEMBL hit: $(awk -F'\t' '$2=="trembl"' "$OUT" | wc -l)"
