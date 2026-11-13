#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_multi_hits
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=mgnifams_multi_hits_%j.log
# classify multi-hit proteins of 1.1 full MSAs (single / split / tandem / mixed / partial_repeat). Read-only.
# Copy classify_multi_hits.py next to this script, then on the cluster: sbatch classify_multi_hits.sh
# CASES="all" runs every family in the full_msa folder (34,419; rows file skipped unless ROWS is set).
set -euo pipefail
: "${V:?set V, see paper/README.md}"
MSA=${MSA:-$V/1.1/output/update_families/full_msa}
# MGnifam 25, 243, 32058; SF_1556: 4542 26748; SF_2317: 8562 14909 26162 29438 30944
CASES=${CASES:-"25 243 32058 4542 26748 8562 14909 26162 29438 30944"}
PY=${PY:-${SLURM_SUBMIT_DIR:-$(dirname "$0")}/classify_multi_hits.py}  # sbatch runs a spool copy, so $0 is not here
ROWS=${ROWS-}
ls -d "$MSA"; md5sum "$PY"
if [ "$CASES" = all ]; then
  python3 "$PY" "$MSA" ${ROWS:+--rows "$ROWS"}
else
  files=(); for id in $CASES; do files+=("$MSA/$id.sto.gz"); done
  ls -la "${files[@]}"
  python3 "$PY" "${files[@]}" --rows "${ROWS:-multi_hits_rows_${SLURM_JOB_ID:-local}.tsv.gz}"
fi
