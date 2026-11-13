#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_architectures
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=mgnifams_architectures_%j.log
# full domain-architecture counts of the case families' member proteins, release 1.1. Read-only.
# Next to this script: architectures.py, classify_multi_hits.py, and
#   wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.1/Pfam-A.clans.tsv.gz -O Pfam38.1.clans.tsv.gz
# Then on the cluster: sbatch architectures.sh. Writes architectures_<job>.tsv.gz in the submit dir.
set -euo pipefail
: "${V:?set V, see paper/README.md}"
MSA=${MSA:-$V/1.1/output/update_families/full_msa}
: "${DB:?set DB, see paper/README.md}"
HERE=${SLURM_SUBMIT_DIR:-$(dirname "$0")}  # sbatch runs a spool copy, so $0 is not here
CLANS=${CLANS:-$HERE/Pfam38.1.clans.tsv.gz}  # the MGnify 2026_07 Pfam parquet is Pfam 38.1
# MGnifam 25, 243, 32058; SF_526: 1107 20804 28802; SF_1556: 4542 26748; SF_2317: 8562 14909 26162 29438 30944
CASES=${CASES:-"25 243 32058 1107 20804 28802 4542 26748 8562 14909 26162 29438 30944"}
: "${CACHE:?set CACHE, see paper/README.md}"
IMG=${IMG:-$(ls "$CACHE"/*e415e86f361373d5bf62a587c79daada9fa797ce076c764f1ccadb32ba2c3861* 2>/dev/null | head -1)}
RUN=${RUN:-singularity exec ${BIND:+-B $BIND} $IMG}   # RUN=env for a local test
echo "== image: ${IMG:-NOT FOUND (set IMG=)}"
md5sum "$HERE/architectures.py" "$HERE/classify_multi_hits.py" "$CLANS"
files=(); for id in $CASES; do files+=("$MSA/$id.sto.gz"); done
ls -la "${files[@]}" "$DB/mgy_proteins_pfam.parquet"
$RUN python3 "$HERE/architectures.py" "${files[@]}" --pfam "$DB/mgy_proteins_pfam.parquet" --clans "$CLANS" \
  --out "$HERE/architectures_${SLURM_JOB_ID:-local}.tsv.gz" --cpus "${SLURM_CPUS_PER_TASK:-1}"
