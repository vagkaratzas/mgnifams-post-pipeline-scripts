#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_neighbourhood
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
#SBATCH --output=mgnifams_neighbourhood_%j.log
# gene neighbourhood (operon check) of one family's member proteins, release 1.1. Read-only.
# Next to this script: neighbourhood.py, architectures.py, classify_multi_hits.py
#, and
#   wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.1/Pfam-A.clans.tsv.gz -O Pfam38.1.clans.tsv.gz
# Then on the cluster: sbatch neighbourhood.sh (FAMILY=32058 sbatch neighbourhood.sh for another family).
# Writes neighbourhood_<family>_<job>_{partners,proteins}.tsv.gz in the submit dir.
set -euo pipefail
: "${V:?set V, see paper/README.md}"
FAMILY=${FAMILY:-25}
MSA=${MSA:-$V/1.1/output/update_families/full_msa/$FAMILY.sto.gz}
: "${DB:?set DB, see paper/README.md}"
HERE=${SLURM_SUBMIT_DIR:-$(dirname "$0")}  # sbatch runs a spool copy, so $0 is not here
CLANS=${CLANS:-$HERE/Pfam38.1.clans.tsv.gz}  # the parquet is Pfam 38.1
: "${CACHE:?set CACHE, see paper/README.md}"
IMG=${IMG:-$(ls "$CACHE"/*e415e86f361373d5bf62a587c79daada9fa797ce076c764f1ccadb32ba2c3861* 2>/dev/null | head -1 || true)}
RUN=${RUN:-singularity exec ${BIND:+-B $BIND} $IMG}   # RUN=env for a local test
OUT=$HERE/neighbourhood_${FAMILY}_${SLURM_JOB_ID:-local}
echo "== image: ${IMG:-NOT FOUND (set IMG=)}"
md5sum "$HERE/neighbourhood.py" "$HERE/architectures.py" "$HERE/classify_multi_hits.py" "$CLANS"
ls -la "$MSA" "$DB/mgy_clusters_metadata.parquet" "$DB/mgy_proteins_metadata.parquet" "$DB/mgy_proteins_pfam.parquet"
$RUN python3 "$HERE/neighbourhood.py" "$MSA" --clusters-metadata "$DB/mgy_clusters_metadata.parquet" \
  --metadata "$DB/mgy_proteins_metadata.parquet" --pfam "$DB/mgy_proteins_pfam.parquet" --clans "$CLANS" \
  --out "${OUT}_partners.tsv.gz" --proteins "${OUT}_proteins.tsv.gz" --cpus "${SLURM_CPUS_PER_TASK:-1}"
