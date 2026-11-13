#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_pfam_per_family
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=32
#SBATCH --output=mgnifams_pfam_per_family_%j.log
# Pfam domains per family and per clan, release 1.1, from the pipeline's query_results/<family>.tsv. Read-only.
# Next to this script: pfam_per_family.py, architectures.py, classify_multi_hits.py, and
#   wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.1/Pfam-A.clans.tsv.gz -O Pfam38.1.clans.tsv.gz
#   wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnifams/releases/1.1/clans.tsv.gz -O mgnifams_1.1_clans.tsv.gz
# Then on the cluster: QR="<work task dir>/query_results ..." sbatch pfam_per_family.sh. Writes pfam_per_family_<job>.tsv.gz.
set -euo pipefail
: "${QR:?set QR to the query_results dir(s), space-separated}"
HERE=${SLURM_SUBMIT_DIR:-$(dirname "$0")}  # sbatch runs a spool copy, so $0 is not here
CLANS=${CLANS:-$HERE/Pfam38.1.clans.tsv.gz}
MCLANS=${MCLANS:-$HERE/mgnifams_1.1_clans.tsv.gz}
: "${CACHE:?set CACHE, see paper/README.md}"
IMG=${IMG:-$(ls "$CACHE"/*e415e86f361373d5bf62a587c79daada9fa797ce076c764f1ccadb32ba2c3861* 2>/dev/null | head -1)}
RUN=${RUN:-singularity exec ${BIND:+-B $BIND} $IMG}   # RUN=env for a local test
echo "== image: ${IMG:-NOT FOUND (set IMG=)}"
md5sum "$HERE/pfam_per_family.py" "$HERE/architectures.py" "$HERE/classify_multi_hits.py" "$CLANS" "$MCLANS"
for d in $QR; do echo "$d $(ls "$d" | wc -l) files"; done
$RUN python3 "$HERE/pfam_per_family.py" $QR --clans "$CLANS" --mgnifam-clans "$MCLANS" \
  --out "$HERE/pfam_per_family_${SLURM_JOB_ID:-local}.tsv.gz" --cpus "${SLURM_CPUS_PER_TASK:-1}"
