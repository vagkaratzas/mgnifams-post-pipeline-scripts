#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_contig_lookup
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --output=mgnifams_contig_lookup_%j.log
# source contigs (longest per member, ERZ, contig name, pipeline version) of the biome-exclusive families
# with a UniProtKB hit, release 1.1, plus every gene call on those contigs. Read-only.
# Next to this script: contig_lookup.py, neighbourhood.py, architectures.py
#, classify_multi_hits.py. Then on the cluster: sbatch contig_lookup.sh
# Writes contig_lookup_<job>_members.tsv and contig_lookup_<job>_genes.tsv.gz in the submit dir.
set -euo pipefail
: "${V:?set V, see paper/README.md}"
# the 21 of the 27 biome-exclusive families with a UniProtKB hit
FAMS=${FAMS:-"250 1613 1925 2553 4142 4840 9300 12231 12410 13101 14322 16382 20748 20910 24217 24643 26779 28004 28224 28897 31302"}
: "${DB:?set DB, see paper/README.md}"
HERE=${SLURM_SUBMIT_DIR:-$(dirname "$0")}  # sbatch runs a spool copy, so $0 is not here
: "${CACHE:?set CACHE, see paper/README.md}"
IMG=${IMG:-$(ls "$CACHE"/*e415e86f361373d5bf62a587c79daada9fa797ce076c764f1ccadb32ba2c3861* 2>/dev/null | head -1 || true)}
RUN=${RUN:-singularity exec ${BIND:+-B $BIND} $IMG}   # RUN=env for a local test
OUT=$HERE/contig_lookup_${SLURM_JOB_ID:-local}
MSAS=$(for f in $FAMS; do echo "$V/1.1/output/update_families/full_msa/$f.sto.gz"; done)
echo "== image: ${IMG:-NOT FOUND (set IMG=)}"
md5sum "$HERE/contig_lookup.py" "$HERE/neighbourhood.py" "$HERE/architectures.py" "$HERE/classify_multi_hits.py"
ls -la $MSAS "$DB/mgy_clusters_metadata.parquet" "$DB/mgy_proteins_metadata.parquet"
$RUN python3 "$HERE/contig_lookup.py" $MSAS --clusters-metadata "$DB/mgy_clusters_metadata.parquet" \
  --metadata "$DB/mgy_proteins_metadata.parquet" --out "${OUT}_members.tsv" --genes "${OUT}_genes.tsv.gz" \
  --cpus "${SLURM_CPUS_PER_TASK:-1}"
