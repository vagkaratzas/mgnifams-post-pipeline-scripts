#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_blastx_uniref
#SBATCH --time=24:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=32
#SBATCH --output=mgnifams_blastx_uniref_%j.log
# DIAMOND blastx of the members' source contigs (fetch_contigs.py) against UniRef90 of the same
# UniProt release as the paper's UniProtKB search (2026_02), so each gene's organism is its best hit in all of UniRef90.
# Range culling keeps the top hits of each contig region (a contig carries many genes). Read-only on the DB.
# Next to this script: contigs.fasta. Container (login node, once):
#   singularity pull docker://quay.io/biocontainers/diamond:2.1.8--h43eeafb_0
# Then: sbatch blastx_uniref.sh    Writes blastx_uniref_<job>.tsv.gz in the submit dir.
set -euo pipefail
HERE=${SLURM_SUBMIT_DIR:-$(dirname "$0")}
: "${D:?set D, see paper/README.md}"
: "${CACHE:?set CACHE, see paper/README.md}"
IMG=${IMG:-$(ls "$CACHE"/diamond_2.1.8--h43eeafb_0.sif "$CACHE"/*diamond*2.1.8* 2>/dev/null | head -1 || true)}
RUN=${RUN:-singularity exec ${BIND:+-B $BIND} $IMG}
Q=${Q:-$HERE/contigs.fasta}
OUT=$HERE/blastx_uniref_${SLURM_JOB_ID:-local}.tsv
echo "== image: ${IMG:-NOT FOUND}"; ls -la "$D"; cat "$D"/README* "$D"/VERSION.txt 2>/dev/null || true
DB=${DB:-$(ls "$D"/*.dmnd | head -1)}
md5sum "$Q"; grep -c '>' "$Q"
$RUN diamond dbinfo -d "$DB"
$RUN diamond blastx -q "$Q" -d "$DB" --sensitive -F 15 --range-culling --top 10 --evalue 0.001 \
  --threads "${SLURM_CPUS_PER_TASK:-1}" -b 4 -c 1 -o "$OUT" \
  --outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qframe slen qlen stitle
wc -l "$OUT"; gzip -f "$OUT"; ls -la "$OUT.gz"
