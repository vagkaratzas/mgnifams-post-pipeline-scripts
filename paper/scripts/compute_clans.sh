#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_compute_clans
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=mgnifams_compute_clans_%j.log
# Release 1.1 clans: `mgnifam compute_clans` 4.1.0 over all 35,459 families, i.e. the 1.1 members of the
# updated families plus members_1.0_not_updated.tsv; bare.fa gives the whole-protein members' lengths. Both come from
# clan_inputs.sh (run first, in the same dir). Default thresholds: clan edge residue Jaccard >= 0.5.
# Sizing (guide "Scale", 2.61 B rows, 35,459 families): pair array 7.5 GB (memory-mapped in scratch) + 512 MB per
# worker; scratch ~89 GB under $W/tmp. Container (login node, once):
#   singularity pull docker://quay.io/biocontainers/mgnifam:4.1.0--pyhdfd78af_0
# Then, in the clan_inputs.sh dir:  sbatch compute_clans.sh    Writes $W/output/1.1_clan*; read-only on the release.
# Release 1.0 alone (diagnosis only, not published), in its own dir:
#   CHUNK=1.0 MEMBERS=$P/1.0/output/generate_families/families/refined_families.tsv FASTA=$P/1.0/output/mgnifams_v2.fa sbatch compute_clans.sh
set -euo pipefail
W=${W:-${SLURM_SUBMIT_DIR:-$PWD}}
: "${P:?set P, see paper/README.md}"
NEW=$P/1.1/output/update_families
: "${CACHE:?set CACHE, see paper/README.md}"
IMG=${IMG:-$CACHE/mgnifam_4.1.0--pyhdfd78af_0.sif}
RUN=${RUN:-singularity exec ${BIND:+-B $BIND} $IMG}   # RUN=env for a local test
CPUS=${SLURM_CPUS_PER_TASK:-4}
CHUNK=${CHUNK:-1.1}
cd "$W"; mkdir -p tmp
if [ -z "${MEMBERS:-}" ]; then
  [ -s members_1.0_not_updated.tsv ] && [ -s bare.fa ] || { echo "run clan_inputs.sh first" >&2; exit 2; }
  MEMBERS="$NEW/refined_families.tsv members_1.0_not_updated.tsv"; FASTA=bare.fa
fi
if [ -s output/${CHUNK}_clans_summary.json ]; then echo "output/${CHUNK}_clans_summary.json exists: done" >&2; exit 0; fi
echo "== image: $IMG"; $RUN mgnifam --version || true
ls -la $MEMBERS "$FASTA"; df -h "$W" | tail -1
$RUN mgnifam compute_clans --members $MEMBERS --fasta_file "$FASTA" \
  -n "$CHUNK" -p "$CPUS" --tmp_dir "$W/tmp" --output_dir output
ls -la output; cat output/${CHUNK}_clans_summary.json
md5sum output/${CHUNK}_clan_edges.csv output/${CHUNK}_clans.csv
