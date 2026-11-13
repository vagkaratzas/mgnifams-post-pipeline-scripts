#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_tm_classes
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=mgnifams_tm_classes_%j.log
# TM / disordered classes for 1.0 and 1.1 (Results catalogue paragraph 2), Methods rules. Read-only.
# Copy tm_classes.py next to this script, then on the cluster: sbatch tm_classes.sh   (or bash, minutes)
# LISTS_OUT=dir also writes the Table 2 category lists (dir/1.0, dir/1.1).
# The 1.0 block must reproduce 5,186 alpha / 372 beta / 989 disordered, +510, 1,777 before 1.1 is used.
set -euo pipefail
: "${V:?set V, see paper/README.md}"
OLD=${OLD:-$V/1.0/output/table_data/mgnifam.csv}
NEW=${NEW:-$V/1.1/output/table_data/mgnifam.csv}
IDS=${IDS:-$V/1.1/output/update_families/families.tsv.gz}
PY=${PY:-${SLURM_SUBMIT_DIR:-$(dirname "$0")}/tm_classes.py}  # sbatch runs a spool copy, so $0 is not here
ls -la "$OLD" "$NEW" "$IDS"
head -1 "$OLD" | cut -c1-400; head -1 "$NEW" | cut -c1-400
python3 "$PY" --old "$OLD" --new "$NEW" --ids "$IDS" ${LISTS_OUT:+--lists-out "$LISTS_OUT"}
