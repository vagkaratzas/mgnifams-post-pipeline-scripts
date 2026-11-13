#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_catalogue_probe
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=mgnifams_catalogue_probe_%j.log
# Unannotated vs annotated families, 1.0 and 1.1 (Results catalogue unannotated and pLDDT paragraphs, Figure 3 bins).
# Read-only. Copy catalogue_probe.py AND annotation_evidence.py (it imports it) next to this script, then:
#   sbatch catalogue_probe.sh   (DUMP=$PWD sbatch … also writes catalogue_1.0/1.1.csv.gz for figure3.py)
set -euo pipefail
: "${V:?set V, see paper/README.md}"
HERE=${SLURM_SUBMIT_DIR:-$(dirname "$0")}  # sbatch runs a spool copy, so $0 is not here
PY=${PY:-$HERE/catalogue_probe.py}
python3 "$PY" ${DUMP:+--dump "$DUMP"} --old "$V/1.0/output" --new "$V/1.1/output" --reann "$V/1.1/output_annotate/annotate_not_updated" \
    --ids "$V/1.1/output/update_families/families.tsv.gz" --prev "$V/1.0/output_db/families.tsv.gz" --pfam_median 1089
