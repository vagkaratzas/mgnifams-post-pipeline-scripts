#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_annotation_evidence
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=mgnifams_annotation_evidence_%j.log
# Annotation evidence counts, 1.0 (must reproduce the published paragraph 6 numbers) and 1.1 (Results catalogue
# paragraph 6, Table S1). Read-only. Needs the Pfam 36.0 and 38.1 family lists next to this script:
#   wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam36.0/Pfam-A.clans.tsv.gz -O Pfam36.0.clans.tsv.gz
#   wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.1/Pfam-A.clans.tsv.gz -O Pfam38.1.clans.tsv.gz
# then: PY=$PWD/annotation_evidence.py sbatch annotation_evidence.sh
# The DUF-only line (2026-11-13, Discussion paragraph 1) takes Pfam names (column 4) from Pfam38.2.clans.tsv.gz, then 38.1:
#   wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/Pfam-A.clans.tsv.gz -O Pfam38.2.clans.tsv.gz
set -euo pipefail
: "${V:?set V, see paper/README.md}"
HERE=${SLURM_SUBMIT_DIR:-$(dirname "$0")}  # sbatch runs a spool copy, so $0 is not here
PY=${PY:-$HERE/annotation_evidence.py}
wc -l "$V"/1.0/output/table_data/mgnifam_model_pfams.csv "$V"/1.1/output/table_data/mgnifam_{pfams,funfams}.csv
python3 "$PY" --old "$V/1.0/output" --new "$V/1.1/output" --reann "$V/1.1/output_annotate/annotate_not_updated" \
    --ids "$V/1.1/output/update_families/families.tsv.gz" --mask10 "$HERE/Pfam36.0.clans.tsv.gz" --mask11 "$HERE/Pfam38.1.clans.tsv.gz" \
    --names "$HERE/Pfam38.2.clans.tsv.gz" "$HERE/Pfam38.1.clans.tsv.gz"
