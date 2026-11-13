#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_biome_breadth
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=mgnifams_biome_breadth_%j.log
# Biome breadth, 1.0 then 1.1 (MGnifam 25 paragraph, Results "Biome analysis"). Read-only.
# Copy biome_breadth.py and post-pipeline bin/biome_analysis.py next to this script, then: sbatch biome_breadth.sh
# The 1.0 block must reproduce 121 leaf paths, median 56, max 85, 24 one-leaf families before 1.1 is used
# (done 2026-11-07; ONLY_NEW=1 skips it). 1.1 = update delta + the not-updated families from the 1.0 DB.
# DUMP=1 also writes leaves_1.1.tsv.gz (id, n_leaves) next to the log, for breadth_size.py (wormhole it back).
set -euo pipefail
: "${V:?set V, see paper/README.md}"
OLD=${OLD:-$V/1.0/output_db/mgnifams_test_copy.sqlite3}
NEW=${NEW:-$V/1.1/output_delta_update/db/mgnifams_1_1_update.sqlite3}
D=${SLURM_SUBMIT_DIR:-$(dirname "$0")}  # sbatch runs a spool copy, so $0 is not here
ls -la "$OLD" "$NEW" "$D/biome_breadth.py" "$D/biome_analysis.py"
python3 -c "import pandas, numpy; print('pandas', pandas.__version__)"
[ -n "${ONLY_NEW:-}" ] || python3 "$D/biome_breadth.py" "$OLD" --family 25
python3 "$D/biome_breadth.py" "$NEW" --base "$OLD" --family 25 ${DUMP:+--dump "$D/leaves_1.1.tsv.gz"}
