#!/usr/bin/env bash
set -euo pipefail

DEFAULT_IMAGE="/home/vangelis/Desktop/Tools/singularity/depot.galaxyproject.org-singularity-hmmer-3.4--hb6cb901_4.img"
DEFAULT_EXTRA_BIND="/home/vangelis/Desktop/Projects/modules/.nf-test/tests/a0341e6c2710429e7dd66d1356b43034/work"

usage() {
    cat <<'EOF'
Usage: bin/annotate_fasta_with_hmm.sh <hmm_file[.gz]> <fasta_file> <domtblout>

Environment overrides:
  HMMER_SINGULARITY_IMAGE  Singularity image path.
  HMMER_EXTRA_BIND         Extra host path to bind if it exists.
  HMMER_EVALUE             Sequence E-value threshold. Default: 0.001
  HMMER_DOM_EVALUE         Domain E-value threshold. Default: 0.001
  HMMER_CPU                Number of HMMER CPU workers. Default: 4
EOF
}

if [[ "$#" -ne 3 ]]; then
    usage
    exit 1
fi

HMM_FILE=$1
FASTA_FILE=$2
OUTPUT_FILE=$3

SINGULARITY_IMAGE=${HMMER_SINGULARITY_IMAGE:-$DEFAULT_IMAGE}
EXTRA_BIND=${HMMER_EXTRA_BIND:-$DEFAULT_EXTRA_BIND}
EVALUE=${HMMER_EVALUE:-0.001}
DOM_EVALUE=${HMMER_DOM_EVALUE:-0.001}
CPU=${HMMER_CPU:-4}

if [[ ! -f "$HMM_FILE" ]]; then
    echo "Error: HMM file does not exist: $HMM_FILE" >&2
    exit 1
fi

if [[ ! -f "$FASTA_FILE" ]]; then
    echo "Error: FASTA file does not exist: $FASTA_FILE" >&2
    exit 1
fi

if [[ ! -f "$SINGULARITY_IMAGE" ]]; then
    echo "Error: Singularity image does not exist: $SINGULARITY_IMAGE" >&2
    exit 1
fi

HMM_ABS=$(realpath "$HMM_FILE")
FASTA_ABS=$(realpath "$FASTA_FILE")
OUTPUT_ABS=$(realpath -m "$OUTPUT_FILE")
OUTPUT_DIR=$(dirname "$OUTPUT_ABS")
mkdir -p "$OUTPUT_DIR"

TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

HMM_FOR_SEARCH=$HMM_ABS
if [[ "$HMM_FILE" == *.gz ]]; then
    HMM_FOR_SEARCH="$TMP_DIR/$(basename "${HMM_FILE%.gz}")"
    gzip -dc "$HMM_ABS" > "$HMM_FOR_SEARCH"
fi

BIND_ARGS=(
    -B "$(dirname "$HMM_ABS")"
    -B "$(dirname "$FASTA_ABS")"
    -B "$OUTPUT_DIR"
    -B "$TMP_DIR"
)

if [[ -d "$EXTRA_BIND" ]]; then
    BIND_ARGS+=(-B "$EXTRA_BIND")
fi

singularity exec --no-home --pid "${BIND_ARGS[@]}" "$SINGULARITY_IMAGE" \
    hmmsearch \
    -E "$EVALUE" \
    --domE "$DOM_EVALUE" \
    --cpu "$CPU" \
    --domtblout "$OUTPUT_ABS" \
    "$HMM_FOR_SEARCH" \
    "$FASTA_ABS" \
    > "$OUTPUT_ABS.log"

echo "Annotation completed. Results saved in $OUTPUT_ABS"
