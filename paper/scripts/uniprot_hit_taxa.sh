#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_uniprot_hit_taxa
#SBATCH --time=4:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=16
#SBATCH --output=mgnifams_uniprot_hit_taxa_%j.log
# Organisms of every UniProtKB hit of chosen families (R5 check of a taxonomy reading, e.g. the 27
# biome-exclusive families of 1.1: Skin -> Malassezia?, marine -> Bathycoccus?). Read-only.
# Same domtbls and rule as uniprot_hits.sh (E-values rescaled to Z = 149,810,139; sequence E and c-E <= 0.001);
# a protein counts once per family. Organism = OS= of the hit description, genus = its first word
# ("Candidatus" keeps the next word).
# Check: proteins per family must equal uniprot_hits_per_family.tsv.
# Output: the log (per family: proteins, genera with counts, top organisms) and uniprot_hit_taxa.tsv
# (family, subset, protein, rescaled E, OS, OX).
set -euo pipefail
: "${U:?set U, see paper/README.md}"
# default: the 27 biome-exclusive families of 1.1
FAMS=${FAMS:-250,1613,1925,2553,4142,4840,5310,5821,9300,12231,12410,13101,13842,14322,14937,16382,20748,20910,21864,22506,24217,24643,26779,28004,28224,28897,31302}
OUT=${OUT:-${SLURM_SUBMIT_DIR:-$PWD}/uniprot_hit_taxa.tsv}
TMPD=$(mktemp -d -p "${SLURM_SUBMIT_DIR:-$PWD}" uniprot_hit_taxa_parts.XXXX)
export TMPD FAMS

one() {  # $1 = domtbl.gz
  case $1 in */swissprot/*) ss=swissprot z=575503 ;; *) ss=trembl z=149234636 ;; esac
  zcat "$1" | awk -v s="$(awk -v z=$z 'BEGIN{printf "%.17g", 149810139/z}')" -v ss=$ss -v fams="$FAMS" 'BEGIN{OFS="\t"; OFMT="%.10g"; split(fams, a, ","); for (i in a) want[a[i]] = 1}
    /^#/ || !($4 in want) {next}
    { e = $7 * s; if (e > 1e-3 || $12 * s > 1e-3) next
      k = $4 SUBSEP $1; if (k in best && best[k] <= e) next
      best[k] = e; d = $23; for (i = 24; i <= NF; i++) d = d " " $i
      os = d; sub(/.* OS=/, "", os); sub(/ OX=.*/, "", os); if (d !~ / OS=/) os = "NA"
      ox = d; sub(/.* OX=/, "", ox); sub(/ .*/, "", ox); if (d !~ / OX=/) ox = "NA"
      row[k] = $4 OFS ss OFS $1 OFS e OFS os OFS ox }
    END { for (k in row) print row[k] }' > "$TMPD/${ss}_$(basename "$1").tsv"
}
export -f one

ls "$U"/swissprot/hmmsearch_mgnifams/*.domtbl.gz "$U"/trembl/hmmsearch_mgnifams/*.domtbl.gz | tee "$TMPD/files.txt" | wc -l
xargs -P "${SLURM_CPUS_PER_TASK:-4}" -n 1 bash -c 'one "$0"' < "$TMPD/files.txt"
ls "$TMPD"/*.domtbl.gz.tsv | wc -l

# a protein in several chunks of one subset cannot happen (chunks split the targets); keep the best E anyway
cat "$TMPD"/*.domtbl.gz.tsv | sort -t$'\t' -k1,1n -k2,2 -k3,3 -k4,4g | awk -F'\t' '!seen[$1 FS $2 FS $3]++' > "$OUT"
rm -r "$TMPD"
echo "rows: $(wc -l < "$OUT")  -> $OUT"

# per family and subset: proteins, genera (count), top 3 organisms
awk -F'\t' 'BEGIN{OFS="\t"}
  { k = $1 OFS $2; n[k]++
    split($5, w, " "); g = (w[1] == "Candidatus" ? w[1] " " w[2] : w[1]); G[k, g]++; if (!((k, g) in seenG)) { seenG[k, g] = 1; gl[k] = gl[k] "|" g }
    O[k, $5]++; if (!((k, $5) in seenO)) { seenO[k, $5] = 1; ol[k] = ol[k] "|" $5 } }
  END { print "family", "subset", "proteins", "genera (proteins)", "organisms (proteins)"
        for (k in n) { gs = ""; m = split(substr(gl[k], 2), a, "|"); for (i = 1; i <= m; i++) gs = gs (i > 1 ? "; " : "") a[i] " " G[k, a[i]]
                       os = ""; m = split(substr(ol[k], 2), b, "|"); for (i = 1; i <= m; i++) os = os (i > 1 ? "; " : "") b[i] " " O[k, b[i]]
                       print k, n[k], gs, os } }' "$OUT" | sort -t$'\t' -k1,1n -k2,2
echo "families asked: $(echo "$FAMS" | tr ',' '\n' | wc -l); with a hit: $(cut -f1 "$OUT" | sort -u | wc -l)"
