#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_survival
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --output=mgnifams_survival_%j.log
# Case-study survival check: do the paper's case-study families keep their subject in 1.1?
# Per family: 1.1 status, full_size 1.0 -> 1.1, representative 1.0 -> 1.1, biome leaves 1.0 -> 1.1.
# Biome breadth as in post-pipeline bin/biome_analysis.py: leaves = `ids` never used as `parents`; exclusive = 1 leaf.
# Read-only. Run on the cluster: sbatch survival.sh (or bash survival.sh; a few minutes).
: "${V:?set V, see paper/README.md}"
OUT=${OUT:-$V/1.1/output}
P10=${P10:-$V/1.0/output_db/families.tsv.gz}
P11=${P11:-$OUT/update_families/families.tsv.gz}
B10=${B10:-$V/1.0/output_db/biome_results}          # 1.0 per-family biome CSVs (override if elsewhere)
B11=${B11:-$OUT/biome_results}
CASES=${CASES:-"25 243 32058 12274 12868 689 1107 20804 28802 4542 26748 8562 14909 26162 29438 30944"}
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
echo "== biome result dirs under 1.0 (in case B10 is wrong)"; find "$V/1.0" -maxdepth 4 -type d -name 'biome*' 2>/dev/null
ls -d "$B10" "$B11"; ls "$B11" | head -2

# family<TAB>number of leaves<TAB>leaf (when exactly one), from a dir of <family>.csv
leaves(){
  find "$1" -maxdepth 1 -name '*.csv' -print0 | xargs -0 awk -F, '
    function flush(){ if (fam == "") return; n = 0; one = ""
      for (i in id) if (!(i in par)) { n++; one = i }
      print fam "\t" n "\t" (n == 1 ? one : "-"); delete id; delete par }
    FNR == 1 { flush(); fam = FILENAME; sub(/.*\//, "", fam); sub(/\.csv$/, "", fam); next }
    NF != 4 { bad++ }
    { id[$1] = 1; if ($3 != "") par[$3] = 1 }
    END { flush(); if (bad) print "WARN rows without 4 fields: " bad > "/dev/stderr" }'
}
[ -d "$B10" ] && leaves "$B10" > "$TMP/b10.tsv"; touch "$TMP/b10.tsv"
leaves "$B11" > "$TMP/b11.tsv"
echo "== families with biome results: 1.0 $(wc -l < "$TMP/b10.tsv"), 1.1 $(wc -l < "$TMP/b11.tsv")"
echo "== exclusive (1 leaf): 1.0 $(awk -F'\t' '$2==1' "$TMP/b10.tsv" | wc -l) (paper: 24), 1.1 $(awk -F'\t' '$2==1' "$TMP/b11.tsv" | wc -l)"
echo "== 1.1 exclusive families by leaf"; awk -F'\t' '$2==1{print $3}' "$TMP/b11.tsv" | sort | uniq -c | sort -rn
EXCL=$(awk -F'\t' '$2==1{print $1}' "$TMP/b10.tsv" | sort -n | tr '\n' ' ')

echo "== per family: id group status reason | size 1.0 1.1 | rep 1.0 -> 1.1 (same MGYP? same region?) | leaves 1.0 1.1 | 1.1 leaf if exclusive"
for f in $CASES; do echo "$f case"; done > "$TMP/list"
for f in $EXCL; do echo "$f excl10"; done >> "$TMP/list"
awk -v P10="$P10" -v P11="$P11" -v B10="$TMP/b10.tsv" -v B11="$TMP/b11.tsv" '
function load(cmd, A,   l, a){ while ((cmd | getline l) > 0) { split(l, a, "\t"); A[a[1]] = l } close(cmd) }
function leafs(file, N, L,   l, a){ while ((getline l < file) > 0) { split(l, a, "\t"); N[a[1]] = a[2]; L[a[1]] = a[3] } close(file) }
BEGIN { load("zcat " P10, F10); load("zcat " P11, F11); leafs(B10, n10, l10); leafs(B11, n11, l11) }
{ f = $1; if (!(f in F10)) { print f, $2, "NOT IN 1.0"; next }
  split(F10[f], a, "\t"); split(F11[f], b, "\t")
  sm = (a[11] == b[11]) ? "same" : "other"; sr = (sm == "same" && a[12] == b[12]) ? "same" : "other"
  printf "%s %s %s %s | %s %s | %s/%s -> %s/%s (%s MGYP, %s region) | %s %s | %s\n", f, $2, b[2], (b[3] == "" ? "-" : b[3]),
    a[9], b[9], a[11], a[12], b[11], b[12], sm, sr, (f in n10 ? n10[f] : "na"), (f in n11 ? n11[f] : "na"), (f in l11 ? l11[f] : "-") }
' "$TMP/list"
