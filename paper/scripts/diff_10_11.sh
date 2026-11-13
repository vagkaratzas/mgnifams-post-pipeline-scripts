#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_diff
#SBATCH --time=24:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --output=mgnifams_diff_%j.log
# 1.0 vs 1.1 diff on MGYPs searched in both releases. Read-only on
# the release outputs; writes only under $W. Submit on the cluster:  W=<scratch dir>/mgnifams_diff sbatch diff_10_11.sh
# Stages skip when their output exists, so a rerun only redoes the cheap report.
set -euo pipefail
: "${W:?set W to a scratch dir}"
: "${V:?set V, see paper/README.md}"
F10=${F10:-$V/1.0/output/generate_families/families}
FAM=${FAM:-$V/1.1/output/update_families}
UR10=${UR10:-$V/1.0/output/mgnifams_v2.fa}          # MGnify-UR 2024_04, the 1.0 searched set
UR11=${UR11:-$FAM/mgnifams_update.faa}              # the 1.1 searched set
P10=${P10:-$V/1.0/output_db/families.tsv.gz}
P11=${P11:-$FAM/families.tsv.gz}
CONV=${CONV:-$W/converged_flag.txt}  # from family_metadata.csv flag; 1.0 converged_families.txt has stale IDs (see SKILL.md)
mkdir -p "$W/tmp"; cd "$W"
[ -s "$CONV" ] || awk -F, '$NF=="True"{print $1}' "$F10/family_metadata.csv" > "$CONV"
export LC_ALL=C TMPDIR=$W/tmp
S="sort -S 80G --parallel=8 -T $W/tmp"
ls -la "$F10/refined_families.tsv" "$FAM/refined_families.tsv" "$UR10" "$UR11" "$P10" "$P11" "$CONV"
head -2 "$F10/refined_families.tsv"

# Stage 1: unique MGYPs of each searched set (headers are mgyp or mgyp_start_end), and their intersection U
ids(){ grep '^>' "$1" | cut -c2- | cut -d_ -f1 | $S -u > "$2.tmp" && mv "$2.tmp" "$2"; }
if [ ! -s u.txt ]; then
  [ -s u10.txt ] || ids "$UR10" u10.txt & a=$!
  [ -s u11.txt ] || ids "$UR11" u11.txt & b=$!
  wait $a; wait $b
  comm -12 u10.txt u11.txt > u.txt.tmp && mv u.txt.tmp u.txt
fi

# Stage 2: per family, member rows (family<TAB>mgyp/start-end) in total (t*) and with the MGYP in U (n*)
mem(){
  awk -F'\t' -v OFS='\t' -v tot="$2.tmp" '{split($2,a,"/"); print a[1],$1; t[$1]++} END{for(f in t) print f,t[f] > tot}' "$1" \
  | $S -t "$(printf '\t')" -k1,1 | join -t "$(printf '\t')" u.txt - \
  | awk -F'\t' '{n[$2]++} END{for(f in n) print f"\t"n[f]}' > "$3.tmp"
  mv "$2.tmp" "$2"; mv "$3.tmp" "$3"
}
[ -s n10.tsv ] || mem "$F10/refined_families.tsv" t10.tsv n10.tsv & a=$!
[ -s n11.tsv ] || mem "$FAM/refined_families.tsv" t11.tsv n11.tsv & b=$!
wait $a; wait $b

# Stage 3: report. Groups: converged in 1.0 vs stopped after round 3 (released HMM has one more round)
echo "== searched MGYPs: 1.0 $(wc -l < u10.txt), 1.1 $(wc -l < u11.txt), both (U) $(wc -l < u.txt)"
awk -v CONV="$CONV" -v P10="$P10" -v P11="$P11" -v RAT="$W/ratio" '
function pct(a, b){ return b ? 100 * a / b : 0 }
function strip(r){ sub(/^MGYP0*/, "", r); return r }
function load(file, arr,   l, a){ while ((getline l < file) > 0) { split(l, a, "\t"); arr[a[1]] = a[2] } close(file) }
BEGIN {
  FS = "\t"
  while ((getline l < CONV) > 0) C[l] = 1
  close(CONV)
  load("n10.tsv", n10); load("n11.tsv", n11); load("t10.tsv", t10); load("t11.tsv", t11)
  c = "zcat " P10
  while ((c | getline l) > 0) { split(l, a, "\t"); if (a[1] == "family_id") continue
    f = a[1]; nf10++; fs10[f] = a[9]; r10[f] = strip(a[11]); g10[f] = a[12]; R[r10[f]] = 1 }
  close(c)
  c = "zcat " P11
  while ((c | getline l) > 0) { split(l, a, "\t"); if (a[1] == "family_id") continue
    f = a[1]; st[f] = a[2]; fs11[f] = a[9]; r11[f] = strip(a[11]); g11[f] = a[12]; R[r11[f]] = 1 }
  close(c)
}
($1 in R) { inU[$1] = 1 }
END {
  for (f in fs10) { if ((f in C)) nconv++; if (t10[f] + 0 != fs10[f] + 0) bad10++ }
  printf "check: 1.0 families %d, converged %d; member rows != full_size: 1.0 %d", nf10, nconv, bad10
  for (f in st) if (st[f] == "updated" && t11[f] + 0 != fs11[f] + 0) bad11++
  printf ", 1.1 %d\n", bad11
  for (f in st) {
    if (st[f] != "updated") continue
    k = (f in C) ? "converged" : "stopped"
    for (i = 1; i <= 2; i++) { g = (i == 1) ? "all" : k
      N[g]++; m0 = n10[f] + 0; m1 = n11[f] + 0; T10[g] += t10[f]; T11[g] += t11[f]; S10[g] += m0; S11[g] += m1
      if (m0 != m1) X[g]++
      if (m1 - m0 > 0.1 * m0 || m0 - m1 > 0.1 * m0) X10[g]++
      if (m0 == 0) Z[g]++
      if (r10[f] == r11[f]) SAME[g]++
      if ((r10[f] in inU) && (r11[f] in inU)) { B[g]++; if (r10[f] != r11[f]) Y[g]++
        else if (g10[f] != "-" && g10[f] != g11[f]) REG[g]++ }
    }
    if (m0 > 0) print m1 / m0, k > RAT
  }
  for (i = 1; i <= 3; i++) { g = (i == 1) ? "all" : (i == 2) ? "converged" : "stopped"
    printf "== %s: updated %d\n", g, N[g]
    printf "  member rows in U: 1.0 %d of %d (%.1f%%), 1.1 %d of %d (%.1f%%)\n", S10[g], T10[g], pct(S10[g],T10[g]), S11[g], T11[g], pct(S11[g],T11[g])
    printf "  x: size in U differs %d (%.1f%%), by >10%% %d (%.1f%%); no 1.0 rows in U %d\n", X[g], pct(X[g],N[g]), X10[g], pct(X10[g],N[g]), Z[g]
    printf "  y: reps both in U %d, rep MGYP differs %d (%.1f%%); same MGYP, other region %d (1.0 '-' regions not counted)\n", B[g], Y[g], pct(Y[g],B[g]), REG[g]
    printf "  same rep MGYP, all updated: %d (%.1f%%)\n", SAME[g], pct(SAME[g],N[g])
  }
}' u.txt
q(){ sort -g | awk '{a[NR]=$1} END{printf "n %d min %s p10 %s p50 %s p90 %s max %s\n",NR,a[1],a[int(NR*.1)+1],a[int(NR*.5)+1],a[int(NR*.9)+1],a[NR]}'; }
for k in converged stopped; do echo "== size in U, 1.1 / 1.0, $k"; awk -v k=$k '$2==k{print $1}' ratio | q; done
