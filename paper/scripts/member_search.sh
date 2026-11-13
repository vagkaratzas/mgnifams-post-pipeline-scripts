#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_member_search
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --output=mgnifams_member_search_%j.log
# do the member regions of a MGnifam match what its representative matched? Read-only on the release.
# Targets: the family's member regions (full MSA rows, ungapped), not whole proteins.
# Queries: every FunFam of one CATH superfamily (FunFam 4.3.0) + one Pfam entry (Pfam-A 38.2).
# Rule of the representative annotation (pipeline ANNOTATE_REPS / FILTER_EXPORT_DOMTBL): hmmsearch -E 0.001
# -Z 35459 (the release's family count, so each region is judged as a representative would be), then a domain
# counts if i-Evalue <= 0.001 and envelope >= 0.8 x model length. Also printed: the same without the length rule.
# Check: the representative region must reproduce its card (MGnifam 25: 29 FunFams, best 2.3e-20, PF00528 3.9e-07);
# E-values here are x 35459/34419 the card values: the 1.1 run searched the 34,419 updated reps in one task, no -Z.
# Expected for MGnifam 25: 29 FunFams, 2.4e-20, 4.1e-07.
#   sbatch member_search.sh            (MGnifam 25; FAM, SF, PF, REP to change)
#   RUN=env MSA=x.sto.gz FUNFAMS=... PFAM=... bash member_search.sh   (local test)
set -euo pipefail
W=${W:-${SLURM_SUBMIT_DIR:-$PWD}}
: "${P:?set P, see paper/README.md}"
: "${REF:?set REF, see paper/README.md}"
FAM=${FAM:-25}
MSA=${MSA:-$P/1.1/output/update_families/full_msa/$FAM.sto.gz}
FUNFAMS=${FUNFAMS:-$REF/funfams/funfam-hmm3-v4_3_0.lib.gz}
PFAM=${PFAM:-$REF/pfam/Pfam-A.38.2.hmm.gz}
SF=${SF:-1.10.3720.10}
PF=${PF:-PF00528}
REP=${REP:-MGYP000395001733/16-290}
REPN=$(sed "s/^MGYP0*//" <<< "$REP")   # full MSA rows name proteins by bare number (395001733/16-290)
Z=${Z:-35459}
E=0.001
: "${IMG:?set IMG, see paper/README.md}"
RUN=${RUN:-singularity exec ${BIND:+-B $BIND} $IMG}   # RUN=env for a local test
CPUS=${SLURM_CPUS_PER_TASK:-4}
NJ=$(( CPUS >= 8 ? CPUS / 4 : 1 ))   # parallel hmmsearch runs, 4 threads each
D=$W/member_search_$FAM
mkdir -p "$D"
$RUN hmmsearch -h | sed -n 2p

# 1. region sequences (Pfam-format Stockholm: one "name aligned" line per row)
zcat "$MSA" | awk '/^#/ || /^\/\// || NF < 2 {next} {s = $2; gsub(/[-.]/, "", s); print ">" $1; print toupper(s)}' > "$D/regions.fa"
echo "regions: $(grep -c '^>' "$D/regions.fa")"

# 2. queries: the superfamily's FunFams (by NAME) and the Pfam entry (by ACC)
pick() { zcat "$1" | awk -v re="$2" -v tag="$3" '/^HMMER3/ {b = ""; k = 0} {b = b $0 "\n"}
  $1 == tag && $2 ~ re {k = 1} /^\/\// {if (k) printf "%s", b}'; }
{ pick "$FUNFAMS" "^${SF//./[.]}-FF-" NAME; pick "$PFAM" "^$PF[.]" ACC; } > "$D/queries.hmm"
echo "queries: $(grep -c '^NAME' "$D/queries.hmm") ($(grep -c "^NAME  ${SF//./[.]}-FF-" "$D/queries.hmm") FunFams of $SF + $PF)"

# 3. search NJ query chunks in parallel; each domtbl is reduced as it streams (the full table would be tens of GB).
# Per region and mode ("rule"; "any" = without the length rule): FunFams passed, best FunFam i-E, best Pfam i-E.
# Distinct FunFams counted against the previous row: a (region, model) pair's domain rows are consecutive.
agg() {
  set -o pipefail
  $RUN hmmsearch --cpu 4 -E $E -Z $Z --noali -o /dev/null --domtblout /dev/stdout "$1" "$D/regions.fa" | awk -v e=$E -v pf="$PF" '
    BEGIN {OFS = "\t"; OFMT = "%.10g"}
    /^#/ {next}
    { r = $1; ff = ($5 !~ "^" pf "[.]"); ie = $13 + 0; if (ie > e) next
      rule = ($21 - $20 + 1 >= 0.8 * $6)
      for (m = 0; m <= rule; m++) {
        k = m SUBSEP r; rs[k]
        if (ff) { if (last[m] != r SUBSEP $4) { nff[k]++; last[m] = r SUBSEP $4 }; if (!(k in bff) || ie < bff[k]) bff[k] = ie }
        else if (!(k in bpf) || ie < bpf[k]) bpf[k] = ie } }
    END { for (k in rs) { split(k, a, SUBSEP)
      print (a[1] ? "rule" : "any"), a[2], nff[k] + 0, (k in bff ? bff[k] : "NA"), (k in bpf ? bpf[k] : "NA") } }' > "$1.tsv"
}
export -f agg; export RUN E Z D PF
rm -f "$D"/chunk_*
awk -v n=$NJ -v d="$D" '/^HMMER3/ {i++} {print > (d "/chunk_" (i % n) ".hmm")}' "$D/queries.hmm"
ls "$D"/chunk_*.hmm | xargs -P $NJ -n 1 bash -c 'agg "$0"'

# 4. merge the chunks: sum the FunFam counts, keep the smallest i-E
cat "$D"/chunk_*.hmm.tsv | awk -F'\t' 'BEGIN {OFS = "\t"; OFMT = "%.10g"}
  { k = $1 OFS $2; n[k] += $3
    if ($4 != "NA" && (!(k in f) || $4 + 0 < f[k] + 0)) f[k] = $4
    if ($5 != "NA" && (!(k in q) || $5 + 0 < q[k] + 0)) q[k] = $5 }
  END { for (k in n) print k, n[k], (k in f ? f[k] : "NA"), (k in q ? q[k] : "NA") }' | sort > "$D/region_hits.tsv"

# 5. fractions of regions and of member proteins (a protein counts if any of its regions does)
echo "== representative $REP (mode, region, FunFams, best FunFam i-E, $PF i-E)"; grep -F "	$REPN	" "$D/region_hits.tsv" || echo "no hit"
awk -F'\t' -v pf="$PF" 'BEGIN {OFMT = "%.10g"}
  FNR == NR { if (/^>/) { r = substr($0, 2); nr++; p = r; sub(/\/[^\/]*$/, "", p); if (!(p in prot)) { prot[p]; np++ } } next }
  { p = $2; sub(/\/[^\/]*$/, "", p); f = ($3 > 0); q = ($5 != "NA")
    R[$1, "ff"] += f; R[$1, "pf"] += q; R[$1, "either"] += (f || q); R[$1, "both"] += (f && q)
    if (f) Pf[$1, p]; if (q) Pq[$1, p]; Pe[$1, p] }
  END { printf "regions %d, member proteins %d\n", nr, np
    for (m = 1; m >= 0; m--) { mode = m ? "rule" : "any"
      for (k in Pe) { split(k, a, SUBSEP); if (a[1] != mode) continue; p = a[2]
        f = ((mode, p) in Pf); q = ((mode, p) in Pq); PF_[mode] += f; PQ[mode] += q; PE[mode] += (f || q); PB[mode] += (f && q) }
      printf "%s: regions FunFam %d (%.10g%%), %s %d (%.10g%%), either %d (%.10g%%), both %d (%.10g%%)\n", mode,
        R[mode, "ff"], 100 * R[mode, "ff"] / nr, pf, R[mode, "pf"], 100 * R[mode, "pf"] / nr,
        R[mode, "either"], 100 * R[mode, "either"] / nr, R[mode, "both"], 100 * R[mode, "both"] / nr
      printf "%s: proteins FunFam %d (%.10g%%), %s %d (%.10g%%), either %d (%.10g%%), both %d (%.10g%%)\n", mode,
        PF_[mode], 100 * PF_[mode] / np, pf, PQ[mode], 100 * PQ[mode] / np, PE[mode], 100 * PE[mode] / np, PB[mode], 100 * PB[mode] / np } }' \
  "$D/regions.fa" "$D/region_hits.tsv"
md5sum "$D/region_hits.tsv"
