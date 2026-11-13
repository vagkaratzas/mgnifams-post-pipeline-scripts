#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_biome_probe
#SBATCH --time=03:00:00
#SBATCH --mem=8G
#SBATCH --output=mgnifams_biome_probe_%j.log
# Why 1.1 has 147 biome-exclusive families vs 24 in 1.0 (survival.sh), though families grew.
# Hypothesis: QUERY_MGNPROTEIN_DB reads the MGnify proteins DB (input/db_config.ini, mgnproteindev); members
# missing from it get no biome. If the DB stops at 2024_04, the 2026_07-only members of 1.1 drop out.
# Also: are 1.0 representatives that moved to another MGYP at the same coordinates still in the 1.1 searched set?
# Read-only. Run from the launch dir (holds input/db_config.ini): sbatch biome_probe.sh
: "${V:?set V, see paper/README.md}"
FAM=${FAM:-$V/1.1/output/update_families}
B10=${B10:-$V/1.0/output_db/biome_results}
B11=${B11:-$V/1.1/output/biome_results}
DBC=${DBC:-input/db_config.ini}
: "${CACHE:?set CACHE, see paper/README.md}"
# pipeline QUERY_MGNPROTEIN_DB image (psycopg2 2.9.10)
IMG=${IMG:-$(ls "$CACHE"/*56121541d6b4f7082058be7bf67ed33c9e6d82dd57015ca431f632b696117179* 2>/dev/null | head -1)}
# 1.1 reps from survival.sh (new IDs first, then two old ones)
IDS=${IDS:-"10197268409 10991932733 9160527382 10264444962 740907732 395001733"}
# 1.1 searched-set MGYPs from diff_10_11.sh, and 1.0 reps that moved to another MGYP at the same coordinates (+ 25 as control)
: "${U11:?set U11, see paper/README.md}"
REPS10=${REPS10:-"2723237289 1011785639 1005985112 1661654766 2712181692 888909820 439634043 2756266471 542704473 741878026"}
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
echo "== image: ${IMG:-NOT FOUND (set IMG=)}"; ls -la "$DBC" "$FAM/refined_families.tsv"

echo "== DB: max mgyp; then per id: mgyp in DB? has biomes?"
singularity exec ${BIND:+-B $BIND} "$IMG" python3 - "$DBC" $IDS <<'PY' | tee "$TMP/db.txt"
import sys, configparser, psycopg2
c = configparser.ConfigParser(); c.read(sys.argv[1])
cur = psycopg2.connect(**dict(c.items("database"))).cursor()
cur.execute("SELECT max(mgyp) FROM sequence_explorer_protein"); print("max", cur.fetchone()[0])
ids = [int(x) for x in sys.argv[2:]]
cur.execute("SELECT mgyp, metadata->'b' IS NOT NULL FROM sequence_explorer_protein WHERE mgyp = ANY(%s::bigint[])", (ids,))
found = dict(cur.fetchall())
for i in ids: print(i, "in_db" if i in found else "MISSING", found.get(i, "-"))
PY
MAX=$(awk '$1=="max"{print $2}' "$TMP/db.txt")
[ -n "$MAX" ] || { echo "no DB max, stop"; exit 1; }

# family<TAB>number of leaves<TAB>leaf (when exactly one); same as survival.sh
leaves(){
  find "$1" -maxdepth 1 -name '*.csv' -print0 | xargs -0 awk -F, '
    function flush(){ if (fam == "") return; n = 0; one = ""
      for (i in id) if (!(i in par)) { n++; one = i }
      print fam "\t" n "\t" (n == 1 ? one : "-"); delete id; delete par }
    FNR == 1 { flush(); fam = FILENAME; sub(/.*\//, "", fam); sub(/\.csv$/, "", fam); next }
    { id[$1] = 1; if ($3 != "") par[$3] = 1 }
    END { flush() }'
}
leaves "$B10" > "$TMP/b10.tsv"; leaves "$B11" > "$TMP/b11.tsv"

echo "== the 1.0 exclusive families: id | 1.0 leaf | 1.1 leaves, leaf"
awk -F'\t' -v L="$TMP/b11.tsv" 'BEGIN { while ((getline l < L) > 0) { split(l, a, "\t"); n[a[1]] = a[2]; lf[a[1]] = a[3] } }
  $2 == 1 { print $1, "|", $3, "|", n[$1], lf[$1], ($3 == lf[$1] ? "same" : "MOVED") }' "$TMP/b10.tsv" | sort -n

echo "== 1.0 representatives: in the 1.1 searched set?"
awk -v W="$REPS10" 'BEGIN { n = split(W, a, " "); for (i = 1; i <= n; i++) w[a[i]] = 1 }
  { x = $1; if (x !~ /^[0-9]/) gsub(/[^0-9]/, "", x); x = x + 0; if (x in w) f[x] = 1 }
  END { for (k in w) print k, (k in f) ? "searched in 1.1" : "NOT searched in 1.1" }' "$U11"

# members are bare numbers (mgyp/start-end); $2 + 0 stops at "/"
echo "== (slow, 2.4B rows) 1.1 members with mgyp > DB max ($MAX), a lower bound on members without biomes: families matched, members, new, new %; by group (exclusive in 1.1 or not): families, mean new %, families all-new"
head -1 "$FAM/refined_families.tsv"
awk -F'\t' -v MAX="$MAX" -v L="$TMP/b11.tsv" '
BEGIN { while ((getline l < L) > 0) { split(l, a, "\t"); n[a[1]] = a[2] } }
{ id = $2 + 0; t[$1]++; T++; if (id + 0 > MAX) { nw[$1]++; N++ } }
END { for (f in t) { if (!(f in n)) continue; m++; g = (n[f] == 1) ? "excl" : "multi"; G[g]++; F[g] += nw[f] / t[f]; if (nw[f] == t[f]) A[g]++ }
  printf "families %d, members %d, new %d (%.1f%%)\n", m, T, N, 100 * N / T
  for (g in G) printf "%s %d families, mean new %.1f%%, all-new %d\n", g, G[g], 100 * F[g] / G[g], A[g] }' "$FAM/refined_families.tsv"
