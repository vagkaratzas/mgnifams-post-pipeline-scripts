#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_clan_inputs
#SBATCH --time=6:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=mgnifams_clan_inputs_%j.log
# Inputs for the 1.1 `mgnifam compute_clans` run: all 35,459 families, i.e. the 1.1 members of the 34,419
# updated families plus the 1.0 members of the 1,040 not updated. Checks that both member files spell rows the same
# way (`<protein>/<start>-<end>`, family id = families.tsv family_id) and counts bare (whole-protein) rows, which need
# a FASTA for their length: bare.fa holds those proteins, from the 1.1 search FASTA first, then the 1.0 one (an MGYP
# sequence never changes between releases). Writes not_updated.txt, members_1.0_not_updated.tsv and bare.fa under $W
# (the run's inputs); read-only on the releases.
#   W=$L/clans_1.1 sbatch clan_inputs.sh
set -euo pipefail
: "${W:?set W to the clans run dir}"
: "${P:?set P, see paper/README.md}"
NEW=$P/1.1/output/update_families
OLD=$P/1.0/output/generate_families/families/refined_families.tsv
FA_NEW=${FA_NEW:-$NEW/mgnifams_update.faa}
FA_OLD=${FA_OLD:-$P/1.0/output/mgnifams_v2.fa}
mkdir -p "$W"; cd "$W"
export LC_ALL=C

echo "== files"; ls -la "$NEW/families.tsv.gz" "$NEW/refined_families.tsv" "$OLD" "$FA_NEW" "$FA_OLD"

echo "== families.tsv.gz status"; : > not_updated.txt
zcat "$NEW/families.tsv.gz" | awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
  { n[$c["status"]]++; if ($c["status"] == "not_updated") print $c["family_id"] > "not_updated.txt" }
  END { for (s in n) print s, n[s] }'
wc -l < not_updated.txt

# shape of a member: R = protein/start-end, B = bare, X = anything else; protein spelling by its first character class
shape='function shape(m) { return m ~ /\/[0-9]+-[0-9]+$/ ? "R" : (m ~ /\// ? "X" : "B") }
  function pclass(m) { sub(/\/.*/, "", m); return m ~ /^[0-9]+$/ ? "numeric" : (m ~ /^MGYP[0-9]+$/ ? "MGYP" : "other") }'

echo "== 1.1 members (updated families): head, rows by shape, protein spelling, families"
head -3 "$NEW/refined_families.tsv"
awk -F'\t' "$shape"'{ s[shape($2)]++; p[pclass($2)]++; f[$1] = 1; if (shape($2) != "R" && x++ < 5) print "  e.g.", $0 }
  END { for (k in s) print "shape", k, s[k]; for (k in p) print "protein", k, p[k]; print "families", length(f) }' \
  "$NEW/refined_families.tsv"

echo "== 1.0 members of the not-updated families -> members_1.0_not_updated.tsv"
head -3 "$OLD"
awk -F'\t' "$shape"'NR == FNR { want[$1] = 1; next } ($1 in want) { print; s[shape($2)]++; p[pclass($2)]++; f[$1] = 1
  if (shape($2) != "R" && x++ < 5) print "  e.g.", $0 > "/dev/stderr" }
  END { for (k in s) print "shape", k, s[k] > "/dev/stderr"; for (k in p) print "protein", k, p[k] > "/dev/stderr"
        print "families matched", length(f), "of", length(want) > "/dev/stderr" }' \
  not_updated.txt "$OLD" > members_1.0_not_updated.tsv.tmp && mv members_1.0_not_updated.tsv.tmp members_1.0_not_updated.tsv
wc -l < members_1.0_not_updated.tsv

echo "== overlap: not-updated ids also in the 1.1 members file (expect 0)"
awk -F'\t' 'NR == FNR { want[$1] = 1; next } ($1 in want) { n++; if (n <= 3) print } END { print n + 0 }' \
  not_updated.txt "$NEW/refined_families.tsv"

echo "== bare proteins -> bare.fa (records named by the bare id, first found wins)"
awk -F'\t' '$2 !~ /\// { print $2 }' "$NEW/refined_families.tsv" members_1.0_not_updated.tsv | sort -u > bare_ids.txt
wc -l < bare_ids.txt
awk 'NR == FNR { want[$1] = 1; n++; next }
  /^>/ { id = substr($1, 2); keep = (id in want) && !(id in got); if (keep) { got[id] = 1; from[FILENAME]++ } }
  keep { print }
  END { for (f in from) print from[f], "from", f > "/dev/stderr"; print length(got), "of", n, "found" > "/dev/stderr"
        for (id in want) if (!(id in got) && m++ < 5) print "  missing", id > "/dev/stderr" }' \
  bare_ids.txt "$FA_NEW" "$FA_OLD" > bare.fa.tmp && mv bare.fa.tmp bare.fa
df -h "$W" | tail -1
