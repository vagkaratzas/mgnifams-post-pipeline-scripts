#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_diff_followup
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --output=mgnifams_diff_followup_%j.log
# Follow-up to diff_10_11.sh: (1) 1.0 converged list vs family_metadata flag (3,733 vs 3,686);
# (2) the 8,876 families whose 1.0 member-TSV rows != full_size. Read-only, minutes. W = the diff_10_11.sh scratch dir (reads its t10.tsv): W=<that dir> sbatch diff_followup.sh
: "${W:?set W to the diff scratch dir}"
: "${V:?set V, see paper/README.md}"
F10=${F10:-$V/1.0/output/generate_families/families}
P10=${P10:-$V/1.0/output_db/families.tsv.gz}
ls "$F10"
echo "== converged_families.txt: lines, unique ids, max id"
awk '{n++; if(!($1 in s)){s[$1]=1; u++}; if($1+0>m) m=$1+0} END{print n, u, m}' "$F10/converged_families.txt"
echo "== family_metadata.csv: rows; first 2 rows (cut); last-field counts"
wc -l < "$F10/family_metadata.csv"; head -2 "$F10/family_metadata.csv" | cut -c1-160
awk -F, '{print $NF}' "$F10/family_metadata.csv" | sort | uniq -c
echo "== flag True? x in list? (count), and up to 5 disagreeing rows: id iteration flag"
awk -F, 'FNR==NR{c[$1]=1; next} {m=($NF=="True"); f=($1 in c); k[m" "f]++; if(m!=f && e<5){e++; print "  ", $1, $2, $NF}} END{for(x in k) print "flag", x, k[x]}' "$F10/converged_families.txt" "$F10/family_metadata.csv"
echo "== 1.0 rows in member TSV vs full_size: mismatched, sum(full_size - rows), min, max, examples id:full_size/rows"
zcat "$P10" | awk -F'\t' -v t="$W/t10.tsv" 'BEGIN{while((getline l < t)>0){split(l,a,"\t"); T[a[1]]=a[2]}} NR>1{d=$9-T[$1]; if(d){n++; s+=d; if(n==1||d<mn) mn=d; if(n==1||d>mx) mx=d; if(n<=3) ex=ex" "$1":"$9"/"T[$1]}} END{print n, s, mn, mx, ex}' | tee "$W/mismatch.txt"
echo "== full MSA of the examples: alignment rows, rows with a repeated name"
for x in $(awk '{for(i=5;i<=NF;i++){split($i,a,":"); print a[1]}}' "$W/mismatch.txt"); do
  f=$(ls "$F10"/full_msa_sto/"$x".sto* 2>/dev/null | head -1)
  [ -n "$f" ] && zcat -f "$f" | awk -v id="$x" '!/^#/ && NF==2{if(!($1 in s)){s[$1]=1; n++} else d++} END{print id, "unique rows", n, "repeated", d+0}'
done
