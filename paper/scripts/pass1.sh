#!/usr/bin/env bash
# Release 1.1 sanity pass 1: counts and completeness. Read-only.
# Run on the cluster from the launch dir (the folder holding mgnifams/, conf/, input/).
set -u
SHEET=${SHEET:-input/samplesheet_update_mgnifams_prod.csv}
: "${OUT:?set OUT, see paper/README.md}"
FAM=$OUT/update_families

echo "== samplesheet"; cat "$SHEET"
echo "== 1.0 inputs md5 (FTP 1.0: hmm f06909ac…, families 7e6ef1d5…)"
tail -n +2 "$SHEET" | cut -d, -f5,7 | tr , '\n' | xargs md5sum
echo "== pipeline commit"
git -C mgnifams rev-parse HEAD; git -C mgnifams log --oneline -3; git -C mgnifams status --short | head
echo "== sqlite files (delta DB)"; find "$OUT/.." -maxdepth 4 -name '*.sqlite3' -exec ls -lh {} \;

echo "== families.tsv.gz: rows, duplicate ids, status, reasons, model release"
zcat "$FAM/families.tsv.gz" | awk -F'\t' 'NR>1{n++; s[$2]++; if($2=="not_updated") r[$3]++; m[$5]++} END{print "families",n; for(k in s) print "status",k,s[k]; for(k in r) print "reason",k,r[k]; for(k in m) print "model_release",k,m[k]}'
echo "duplicate ids: $(zcat "$FAM/families.tsv.gz" | cut -f1 | sort | uniq -d | wc -l)"
echo "== updated_delta.csv outcomes (must equal the status/reason counts)"
awk -F, 'NR>1{n++; o[$NF]++} END{print "rows",n; for(k in o) print k,o[k]}' "$FAM/updated_delta.csv"

echo "== line counts (header included; per-family files must equal 'updated' + 1)"
wc -l "$FAM"/*.csv "$FAM"/*.txt "$FAM"/*.tsv "$OUT"/table_data/*.csv 2>/dev/null
echo "full_msa $(ls "$FAM/full_msa" | wc -l)  biome $(ls "$OUT/biome_results" | wc -l)  domain $(ls "$OUT/domain_results" | wc -l)"
echo "== update_info"; cat "$FAM/update_info.csv"
echo "== searched set (cluster reps, not all 5.7B proteins)"; cat "$OUT"/qc/*.tsv
echo "== params"; grep -h -E 'parquet|input|release|skip|recruit|outdir' "$OUT"/pipeline_info/params_*.json | sort -u
echo "== versions"; grep -h -i -E 'mgnifam|pyhmmer|nextflow' "$OUT"/pipeline_info/*versions*.yml | sort -u

echo "== members (sum of full_size): 1.0, then 1.1 by status"
zcat "$(tail -n +2 "$SHEET" | cut -d, -f7)" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next}{s+=$c["full_size"]} END{print "1.0",s}'
zcat "$FAM/families.tsv.gz" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next}{s[$c["status"]]+=$c["full_size"]} END{for(k in s) print "1.1",k,s[k]}'
