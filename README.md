# mgnifams-post-pipeline-scripts
Contains post-processing scripts for various stats, after the main MGnifams pipeline finishes execution

## bin/extract_true_novel_superfamilies.py
Writes a TXT list of true-novel superfamily ids. A superfamily is included only
when every member family listed in the `Family Ids` column of
`superfamily_statistics.csv` is present in `novel_ids.txt`.

Default inputs:
- `assets/mgnifams_v2_results/generate_families/novel/novel_ids.txt`
- `assets/mgnifams_v2_results/generate_families/network/superfamily_statistics.csv`

Default output:
`assets/mgnifams_v2_results/generate_families/network/true_novel_superfamilies.txt`

```
python bin/extract_true_novel_superfamilies.py
```

Custom paths can be supplied with:

```
python bin/extract_true_novel_superfamilies.py <path/to/novel_ids.txt> \
  <path/to/superfamily_statistics.csv> \
  --output <path/to/true_novel_superfamilies.txt>
```

## bin/calculate_true_novel_superfamily_novelty.py
Filters `superfamily_statistics.csv` down to the true-novel superfamilies listed
in `true_novel_superfamilies.txt`, computes a percentage novelty score for each
row from the family IDs present in `input/mgnifams_l100_plddt70_novel.csv`,
drops rows with no matching family IDs, and sorts the output by novelty score
descending, then family size descending.

Default inputs:
- `input/mgnifams_l100_plddt70_novel.csv`
- `assets/mgnifams_v2_results/generate_families/network/superfamily_statistics.csv`
- `assets/mgnifams_v2_results/generate_families/novel/true_novel_superfamilies.txt`

Default output:
`assets/mgnifams_v2_results/generate_families/network/true_novel_superfamily_novelty_scores.csv`

```
python bin/calculate_true_novel_superfamily_novelty.py
```

Custom paths can be supplied with:

```
python bin/calculate_true_novel_superfamily_novelty.py \
  <path/to/mgnifams_l100_plddt70_novel.csv> \
  <path/to/superfamily_statistics.csv> \
  <path/to/true_novel_superfamilies.txt> \
  --output <path/to/true_novel_superfamily_novelty_scores.csv>
```

## bin/build_superfamilies.py
Builds weighted superfamily clusters from the MGnifams family similarity MultiQC CSV.
The script removes the display-only `Row` column, filters out edges where either
family id contains `_` (pre-redundant leftovers), writes the filtered edge list, identifies connected
components as superfamilies, selects the representative family by highest
within-cluster summed Jaccard weight, exports canonical family ids in `1..35459`
that are absent from both edge columns as singletons, writes a size-descending
statistics CSV, writes a total family count sanity check, and renders two size
distribution barplots: non-singleton exact sizes `2..10`, then 10-wide ranges
from `11` onward.

Default input:
`assets/mgnifams_v2_results/generate_families/similarity_mqc.csv`

Default outputs are written under
`assets/mgnifams_v2_results/generate_families/network/`:
- `similarity_filtered_superfamilies.csv`
- `singleton_ids.txt`
- `superfamily_statistics.csv`
- `total_families.txt`
- `superfamily_size_distribution_2_to_10.png`
- `superfamily_size_distribution_11_plus.png`

```
python bin/build_superfamilies.py
```

Custom paths can be supplied with:

```
python bin/build_superfamilies.py <path/to/similarity_mqc.csv> \
  --filtered-csv <path/to/filtered_edges.csv> \
  --singleton-ids <path/to/singleton_ids.txt> \
  --report-csv <path/to/superfamily_statistics.csv> \
  --total-families <path/to/total_families.txt> \
  --plot-small-png <path/to/superfamily_size_distribution_2_to_10.png> \
  --plot-large-png <path/to/superfamily_size_distribution_11_plus.png>
```

## bin/rank_hmmstat.py
ank profile HMMs from a HMMER `hmmstat` report by a composite quality score

## bin/update_seed_msa_blobs.py
Updates seed MSA blobs (used after trimming off env coords)

## bin/trim_seed_msa_envelopes.py
Trims off any start/end parts of seed MSA sequences that belong to the envelope and re-calculates the sequence coordinates appropriately.

## bin/plot_family_length_distribution.py
Produces 3 stacked bar chart PNGs of family distribution by annotation status, split into short,
medium, and long groups. Binning and splitting can be done by HMM consensus length (aa) or family
size (number of sequences) via `--mode length|size`.
```
python bin/plot_family_length_distribution.py \
  --metadata <path/to/metadata_mqc.csv> \
  --domtbl <path/to/funfams.domtbl.gz> <path/to/pfam.domtbl.gz> \
  --output-prefix <path/to/family_length> \
  --mode length \
  --short-max 300 --med-max 1000 \
  --short-bin 10 --med-bin 50 --long-bin 100
```

## calculate_annotation_stats.py
Calculates total annotated sequences and total annotated amino acids from a
protein CSV or CSV.GZ metadata column. Use `--annotation-keys p` for Pfam-only
statistics and `--annotation-keys p,m` after appending MGnifam annotations.
Lives under `pipelines/annotation_percentages/bin/` (used by that pipeline).

```
python pipelines/annotation_percentages/bin/calculate_annotation_stats.py \
  --input <proteins.csv.gz> \
  --annotation-keys p,m \
  --label pfam_mgnifam \
  --output <annotation_stats.csv>
```

## calculate_stars.py
Through sqlite queries, calculates each family stars for both quality and novelty.

## calculate_stars_from_folders.py - not updated with extra logic
Through parsing the MGnifams output folders, calculates each family stars for both quality and novelty.

## extract_tm_families.py
From the mgnifams workdir, find and extract TM families into a dedicated output folder.

## extract_fasta_from_proteins_csv.py
From the initial MGnifams input (PLP output CSV or CSV.GZ), extract a FASTA file.
Lives under `pipelines/annotation_percentages/bin/` (used by that pipeline).

## extract_reps_fasta_from_msa.py
From a .fas format MSA folder, extract all sequence representatives in a single fasta file.

## annotate_fasta_with_hmm.sh
Annotates FASTA sequences with `hmmsearch` inside the configured HMMER
Singularity image. Gzipped HMM libraries are decompressed into a temporary
work directory before searching.

```
bash bin/annotate_fasta_with_hmm.sh \
  assets/mgnifams_v2_results/generate_families/families/mgnifams_hmm.lib.gz \
  <proteins.fasta> \
  <mgnify_proteins_mgnifams.domtbl>
```

## append_mgnifams_annot.py
Appends MGnifam `hmmsearch --domtblout` annotations back onto the initial CSV or
CSV.GZ sequences by replacing `metadata["m"]`. Reads plain or gzipped domtbl.
Lives under `pipelines/annotation_percentages/bin/` (used by that pipeline).

## compare_annotation_stats.py
Compares two outputs from `calculate_annotation_stats.py` and writes the
percentage-point and relative annotation increase.
Lives under `pipelines/annotation_percentages/bin/` (used by that pipeline).

## pipelines/annotation_percentages
Mini local Nextflow workflow (loose nf-core layout: one module per process under
`modules/`, Conda env + container per module). Takes a local protein CSV/CSV.GZ,
extracts FASTA, runs MGnifam `hmmsearch` (nf-core `hmmer/hmmsearch` module),
appends MGnifam annotations, computes Pfam-only and Pfam+MGnifam annotation
statistics, and compares the increase. Its python scripts live in that pipeline's
own `bin/`. See `pipelines/annotation_percentages/README.md`.

```
nextflow run pipelines/annotation_percentages/main.nf -profile singularity
```

## extract_pdb_scores.py
From the mgnifams workdir, find and extract esmfold predicted plddt and ptm scores along with name and length in csv format.

## calculate_plddt_from_cif.py
Fallback measure for those 10 sequences that pdb predictions were lost, and we average the predicted scores from their cif file. ptm 0.0
- didn't use in the end

## extract_pdb_scores.py
From the mgnifams workdir, find, concat and extract foldseek output m8s.

## extract_parsed_domains.py
From the mgnifams workdir, find and output all parsed domain data files in an output folder.

## subset_mgnifams_sqlite.py
Given a file with family ids, subset a new sqlite db.

## calculate_fam_similarities.py
Given two of the mgnifams output tables csvs as input (mgnifam_proteins.csv and mgnifam.csv), calculates similarity output csv based on family and aa jackard indices.

## report_redundant_fams.py
Given a similarities.csv, produces a list of redundant fam ids as well as the similarities edgelist from 0.5 <= x < 0.95.
If similarity >= 0.95, keep the bigger family. If same, keep family with smaller id.

## assign_correct_domain_ids.py
After removing redundant, filenames do not correspond to family ids inside the domain architecture.
This script map it properly to the basename of the file.

## update_domain_blobs.py
From a domain_results folder, update blobs of sqlite.

## update_stars_sqlite.py
Checks if column exists. Create if not. Then append data.

## identify_unmgnifamed_space.py
This scripts receives two inputs; the original linclust clusters produced by the MGnifams pipeline
and the final mgnifam_proteins.csv table. All initial clusters that have at least one of their MGYP
members in the final proteins are removed. The remaining cluster reps along with their cluster sizes are reported.

## calc_annot_diffs.py
This script receives two annotated files with three columns (id,protein,region).
The id is relative to the HMM lib we annotated against.
The annotated proteins between these two files must much.
The script calculates the annotation differences and reports them.

## concatenate_panther_subfolder_hmms.sh
Bash script to aggregate PANTHER extracted hmms into one hmm file.

## extract_hmmsearch_mgnifams_exclusive.py
Write out exclusive annotated protein regions of second given file, compared to the first given file.
CSV File format: mgnifam_id,protein,region
