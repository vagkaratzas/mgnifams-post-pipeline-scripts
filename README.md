# mgnifams-post-pipeline-scripts
Contains post-processing scripts for various stats, after the main MGnifams pipeline finishes execution

## bin/biome_analysis.py
Streams MGnifam biome blobs from the `mgnifam` SQLite table, counts leaf-level
biome paths per family, writes a text report, and optionally renders a PNG bar
plot. Leaf biomes are keyed by their full `root:...` path, so duplicate terminal
labels such as `Sediment` or `Fecal` remain separate when they come from
different branches. The report also lists duplicate terminal leaf labels to show
where last-label-only grouping would collapse distinct paths.

Default outputs are written under `output/`:
- `mgnifams_biome_distribution.txt`
- `mgnifams_biome_distribution.png`

```
python bin/biome_analysis.py <path/to/mgnifams.sqlite3>
```

Useful options:

```
python bin/biome_analysis.py <path/to/mgnifams.sqlite3> \
  --output-dir <path/to/output> \
  --no-plot \
  --log-every 1000
```

## bin/extract_true_novel_superfamilies.py
Writes a TXT list of true-novel superfamily ids. A superfamily is included only
when every member family listed in the `Family Ids` column of
`clan_membership.csv` is present in `mgnifams_no_annotation_ids.txt`.

Default inputs:
- `assets/mgnifams_v2_results/generate_families/novel/mgnifams_no_annotation_ids.txt`
- `assets/mgnifams_v2_results/generate_families/network/clan_membership.csv`

Default output:
`assets/mgnifams_v2_results/generate_families/network/true_novel_clans.txt`

```
python bin/extract_true_novel_superfamilies.py
```

Custom paths can be supplied with:

```
python bin/extract_true_novel_superfamilies.py <path/to/mgnifams_no_annotation_ids.txt> \
  <path/to/clan_membership.csv> \
  --output <path/to/true_novel_clans.txt>
```

## bin/annotate_novel_through_domain_architecture.py
Parses `domain_blob` JSON from the `mgnifam` SQLite table for the novel family
ids listed in `mgnifams_no_annotation_ids.txt`, extracts unique Pfam ids from
each domain architecture, and reports true-novel clans/singletons whose member
families gained Pfam annotations through those architectures.

Default inputs:
- `assets/mgnifams_v2_results/mgnifams_test.sqlite3`
- `assets/mgnifams_v2_results/generate_families/novel/mgnifams_no_annotation_ids.txt`
- `assets/mgnifams_v2_results/generate_families/network/clan_membership.csv`

Default outputs:
- `assets/mgnifams_v2_results/generate_families/novel/mgnifams_no_annotation_annotate_novel_through_domain_architecture.csv`
- `assets/mgnifams_v2_results/generate_families/network/transiently_annotated_clans.csv`
- `assets/mgnifams_v2_results/generate_families/novel/true_novel_without_annotate_novel_through_domain_architecture.txt`

```
python bin/annotate_novel_through_domain_architecture.py
```

Custom paths can be supplied with:

```
python bin/annotate_novel_through_domain_architecture.py \
  --sqlite <path/to/mgnifams.sqlite3> \
  --novel-ids <path/to/mgnifams_no_annotation_ids.txt> \
  --clan-membership <path/to/clan_membership.csv> \
  --family-pfams-output <path/to/family_pfams.csv> \
  --clan-output <path/to/transiently_annotated_clans.csv> \
  --true-novel-output <path/to/true_novel_without_pfams.txt>
```

## bin/extract_family_sequences_fasta.py
Writes FASTA records from `metadata_mqc.csv` for families listed in a TXT file.
The FASTA header is the `Family Id` value, and the sequence is taken from the
`Sequence` column. Input selector IDs such as `Singleton_243` are matched by the
second underscore-delimited field, so `Singleton_243` selects `Family Id` `243`.

```
python bin/extract_family_sequences_fasta.py \
  assets/mgnifams_v2_results/generate_families/metadata_mqc.csv \
  assets/mgnifams_v2_results/generate_families/novel/true_novel_without_annotate_novel_through_domain_architecture.txt \
  --output assets/mgnifams_v2_results/generate_families/novel/true_novel_without_annotate_novel_through_domain_architecture.fasta
```

## bin/calculate_true_novel_superfamily_novelty.py
Filters `clan_membership.csv` down to the true-novel superfamilies listed
in `true_novel_clans.txt`, computes a percentage novelty score for each
row from the family IDs present in `input/mgnifams_l100_plddt70_novel.csv`,
drops rows with no matching family IDs, and sorts the output by novelty score
descending, then family size descending.

Default inputs:
- `input/mgnifams_l100_plddt70_novel.csv`
- `assets/mgnifams_v2_results/generate_families/network/clan_membership.csv`
- `assets/mgnifams_v2_results/generate_families/novel/true_novel_clans.txt`

Default output:
`assets/mgnifams_v2_results/generate_families/network/true_novel_clan_l100_plddt70_scores.csv`

```
python bin/calculate_true_novel_superfamily_novelty.py
```

Custom paths can be supplied with:

```
python bin/calculate_true_novel_superfamily_novelty.py \
  <path/to/mgnifams_l100_plddt70_novel.csv> \
  <path/to/clan_membership.csv> \
  <path/to/true_novel_clans.txt> \
  --output <path/to/true_novel_clan_l100_plddt70_scores.csv>
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
- `similarity_mapped_mgnifams.csv`
- `unclustered_mgnifam_ids.txt`
- `clan_membership.csv`
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
  --singleton-ids <path/to/unclustered_mgnifam_ids.txt> \
  --report-csv <path/to/clan_membership.csv> \
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

## bin/plot_mgnifam_metadata_distributions.py
Family-level distribution figures for the catalogue (paper Figure X). Every panel is a stacked
barplot of annotated vs unannotated (novel) families, where novel ids are read from
`mgnifams_no_annotation_ids.txt`; the percentage above each bar is the share of unannotated
families in that bin — printed on **every** bin, rotated upright where the bar is too short to
carry it horizontally.

Written to the repo's `scientific-plotting-skill` standard: **plotnine** only, vector **PDF**
output (with a 300-dpi PNG preview alongside), one **Times** text size for every text element
(`TEXT_SIZE_PT`), no plot titles — context belongs in the manuscript caption — and the
colorblind-safe **Wong** palette (orange = annotated, sky blue = novel). Everything tunable about
how the figures look sits in a single parameter block at the top of the script.

Output is split by where the plot is meant to end up (`--mode`, default `all`):

**`figure` → `<output-dir>/figures/`** — the main-text Figure X, one panel per claim in the
manuscript paragraph:
- **A** `figure_Xa_size.pdf` — family size, log₁₀ bins, full range
- **B** `figure_Xb_length.pdf` — representative sequence length, 100-aa bins, full range (75–2,000 aa)
- **C** `figure_Xc_plddt.pdf` — mean pLDDT of the representative structure, 5-unit bins
- `figure_X_family_metadata.pdf` — the three panels stacked and tagged A/B/C, 180 mm wide

Size uses log₁₀ bins here because its range (29–1,515,677) cannot be binned legibly on a linear
axis; a single panel is what lets the reader see the full span *and* the annotated/novel shift.
Length needs no such trick — 75–2,000 aa fits one panel at 100-aa bins. Length bins are
right-closed (`1,901–2,000`, not `2,000–2,099`) so the axis closes exactly on the longest
representative sequence.

**`supplementary` → `<output-dir>/supplementary_figures/`** — the fine-grained view behind panels
A and B, plus the companion structural metric:
- `size_{small,medium,large}.pdf` — linear bins, half-open cuts so each cut lands on a bin edge
- `length_{short,medium,long}.pdf` — same
- `ptm.pdf` — pTM of the representative structure

Also writes `<output-dir>/mgnifam_metadata_stats.txt` with min/Q1/median/Q3/max/mean per metric,
reported overall and split by annotated vs novel — the numbers quoted in the manuscript.

Input is the `mgnifam` table CSV (`id,full_size,rep_length,...,plddt,ptm`): the full-catalogue
`mgnifam_codon.csv`, or `table_data/mgnifam.csv` for a quick test. Note this is *not* the
MultiQC `metadata_mqc*.csv`, which carries no pLDDT/pTM columns.

Needs `plotnine` (`pip install plotnine`) and a Times-metric serif font. Times New Roman is used
when present; otherwise the script falls back to Nimbus Roman, then Liberation Serif (both ship
with most Linux distributions), and reports which one it used.

```
python bin/plot_mgnifam_metadata_distributions.py \
  --metadata mgnifam_codon.csv \
  --novel-ids assets/mgnifams_v2_results/generate_families/novel/mgnifams_no_annotation_ids.txt \
  --output-dir output
```

Bin sizes and split thresholds are tunable; defaults are tuned to the full 35,459-family
catalogue, so a small test dataset needs smaller ones:

```
python bin/plot_mgnifam_metadata_distributions.py \
  --metadata assets/mgnifams_v2_results/table_data/mgnifam.csv \
  --novel-ids assets/mgnifams_v2_results/generate_families/novel/mgnifams_no_annotation_ids.txt \
  --output-dir output/test \
  --mode all \
  --size-small-max 100 --size-medium-max 1000 \
  --size-small-bin 20 --size-medium-bin 200 --size-large-bin 1000 \
  --length-short-max 120 --length-medium-max 200 \
  --length-short-bin 10 --length-medium-bin 20 --length-long-bin 50 \
  --plddt-bin 5 --ptm-bin 0.05
```

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
