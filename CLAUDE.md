# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Post-processing scripts for the MGnifams (MGnify protein families) pipeline — part of the MGnify platform. Scripts run after the main MGnifams pipeline completes to compute quality metrics, manage redundancy, annotate sequences, and populate the SQLite database for the web portal.

## Running Scripts

All Python scripts use `argparse` and accept `--help`. Most follow this pattern:

```bash
python <script>.py --input <file> --output <file> [options]
```

Bash scripts are run directly:
```bash
bash annotate_fasta_with_hmm.sh <hmm_file[.gz]> <fasta> <domtblout>
```

Test data is available in `assets/test_data/` for manual testing.

## Pipeline Stages & Script Relationships

Scripts are meant to be run in this logical order:

1. **Extraction** — Pull data from the MGnifams work directory output:
   - `extract_pdb_scores.py`, `extract_foldseek_m8s.py`, `extract_parsed_domains.py`
   - `extract_reps_fasta_from_msa.py`, `extract_tm_families.py`

2. **Metrics** — Compute similarity and quality scores:
   - `calculate_fam_similarities.py` → `report_redundant_fams.py` (Jaccard-based redundancy detection)
   - `calculate_stars.py` → `update_stars_sqlite.py` (quality/novelty ranks into DB)
   - `calculate_annotation_stats.py` and `compare_annotation_stats.py`

3. **Redundancy removal** — After filtering redundant families:
   - `assign_correct_domain_ids.py` (remaps MGnifam IDs in domain JSON files)
   - `remove_ids_from_csv.py`, `map_ids.py`, and related `map_*` / `remove_*` scripts

4. **Annotation** — Self-annotation with MGnifam HMMs:
   - `pipelines/annotation_percentages/main.nf` runs the local CSV/CSV.GZ → FASTA → hmmsearch (nf-core `hmmer/hmmsearch` module) → append → stats comparison workflow
   - That pipeline's python scripts (`extract_fasta_from_proteins_csv.py`, `append_mgnifams_annot.py`, `calculate_annotation_stats.py`, `compare_annotation_stats.py`) now live in `pipelines/annotation_percentages/bin/`, not the repo-root `bin/`
   - `calc_annot_diffs.py`, `extract_hmmsearch_mgnifams_exclusive.py`

5. **Database population** — Finalize the SQLite database:
   - `update_domain_blobs.py`, `update_secondary_structure_blobs.py`
   - `subset_mgnifams_sqlite.py`, `parse_s4pred_to_feature_viewer.py`

## SQLite Schema (key tables)

| Table | Key Columns |
|-------|-------------|
| `mgnifam` | `id`, `plddt`, `family_size`, `rep_length`, `protein_rep`, `rep_region`, `domain_architecture_file`, `domain_architecture_blob`, `pred_secondary_structure_file`, `pred_secondary_structure_blob`, `quality_rank`, `novelty_rank` |
| `mgnifam_proteins` | `mgnifam_id`, `protein`, `region` |
| `mgnifam_pfams` | `mgnifam_id`, `template_hmm_range` |
| `mgnifam_folds` | `mgnifam_id`, `target_structure`, `aligned_length`, `rep_length` |

## External Tool Dependencies

- **HMMER** (`hmmsearch` via Singularity for the local annotation workflow)
- **Python**: `pandas`, `biopython`, standard library only for most scripts

## Key Data Formats

- **region** fields use `start-end` string format (e.g., `"1-150"`)
- **metadata** in CSV files is stored as JSON strings in a `metadata` column; MGnifam annotations are stored under the `"m"` key as `[family_id, domain_i_evalue, domain_score, ali_from, ali_to]`
- **Domain architecture files** are JSON; after redundancy removal, internal MGnifam ID references must be updated via `assign_correct_domain_ids.py`
- **S4Pred output** (3-line format: Conf/Pred/AA) is parsed by `parse_s4pred_to_feature_viewer.py` into FeatureViewer-compatible JSON
- Similarity files use CSV with `fam1,fam2,family_jaccard,aa_jaccard` columns; redundancy threshold is `>= 0.95`
