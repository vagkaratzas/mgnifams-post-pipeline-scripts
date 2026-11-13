# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Post-processing scripts for the MGnifams (MGnify protein families) pipeline — part of the MGnify platform. Scripts run after the main MGnifams pipeline completes to compute quality metrics, manage redundancy, annotate sequences, and populate the SQLite database for the web portal.

## Working Conventions

- **Git**: Never commit or push except if explicitly asked to.
- **Tooling**: Always use the `token-saviour` skill.
- **Response style**: Keep responses concise; avoid exceeding output token limits by summarizing rather than dumping full files.

## Running Scripts

All Python scripts use `argparse` and accept `--help`. Most follow this pattern:

```bash
python <script>.py --input <file> --output <file> [options]
```

Bash scripts are run directly:
```bash
bash annotate_fasta_with_hmm.sh <hmm_file[.gz]> <fasta> <domtblout>
```

## Layout

- **Publication scripts** (every number, figure and table of the MGnifams paper): `pipelines/annotation_percentages`, `pipelines/uniprot_annotation_percentages`, `pipelines/mgnifam_uniprot_coverage`, `bin/biome_analysis.py`, `bin/pfam_version_checks.py`, and `paper/` (`scripts/` side by side because they import each other, `data/`, `figures/`, `tests/`). The root `README.md` maps each paper item to its script; keep it to publication scripts only.
- **Playground** (`playground/`): superseded analyses, exploration and the website-database tools (`update_*_blobs.py`, `update_stars_sqlite.py`, `subset_mgnifams_sqlite.py`, `calculate_stars*.py`), with their tests and fixtures; described in `playground/README.md`.
- Tests: `python -m pytest tests`, `python -m pytest playground/tests`, and `python paper/tests/test_<name>.py`.
- No local or cluster paths in committed files: cluster locations are environment variables the `.sh` wrappers require.

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
- **S4Pred output** (3-line format: Conf/Pred/AA) is parsed by `parse_s4pred_to_feature_viewer.py` into FeatureViewer-compatible JSON
- Similarity files use CSV with `fam1,fam2,family_jaccard,aa_jaccard` columns; redundancy threshold is `>= 0.95`
