# mgnifams-post-pipeline-scripts

Scripts behind the numbers, figures and tables of the MGnifams paper (releases 1.0 and 1.1, https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnifams/). They run on the outputs of the [MGnifams pipeline](https://github.com/EBI-Metagenomics/mgnifams) and on public MGnify, UniProt and Pfam data.

Scripts that no paper number depends on (superseded analyses, exploration, website database tools) are in [`playground/`](playground/README.md).

## Layout

| Path | What |
|---|---|
| `pipelines/annotation_percentages/` | Coverage of MGnify90 by Pfam and MGnifams (Nextflow) |
| `pipelines/uniprot_annotation_percentages/` | Pfam and MGnifams searches of UniProtKB (Nextflow) |
| `pipelines/mgnifam_uniprot_coverage/` | UniProtKB coverage statistics, Tables 1 and 2, from those searches (Nextflow) |
| `bin/` | `biome_analysis.py` (biome leaves per family), `pfam_version_checks.py` (Pfam release of a set of annotations) |
| `paper/scripts/` | Analyses and figure / table generators; Python files import each other, so they stay side by side |
| `paper/data/` | Small inputs and intermediate tables the generators read |
| `paper/figures/` | Figure 1 (nf-metro spec) and Figure 2 (Excalidraw scene) |
| `paper/tests/`, `tests/` | Tests |

## Paper item → script

| Paper item | Scripts | Inputs |
|---|---|---|
| Release 1.1 counts (updated and not-updated families, by reason) | `paper/scripts/pass1.sh` | release 1.1 outputs |
| MGnify90 coverage (Results, Abstract) | `pipelines/annotation_percentages` | MGnify Proteins 2024_04 MGnify90, `mgnifams_hmm.lib.gz` |
| UniProtKB coverage, Tables 1 and 2 | `pipelines/uniprot_annotation_percentages` (searches), `pipelines/mgnifam_uniprot_coverage` (statistics, one search space), `paper/scripts/tm_classes.py` and `catalogue_lists.py` (Table 2 categories), `tables_1_2.py` (tables) | UniProtKB 2026_02, Pfam-A 38.2 |
| Annotation evidence (Results, Supplementary Table S1, Discussion) | `annotation_evidence.py` / `.sh`, `table_s1.py` | release 1.0 and 1.1 `table_data/`, Pfam clan lists 36.0, 38.1, 38.2 |
| Figure 3, Supplementary Figures S1–S3 | `catalogue_probe.py` / `.sh` (writes `data/catalogue_*.csv.gz`), `figure3.py` | release 1.0 and 1.1 outputs |
| Transmembrane and disordered families | `tm_classes.py` / `.sh` | DeepTMHMM, s4pred and ESMFold outputs of the release |
| Pfam family sizes and model lengths | `pfam_recount.sh`, `pfam_lengths.py` | MGnify Proteins 2026_07 Pfam parquet, Pfam-A 38.1 HMMs |
| Pfam release behind MGnify Proteins 2024_04 (Methods) | `bin/pfam_version_checks.py` | dated annotation samples, Pfam family lists |
| Biome breadth and biome-exclusive families, Supplementary Table S2 | `biome_probe.sh`, `biome_breadth.py` / `.sh` (imports `bin/biome_analysis.py`), `breadth_size.py`, `survival.sh`, `uniprot_hits.sh`, `uniprot_hit_taxa.sh`, `table_s2.py` | release biome tables, UniProtKB search outputs |
| Habitat check of the biome-exclusive families (BLASTx) | `contig_lookup.py` / `.sh`, `fetch_contigs.py`, `fetch_hits.py`, `blastx_uniref.sh`, `blastx_summary.py`, `uniref_summary.py` | MGnify API, UniProt API, UniRef90 2026_02 |
| Clans | `clan_inputs.sh`, `compute_clans.sh` (`mgnifam compute_clans` 4.1.0), `clan_stats.py` | release 1.1 full MSAs |
| Near-identical family pairs (Limitations) | `near_identical.py` | `clan_edges.tsv.gz` of two runs |
| Change between releases 1.0 and 1.1 (Discussion) | `diff_10_11.sh`, `diff_followup.sh` | release 1.0 and 1.1 full MSAs |
| Multiple hits per protein (Methods, case studies) | `classify_multi_hits.py` / `.sh`, `rows_stats.py` | release 1.1 full MSAs |
| Domain architectures, Figure 4; Pfam domains per family | `architectures.py` / `.sh`, `figure_architectures.py`, `pfam_per_family.py` / `.sh` | full MSAs, MGnify Proteins 2026_07 Pfam parquet |
| Member regions of MGYF0000000025 against Pfam and FunFam | `member_search.sh` | full MSA, Pfam-A 38.2, CATH-FunFams 4.3.0 |
| Gene neighbourhoods (MGYF0000000025, MGYF0000000243, MGYF0000032058), Figure 6 | `neighbourhood.py` / `.sh`, `independence.py`, `figure6.py` | full MSAs, MGnify Proteins 2026_07 parquets, MGnify API |
| Structure of the 1.1_clan_3863 representative, Figure 5 | `figure5.py`, `render_structure.py` (PyMOL) | `data/esmfold_*`, Pfam 38.2 domtbl |
| Novel families against NMPFamsDB (Supplementary Table S3) | `novel_reps.py` | `data/novel_*`, NMPFamsDB v1 and NMPFamsDB2 HMMs |
| Representatives no longer in MGnify Proteins (Limitations) | `rep_resolve.py` | release 1.1 `families.tsv.gz`, MGnify Proteins API |
| Figure 1 | `paper/figures/figure-1.mmd` | rendered with nf-metro |
| Figure 2 | `figure2.gen.py` | writes `paper/figures/figure-2.excalidraw` |

## Running

- Python 3.10 or later. Most scripts use the standard library only; `figure3.py` needs plotnine, `figure_architectures.py`, `figure5.py` and `figure6.py` need matplotlib, `render_structure.py` needs PyMOL, and the parquet readers need pyarrow or polars. External tools: HMMER 3.4, DIAMOND 2.1.8, mgnifam 4.1.0 (all as Singularity images).
- Run generators from the repository root (`python paper/scripts/figure6.py --out figures/figure-6`); their defaults read `paper/data/`.
- The `.sh` files are Slurm jobs for a cluster that holds the release outputs. Each stops with a message naming the variables to set (`V`: the directory holding the `1.0/` and `1.1/` release outputs; `DB`: an MGnify Proteins release; `CACHE` / `IMG`: Singularity images; `BIND`: extra bind paths; `W`: a scratch directory; others are named in the script header). Copy the Python files a wrapper names next to it before submitting.
- Tests: `python -m pytest tests` and `for t in paper/tests/test_*.py; do python $t; done`; many scripts also take `--selftest`.

