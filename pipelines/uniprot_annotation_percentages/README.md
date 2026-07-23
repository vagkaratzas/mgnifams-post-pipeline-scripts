# UniProt Annotation Percentages Pipeline

Sister pipeline to `annotation_percentages`. Annotates a **UniProt FASTA** with **both Pfam and
Mgnifams HMMs** (two real `hmmsearch` runs) and compares sequence-level and residue-level
annotation coverage, to quantify whether Mgnifams enlarges the annotation space beyond Pfam.

Unlike `annotation_percentages` (which takes an MGnify CSV with Pfam already baked in and runs only
the Mgnifams search), this pipeline starts from a plain FASTA with no pre-computed annotations and
derives all statistics directly from the two searches' `--domtblout` outputs.

## Scoring

- **Pfam → `--cut_ga`**: Pfam models carry curated per-model gathering thresholds. Bit-score based,
  so chunk-safe and independent of `-Z`/effective DB size.
- **Mgnifams → E-value** (`-E`/`--domE`, default `0.001`): the Mgnifams library has no GA/TC/NC
  lines, so `--cut_ga` is not possible. `-Z`/`--domZ` pin the global effective DB size across chunks
  (set `--effective_db_size` to the total input sequence count; `0` omits `-Z`).

Search-time filters are authoritative, so the stats step applies no further E-value threshold.

## Layout

```
uniprot_annotation_percentages/
├── main.nf
├── nextflow.config
├── conf/{base,modules,test}.config
├── bin/
│   ├── calculate_stats_from_domtbl.py   # FASTA + domtblouts -> stats CSV
│   ├── compare_annotation_stats.py      # before/after comparison (salvaged)
│   ├── provenance_report.py             # DB md5s + tool/run versions
│   └── tests/test_calculate_stats_from_domtbl.py
└── modules/
    ├── nf-core/hmmer/hmmsearch/         # salvaged
    └── local/{calculate_stats_from_domtbl,compare_annotation_stats,provenance_report}/
```

## Run

```bash
export NXF_SINGULARITY_CACHEDIR="/home/vangelis/Desktop/Tools/singularity"
nextflow run main.nf -profile singularity \
  --input_fasta   /home/vangelis/Desktop/Databases/uniprot/uniprot_sprot.2026_02.fasta.gz \
  --pfam_hmm      /home/vangelis/Desktop/Databases/pfam/Pfam-A.38.2.hmm.gz \
  --mgnifams_hmm  /home/vangelis/Desktop/Databases/mgnifams/mgnifams_hmm.lib.gz \
  --effective_db_size <total SwissProt seq count> \
  --outdir output/uniprot_annotation
```

Test run (100 seqs vs 14 Pfam + 108 Mgnifams models, output to `tmp/`):

```bash
export NXF_SINGULARITY_CACHEDIR="/home/vangelis/Desktop/Tools/singularity"
nextflow run main.nf -profile test,singularity
```

## Outputs

- `pfam_annotation_stats.csv` — Pfam-only, plus Pfam-**exclusive** residues/sequences (space
  Mgnifams misses) in `exclusive_*` columns.
- `mgnifam_annotation_stats.csv` — Mgnifams-only, plus Mgnifams-**exclusive** residues/sequences
  (space Pfam misses) in `exclusive_*` columns.
- `pfam_mgnifam_annotation_stats.csv` — Pfam + Mgnifams union.
- `annotation_percentage_increase.csv` — Pfam vs Pfam+Mgnifams (percentage-point + relative increase).
- `provenance.txt` — DB basenames/sizes/md5sums + collected tool versions + run metadata.
- `versions.yml` — collected tool versions (`versions` topic channel).
