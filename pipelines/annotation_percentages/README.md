# Annotation Percentages Pipeline

Mini local Nextflow workflow for estimating how much residue-level annotation coverage increases when MGnifam `hmmsearch` hits are appended to a Pfam-annotated MGnify protein CSV.

It loosely follows nf-core conventions: one module per process, a Conda `environment.yml` and a Docker/Singularity container per module, resource `label`s, the nf-core `hmmer/hmmsearch` module under `modules/nf-core/`, and centralised publishing/arguments in `conf/modules.config`.

## Layout

```
annotation_percentages/
├── main.nf                       # workflow: wires the modules together
├── nextflow.config               # params + container/conda profiles
├── conf/
│   ├── base.config               # resource labels + resourceLimits
│   └── modules.config            # publishDir + per-module ext.args / ext.prefix
├── bin/                          # python scripts staged onto the task PATH
├── modules/local/
│   ├── extract_fasta/            # CSV/CSV.GZ -> FASTA
│   ├── append_mgnifams/          # append (gzipped) domtbl hits into metadata["m"]
│   ├── calculate_annotation_stats/  # used twice (pfam, pfam+mgnifam)
│   └── compare_annotation_stats/    # before/after comparison
└── modules/nf-core/
    └── hmmer/hmmsearch/          # nf-core HMMER_HMMSEARCH module
```

## Run

```bash
nextflow run pipelines/annotation_percentages/main.nf -profile singularity
```

Swap `-profile singularity` for `docker`, `conda`, or `mamba` as available.

Useful overrides:

```bash
nextflow run pipelines/annotation_percentages/main.nf -profile singularity \
  --input_csv /path/to/sequence_explorer_protein.csv.gz \
  --hmm_lib assets/mgnifams_v2_results/generate_families/families/mgnifams_hmm.lib.gz \
  --outdir output/annotation_percentages \
  --evalue_cutoff 0.001 \
  --fasta_records_per_chunk 100000 \
  --effective_db_size 717738164
```

`--input_csv` takes a local `.csv` or `.csv.gz` path; the python scripts read either
transparently. `hmmsearch` (HMMER 3.4) reads the gzipped HMM library directly, so no
decompression step is needed.

The extracted FASTA is split with Nextflow's `splitFasta` operator before `hmmsearch`, so
the HMMER searches can run in parallel. `--fasta_records_per_chunk` controls the chunk size.
`--effective_db_size` is passed to `hmmsearch` `-Z` and `--domZ` so chunked searches use the
global effective database size for E-value calculation; set it to `0` to omit `-Z`.

## Outputs

- `mgnify_proteins.fasta`: FASTA extracted from the input CSV.
- `mgnify_proteins_chunk_*_mgnifams.domtbl.gz`: chunked `hmmsearch --domtblout` results (gzipped by the nf-core module).
- `mgnify_proteins_chunk_*_mgnifams.txt.gz`: chunked full `hmmsearch` human-readable outputs.
- `proteins_mgnifams.csv.gz`: input CSV with `metadata["m"]` MGnifam annotations replaced from the domtblout.
- `pfam_annotation_stats.csv`: Pfam-only annotation statistics.
- `pfam_mgnifam_annotation_stats.csv`: Pfam + MGnifam annotation statistics.
- `annotation_percentage_increase.csv`: percentage-point and relative annotation increases.
