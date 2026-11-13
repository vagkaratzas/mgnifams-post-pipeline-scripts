# Run Prediction Annotations Pipeline

Mini local Nextflow workflow that runs two deep-learning protein annotation predictors —
[FUNGTION](https://github.com/nf-core/modules) and [CAALM](https://huggingface.co/lczong/CAALM) —
over a single amino acid FASTA file.

It follows nf-core conventions: the `fungtion/fungtion`, `fungtion/downloadmodels`,
`caalm/caalm` and `caalm/downloadmodels` modules live under `modules/nf-core/`, resources come
from `label`s in `conf/base.config`, and publishing/arguments are centralised in
`conf/modules.config`. Both predictor modules carry a `process_gpu` label, so an
`accelerator` is requested (and the CUDA container/conda environment is selected) whenever the
pipeline is run with `-profile gpu`.

## Layout

```
run_prediction_annotations/
├── main.nf                       # workflow: wires the modules together
├── nextflow.config               # params + container/conda/gpu profiles
├── conf/
│   ├── base.config               # resource labels (incl. process_gpu) + resourceLimits
│   └── modules.config            # publishDir + per-module ext.args
└── modules/nf-core/
    ├── fungtion/downloadmodels/  # downloads the ESM-1b pretrain checkpoint
    ├── fungtion/fungtion/        # function prediction
    ├── caalm/downloadmodels/     # downloads the 3-level CAALM models from HuggingFace
    └── caalm/caalm/              # hierarchical annotation prediction
```

## Run

```bash
nextflow run pipelines/run_prediction_annotations/main.nf \
  -profile singularity,gpu,slurm \
  --input_fasta /path/to/proteins.fasta
```

Drop `gpu` from the profile list to run on CPU (CPU containers/environments are used, and no
accelerator is requested). Swap `singularity` for `docker`, `conda`, or `mamba` as available.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input_fasta` | `null` | **Required.** Amino acid FASTA. |
| `--prefix` | FASTA basename | Output file prefix / `meta.id`. |
| `--fungtion_models` | `null` | Directory holding `esm1b_t33_650M_UR50S.pt`. If set, `FUNGTION_DOWNLOADMODELS` is skipped. |
| `--caalm_models` | `null` | Directory holding `level0/`, `level1/`, `level2/`. If set, `CAALM_DOWNLOADMODELS` is skipped. |
| `--fungtion_args` | `''` | Extra `fungtion` CLI args (e.g. `--html-report`, `--skip-visualization`). |
| `--caalm_args` | `''` | Extra `caalm` CLI args (e.g. `--save-level0-embeddings`). |
| `--outdir` | `output/run_prediction_annotations` | Publish directory. |

Model paths always take priority: when `--fungtion_models` / `--caalm_models` are given, the
corresponding download process is not executed at all. On a first run, leave them unset — the
downloaded models are published to `<outdir>/models/{fungtion,caalm}` and can be passed back in
on later runs:

```bash
nextflow run pipelines/run_prediction_annotations/main.nf \
  -profile singularity,gpu,slurm \
  --input_fasta /path/to/proteins.fasta \
  --fungtion_models output/run_prediction_annotations/models/fungtion/models \
  --caalm_models output/run_prediction_annotations/models/caalm/models
```

## Outputs

Under `<outdir>/fungtion/`:

- `<prefix>/<prefix>.csv`: FUNGTION predictions.
- `<prefix>/<prefix>_analysis`, `<prefix>/<prefix>.html`, `<prefix>/<prefix>_assets`: optional visualization/report outputs (`--fungtion_args`).
- `<prefix>.log`: FUNGTION log.

Under `<outdir>/caalm/`:

- `<prefix>_predictions.tsv`: CAALM predictions.
- `<prefix>_probabilities.jsonl`: per-label probabilities.
- `<prefix>_statistics.tsv`: run statistics.
- `<prefix>_level{0,1,2}_embeddings.npy`: optional embeddings (`--caalm_args`).
- `<prefix>.log`: CAALM log.
