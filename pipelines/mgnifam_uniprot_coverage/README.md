# MGnifam × UniProt Coverage Pipeline

How much of UniProtKB do the MGnifam models cover that **Pfam does not already explain**, at
residue and sequence level, broken down by family category (novel, membrane-α/β, disordered)?

Consumes `hmmsearch --domtblout` chunks that already exist — it runs no search of its own. A new
UniProt release means running [`uniprot_annotation_percentages`](../uniprot_annotation_percentages)
first to produce the domtbls, then this pipeline over them.

## Why a pipeline and not a shell loop

The analysis is one `map` per chunk plus a `reduce`, which is easy to write as a SLURM array and
easy to get quietly wrong. The exclusive pass needs each MGnifams chunk paired with the Pfam
per-target file **for the same chunk**; a hand-built mask path that does not resolve produces
*total* coverage with every column still labelled *exclusive*. Here the two searches are joined on
the chunk key with `failOnMismatch`, so a missing chunk aborts the run instead. Same story for a
partial chunk set, which reduces perfectly cleanly and merely understates every whole-database
percentage — `--expect-chunks` turns that into a failure too.

## Run

```bash
export NXF_SINGULARITY_CACHEDIR="/path/to/singularity"

nextflow run main.nf -profile singularity \
  --input        samplesheet.csv \
  --lists_dir    /path/to/lists \
  --mgnifams_hmm /path/to/mgnifams_hmm.lib.gz \
  --reference_csv /path/to/annotation_percentage_increase.csv \
  --outdir       results
```

Test profile (real hmmsearch output over 100 SwissProt sequences, ~1 minute):

```bash
nextflow run main.nf -profile test,singularity
```

### Samplesheet

One row per database subset. Directory columns may be relative — they resolve against the
samplesheet's own location, so a sheet travels with its data.

```csv
subset,mgnifams_domtbl_dir,pfam_domtbl_dir,total_sequences,total_residues,n_chunks
swissprot,/nfs/.../swissprot/hmmsearch_mgnifams,/nfs/.../swissprot/hmmsearch_pfams,575503,208906902,1
trembl,/nfs/.../trembl/hmmsearch_mgnifams,/nfs/.../trembl/hmmsearch_pfams,149234636,58049358499,150
```

`total_sequences` / `total_residues` are the **whole** database, sequences with no hit included —
neither search sees those, and without them there is no denominator for a percentage-point gain.
`n_chunks` is the `--expect-chunks` guard, not decoration.

### Parameters

| Param | Default | Meaning |
|---|---|---|
| `--input` | — | samplesheet above |
| `--lists_dir` | — | directory of curated `*.txt` MGnifam category lists |
| `--mgnifams_hmm` | — | HMM library; `hmmstat` gives the library-size denominator |
| `--reference_csv` | none | `annotation_percentage_increase.csv` to reconcile against |
| `--coords` | `['ali','env']` | alignment coordinates reproduce the published figures, envelope is the sensitivity check |

## What it does

```
BUILD_CATEGORY_LISTS   5 curated lists + HMM library -> 8 categories + library size
  |
  |  per chunk, per coordinate system:
MAP_PFAM               Pfam coverage; its per-target spans become the mask
MAP_TOTAL              total MGnifam coverage, by category
MAP_EXCLUSIVE          MGnifam coverage minus the Pfam spans of the SAME chunk
  |
COVERAGE_REDUCE        3 views (swissprot, trembl, pooled uniprotkb) x 3 passes x 2 coords
VALIDATE_COVERAGE      hard assertions; the report is gated on it
COVERAGE_REPORT        tables + filled prose, Markdown and self-contained HTML
COVERAGE_FIGURES       four vector PDFs
PROVENANCE_REPORT      input checksums, tool versions, command line
```

Pooling SwissProt and TrEMBL into UniProtKB is arithmetically sound because hmmsearch chunks
partition the target database — no sequence is counted in two subsets.

### The three derived categories

`BUILD_CATEGORY_LISTS` adds `tm_any`, `tm_or_disorder` and `not_tm_disorder` to the five curated
lists, because two questions cannot be answered by adding table rows:

- **Enrichment** needs `tm_or_disorder` as a real **union**. Summing the `membrane_a`,
  `membrane_b` and `disorder` rows double counts residues wherever two such families overlap on
  the same protein.
- **"Excluding these families"** needs `not_tm_disorder`, the complement against the whole
  library. There is no subtraction of residue counts that gives this, for the same reason.

## Outputs

| Path | Content |
|---|---|
| `report/coverage_report.{md,html}` | Tables 1–4 plus the generated prose block |
| `report/validation.txt` | every check, PASS/FAIL |
| `tables/<view>_<pass>_<coords>.reduced.tsv` | the numbers behind every table |
| `tables/*.list_overlaps.tsv` | pairwise category overlap; the lists are not disjoint |
| `figures/fig{1..4}_*.pdf` | vector, 85 mm single column, Times, Wong palette |
| `lists/`, `library_sizes.tsv` | the eight categories and the library denominator |
| `chunks/<subset>/<coords>/<pass>/` | per-chunk summaries and family lists |
| `provenance.txt` | input md5s, tool versions, command line |

Per-target interval files are not published — hundreds of MB at UniProt scale, and they exist only
to be the mask.

## Validation

`VALIDATE_COVERAGE` fails the run on any of:

- exclusive residues/targets exceeding the total pass, per category
- any category exceeding its `any` row
- newly annotated sequences exceeding hit sequences
- newly annotated counts differing between `ali` and `env` (it is a presence/absence quantity, so
  a difference is a bug, not a modelling choice)
- the pooled UniProtKB view not equalling the sum of its subsets
- with `--reference_csv`, the SwissProt `ali` numbers not reproducing the independently computed
  `annotation_percentage_increase.csv` exactly

That last one is the check worth having: it exercises the whole DAG end to end. The engine
reproduces the reference to the residue — 144,906,700 Pfam aa, 14,514,544 exclusive aa, 6.9479 pp,
2,616 newly annotated sequences, 0.4546 pp.

It also warns, without failing, when the Pfam share of the database is far from the reference's —
the signature of a partial chunk set.

## Cost

Measured, per chunk, single core:

    peak RSS ≈ 0.09 KB × kept rows + 0.55 KB × distinct targets   (≈ 1 GB per 8 M rows)
    wall     ≈ 4 s per 1 M rows (+ gzip decode)

`--mask` adds ~0.25 GB, scaling with the *mask's* target count, and `not_tm_disorder` roughly
doubles the interval set. Largest real chunk seen: 8.2 M rows, 374 k targets → 1.28 GB, 47 s.
`process_low` (12 GB) is ample. A full run is ~900 short tasks — 151 chunks × 3 passes × 2
coordinate systems — a few CPU-hours.

## Layout

```
mgnifam_uniprot_coverage/
├── main.nf
├── nextflow.config
├── conf/{base,modules,test}.config
├── assets/test_data/            # real hmmsearch output, 3 chunks, toy lists
├── bin/
│   ├── mgnifam_uniprot_coverage_stats_from_domtbl.py   # map/reduce engine
│   ├── validate_coverage.py
│   ├── coverage_report.py
│   ├── coverage_figures.py
│   ├── provenance_report.py
│   └── tests/
└── modules/local/{build_category_lists,coverage_map,coverage_reduce,
                   validate_coverage,coverage_report,coverage_figures,provenance_report}/
```

Unit tests: `python3 -m pytest bin/tests/ -q`.
