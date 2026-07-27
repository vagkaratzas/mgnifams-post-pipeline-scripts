# PLAN — Domain architectures from the MGnify re-annotation

Behaviour is frozen in [SPEC.md](SPEC.md). This file tracks execution.

**Iron law: no production code without a failing test first.** Each box below is one
red → green → refactor cycle. Watch the test fail for the right reason, write the minimum to pass,
then move on. After each box:

```bash
python -m pytest tests/test_parse_domain_architectures.py -q
```

## Phase 0 — docs and fixtures

- [x] Write `SPEC.md` and `PLAN.md` at the repo root
- [ ] `assets/test_data/domain_architecture/clan_membership_dummy.csv` — a multi-family clan, a
      second clan, a singleton clan, and one family that gets zero hits
- [ ] `assets/test_data/domain_architecture/pfam_mapping_dummy.tsv` — deliberately omitting one
      accession used by the dummy proteins
- [ ] `assets/test_data/domain_architecture/proteins_dummy.csv.gz` — covering: a row with no `"m"` ·
      three same-clan same-region hits · two same-clan disjoint hits · two different-clan
      overlapping hits · a Pfam whose `ali_from` precedes the MGnifam while its `hmm_from` is 1 ·
      a repeated Pfam · an accession missing from the mapping · a 50%-vs-51% overlap pair ·
      equal-start different-length hits · equal-start equal-length hits

## Phase 1 — pure functions (`bin/parse_domain_architectures.py`)

- [ ] `load_clan_membership` returns family → clan and clan → rep
- [ ] `overlaps` is True at 51%, **False at exactly 50%**
- [ ] `cluster_hits` merges three same-clan same-region hits into one group
- [ ] `cluster_hits` keeps same-clan disjoint hits as two groups
- [ ] `cluster_hits` never merges different-clan hits, however much they overlap
- [ ] `cluster_hits` chains A–B and B–C into one group when A and C do not overlap (single linkage)
- [ ] `mgnifam_chip` with >1 distinct family → `SF_233` / `MGnifam clan 233` / rep link
- [ ] `mgnifam_chip` with 1 distinct family → `470` / `MGnifam470` / own link
- [ ] `pfam_chip` falls back to the raw accession when it is missing from the mapping (bug 6)
- [ ] `string_to_hex_color` parity with the old output: `"ATP synthase alpha/beta family,
      nucleotide-binding domain"` → `#8b5115`, asserted against the existing asset JSON
- [ ] `decide_font_color` returns `white` below the 0.2 luminosity cut and `black` above it

## Phase 2 — per-row assembly

- [ ] `build_chips` puts a Pfam before the MGnifam when its `ali_from` is smaller — **bug 1
      regression; uses `p[5]`, not `p[3]`**
- [ ] `build_chips` orders equal-start chips shortest-first
- [ ] `build_chips` orders equal-start equal-length chips alphabetically by id
- [ ] `architecture_key` preserves repeats: the same Pfam twice → two entries
- [ ] a row at 50% overlap and a row at 51% overlap yield two different keys

## Phase 3 — streaming and counting

- [ ] rows without `"m"` are skipped
- [ ] a family is credited for a sequence even when its own hit merged into a clan chip
- [ ] every distinct family id on a row is credited exactly once, even when it hits twice
- [ ] `--no-prefilter` and the `zcat | grep` path produce identical output

## Phase 4 — output and validation

- [ ] `write_family_json` emits the exact schema, top-N, count-descending
- [ ] a family in the clan file with zero hits gets `{"architecture_containers": []}`
- [ ] that family is listed in `missing_families.txt` and a WARNING is logged
- [ ] end-to-end CLI run on the dummy gz produces the expected JSON set

## Phase 5 — `bin/update_domain_blobs.py`

- [ ] updates `domain_blob` for every db id with a matching JSON
- [ ] leaves rows with no JSON untouched and reports the count
- [ ] reports JSON files with no matching db row
- [ ] rejects a `--column` absent from `PRAGMA table_info` instead of interpolating it

## Phase 6 — wrap-up

- [ ] README sections for both scripts, in the existing style, under stage 5 — short description
      plus a CLI call each; include the "do not run `assign_correct_domain_ids.py` on these JSONs"
      warning
- [ ] full `python -m pytest -q` green
- [ ] every box above ticked

## Verification

```bash
# 1. unit tests
python -m pytest tests/test_parse_domain_architectures.py tests/test_update_domain_blobs.py -q

# 2. end to end on dummy data
python bin/parse_domain_architectures.py \
  --proteins assets/test_data/domain_architecture/proteins_dummy.csv.gz \
  --clan-membership assets/test_data/domain_architecture/clan_membership_dummy.csv \
  --pfam-mapping assets/test_data/domain_architecture/pfam_mapping_dummy.tsv \
  --output-dir /tmp/dom_out
cat /tmp/dom_out/missing_families.txt

# 3. load into a scratch copy of the test db, confirm the blob round-trips
cp assets/mgnifams_v2_results/mgnifams_test.sqlite3 /tmp/t.sqlite3
python bin/update_domain_blobs.py --db /tmp/t.sqlite3 --json-dir /tmp/dom_out
python -c "import sqlite3, json; print(json.loads(sqlite3.connect('/tmp/t.sqlite3').execute(
  'select domain_blob from mgnifam where domain_blob is not null limit 1').fetchone()[0]))"

# 4. visual check against the real viewer
#    point mgnifams-site at /tmp/t.sqlite3, open a details page, confirm the
#    Domain architecture card renders clan chips bold and in ali_from order
```

Real run:

```bash
python bin/parse_domain_architectures.py \
  --proteins /nfs/production/rdf/metagenomics/users/vangelis/mgnifams_annotation_coverage/proteins_mgnifams.csv.gz \
  --clan-membership assets/mgnifams_v2_results/generate_families/network/clan_membership.csv \
  --pfam-mapping /path/to/pfam_mapping.tsv \
  --output-dir output/domain_results
```

The validation line must report 35459 families written and an empty `missing_families.txt`.
