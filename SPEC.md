# SPEC — Domain architectures from the MGnify re-annotation

Status: frozen. Implementation plan lives in [PLAN.md](PLAN.md).

## 1. Context

Every MGnifam page carries a "Domain architecture" card (`details.html:415-428`, rendered by
`renderArchitecture` in `explorer/static/explorer/js/details.js:214`). Its JSON used to be produced
during the pipeline by `mgnifams/bin/parse_domains.py` from per-family TSVs queried out of the old
MGnify protein database. That input no longer exists.

The replacement is a single re-annotation over all ~718M MGnify sequences:

```
mgyp,sequence,full_length,cluster_size,metadata

metadata = {"s": …,
            "b": …,
            "p": [[pfam_acc, evalue, score, hmm_from, hmm_to, ali_from, ali_to], …],
            "m": [[family_id, i_evalue, score, ali_from, ali_to], …]}
```

`"m"` is written by `pipelines/annotation_percentages/bin/append_mgnifams_annot.py` and is already
filtered at independent domain E-value ≤ 0.001. It now carries **many** MGnifam hits per sequence —
the result of hmmsearching all ~35K MGnifam HMMs against every sequence — because highly similar
families all land on the same region. One observed row has 22 hits covering ~25-105. The old
"one MGnifam plus some Pfams" model cannot render that.

The fix is to collapse co-located MGnifam hits by **clan**, using
`assets/mgnifams_v2_results/generate_families/network/clan_membership.csv`, which covers all 35459
families across 19988 clans (15979 of them singletons).

Not all ~718M sequences made it into a family; only rows carrying an `"m"` annotation are relevant.

## 2. Field layout, and why `p[5]` is the start

`"p"` has seven fields, carrying **both** HMM and alignment coordinates:

```
["PF02954", 0.0019, 24.5, 8, 37, 126, 160]
                          ^^^^^  ^^^^^^^^
                          hmm    ali
```

PF02954 / HTH_8 is a ~42-aa model and the sequence is 164 aa, so `8,37` can only be HMM coordinates
and `126,160` can only be alignment coordinates. Therefore `p[3:5]` is `hmm_from,hmm_to` and
`p[5:7]` is `ali_from,ali_to`.

`"m"` has five fields and only alignment coordinates: `m[3:5]` is `ali_from,ali_to`.

The old script sorted Pfams by `p[3]`. See §7 bug 1.

## 3. Chip model

One sequence produces an ordered list of **chips**, which is reduced to an **architecture key**.

### 3.1 Pfam chips

Every `"p"` hit becomes exactly one chip spanning `p[5]`–`p[6]`. No E-value filtering, no merging,
no overlap resolution. Repeats are preserved; a Pfam hitting twice yields two chips.

### 3.2 MGnifam chips

`"m"` hits are partitioned by clan. Within a clan, hits are clustered by **single linkage** using:

```
overlap(a, b) > overlap_fraction * min(len(a), len(b))
```

with `overlap_fraction` defaulting to `0.5` and the comparison **strictly greater**. Hits in
different clans never merge, however much they overlap — both chips show.

Single linkage means A–B and B–C merging pulls A and C into the same group even when A and C do not
overlap each other.

Each resulting group becomes one chip:

| Group | id | name | link |
|---|---|---|---|
| >1 distinct family id | `SF_233` | `MGnifam clan 233` | `<base>/MGYF0000021329` (clan's `Family Rep Id`) |
| 1 distinct family id | `470` | `MGnifam470` | `<base>/MGYF0000000470` |

Labelling is **per merged group**, not per sequence, so one row can mix both forms:

```
m: 470(10-100), 22289(150-250), 21329(155-255)     # all clan SF_233
   22289 + 21329 overlap; 470 is disjoint
→ [MGnifam470][MGnifam clan 233]
```

A family in a singleton clan can never merge, so it always renders as `MGnifam<id>`.

The strict `>` boundary is deliberate and load-bearing: two otherwise identical sequences, one at
50% overlap and one at 51%, must yield two *different* architectures.

### 3.3 Ordering

All chips — Pfam and MGnifam alike — are sorted by `(ali_from, ali_to, id)`:

1. `ali_from` ascending: whichever appears at the earliest amino acid shows first;
2. then `ali_to` ascending: with equal starts this is shortest-first;
3. then `id` ascending as a plain string: deterministic last resort, may never fire in practice.

### 3.4 Architecture key

Chip ids joined with `\t` and `sys.intern`ed. Repeats are significant: `MGnifam1\tMGnifam1` is a
different architecture from `MGnifam1`.

### 3.5 Tally

A family's page counts **every sequence carrying an `"m"` hit to that family**, whether or not that
family's own hit ended up folded into a clan chip. A sequence hit twice by the same family still
counts once for it.

## 4. Output

`<output-dir>/<family_id>.json`, schema unchanged from the old pipeline so the viewer needs no
edits:

```json
{
  "architecture_containers": [
    {
      "architecture_text": "<count>",
      "domains": [
        {
          "id": "SF_233",
          "color": "#…",
          "link": "http://mgnifams-demo.mgnify.org/details/MGYF0000021329",
          "name": "MGnifam clan 233",
          "font_color": "white"
        }
      ]
    }
  ]
}
```

Top `--top` (default 15) architectures by count, descending. `color = string_to_hex_color(name)` for
every chip type — one rule, unlike the old script (§7 bug 8). `font_color` from `decide_font_color`.

## 5. Compatibility with the site

No site changes required.

- `MGnifam clan 233` contains the substring `MGnifam`, so the bold-and-pad branch at
  `details.js:237` still fires for clan chips. It is 18 characters, under the 30-character
  truncation limit at `details.js:228`, so no tooltip is needed.
- `SF_233` does not match `PFAM_RE` in `bin/annotate_novel_through_domain_architecture.py:81`, so
  that downstream consumer keeps extracting exactly the Pfam accessions it did before.
- `bin/assign_correct_domain_ids.py` rewrites any domain whose `name` contains `MGnifam`. It is an
  old-pipeline redundancy-removal fixup and **must not be run against these JSONs** — it would
  flatten clan chips back into family chips.

## 6. Scripts

### 6.1 `bin/parse_domain_architectures.py` (new)

```
--proteins           proteins_mgnifams.csv[.gz]
--clan-membership    clan_membership.csv
--pfam-mapping       pfam_mapping.tsv
--output-dir         DIR
--overlap-fraction   0.5
--top                15
--base-url           http://mgnifams-demo.mgnify.org/details/
--log-every          1000000
--no-prefilter
```

1. Load `clan_membership.csv` (`Cluster Id`, `Family Rep Id`, `Family Ids`) into `family → clan`,
   `clan → rep`, and the **expected-family set**. Load `pfam_mapping.tsv` into `acc → name`, with
   a `.get(acc, acc)` fallback for accessions missing after a Pfam version bump.
2. Stream the proteins file. The header is read from the gzip directly; the body is piped through
   `zcat … | LC_ALL=C grep -F '""m"":'` so decompression and the discard of every `"m"`-less row
   happen in C, and only survivors reach `json.loads`. `--no-prefilter`, and any non-`.gz` input,
   falls back to an in-Python substring check on the raw line. Both paths must produce identical
   output. `csv.field_size_limit` is raised as in `append_mgnifams_annot.py:10`.
3. Build chips per row, sort, intern the key, and increment `counters[family][key]` for every
   distinct family id in `"m"`.
4. Write one JSON per family.
5. **Validate.** Diff the families written against the expected set from `clan_membership.csv`.
   Every missing family still gets a valid `{"architecture_containers": []}` file so the database
   update cannot fail, its id is listed in `<output-dir>/missing_families.txt`, and a WARNING is
   logged with the count.

Logging (INFO): clans and Pfams loaded, rows read / rows kept every `--log-every`, families written,
families missing, elapsed time.

Memory ceiling: `dict[int, Counter[str]]` over ~35K families is held in RAM for the whole pass.
Marked in the source with a `ponytail:` comment naming the upgrade path (emit `family\tarch` to
disk, `sort | uniq -c`) if it ever stops fitting.

### 6.2 `bin/update_domain_blobs.py` (rewritten in place)

The existing file targets a `domain_architecture_blob` column and reads a `domain_architecture_file`
column. Neither exists in the current schema — both `assets/mgnifams_v2_results/mgnifams_test.sqlite3`
and the site database have a plain `domain_blob`. The script is dead against today's database, so it
is rewritten rather than duplicated.

```
--db DB --json-dir DIR [--table mgnifam] [--column domain_blob]
```

`SELECT id FROM <table>`, read `<json-dir>/<id>.json`, one connection, parameterised
`UPDATE … SET <column> = ?`, single commit. Table and column names are validated against
`PRAGMA table_info` before being interpolated into the statement. Logs rows updated, database rows
with no JSON file, and JSON files with no database row.

## 7. Bugs in the superseded scripts

Recorded here because the rewrite has to not reproduce them. The old repo is not modified.

1. **`parse_domains.py:56` sorts Pfams by `pfam[3]` = `hmm_from`, not `ali_from` (`pfam[5]`).** This
   is why the MGnifam nearly always rendered first: `hmm_from` is almost always 1-3, while
   `calculate_mgnifam_start` returns `1` for a whole-protein member, so nearly every chip tied at
   ~1 and `sorted(zip(start_points, fam_names))` broke the tie on the *name* string, where
   `"1" < "PF…"`. A Pfam genuinely upstream of the MGnifam rendered after it.
2. `calculate_mgnifam_start` (`:29`) leaves `start` unbound when the protein name has more than
   three underscores → `UnboundLocalError`; it also returns `str` from one branch and `int` from
   the others.
3. `construct_domain_architecture:60` — `zip(*sorted_data)` raises `ValueError` when both the Pfam
   list and the MGnifam-start list are empty.
4. `construct_solo_domain_architecture:66` with empty `mgnifam_starts` yields `""`, counted as a
   legitimate and typically top-ranking architecture.
5. `count_domain_architectures:79` — a sequence whose MGYP is absent from `mgyp_lookup` still
   contributes a Pfam-only architecture to that family's tally.
6. `translate_architecture:175` — `pfam_mapping[domain["id"]]` raises `KeyError` on any accession
   missing from the mapping file, which a Pfam version bump guarantees.
7. `subset_json:119` — the top-15 cut is a hardcoded default, never exposed as a CLI argument.
8. Colour is computed from `id` and then silently overwritten from `name` for Pfams, but from `id`
   only for MGnifams — two rules for one thing.
9. `update_sqlite_blobs.py:54` opens, commits and closes a fresh SQLite connection per column per
   row: 8 × N connections.
10. `assign_correct_domain_ids.py:32` writes `details/?id=MGYF…`, but the live route is
    `details/<str:pk>/` (`mgnifams-site/…/explorer/urls.py:11`). The `parse_domains.py:178` form is
    the correct one.

## 8. Test fixtures

`assets/test_data/domain_architecture/` — `proteins_dummy.csv.gz`, `clan_membership_dummy.csv`,
`pfam_mapping_dummy.tsv`, hand-built to exercise every branch above: a row with no `"m"`, three
same-clan same-region hits, two same-clan disjoint hits, two different-clan overlapping hits, a Pfam
whose `ali_from` precedes the MGnifam while its `hmm_from` is 1, a repeated Pfam, an accession
missing from the mapping, a 50%-vs-51% overlap pair, equal-start different-length hits, equal-start
equal-length hits, and a family with zero hits.
