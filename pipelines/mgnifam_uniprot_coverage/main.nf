//
// MODULE: Local modules, one per process
//
include { BUILD_CATEGORY_LISTS                        } from './modules/local/build_category_lists/main.nf'
include { COVERAGE_MAP as MAP_PFAM                    } from './modules/local/coverage_map/main.nf'
include { COVERAGE_MAP as MAP_TOTAL                   } from './modules/local/coverage_map/main.nf'
include { COVERAGE_MAP as MAP_EXCLUSIVE               } from './modules/local/coverage_map/main.nf'
include { COVERAGE_REDUCE                             } from './modules/local/coverage_reduce/main.nf'
include { VALIDATE_COVERAGE                           } from './modules/local/validate_coverage/main.nf'
include { COVERAGE_REPORT                             } from './modules/local/coverage_report/main.nf'
include { COVERAGE_FIGURES                            } from './modules/local/coverage_figures/main.nf'
include { PROVENANCE_REPORT                           } from './modules/local/provenance_report/main.nf'

// Directory columns in the samplesheet may be relative; they resolve against the
// sheet's own location so a samplesheet can travel with its data.
def resolveDir(String base, String path) {
    return path.startsWith('/') ? file(path) : file("${base}/${path}")
}

workflow {

    // ---------------------------------------------------------------- inputs

    // One row per database subset. n_chunks is not decoration: it is the
    // --expect-chunks guard, the thing that stops a partial chunk set from
    // reducing cleanly into an understated whole-database percentage.
    def sheet_dir = file(params.input).parent.toString()

    ch_subsets = channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            [
                subset          : row.subset,
                mgnifams_dir    : resolveDir(sheet_dir, row.mgnifams_domtbl_dir),
                pfam_dir        : resolveDir(sheet_dir, row.pfam_domtbl_dir),
                total_sequences : row.total_sequences as long,
                total_residues  : row.total_residues as long,
                n_chunks        : row.n_chunks as int,
            ]
        }

    ch_lists_dir = channel.value( file(params.lists_dir,    checkIfExists: true) )
    ch_hmm_lib   = channel.value( file(params.mgnifams_hmm, checkIfExists: true) )

    // The five curated lists plus the three derived ones (tm_any, tm_or_disorder,
    // not_tm_disorder) and the library size that the enrichment figure divides by.
    BUILD_CATEGORY_LISTS(
        channel.value( [ id: 'mgnifams' ] ).combine(ch_lists_dir).combine(ch_hmm_lib)
    )
    ch_lists = BUILD_CATEGORY_LISTS.out.lists.map { _m, d -> d }.first()

    // ------------------------------------------------------- pair the chunks

    // Key on the chunk stem so the two searches are joined by identity rather
    // than by rebuilding one path from the other. A chunk missing from either
    // side then fails the join instead of silently dropping the mask.
    ch_pfam_chunks = ch_subsets
        .flatMap { s -> files("${s.pfam_dir}/*_pfam.domtbl.gz").collect { f ->
            [ [ s.subset, f.name - ~/_pfam\.domtbl\.gz$/ ], s, f ]
        } }
    ch_mg_chunks = ch_subsets
        .flatMap { s -> files("${s.mgnifams_dir}/*_mgnifams.domtbl.gz").collect { f ->
            [ [ s.subset, f.name - ~/_mgnifams\.domtbl\.gz$/ ], s, f ]
        } }

    ch_pairs = ch_mg_chunks
        .map { key, s, f -> [ key, s, f ] }
        .join( ch_pfam_chunks.map { key, _s, f -> [ key, f ] }, failOnMismatch: true, failOnDuplicate: true )
        .map { key, s, mg, pfam -> [ key, s, mg, pfam ] }

    // ali reproduces the published figures; env is the sensitivity check.
    ch_coords = channel.fromList(params.coords instanceof List ? params.coords : params.coords.tokenize(','))

    ch_work = ch_pairs.combine(ch_coords)
        .map { key, s, mg, pfam, coords ->
            def meta = [
                id             : key[0],
                chunk          : key[1],
                coords         : coords,
                total_sequences: s.total_sequences,
                total_residues : s.total_residues,
                n_chunks       : s.n_chunks,
            ]
            [ meta, mg, pfam ]
        }

    // ------------------------------------------------------------ map stages

    // Pfam first: its per-target spans are the mask everything else subtracts.
    MAP_PFAM(
        ch_work.map { meta, _mg, pfam -> [ meta + [ pass: 'pfam' ], pfam, [], [] ] }
    )
    // Total MGnifam coverage, categorised.
    MAP_TOTAL(
        ch_work.combine(ch_lists).map { meta, mg, _pfam, lists -> [ meta + [ pass: 'total' ], mg, lists, [] ] }
    )
    // MGnifam-exclusive: same chunk, same coords, Pfam spans removed.
    ch_masks = MAP_PFAM.out.pertarget.map { meta, pt -> [ [ meta.id, meta.chunk, meta.coords ], pt ] }
    MAP_EXCLUSIVE(
        ch_work.map { meta, mg, _pfam -> [ [ meta.id, meta.chunk, meta.coords ], meta, mg ] }
            .join( ch_masks, failOnMismatch: true, failOnDuplicate: true )
            .map { _key, meta, mg, pt -> [ meta + [ pass: 'exclusive' ], mg, pt ] }
            .combine(ch_lists)
            .map { meta, mg, pt, lists -> [ meta, mg, lists, pt ] }
    )

    // --------------------------------------------------------- reduce stages

    ch_mapped = MAP_PFAM.out.summary.join(MAP_PFAM.out.families)
        .mix( MAP_TOTAL.out.summary.join(MAP_TOTAL.out.families) )
        .mix( MAP_EXCLUSIVE.out.summary.join(MAP_EXCLUSIVE.out.families) )

    // Per subset, and again with the subsets pooled into the whole of UniProtKB.
    // Pooling is arithmetically sound because hmmsearch chunks partition the
    // target database, so no sequence is counted in two subsets.
    ch_per_subset = ch_mapped
        .map { meta, sm, fam -> [ [ meta.id, meta.pass, meta.coords ], meta, sm, fam ] }
        .groupTuple()
        .map { key, metas, sms, fams ->
            def m = metas[0]
            [ [ id: key[0], pass: key[1], coords: key[2],
                total_sequences: m.total_sequences, total_residues: m.total_residues,
                n_chunks: m.n_chunks ], sms, fams ]
        }

    ch_pooled = ch_mapped
        .map { meta, sm, fam -> [ [ meta.pass, meta.coords ], meta, sm, fam ] }
        .groupTuple()
        .map { key, metas, sms, fams ->
            // one meta per subset, deduplicated, then summed
            def by_subset = metas.collectEntries { [ (it.id): it ] }.values()
            [ [ id: 'uniprotkb', pass: key[0], coords: key[1],
                total_sequences: by_subset.sum { it.total_sequences },
                total_residues : by_subset.sum { it.total_residues },
                n_chunks       : by_subset.sum { it.n_chunks } ], sms, fams ]
        }

    COVERAGE_REDUCE(
        ch_per_subset.mix(ch_pooled)
            .combine(ch_lists)
            .map { meta, sms, fams, lists ->
                // the Pfam pass has no MGnifam categories to report on
                [ meta, sms, fams, meta.pass == 'pfam' ? [] : lists ]
            }
    )

    // ----------------------------------------------- validate, then report

    // The reduced tables are named <id>_<pass>_<coords>.reduced.tsv, so the view
    // identity travels with the file and needs no parallel key channel.
    ch_reduced = COVERAGE_REDUCE.out.reduced
        .map { _meta, tsv -> tsv }
        .collect()
        .map { tsvs -> [ [ id: 'coverage' ], tsvs ] }

    // An optional path has to be appended explicitly: combining with a channel
    // carrying an empty list contributes nothing and silently shortens the tuple.
    def reference = params.reference_csv ? file(params.reference_csv, checkIfExists: true) : []

    // Hard assertions first: nothing downstream should render numbers that have
    // not been checked for internal consistency.
    VALIDATE_COVERAGE( ch_reduced.map { m, tsvs -> [ m, tsvs, reference ] } )

    ch_sizes      = BUILD_CATEGORY_LISTS.out.sizes.map { _m, s -> s }.first()
    ch_validation = VALIDATE_COVERAGE.out.ok.map { _m, f -> f }
    ch_sheet      = channel.value( file(params.input) )

    // Gated on the validation output, so a report can never be rendered from
    // numbers that failed a check.
    COVERAGE_REPORT(
        ch_reduced.combine(ch_sizes).combine(ch_validation),
        ch_sheet
    )
    COVERAGE_FIGURES( ch_reduced.combine(ch_sizes), ch_sheet )

    // ------------------------------------------------------------ provenance

    ch_versions = channel.topic('versions')
        .map { row -> "${row[1]}: ${row[2]}" }
        .unique()
        .collectFile(name: 'versions.yml', newLine: true, sort: true)

    def workflow_json = groovy.json.JsonOutput.toJson([
        nextflow_version : workflow.nextflow.version.toString(),
        pipeline_revision: (workflow.revision ?: 'N/A').toString(),
        commit_id        : (workflow.commitId ?: 'N/A').toString(),
        command_line     : workflow.commandLine.toString(),
        start            : workflow.start.toString(),
        container_engine : (workflow.containerEngine ?: 'none').toString(),
        profile          : workflow.profile.toString(),
        coords           : params.coords.toString(),
        samplesheet      : params.input.toString(),
    ])

    // Checksum the individual list files, not the directory holding them: the
    // provenance script hashes file contents.
    def prov_inputs = [ file(params.input), file(params.mgnifams_hmm) ] + files("${params.lists_dir}/*.txt")
    ch_provenance = channel.value( [ [ id: 'coverage' ], prov_inputs ] )
        .combine( ch_versions )
    PROVENANCE_REPORT( ch_provenance, workflow_json )
}
