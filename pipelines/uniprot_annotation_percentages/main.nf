//
// MODULE: Local modules, one per process
//
include { CALCULATE_STATS_FROM_DOMTBL as STATS_PFAM         } from './modules/local/calculate_stats_from_domtbl/main.nf'
include { CALCULATE_STATS_FROM_DOMTBL as STATS_MGNIFAM      } from './modules/local/calculate_stats_from_domtbl/main.nf'
include { CALCULATE_STATS_FROM_DOMTBL as STATS_PFAM_MGNIFAM } from './modules/local/calculate_stats_from_domtbl/main.nf'
include { COMPARE_ANNOTATION_STATS                          } from './modules/local/compare_annotation_stats/main.nf'
include { PROVENANCE_REPORT                                 } from './modules/local/provenance_report/main.nf'

//
// MODULE: Installed directly from nf-core/modules
//
include { HMMER_HMMSEARCH as HMMSEARCH_PFAM     } from './modules/nf-core/hmmer/hmmsearch/main.nf'
include { HMMER_HMMSEARCH as HMMSEARCH_MGNIFAMS } from './modules/nf-core/hmmer/hmmsearch/main.nf'

workflow {
    def meta = [ id: file(params.input_fasta).simpleName ]

    ch_fasta        = channel.value( file(params.input_fasta,   checkIfExists: true) )
    ch_pfam_hmm     = channel.value( file(params.pfam_hmm,      checkIfExists: true) )
    ch_mgnifams_hmm = channel.value( file(params.mgnifams_hmm,  checkIfExists: true) )

    // Split the input FASTA into chunks so the hmmsearches parallelise; parent_id regroups them.
    ch_chunks = ch_fasta.flatMap { fasta ->
        fasta.splitFasta(by: params.fasta_records_per_chunk.toInteger(), file: true).withIndex().collect { chunk_file, chunk_index ->
            def chunk_meta = meta + [ id: "${meta.id}_chunk_${String.format('%06d', chunk_index + 1)}", parent_id: meta.id ]
            [ chunk_meta, chunk_file ]
        }
    }

    // Two real searches: Pfam (--cut_ga, in modules.config) and Mgnifams (E-value). domtbl only.
    HMMSEARCH_PFAM(     ch_chunks.combine(ch_pfam_hmm    ).map { m, fa, hmm -> [ m, hmm, fa, false, false, true ] } )
    HMMSEARCH_MGNIFAMS( ch_chunks.combine(ch_mgnifams_hmm).map { m, fa, hmm -> [ m, hmm, fa, false, false, true ] } )

    // Regroup each search's chunked domtblouts back to the parent input.
    ch_pfam_domtbls = HMMSEARCH_PFAM.out.domain_summary
        .map { cm, d -> [ [ id: cm.parent_id ?: cm.id ], d ] }.groupTuple()
    ch_mgnifams_domtbls = HMMSEARCH_MGNIFAMS.out.domain_summary
        .map { cm, d -> [ [ id: cm.parent_id ?: cm.id ], d ] }.groupTuple()

    // Search-time filters are authoritative (Pfam --cut_ga, Mgnifams -E/--domE), so the stats
    // script applies no extra E-value threshold. Full input FASTA gives the total-sequence denominator.
    // Pfam-only, subtracting Mgnifams coverage -> also reports Pfam-exclusive residues.
    STATS_PFAM(
        ch_pfam_domtbls.join(ch_mgnifams_domtbls).map { m, p, g -> [ m, 'pfam', file(params.input_fasta), p, g ] }
    )
    // Mgnifams-only, subtracting Pfam coverage -> also reports Mgnifams-exclusive residues.
    STATS_MGNIFAM(
        ch_mgnifams_domtbls.join(ch_pfam_domtbls).map { m, g, p -> [ m, 'mgnifam', file(params.input_fasta), g, p ] }
    )
    // Union of Pfam + Mgnifams coverage.
    STATS_PFAM_MGNIFAM(
        ch_pfam_domtbls.join(ch_mgnifams_domtbls).map { m, p, g -> [ m, 'pfam_mgnifam', file(params.input_fasta), p + g, [] ] }
    )

    COMPARE_ANNOTATION_STATS(
        STATS_PFAM.out.stats.join( STATS_PFAM_MGNIFAM.out.stats )
    )

    // Provenance: DB checksums + collected tool versions + run metadata.
    ch_versions = channel.topic('versions')
        .map { row -> "${row[1]}: ${row[2]}" }
        .unique()
        .collectFile(name: 'versions.yml', newLine: true, sort: true)

    def workflow_json = groovy.json.JsonOutput.toJson([
        nextflow_version       : workflow.nextflow.version.toString(),
        pipeline_revision      : (workflow.revision ?: 'N/A').toString(),
        commit_id              : (workflow.commitId ?: 'N/A').toString(),
        command_line           : workflow.commandLine.toString(),
        start                  : workflow.start.toString(),
        container_engine       : (workflow.containerEngine ?: 'none').toString(),
        profile                : workflow.profile.toString(),
        pfam_threshold         : 'cut_ga',
        evalue_cutoff          : params.evalue_cutoff,
        effective_db_size      : params.effective_db_size,
        fasta_records_per_chunk: params.fasta_records_per_chunk,
    ])

    ch_provenance = channel.value( [ meta, [ file(params.input_fasta), file(params.pfam_hmm), file(params.mgnifams_hmm) ] ] )
        .combine( ch_versions )
    PROVENANCE_REPORT( ch_provenance, workflow_json )
}
