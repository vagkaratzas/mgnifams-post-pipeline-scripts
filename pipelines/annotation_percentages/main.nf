//
// MODULE: Local modules, one per process
//
include { EXTRACT_FASTA                                              } from './modules/local/extract_fasta/main.nf'
include { APPEND_MGNIFAMS                                            } from './modules/local/append_mgnifams/main.nf'
include { CALCULATE_ANNOTATION_STATS as CALCULATE_PFAM_STATS         } from './modules/local/calculate_annotation_stats/main.nf'
include { CALCULATE_ANNOTATION_STATS as CALCULATE_PFAM_MGNIFAM_STATS } from './modules/local/calculate_annotation_stats/main.nf'
include { COMPARE_ANNOTATION_STATS                                   } from './modules/local/compare_annotation_stats/main.nf'

//
// MODULE: Installed directly from nf-core/modules
//
include { HMMER_HMMSEARCH                                            } from './modules/nf-core/hmmer/hmmsearch/main.nf'

workflow {
    def meta = [ id: 'mgnify_proteins' ]

    // Local CSV / CSV.GZ input; the python scripts read either transparently
    ch_input = channel.value( file(params.input_csv, checkIfExists: true) ).map { csv -> [ meta, csv ] }
    ch_hmm   = channel.value( file(params.hmm_lib,   checkIfExists: true) )

    EXTRACT_FASTA( ch_input )

    // Run hmmsearch independently for each FASTA chunk. The parent_id keeps chunks grouped
    // back to the original input CSV before annotations are appended.
    ch_fasta_chunks = EXTRACT_FASTA.out.fasta
        .flatMap { parent_meta, fasta ->
            fasta.splitFasta(by: params.fasta_records_per_chunk.toInteger(), file: true).withIndex().collect { chunk_file, chunk_index ->
                def chunk_id = "${parent_meta.id}_chunk_${String.format('%06d', chunk_index + 1)}"
                def chunk_meta = parent_meta + [ id: chunk_id, parent_id: parent_meta.id ]
                [ chunk_meta, chunk_file ]
            }
        }

    // hmmsearch reads the gzipped HMM library directly; only the domain table is requested
    ch_hmmsearch = ch_fasta_chunks
        .combine( ch_hmm )
        .map { hmeta, fasta, hmm -> [ hmeta, hmm, fasta, false, false, true ] }
    HMMER_HMMSEARCH( ch_hmmsearch )

    // Recombine all chunked domtbl outputs before replacing metadata["m"] in the input CSV
    ch_domain_summaries = HMMER_HMMSEARCH.out.domain_summary
        .map { chunk_meta, domtbl -> [ [ id: chunk_meta.parent_id ?: chunk_meta.id ], domtbl ] }
        .groupTuple()
    APPEND_MGNIFAMS( ch_input.join( ch_domain_summaries ) )

    // Pfam-only (before) vs Pfam+MGnifam (after) annotation statistics
    CALCULATE_PFAM_STATS( ch_input, 'p', 'pfam' )
    CALCULATE_PFAM_MGNIFAM_STATS( APPEND_MGNIFAMS.out.csv, 'p,m', 'pfam_mgnifam' )

    COMPARE_ANNOTATION_STATS(
        CALCULATE_PFAM_STATS.out.stats.join( CALCULATE_PFAM_MGNIFAM_STATS.out.stats )
    )
}
