//
// MODULE: Installed directly from nf-core/modules
//
include { FUNGTION_DOWNLOADMODELS } from './modules/nf-core/fungtion/downloadmodels/main.nf'
include { FUNGTION_FUNGTION       } from './modules/nf-core/fungtion/fungtion/main.nf'
include { CAALM_DOWNLOADMODELS    } from './modules/nf-core/caalm/downloadmodels/main.nf'
include { CAALM_CAALM             } from './modules/nf-core/caalm/caalm/main.nf'

workflow {
    def fasta = file(params.input_fasta, checkIfExists: true)
    def meta  = [ id: params.prefix ?: fasta.simpleName ]

    ch_fasta = channel.value( [ meta, fasta ] )

    // Pre-downloaded model paths win over re-running the download modules
    ch_fungtion_models = params.fungtion_models
        ? channel.value( file(params.fungtion_models, checkIfExists: true) )
        : FUNGTION_DOWNLOADMODELS().models

    // CAALM expects the three model levels as a tuple; snapshot_download lays them
    // out as <models>/level0, <models>/level1, <models>/level2
    ch_caalm_models = params.caalm_models
        ? channel.value( [ 'level0', 'level1', 'level2' ].collect { lvl ->
              file("${params.caalm_models}/${lvl}", checkIfExists: true)
          } )
        : CAALM_DOWNLOADMODELS().models

    FUNGTION_FUNGTION( ch_fasta, ch_fungtion_models )
    CAALM_CAALM( ch_fasta, ch_caalm_models )
}
