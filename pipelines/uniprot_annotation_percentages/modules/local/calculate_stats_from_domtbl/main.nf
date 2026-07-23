process CALCULATE_STATS_FROM_DOMTBL {
    tag "${meta.id}:${label}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), val(label), path(fasta), path(domtbls), path(subtract)

    output:
    tuple val(meta), path("*_annotation_stats.csv"), emit: stats

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def domtbl_list   = domtbls instanceof List ? domtbls : [domtbls]
    def domtbl_arg    = "--domtbl ${domtbl_list.join(' ')}"
    def subtract_list = subtract instanceof List ? subtract : (subtract ? [subtract] : [])
    def subtract_arg  = subtract_list ? "--subtract ${subtract_list.join(' ')}" : ''
    """
    calculate_stats_from_domtbl.py \\
        --fasta ${fasta} \\
        --label ${label} \\
        ${domtbl_arg} \\
        ${subtract_arg} \\
        ${args} \\
        --output ${label}_annotation_stats.csv
    """

    stub:
    """
    touch ${label}_annotation_stats.csv
    """
}
