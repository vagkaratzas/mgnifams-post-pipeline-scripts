process APPEND_MGNIFAMS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(input_csv), path(domtbl)

    output:
    tuple val(meta), path("*_mgnifams.csv.gz"), emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: 'proteins'
    def domtbl_args = domtbl instanceof List ? domtbl.sort { domtbl_file -> domtbl_file.toString() }.join(' ') : domtbl
    """
    append_mgnifams_annot.py \\
        ${domtbl_args} \\
        ${input_csv} \\
        ${prefix}_mgnifams.csv.gz \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: 'proteins'
    """
    echo "" | gzip > ${prefix}_mgnifams.csv.gz
    """
}
