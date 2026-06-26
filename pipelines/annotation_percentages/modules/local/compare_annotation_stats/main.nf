process COMPARE_ANNOTATION_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(before_stats), path(after_stats)

    output:
    tuple val(meta), path("*increase.csv"), emit: comparison

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: 'annotation_percentage_increase'
    """
    compare_annotation_stats.py \\
        --before ${before_stats} \\
        --after ${after_stats} \\
        --output ${prefix}.csv
    """

    stub:
    def prefix = task.ext.prefix ?: 'annotation_percentage_increase'
    """
    touch ${prefix}.csv
    """
}
