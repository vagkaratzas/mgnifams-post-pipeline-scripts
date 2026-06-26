process CALCULATE_ANNOTATION_STATS {
    tag "${meta.id}:${label}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(input_csv)
    val annotation_keys
    val label

    output:
    tuple val(meta), path("*_annotation_stats.csv"), emit: stats

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    calculate_annotation_stats.py \\
        --input ${input_csv} \\
        --annotation-keys ${annotation_keys} \\
        --label ${label} \\
        --output ${label}_annotation_stats.csv
    """

    stub:
    """
    touch ${label}_annotation_stats.csv
    """
}
