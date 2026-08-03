process COVERAGE_REPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(reduced), path(library_sizes), path(validation)
    path samplesheet

    output:
    tuple val(meta), path("${prefix}.md")  , emit: markdown
    tuple val(meta), path("${prefix}.html"), emit: html
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix   = task.ext.prefix ?: 'coverage_report'
    def args = task.ext.args ?: ''
    """
    coverage_report.py \\
        --reduced ${reduced} \\
        --samplesheet ${samplesheet} \\
        --library-sizes ${library_sizes} \\
        --validation ${validation} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: 'coverage_report'
    """
    touch ${prefix}.md ${prefix}.html
    """
}
