process VALIDATE_COVERAGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(reduced), path(reference)

    output:
    tuple val(meta), path("validation.txt"), emit: ok
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def ref_arg = reference ? "--reference ${reference}" : ''
    """
    validate_coverage.py \\
        --reduced ${reduced} \\
        ${ref_arg} \\
        --output validation.txt \\
        ${args}
    """

    stub:
    """
    touch validation.txt
    """
}
