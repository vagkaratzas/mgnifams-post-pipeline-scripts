process VALIDATE_COVERAGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(reduced), path(references), val(reference_subsets)

    output:
    tuple val(meta), path("validation.txt"), emit: ok
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Pair each staged CSV with the subset it belongs to; both lists keep the
    // order they were built in, so a subset can never be checked against
    // another subset's reference.
    def ref_list = references instanceof List ? references : (references ? [references] : [])
    def sub_list = reference_subsets instanceof List ? reference_subsets : (reference_subsets ? [reference_subsets] : [])
    def ref_arg  = ref_list ? "--references " + [sub_list, ref_list].transpose().collect { s, f -> "${s}=${f}" }.join(' ') : ''
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
