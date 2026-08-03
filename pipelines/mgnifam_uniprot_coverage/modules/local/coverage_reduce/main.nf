process COVERAGE_REDUCE {
    tag "${meta.id}:${meta.pass}:${meta.coords}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(summaries), path(families), path(lists)

    output:
    tuple val(meta), path("${prefix}.reduced.tsv")       , emit: reduced
    tuple val(meta), path("${prefix}.list_overlaps.tsv") , emit: overlaps, optional: true
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix        = task.ext.prefix ?: "${meta.id}_${meta.pass}_${meta.coords}"
    def args      = task.ext.args ?: ''
    def lists_arg = lists ? "--lists ${lists}" : ''
    // --expect-chunks is the guard that turns "some chunks are missing" from a
    // plausible-looking smaller number into a failed run.
    """
    mkdir -p reduce
    cp -L *.summary.tsv *.families.tsv.gz reduce/

    mgnifam_uniprot_coverage_stats_from_domtbl.py reduce \\
        --outdir reduce \\
        ${lists_arg} \\
        --total-sequences ${meta.total_sequences} \\
        --total-residues ${meta.total_residues} \\
        --expect-chunks ${meta.n_chunks} \\
        ${args}

    mv reduce/reduced.tsv ${prefix}.reduced.tsv
    [ -f reduce/list_overlaps.tsv ] && mv reduce/list_overlaps.tsv ${prefix}.list_overlaps.tsv || true
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.pass}_${meta.coords}"
    """
    touch ${prefix}.reduced.tsv ${prefix}.list_overlaps.tsv
    """
}
