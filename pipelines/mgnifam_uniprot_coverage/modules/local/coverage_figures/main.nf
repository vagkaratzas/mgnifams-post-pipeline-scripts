process COVERAGE_FIGURES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // plotnine is conda-forge, so there is no biocontainer for it and nothing on
    // depot.galaxyproject.org; these are Seqera community builds of environment.yml.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_plotnine_pandas:2b0b556e61a10ef2' :
        'community.wave.seqera.io/library/python_plotnine_pandas:dae7bdcea75615d2' }"

    input:
    tuple val(meta), path(reduced), path(library_sizes)
    path samplesheet

    output:
    tuple val(meta), path("*.pdf"), emit: figures
    tuple val("${task.process}"), val('plotnine'), eval("python3 -c \"import importlib.metadata; print(importlib.metadata.version('plotnine'))\""), emit: versions_plotnine, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    coverage_figures.py \\
        --reduced ${reduced} \\
        --samplesheet ${samplesheet} \\
        --library-sizes ${library_sizes} \\
        --outdir . \\
        ${args}
    """

    stub:
    """
    touch fig1_coverage_by_database.pdf
    """
}
