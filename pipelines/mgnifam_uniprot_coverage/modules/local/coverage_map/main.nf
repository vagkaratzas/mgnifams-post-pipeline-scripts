process COVERAGE_MAP {
    tag "${meta.id}:${meta.pass}:${meta.coords}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(domtbl), path(lists), path(mask)

    output:
    tuple val(meta), path("out/*.summary.tsv")      , emit: summary
    tuple val(meta), path("out/*.families.tsv.gz")  , emit: families
    tuple val(meta), path("out/*.pertarget.tsv.gz") , emit: pertarget
    tuple val(meta), path("out/logs/*.log")         , emit: log
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    // meta.coords is authoritative: 'ali' reproduces the published figures,
    // 'env' is the sensitivity check. Never silently defaulted.
    def coords    = meta.coords == 'env' ? '--env' : '--no-env'
    def lists_arg = lists ? "--lists ${lists}" : ''
    def mask_arg  = mask  ? "--mask ${mask}"   : ''
    """
    mgnifam_uniprot_coverage_stats_from_domtbl.py map \\
        --domtbl ${domtbl} \\
        --outdir out \\
        ${coords} \\
        ${lists_arg} \\
        ${mask_arg} \\
        ${args}
    """

    stub:
    def stem = domtbl.name.replaceAll(/\\.domtbl(\\.gz)?$/, '')
    """
    mkdir -p out/logs
    touch out/${stem}.summary.tsv out/logs/${stem}.log
    echo "" | gzip > out/${stem}.families.tsv.gz
    echo "" | gzip > out/${stem}.pertarget.tsv.gz
    """
}
