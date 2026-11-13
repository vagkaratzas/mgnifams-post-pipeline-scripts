process PROVENANCE_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(databases), path(versions)
    val workflow_json

    output:
    tuple val(meta), path("provenance.txt"), emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    def db_list  = databases instanceof List ? databases : [databases]
    def db_arg   = "--db ${db_list.join(' ')}"
    def ver_arg  = versions ? "--versions ${versions}" : ''
    """
    cat > workflow.json <<'JSON'
${workflow_json}
JSON

    provenance_report.py \\
        ${db_arg} \\
        ${ver_arg} \\
        --workflow-json workflow.json \\
        --output provenance.txt
    """

    stub:
    """
    touch provenance.txt
    """
}
