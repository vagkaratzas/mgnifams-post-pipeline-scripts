process BUILD_CATEGORY_LISTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hb6cb901_4' :
        'quay.io/biocontainers/hmmer:3.4--hb6cb901_4' }"

    input:
    tuple val(meta), path(lists_dir, stageAs: 'base_lists'), path(hmm_lib)

    output:
    tuple val(meta), path("lists")            , emit: lists
    tuple val(meta), path("library_sizes.tsv"), emit: sizes
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | sed '2!d;s/^# HMMER *//;s/ .*//'"), emit: versions_hmmer, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The derived categories are pure set algebra, so sort/comm do it directly --
    // no script to maintain, and the ID form is the one the engine already accepts.
    def decompress = hmm_lib.name.endsWith('.gz') ? "gunzip -c ${hmm_lib} > lib.hmm" : "cp ${hmm_lib} lib.hmm"
    """
    mkdir -p lists
    cp base_lists/*.txt lists/

    ${decompress}

    # library_ids: every model in the HMM library. Names are bare integers there
    # and padded MGYF in the lists, so pad here and both sides compare directly.
    hmmstat lib.hmm | awk '!/^#/ && NF { printf "MGYF%010d\\n", \$2 }' | sort -u > lists/library_all.txt

    # transmembrane of either flavour
    sort -u lists/membrane_a.txt lists/membrane_b.txt > lists/tm_any.txt
    # TM or disordered as a UNION: summing the per-category residue rows would
    # double count spans where two such families overlap on the same protein
    sort -u lists/tm_any.txt lists/disorder.txt > lists/tm_or_disorder.txt
    # the complement, which is what "excluding these families" has to mean
    comm -23 lists/library_all.txt lists/tm_or_disorder.txt > lists/not_tm_disorder.txt

    # library_all is the denominator, not a category to report on
    mv lists/library_all.txt library_all.txt

    printf 'category\\tn_families\\n' > library_sizes.tsv
    printf '__library__\\t%s\\n' "\$(wc -l < library_all.txt)" >> library_sizes.tsv
    for f in lists/*.txt; do
        printf '%s\\t%s\\n' "\$(basename \$f .txt)" "\$(wc -l < \$f)" >> library_sizes.tsv
    done
    """

    stub:
    """
    mkdir -p lists
    touch lists/novel.txt
    printf 'category\\tn_families\\n__library__\\t0\\n' > library_sizes.tsv
    """
}
