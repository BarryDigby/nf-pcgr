process REFORMAT_VCF {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/barryd237/pysam-xcmds:latest' :
        'docker.io/barryd237/pysam-xcmds:latest' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    reformat_vcf.py \\
        reformat_vcf \\
        -vcf_file $vcf \\
        -out ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}

process REFORMAT_CNA {
    tag "${meta.patient}:${meta.sample}:${cna}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/barryd237/pysam-xcmds:latest' :
        'docker.io/barryd237/pysam-xcmds:latest' }"

    input:
    tuple val(meta), path(cna)

    output:
    tuple val(meta), path("${meta.id}.*.tsv"), emit: cna
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    reformat_cna.py \
        reformat_cna \
        -cna_file $cna \
        -sample $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
