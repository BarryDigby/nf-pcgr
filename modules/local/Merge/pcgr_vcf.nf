process PCGR_VCF {
    tag "$meta"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/barryd237/pysam-xcmds:latest' :
        'docker.io/barryd237/pysam-xcmds:latest' }"

    input:
    tuple val(meta), path(keys), path(files)
    path(pcgr_header)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi") emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta}"
    """
    python3.6 "${projectDir}/bin/pcgr_vcf.py" \
        -sample ${prefix}
    """
}
