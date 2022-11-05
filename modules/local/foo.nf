process VCF_MERGE {
    tag "$meta"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}.merged.vcf.gz*"), emit: merge

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.tool}"
    """

    """
}
