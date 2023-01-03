process PICARD_SPLITVCF {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::picard=2.27.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard%3A2.27.4--hdfd78af_0':
        'quay.io/biocontainers/picard:2.27.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*indels*"), emit: indels
    tuple val(meta), path("*snps*"),   emit: snps
    path "versions.yml"           ,    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard SortVcf] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    picard \\
        SplitVcfs \\
        -Xmx${avail_mem}g \\
        --INPUT $vcf \\
        $args \\
        --INDEL_OUTPUT ${prefix}_indels.vcf.gz \\
        --SNP_OUTPUT ${prefix}_snps.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard SplitVcfs --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_indels.vcf.gz
    touch ${prefix}_snps.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard SplitVcfs --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
