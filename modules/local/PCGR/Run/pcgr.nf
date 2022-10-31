process PCGR {
    tag "${meta.id}:${meta.tool}"
    label 'process_medium'

    conda (params.enable_conda ? "pcgr::pcgr=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/sigven/pcgr:1.1.0':
        'docker.io/sigven/pcgr:1.1.0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(cna)

    output:
    path "versions.yml"           , emit: versions
    tuple val(meta), path("${prefix}"), emit: pcgr_reports

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome   = task.ext.genome ?: ''
    def database = task.ext.database ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.tool}"
    def cna      = params.cna_analysis ? "--input_cna $cna" : ''
    def args     = task.ext.args ?: ''
    """
    mkdir -p $prefix

    pcgr \\
        --input_vcf $vcf \\
        --pcgr_dir $database \\
        --output_dir $prefix \\
        --genome_assembly $genome \\
        --sample_id $prefix \\
        $cna \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pcgr: \$(echo \$( pcgr --version | sed 's/pcgr//g' ))
    END_VERSIONS
    """
}
