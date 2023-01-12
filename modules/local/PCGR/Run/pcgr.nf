process PCGR {
    tag "${meta.patient}:${meta.sample}"
    label 'process_high'

    conda (params.enable_conda ? "pcgr::pcgr=1.2.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/sigven/pcgr:1.2.0':
        'docker.io/sigven/pcgr:1.2.0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(cna)
    path(pcgr_dir), stageAs: "PCGR/data/${params.genome.toLowerCase()}"
    path pon

    output:
    tuple val(meta), path("${prefix}"), emit: pcgr_reports
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome   = task.ext.genome ?: ''
    def database = './PCGR'
    def args     = task.ext.args ?: ''
    prefix       = task.ext.prefix ?: "${meta.id}"
    def cna      = params.cna_analysis ? "--input_cna $cna" : ''
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
        pcgr: \$(echo \$( pcgr --version | sed 's/pcgr //g' ))
    END_VERSIONS
    """
}
