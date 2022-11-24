process CPSR {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "pcgr::pcgr=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/sigven/pcgr:1.1.0':
        'docker.io/sigven/pcgr:1.1.0' }"

    input:
    // tuple [ meta, [vcf] , [vcf.tbi] ]
    tuple val(meta), path(vcf), path(tbi), path(cna)
    path(pcgr_dir), stageAs: "PCGR/data/${params.genome}"

    output:
    path "versions.yml"           , emit: versions
    tuple val(meta), path("${meta.id}"), emit: cpsr_reports

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome   = task.ext.genome ?: ''
    def database = "./PCGR"
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def args     = task.ext.args ?: ''
    """
    mkdir -p $prefix

    cpsr \\
        --input_vcf $vcf \\
        --pcgr_dir $database \\
        --output_dir $prefix \\
        --genome_assembly $genome \\
        --sample_id $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pcgr: \$(echo \$( pcgr --version | sed 's/pcgr//g' ))
    END_VERSIONS
    """
}
