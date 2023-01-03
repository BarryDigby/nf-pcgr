process CPSR_VALIDATE_INPUT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "pcgr::pcgr=1.2.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/sigven/pcgr:1.2.0':
        'docker.io/sigven/pcgr:1.2.0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path(pcgr_dir), stageAs: "PCGR/data/${params.genome.toLowerCase()}"

    output:
    tuple val(meta), path("${prefix}/*cpsr_ready.vcf"), emit: validated_vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def database = "./PCGR"
    prefix   = task.ext.prefix ?: "${meta.id}"
    def args     = task.ext.args ?: ''
    """
    mkdir $prefix

    cpsr_validate_input.py \\
        $database \\
        $vcf \\
        "None" \\
        "None" \\
        0 \\
        ${params.genome.toLowerCase()} \\
        $prefix \\
        $params.panel_id \\
        0 \\
        --output_dir $prefix \\
        --debug
    """
}
