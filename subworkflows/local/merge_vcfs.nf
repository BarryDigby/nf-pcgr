include { CPSR_VALIDATE_INPUT } from '../../modules/local/cpsr/validate_input'
include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX } from '../../modules/nf-core/tabix/tabix/main'
include { ISEC_SOMATIC_VCFS as INTERSECT_SOMATIC_VARIANTS } from '../../modules/local/Merge/isec_vcfs'
include { PCGR_VCF as PCGR_READY_VCF } from '../../modules/local/Merge/pcgr_vcf'

workflow MERGE_VCFS {
    take:
    somatic_files
    germline_files
    fasta
    pcgr_header
    pcgr_dir

    main:
    // Use CPSR scripts to simplify germline VCF files.
    // Strategy is to consolidate variants from all germline tools.
    CPSR_VALIDATE_INPUT( germline_files, pcgr_dir.collect() )
    TABIX_BGZIPTABIX( CPSR_VALIDATE_INPUT.out.validated_vcf )




    // create master TSV file with variant <-> tool mapping
    // Extract VCF and TBI from channel, choose suitable meta info for merging samples (pop meta.tool, meta.status)
    // < [[ meta.patient, meta.sample], all tool vcfs, all tool tbi ]
    per_sample_somatic_vcfs = somatic_files.map{ it -> return it[0..2] }.map{ meta, vcf, tbi -> var = [:]; var.patient = meta.patient; var.sample = meta.sample; return [ var, vcf, tbi ] }.groupTuple()
    INTERSECT_SOMATIC_VARIANTS( per_sample_somatic_vcfs )

    // merge mapping key back with sample VCFs, produce PCGR ready VCFs. -- massaging to remove lists introduced by grouptuple
    sample_vcfs_keys = INTERSECT_SOMATIC_VARIANTS.out.variant_tool_map.join(per_sample_somatic_vcfs)

    PCGR_READY_VCF( sample_vcfs_keys, pcgr_header )

    // Add the CNVkit file back to the PCGR ready VCFs
    emit:
    pcgr_ready_vcf = params.cna_analysis ? PCGR_READY_VCF.out.vcf.join( somatic_files.map{ it -> return it[3] }.flatten().take(1).map{ it -> meta = [:]; meta.id = it.simpleName; return [ meta, it ] } ) : PCGR_READY_VCF.out.vcf.map{ meta, vcf, tbi -> return [ meta, vcf, tbi, [] ] }
    cpsr_ready_vcf = CPSR_VALIDATE_INPUT.out.validated_vcf
}
