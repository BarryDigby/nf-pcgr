include { CPSR_VALIDATE_INPUT } from '../../modules/local/cpsr/validate_input'
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
    CPSR_VALIDATE_INPUT( germline_files, pcgr_dir.collect() )

    // create master TSV file with variant <-> tool mapping
    // Extract VCF and TBI from channel, choose suitable meta info for merging samples (pop meta.tool) -- message me if you know how to iterate LinkedHashmaps like meta[1..3]
    // < [[ meta.patient, meta.status, meta.sample], all tool vcfs, all tool tbi ]
    per_sample_somatic_vcfs = somatic_files.map{ it -> return it[0..2] }.map{ meta, vcf, tbi -> var = [:]; var.patient = meta.patient; var.status = meta.status; var.sample = meta.sample; return [ var, vcf, tbi ] }.groupTuple()
    INTERSECT_SOMATIC_VARIANTS( per_sample_somatic_vcfs )

    // merge mapping key back with sample VCFs, produce PCGR ready VCFs.
    sample_vcfs_keys = INTERSECT_SOMATIC_VARIANTS.out.variant_tool_map.join(per_sample_somatic_vcfs)

    PCGR_READY_VCF( sample_vcfs_keys, pcgr_header )

    // Add the CNVkit file back to the PCGR ready VCFs
    emit:
    pcgr_ready_vcf = params.cna_analysis ? PCGR_READY_VCF.out.vcf.join( somatic_files.map{ it -> return it[3] }.flatten().take(1).map{ it -> meta = [:]; meta.id = it.simpleName; return [ meta, it ] } ) : PCGR_READY_VCF.out.vcf.map{ meta, vcf, tbi -> return [ meta, vcf, tbi, [] ] }

}
