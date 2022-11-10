
include { ISEC_VCFS } from '../../modules/local/Merge/isec_vcfs'
include { PCGR_VCF  } from '../../modules/local/Merge/pcgr_vcf'

workflow MERGE_VCFS {
    take:
    files
    fasta

    main:
    println('Viewing files channel')
    files.view()
    // create master TSV file with variant <-> tool mapping
    // Extract VCF and TBI from channel, keeping the meta information.
    sample_vcfs = files.map{ it -> return it[1..2] }.flatten().map{ it -> meta = [:]; meta.id = it.simpleName; return [ meta, it ] }.groupTuple()
    println('Viewing sample_vcfs..')
    sample_vcfs.view()
    ISEC_VCFS( sample_vcfs )

    // merge back with sample VCFs, produce PCGR ready VCFs.
    sample_vcfs_keys = ISEC_VCFS.out.variant_tool_map.join(sample_vcfs)

    PCGR_VCF( sample_vcfs_keys, "${projectDir}/bin/pcgr_header.txt")

    // rework to make meta a tuple
    //println('viewing foo channel..')
    //foo = PCGR_VCF.out.vcf.map{ meta, vcf, tbi -> var = [:]; var.id = meta; return [ var, vcf, tbi ]}.view()
    // TODO: view the files channel and figure out how merge cna file back by meta
    emit:
    // flattencollect to satisfy input cardinality. nested tuple [ meta[vcf, tbi], cna] otherwise
    pcgr_ready_vcf = params.cna_analysis ? PCGR_VCF.out.vcf.join( files.map{ it -> return [ it[0], it[3] ]} ) : PCGR_VCF.out.vcf.map{ meta, vcf, tbi -> return [ meta, vcf, tbi, [] ] }
    //sample_vcfs
}
