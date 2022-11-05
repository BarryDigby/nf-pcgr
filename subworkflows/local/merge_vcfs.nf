
include { ISEC_VCFS } from '../../modules/local/Merge/isec_vcfs'

workflow MERGE_VCFS {
    take:
    files
    fasta

    main:
    // create master TSV file with variant <-> tool mapping
    sample_vcfs = files.map{ it -> return it[1..2] }.flatten().map{ it -> meta = it.simpleName; return [ meta, it ] }.groupTuple()
    ISEC_VCFS( sample_vcfs )

    // merge back with sample VCFs, produce PCGR ready VCFs.
    ISEC_VCFS.out.variant_tool_map.join(sample_vcfs).view()

    emit:
    foo = ISEC_VCFS.out.sample_keys

}
