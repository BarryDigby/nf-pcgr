
include { ISEC_VCFS } from '../../modules/local/merge_vcf/isec_vcfs'

workflow MERGE_VCFS {
    take:
    files
    fasta

    main:
    // create master TSV file with variant <-> tool mapping
    sample_vcfs = files.map{ it -> return it[1..2] }.flatten().map{ it -> meta = it.simpleName; return [ meta, it ] }.groupTuple()
    ISEC_VCFS( sample_vcfs )

    emit:
    foo = ISEC_VCFS.out.sample_keys

}
