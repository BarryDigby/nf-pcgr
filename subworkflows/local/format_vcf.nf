include { ADD_HEADER_INFO } from '../../modules/local/PCGR/Format/pcgr_reformat'

workflow FORMAT_VCF {
    take:
    fasta  // channel fasta [<>.{fa,fasta}]
    files  // channel [ val(meta), [ vcf.gz] , [ vcf.gz.tbi ] ] OR [ val(meta), [vcf.gz], [vcf.gz.tbi], [CNA] ]

    main:
    ADD_HEADER_INFO( fasta, files )

    emit:
    files = ADD_HEADER_INFO.out.files // channel: same as above, but VCF now has correct header and info fields for PCGR. Fasta discarded.
}
