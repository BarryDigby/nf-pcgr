include { BCFTOOLS_NORM as NORMALISE_VARIANTS } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_FILTER as FILTER_VARIANTS  } from '../../modules/nf-core/bcftools/filter/main'
include { TABIX_TABIX as TABIX_VARIANTS       } from '../../modules/nf-core/tabix/tabix/main'
include { REFORMAT_VCF } from '../../modules/local/PCGR/Format/pcgr_reformat'
include { REFORMAT_CNA } from '../../modules/local/PCGR/Format/pcgr_reformat'

workflow FORMAT_FILES {
    take:
    fasta  // channel fasta [<>.{fa,fasta}]
    files  // channel [ val(meta), [ vcf.gz] , [ vcf.gz.tbi ] ] OR [ val(meta), [vcf.gz], [vcf.gz.tbi], [CNA] ]

    main:
    NORMALISE_VARIANTS( files.map{ it -> return [ it[0], it[1], it[2] ] }, fasta )
    FILTER_VARIANTS( NORMALISE_VARIANTS.out.vcf )
    TABIX_VARIANTS( FILTER_VARIANTS.out.vcf )

    REFORMAT_VCF( FILTER_VARIANTS.out.vcf.join( TABIX_VARIANTS.out.tbi ) )
    REFORMAT_CNA( files.map{ it -> return [ it[0], it[3] ]} )

    REFORMAT_CNA.out.cna.view()

    files = params.cna_analysis ? REFORMAT_VCF.out.vcf.join( REFORMAT_CNA.out.cna ) : REFORMAT_VCF.out.vcf.map{ meta, vcf, tbi -> return [ meta, vcf, tbi, [] ] }

    emit:
    files
}
