include { BCFTOOLS_NORM   } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_FILTER } from '../../modules/nf-core/bcftools/filter/main'
include { TABIX_TABIX as TBX } from '../../modules/nf-core/tabix/tabix/main'
include { FORMAT_VCF      } from '../../modules/local/PCGR/Format/pcgr_reformat'
include { FORMAT_CNA      } from '../../modules/local/PCGR/Format/pcgr_reformat'

workflow FORMAT_FILES {
    take:
    fasta  // channel fasta [<>.{fa,fasta}]
    files  // channel [ val(meta), [ vcf.gz] , [ vcf.gz.tbi ] ] OR [ val(meta), [vcf.gz], [vcf.gz.tbi], [CNA] ]

    main:
    BCFTOOLS_NORM( files.map{ it -> return [ it[0], it[1], it[2] ] }, fasta )
    BCFTOOLS_FILTER( BCFTOOLS_NORM.out.vcf )
    TBX( BCFTOOLS_FILTER.out.vcf )

    FORMAT_VCF( BCFTOOLS_FILTER.out.vcf.join( TBX.out.tbi ) )
    FORMAT_CNA( files.map{ it -> return [ it[0], it[3] ]} )

    FORMAT_CNA.out.cna.view()

    files = params.cna_analysis ? FORMAT_VCF.out.vcf.join( FORMAT_CNA.out.cna ) : FORMAT_VCF.out.vcf.map{ meta, vcf, tbi -> return [ meta, vcf, tbi, [] ] }

    emit:
    files
}
