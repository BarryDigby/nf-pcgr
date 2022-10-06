include { FORMAT_VCF    } from '../../modules/local/PCGR/Format/pcgr_reformat'
include { FORMAT_CNA    } from '../../modules/local/PCGR/Format/pcgr_reformat'

workflow FORMAT_FILES {
    take:
    fasta  // channel fasta [<>.{fa,fasta}]
    files  // channel [ val(meta), [ vcf.gz] , [ vcf.gz.tbi ] ] OR [ val(meta), [vcf.gz], [vcf.gz.tbi], [CNA] ]

    main:
    FORMAT_VCF( fasta, files )
    FORMAT_CNA( FORMAT_VCF.out.files )
    FORMAT_CNA.out.files.view()

    emit:
    files = params.cna_analysis ? FORMAT_CNA.out.files : FORMAT_VCF.out.files

}
