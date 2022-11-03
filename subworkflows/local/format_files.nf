include { FORMAT_VCF      } from '../../modules/local/PCGR/Format/pcgr_reformat'
include { FORMAT_CNA      } from '../../modules/local/PCGR/Format/pcgr_reformat'
include { FOO             } from '../../modules/local/foo'

workflow FORMAT_FILES {
    take:
    fasta  // channel fasta [<>.{fa,fasta}]
    files  // channel [ val(meta), [ vcf.gz] , [ vcf.gz.tbi ] ] OR [ val(meta), [vcf.gz], [vcf.gz.tbi], [CNA] ]

    main:
    FORMAT_VCF( fasta, files )
    FORMAT_CNA( FORMAT_VCF.out.files )
    FORMAT_CNA.out.files.view()

    // grouptuple and dump VCF files into pcgr merge/simplify process - even if only one tool used. at min, we will strip genotype data.
    files = params.cna_analysis ? FORMAT_CNA.out.files : FORMAT_VCF.out.files
    sample_files = files.map{ it -> return it[1..-1] }.flatten().unique().map{ it -> meta = it.simpleName; return [ meta, it ] }.groupTuple().view()

    FOO( sample_files )

    emit:
    files

}
