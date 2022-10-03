#!/usr/bin/env nextflow

files = collect_files(file(params.input))
files.view()

def collect_files(input){
    vcf_files = []
    cna_files = []
    input.eachFileRecurse{ it ->
        vcf = it.name.contains('_vs_') && ( it.name.endsWith('.vcf') || it.name.endsWith('.vcf.gz') ) && !it.name.endsWith('.tbi') ? file(it) : []
        cna = it.name.endsWith('.cns') ? file(it) : []
        ids = cna == [] ? it.simpleName.tokenize('_')[0] : cna.simpleName
        vcf_files << [ ids, vcf ]
        cna_files << [ ids, cna ]
        }
    Channel.fromList( vcf_files ).filter{ ids, vcf -> vcf.toString().contains('.vcf') }.set{ collect_vcf }
    Channel.fromList( cna_files ).filter{ ids, cna -> cna.toString().contains('.cns') }.set{ collect_cna }
    sarek_files = collect_vcf.combine(collect_cna, by:0)
    return sarek_files
}
