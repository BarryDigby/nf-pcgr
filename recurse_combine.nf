#!/usr/bin/env nextflow

files = collect_sarek_files(file(params.input))
files.view()
files.collectFile( name: 'constructed_samplesheet.csv', newLine:false, storeDir: "${params.outdir}", keepHeader: true ){ ids, vcf, cna -> "sample,vcf,cna" + "\n" + "$ids,$vcf,$cna" + "\n"}.set{ constructed_samplesheet }
constructed_samplesheet.view()
def collect_sarek_files(input){
    vcf_files = []
    cna_files = []
    input.eachFileRecurse{ it ->
        vcf = it.name.contains('_vs_') && ( it.name.endsWith('.vcf') || it.name.endsWith('.vcf.gz') ) && !it.name.endsWith('.tbi') ? file(it) : []
        cna = it.name.endsWith('.cns') ? file(it) : params.cna_analysis ? [] : 'NA' // catch for sarek run without CNVkit results
        ids = it.simpleName.tokenize('_')[0]
        vcf_files << [ ids, vcf ]
        cna_files << [ ids, cna ]
        }
    Channel.fromList( vcf_files ).filter{ ids, vcf -> vcf.toString().contains('.vcf') }.set{ collect_vcf }
    if(params.cna_analysis){
        Channel.fromList( cna_files ).filter{ ids, cna -> cna.toString().contains('.cns') }.set{ collect_cna }
    }else{
        Channel.fromList( cna_files ).set{ collect_cna }
    }

    collect_cna.view()
    sarek_files = collect_vcf.combine(collect_cna, by:0)
    return sarek_files
}
