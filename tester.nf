#!/usr/bin/env nextflow

// mimick sarek file structure:
// touch 123_vs_456.{freebayes,strelka,mutect2}.vcf.gz
// touch ABC_vs_DEF.{freebayes,strelka,mutect2}.vcf.gz
// nextflow run tester.nf --input "path to touch files"

// input can be path or file
if (params.input) { ch_input = file( params.input ) } else { exit 1, "Input does not exist" }

if(ch_input.toString().endsWith('.csv')){
    println "samplesheet.csv as input.."
}else{
    vcf = collect_vcf(ch_input)
    cna = collect_cna(ch_input)
    files = vcf.combine(cna, by:0)
    files.view()
}
// scan for vcf files and cna files using different functions baz. asking too much to use 'return twice' when combine is not behaving as expected.
// works. commit and clean in the morning.
def collect_vcf(input){

    // stage vars outside conditional scope

    // not sure where to stage these, inside or outside of eachFileRecurse...

    //collect_vcf = Channel.empty()
    vcf_files = []

    input.eachFileRecurse{ it ->

        if( ( it.name.contains('_vs_') ) && !it.name.endsWith('.tbi') ){

            // Sarek output naming convention allows use of tokenize()
            sample = it.simpleName.tokenize('_')[0]
            vcf    = it.name.endsWith('.vcf') || it.name.endsWith('.vcf.gz') ? file(it) : []
            //cna    = it.name.endsWith('.cns') ? file(it) : []

            vcf_files << [ sample, vcf ]
            //cna_files << [ sample, cna ]

            // remove empty slots '[]'
            //Channel.of( [ sample, vcf ] ).filter{ sample, vcf -> vcf.toString().contains("_") }.set{ collect_vcf }
            //Channel.fromList( [ sample, vcf ] ).collate(2).filter{ sample, vcf ->  vcf.toString().contains("_") }.set{ collect_vcf }
            //Channel.fromList( [ sample, cna ] ).collate(2).filter{ sample, cna -> cna.toString().contains(".cns") }.set{ collect_cna }

        }

    }

    Channel.fromList(vcf_files).flatten().collate(2).filter{ sample, vcf -> vcf.toString().contains("_") }.map{ sample, vcf -> [ sample, vcf ] }.set{ collect_vcf_files }
    //Channel.fromList(cna_files).flatten().collate(2).filter{ sample, cna -> cna.toString().contains(".cns") }.map{ sample, cna -> [ sample, cna ] }.set{ collect_cna }
    return collect_vcf_files
}


def collect_cna(input){

    // stage vars outside conditional scope

    // not sure where to stage these, inside or outside of eachFileRecurse...

    //collect_vcf = Channel.empty()
    cna_files = []

    input.eachFileRecurse{ it ->

        if( it.name.contains('.cns') ){

            // Sarek output naming convention allows use of tokenize()
            sample = it.simpleName.tokenize('_')[0]
            //vcf    = it.name.endsWith('.vcf') || it.name.endsWith('.vcf.gz') ? file(it) : []
            cna    = it.name.endsWith('.cns') ? file(it) : []

            //vcf_files << [ sample, vcf ]
            cna_files << [ sample, cna ]

            // remove empty slots '[]'
            //Channel.of( [ sample, vcf ] ).filter{ sample, vcf -> vcf.toString().contains("_") }.set{ collect_vcf }
            //Channel.fromList( [ sample, vcf ] ).collate(2).filter{ sample, vcf ->  vcf.toString().contains("_") }.set{ collect_vcf }
            //Channel.fromList( [ sample, cna ] ).collate(2).filter{ sample, cna -> cna.toString().contains(".cns") }.set{ collect_cna }

        }

    }

    //Channel.fromList(vcf_files).flatten().collate(2).filter{ sample, vcf -> vcf.toString().contains("_") }.map{ sample, vcf -> [ sample, vcf ] }.set{ collect_vcf_files }
    Channel.fromList(cna_files).flatten().collate(2).filter{ sample, cna -> cna.toString().contains(".cns") }.map{ sample, cna -> [ sample, cna ] }.set{ collect_cna_files }
    return collect_cna_files
}
