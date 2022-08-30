
include { TABIX_TABIX      } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_BGZIPTABIX } from '../../modules/nf-core/modules/tabix/bgziptabix/main'

workflow INPUT_CHECK {
    take:
    samplesheet // samplesheet file or path to VCF files

    main:
    check_input(samplesheet)
    TABIX_BGZIPTABIX( files.map{ meta, vcf, tbi, cna -> [meta, vcf]} )


    emit:
    files  // channel: [ val(meta), [ vcf.gz] , [ vcf.gz.tbi ] ] OR [ val(meta), [vcf.gz], [vcf.gz.tbi], [CNA] ]
}

def check_input(input){

    // Function performs the following checks:
    // 1. VCF file:
    //    a. Check if the VCF column exists in samplesheet
    //    b. Check if the VCF file exists
    //    c. Check if the VCF file is bgzipped (meta.gzip_vcf = true)
    //    d. Check if the VCF file is tabixed (meta.tabix_vcf = true)
    //
    //    when points c & d are true, bgzip/tabix the file for user
    //
    // 2. CNA file
    //    < when mode == 'pcgr' >
    //    a. Check CNA column exists in samplesheet
    //    b. Check the CNA file exists
    //    < when mode == 'cpsr' >
    //    a. Output empty channel for CNA

    Channel.from(input).splitCsv(header:true, sep:',')
        .map{ row ->

            // Check if the VCF column exists in samplesheet
            if (row.vcf) vcf = file(row.vcf)
            else vcf = 'NA'

            // Exit and tell user to add VCF column to samplesheet
            if(vcf == 'NA'){
                log.error("ERROR: Your input file '(${input})'' does not have a 'vcf' column.\n\nYou must add a 'vcf' column in the samplesheet specifying paths to input VCF files.")
                System.exit(1)
            }

            // Check if the VCF file exists
            if(!file(vcf).exists()){
                log.error("ERROR: Check input file (${input}). VCF file does not exist at path: ${vcf}")
                System.exit(1)
            }else{

                // Capture sample ID
                def meta = [:]
                vcf      = file(row.vcf)
                meta.id  = vcf.simpleName

                // Check if the VCF file is bgzipped
                if(!vcf.toString().endsWith('.gz') && vcf.toString().endsWith('.vcf')){
                    log.warn("The input VCF file '${vcf}' is not bgzipped.")
                    meta.bgzip_vcf = true
                }

                // Check existence of TBI indexed VCF file (!presumed to be in the same directory!)
                // Unsure how this behaves on a cloud instance.
                tbi  = vcf.toString() + '.tbi'
                if(!file(tbi).exists()){
                    log.warn("The input VCF file '${vcf}' is not tabix indexed.")
                    meta.tbi_vcf = true
                    tbi = []
                }else{
                    tbi = [ file(tbi) ]
                }

                // CNA only available in PCGR mode
                if(params.mode.toLowerCase() == 'pcgr'){

                    // Stage CNA (NA in samplesheet evals as null. Explicitly set as 'NA' here)
                    if (row.cna) cna = file(row.cna)
                    else cna = 'NA'

                    // If user does not select CNA_analysis, output empty slot in channel.
                    if(!params.cna_analysis){
                        // Output PCGR channel with empty slot for CNA (so process does not complain about input cardinality)
                        return [ meta, [ file(vcf) ], tbi, [] ]
                    }

                    // If user selects params.cna_analysis but the entries are NA or not valid, exit.
                    if(params.cna_analysis && cna == 'NA'){
                        // Produce Error message, user wants CNA analysis but did not provide file
                        log.error('ERROR: CNA analysis selected but no copy number alteration files provided in samplesheet.')
                        System.exit(1)
                    }else if(params.cna_analysis && !file(cna).exists()){
                        // Produce Error message, user wants CNA analysis but did not provide valid file
                        log.error('ERROR: CNA analysis selected but copy numer alteration file ' + row.cna.toString() + ' does not exist.')
                        System.exit(1)
                    }else{
                        // Valid file for CNA? Output the final channel
                        return [ meta, [ file(vcf) ], tbi, [ file(cna) ] ]
                    }
                }
                // File channel for CPSR
                return  [ meta, [ file(vcf) ], tbi   ]
            }
        }
        .set { files }
}
