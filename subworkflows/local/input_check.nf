//
// Check input and get vcf channels
// Sarek outputs VCF files in the form of:
// ${meta.id}.{tool}.vcf.gz..
// 'simpleName' should suffice here.
// Allow file or path as input, automatically check if VCF files are indexed
//

workflow INPUT_CHECK {
    take:
    samplesheet // samplesheet file or path to VCF files

    main:
    check_input(samplesheet)

    emit:
    files  // channel: [ val(meta), [ vcf.gz] , [ vcf.gz.tbi ] ] OR [ val(meta), [vcf.gz], [vcf.gz.tbi], [CNA] ]
}

def check_input(input){

    Channel.from(input).splitCsv(header:true, sep:',')
        .map{ row ->

            // Stage the VCF files. NA not permitted
            vcf = file(row.vcf)

            // Check if they exist
            if(!file(vcf).exists()){
                log.error("ERROR: Check input file (${input}). VCF file does not exist at path: ${vcf}")
                System.exit(1)
            }else{

                // Capture sample ID
                def meta = [:]
                vcf      = file(row.vcf)
                meta.id  = vcf.simpleName

                // Check existence of TBI indexed VCF file (!presumed to be in the same directory!)
                tbi      = vcf.toString() + '.tbi'
                if(!file(tbi).exists()){
                    log.error("ERROR: Please make sure the VCF files are indexed with tabix. Offending file: ${tbi}")
                    System.exit(1)
                }

                // CNA only available in PCGR mode
                if(params.mode.toLowerCase() == 'pcgr'){

                    // Stage CNA (NA in samplesheet evals as null. Explicitly set as 'NA' here)
                    if (row.cna) cna = file(row.cna)
                    else cna = 'NA'

                    // If user does not select CNA_analysis, output empty slot in channel.
                    if(!params.cna_analysis){
                        // Output PCGR channel with empty slot for CNA (so process does not complain about input cardinality)
                        return [ meta, [ file(vcf) ], [ file(tbi) ], [] ]
                    }

                    // If user selects params.cna_analysis but the entries are NA or not valid, exit.
                    if(params.cna_analysis && cna == 'NA'){
                        // Produce Error message, user wants CNA analysis but did not provide file
                        log.error('ERROR: CNA analysis selected but copy number alteration column contains NA values.')
                        System.exit(1)
                    }else if(params.cna_analysis && !file(cna).exists()){
                        // Produce Error message, user wants CNA analysis but did not provide valid file
                        log.error('ERROR: CNA analysis selected but copy numer alteration file ' + row.cna.toString() + ' does not exist.')
                        System.exit(1)
                    }else{
                        // Valid file for CNA? Output the final channel
                        return [ meta, [ file(vcf) ], [ file(tbi) ], [ file(cna) ] ]
                    }
                }
                // File channel for CPSR
                return  [ meta, [ file(vcf) ], [ file(tbi) ]  ]
            }
        }
        .set { files }
}
