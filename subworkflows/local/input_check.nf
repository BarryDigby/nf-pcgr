
include { TABIX_TABIX      } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_BGZIPTABIX } from '../../modules/nf-core/modules/tabix/bgziptabix/main'

workflow INPUT_CHECK {
    take:
     // samplesheet file or path to files accpeted
    ch_input

    main:
    // Step 1.
    // allow the user to provide a samplesheet, or path to sarek directory.
    // Functions at end of file.
    if( ch_input.toString().endsWith('.csv') ){
        check_input(ch_input)
    }else{
        sarek_files = collect_sarek_files(ch_input)
        sarek_files.collectFile( name: 'constructed_samplesheet.csv', newLine:false, storeDir: "${params.outdir}", keepHeader: true ){ ids, vcf, cna -> "sample,vcf,cna" + "\n" + "$ids,$vcf,$cna" + "\n"}.set{ constructed_samplesheet }
        check_input(constructed_samplesheet)
    }

    // Step 2.
    // Determine if vcf files need to be bgzipped and or tabixed.
    // 0: meta, 1: vcf, 2: tbi, 3: cna.
    // must use it instead of names here due to different input tuple len for modes.
    TABIX_BGZIPTABIX( files.map{ it -> [it[0], it[1]]} )
    TABIX_TABIX( files.map{ it -> [it[0], it[1]]} )

    ch_tabix_bgzip = TABIX_BGZIPTABIX.out.gz_tbi.ifEmpty([])
    ch_tabix_tabix = TABIX_TABIX.out.tbi.ifEmpty([])

    // Step 3.
    // Create channels for CPSR/PCGR. Requires some work as user may choose to omit cna analysis.
    // PCGR:
    // Mix channels, group by meta.id, flatten to remove [] introduced by ifEmpty
    // collate 4 (meta, vcf, tbi, cna), remove metadata associated with tabix logic.
    if(params.mode.toLowerCase() == 'pcgr'){
        // CNA analysis dictates collate #
        if (params.cna_analysis){
            ch_files = files.mix(ch_tabix_bgzip, ch_tabix_tabix)
                            .groupTuple(by: [0,0])
                            .flatten()
                            .collate( 4, false)
                            .map{ meta, vcf, tbi, cna ->
                                    var = [:]
                                    var.id = meta.id
                                    return [var, vcf, tbi, cna]
                            }
        }else{
            // No cna analysis, leave empty slot in tuple.
            ch_files = files.mix(ch_tabix_bgzip, ch_tabix_tabix)
                            .groupTuple(by: [0,0])
                            .flatten()
                            .collate( 3, false)
                            .map{ meta, vcf, tbi ->
                                    var = [:]
                                    var.id = meta.id
                                    return [var, vcf, tbi, [] ]
                            }
        }
    }else{
        // CPSR mode, input tuple len = 3
        ch_files = files.mix(ch_tabix_bgzip, ch_tabix_tabix)
                        .groupTuple(by: [0,0])
                        .flatten()
                        .collate( 3, false)
                        .map{ meta, vcf, tbi ->
                                var = [:]
                                var.id = meta.id
                                return [var, vcf, tbi ]
                        }
    } // open to more elgant solutions to the if else statements above

    emit:
    ch_files  // channel: [ [meta:id], vcf.gz, vcf.gz.tbi, [] ] OR [ [meta:id], vcf.gz, vcf.gz.tbi, CNA ]
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

            // Check that sample column exists in samplesheet
            if (row.sample) sample = row.sample
            else sample = 'NA'

            // Exit and tell user to add Sample column to samplesheet
            if(sample == 'NA'){
                log.error("ERROR: Your input file '(${input})'' does not have a 'sample' column.\n\nYou must add a 'sample' column in the samplesheet denoting the sample ID.")
                System.exit(1)
            }

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
                meta.id  = sample

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
                    meta.tabix_vcf = true
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
                return  [ meta, [ file(vcf) ], tbi ]
            }
        }
        .set{ files }
}

// Important the user sets params.cna_analysis to false if CNVkit was not used in Sarek.
def collect_sarek_files(input){
    vcf_files = []
    cna_files = []
    input.eachFileRecurse{ it ->
        // match tumor vs normal VCF files
        vcf = it.name.contains('_vs_') && ( it.name.endsWith('.vcf') || it.name.endsWith('.vcf.gz') ) && !it.name.endsWith('.tbi') ? file(it) : []
        // Match CNVkit output file OR produce NA for samplesheet
        if(params.cna_analysis){ cna = it.name.contains('.cns') ? file(it) : [] }else{ cna = 'NA' }
        // Only grab IDs for VCF or CNVkit file, else NA
        ids = ( it.name.contains('.cns') || it.name.contains('_vs_') && ( it.name.endsWith('.vcf') || it.name.endsWith('.vcf.gz') ) && !it.name.endsWith('.tbi') ) ? it.simpleName.tokenize('_')[0] : 'NA'

        vcf_files << [ ids, vcf ]
        cna_files << [ ids, cna ]
        }
    // Filter out the empty tuple slots '[]' in VCF array
    collect_vcf = Channel.fromList( vcf_files ).filter{ ids, vcf -> vcf.toString().contains('.vcf') }
    // As above, with catch for no CNVkit files.
    collect_cna = params.cna_analysis ? Channel.fromList( cna_files ).filter{ ids, cna -> cna.toString().contains('.cns') } : Channel.fromList( cna_files )
    sarek_files = collect_vcf.combine(collect_cna, by:0).unique()
    return sarek_files
}

