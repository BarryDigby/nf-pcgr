
include { TABIX_TABIX as TABIX_INPUT_VCF } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIPTABIX as BGZIP_INPUT_VCF } from '../../modules/nf-core/tabix/bgziptabix/main'

workflow INPUT_CHECK {
    take:
    ch_input // staged as file object using file(params.input)

    main:
    if( ch_input.toString().endsWith('.csv') ){
        samplesheet = Channel.of(ch_input)
        check_input(samplesheet)
    }else{
        sarek_files = collect_sarek_files(ch_input)
        sarek_files.collectFile( name: 'constructed_samplesheet.csv', newLine:false, storeDir: "${params.outdir}/pipeline_info", keepHeader: true ){ ids, vcf, cna -> "sample,vcf,cna" + "\n" + "$ids,$vcf,$cna" + "\n"}.set{ constructed_samplesheet }
        samplesheet = constructed_samplesheet.map{ it ->
                                                    samp_file = file(it)
                                                    return samp_file }
        check_input(samplesheet)
    }

    // Determine if vcf files need to be bgzipped and or tabixed.
    // 0: meta, 1: vcf, 2: tbi, 3: cna.
    // must use it instead of names here due to different input tuple len for modes.
    BGZIP_INPUT_VCF( files.map{ it -> [it[0], it[1]]} )
    TABIX_INPUT_VCF( files.map{ it -> [it[0], it[1]]} )

    // If procs were not run, flatten takes care of ifEmpty([]) in step below.
    ch_tabix_bgzip = BGZIP_INPUT_VCF.out.gz_tbi.ifEmpty([])
    ch_tabix_tabix = TABIX_INPUT_VCF.out.tbi.ifEmpty([])

    // Step 3.
    // Using meta as the grouping key, combine the newly bgzipped/tabixed VCF files appropriately.
    // Catch here is to remove the original raw VCF using flatten + filter. Restore tuple using collate.
    // Made the decision to use the same input tuple for (!params.cna_analysis && params.mode = 'cpsr')
    if(params.cna_analysis){
        files.mix(ch_tabix_bgzip, ch_tabix_tabix)
                .groupTuple(by: 0)
                .flatten()
                .filter{ it -> !it.toString().endsWith('.vcf')}
                .collate(4, false)
                .map{ meta, vcf, tbi, cna ->
                        var = [:]
                        var.id      = meta.id
                        var.patient = meta.patient
                        var.status  = meta.status
                        var.sample  = meta.sample
                        var.tool    = meta.tool
                        return [var, vcf, tbi, cna]}
                .set{ ch_files }
    }else{
        files.mix(ch_tabix_bgzip, ch_tabix_tabix)
                .groupTuple(by: 0)
                .flatten()
                .filter{ it -> !it.toString().endsWith('.vcf')}
                .collate(3, false)
                .map{ meta, vcf, tbi ->
                        var = [:]
                        var.id      = meta.id
                        var.patient = meta.patient
                        var.status  = meta.status
                        var.sample  = meta.sample
                        var.tool    = meta.tool
                        return [var, vcf, tbi, [] ]}
                .set{ ch_files }
    }

    ch_files.view()

    emit:
    ch_files  // channel: [ [meta], vcf.gz, vcf.gz.tbi, [] ] OR [ [meta], vcf.gz, vcf.gz.tbi, CNA ]
}

def check_input(input){

    input.splitCsv(header:true, sep:',')
        .map{ row ->

            // Check that the patient column exists in samplesheet
            if (row.patient) patient = row.patient
            else patient = 'NA'

            // Exit and tell user to add Patient column to samplesheet
            if(patient == 'NA'){
                log.error("ERROR: The input file '(${input})'' does not have a 'patient' column.\n\nYou must add a 'patient' column in the samplesheet denoting the patient/subject ID.")
                System.exit(1)
            }

            // Check status column exists in samplesheet
            if (row.status) status = row.status
            else status = 'NA'

            // Exit and tell user to add Status column to samplesheet
            if(status == 'NA'){
                log.error("ERROR: The input file '(${input})'' does not have a 'status' column.\n\nYou must add a 'status' column in the samplesheet denoting the status of the sample: Normal (0), Tumor (1).")
                System.exit(1)
            }

            // Check Status is valid
            if (!status in ['0', '1']){
                log.error("ERROR: The input file '(${input})'' must contain 0 or 1 values for the colum 'status'. Offending entry: ${status}.")
                System.exit(1)
            }

            // Check that sample column exists in samplesheet
            if (row.sample) sample = row.sample
            else sample = 'NA'

            // Exit and tell user to add Sample column to samplesheet
            if(sample == 'NA'){
                log.error("ERROR: The input file '(${input})'' does not have a 'sample' column.\n\nYou must add a 'sample' column in the samplesheet denoting the sample type belonging to the patient.")
                System.exit(1)
            }

            // Check if the VCF column exists in samplesheet
            if (row.vcf) vcf = file(row.vcf)
            else vcf = 'NA'

            // Exit and tell user to add VCF column to samplesheet
            if(vcf == 'NA'){
                log.error("ERROR: The input file '(${input})'' does not have a 'vcf' column.\n\nYou must add a 'vcf' column in the samplesheet specifying paths to input VCF files.")
                System.exit(1)
            }

            // Check if the VCF file exists
            if(!file(vcf).exists()){
                log.error("ERROR: Check input file (${input}). VCF file does not exist at path: ${vcf}")
                System.exit(1)
            }else{

                // Capture metadata
                def meta     = [:]
                meta.patient = patient
                meta.status  = status == "1" ? 'somatic' : 'germline'
                meta.sample  = sample

                // Capture tool name (users must follow sarek naming conventions)
                // This is crucial for properly combining the outputs of BGZIP/TABIX
                vcf = file(row.vcf)
                filename = vcf.getName()
                meta.tool = filename.toString().tokenize('.')[1]
                if( meta.tool == 'strelka' ) { meta.tool = filename.tokenize('.')[1,2].join('.') }

                // meta.id for process tags
                meta.id = "${meta.patient}:${meta.sample}:${meta.tool}"

                // Check if the VCF file is bgzipped
                if(!vcf.toString().endsWith('.gz') && vcf.toString().endsWith('.vcf')){
                    meta.bgzip_vcf = true
                }else{
                    meta.bgzip_vcf = false
                }

                // Check existence of TBI indexed VCF file (presumed to be in the same directory)
                // Unsure how this behaves on a cloud instance.
                tbi  = vcf.toString() + '.tbi'
                if(!file(tbi).exists()){
                    meta.tabix_vcf = true
                    tbi = []
                }else{
                    meta.tabix_vcf = false
                    tbi = [ file(tbi) ]
                }

                // CNA only available in PCGR
                // Error logging for invalid CNA files etc.
                if(meta.status == 'somatic' && params.cna_analysis){

                    // Stage CNA (NA in samplesheet evals as null. Explicitly set as 'NA' here)
                    if (row.cna) cna = file(row.cna)
                    else cna = 'NA'

                    // If user selects params.cna_analysis but the entries are NA or not valid, exit.
                    if(cna == 'NA'){
                        // Produce Error message, user wants CNA analysis but did not provide file
                        log.error('ERROR: CNA analysis selected but no copy number alteration files provided with somatic VCF files in samplesheet.')
                        System.exit(1)
                    }else if(!file(cna).exists()){
                        // Produce Error message, user wants CNA analysis but did not provide valid file
                        log.error('ERROR: CNA analysis selected but copy numer alteration file ' + row.cna.toString() + ' does not exist.')
                        System.exit(1)
                    }
                }
                // File channel for CPSR
                // cna file checks have been run above
                cna_handling = ( meta.status == 'somatic' && file(cna).exists() ) ? [ file(cna) ] : []
                return  [ meta, [ file(vcf) ], tbi, cna_handling ]
            }
        }
        .set{ files }
}

// Important the user sets params.cna_analysis to false if CNVkit was not used in Sarek.
def collect_sarek_files(input){
    // Init empty array outside eachFileRecurse scope
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
    // Filter out the empty tuple slots '[]' in VCF array i.e select appropriate files
    collect_vcf = Channel.fromList( vcf_files ).filter{ ids, vcf -> vcf.toString().contains('.vcf') }
    // As above, with catch for no CNVkit files.
    collect_cna = params.cna_analysis ? Channel.fromList( cna_files ).filter{ ids, cna -> cna.toString().contains('.cns') } : Channel.fromList( cna_files )
    // looking here because workflow exists with no CNS files? Check that --cna_analysis is false!
    sarek_files = collect_vcf.combine(collect_cna, by:0).unique()
    return sarek_files
}

