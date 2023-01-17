
include { TABIX_TABIX as TABIX_INPUT_VCF } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIPTABIX as BGZIP_INPUT_VCF } from '../../modules/nf-core/tabix/bgziptabix/main'

workflow INPUT_CHECK {
    take:
    ch_input

    main:
    ch_versions = Channel.empty()
    samplesheet = Channel.of(ch_input)
    check_input(samplesheet)

    samplesheet_files.multiMap{ meta, vcf, tbi, cna ->
                                vcf_files : [ meta, vcf, tbi ]
                                cna_files : [ meta, cna ]
                                }.set{ files }

    // simplify metadata in CNA channels for re-merging to processed VCF files.
    vcf_files = files.vcf_files
    ch_cna_files = files.cna_files.map{ meta, cna -> var = [:]; var.patient = meta.patient; var.status = meta.status; var.sample = meta.sample; return [ meta, cna ] }

    // BGZIP/TABIX based on metadata generated from samplesheet scan.
    BGZIP_INPUT_VCF( vcf_files.map{ meta, vcf, tbi -> [ meta, vcf ] } )
    TABIX_INPUT_VCF( vcf_files.map{ meta, vcf, tbi -> [ meta, vcf ] } )

    // Using meta as the grouping key, combine the newly bgzipped/tabixed VCF files appropriately.
    // Catch here is to remove the original raw VCF using flatten + filter. Restore tuple using collate.
    vcf_files.mix( BGZIP_INPUT_VCF.out.gz_tbi.ifEmpty([]), TABIX_INPUT_VCF.out.tbi.ifEmpty([]) )
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
                        return [var, vcf, tbi ]}
                .set{ ch_vcf_files }

    ch_versions = ch_versions.mix( BGZIP_INPUT_VCF.out.versions )
    ch_versions = ch_versions.mix( TABIX_INPUT_VCF.out.versions )

    emit:
    ch_vcf_files // channel: [ [meta], vcf.gz, vcf.gz.tbi ]
    ch_cna_files //  channel: [meta, cna]
    versions = ch_versions
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
                meta.id = "${meta.patient}.${meta.sample}.${meta.tool}"

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
                cna_handling = ( meta.status == 'somatic' && params.cna_analysis ) ? file(cna) : []
                return [ meta, [ file(vcf) ], tbi, cna_handling ]
            }
        }
        .set{ samplesheet_files }
}
