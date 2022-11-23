/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Stage
//println(params.mode.toLowerCase())
if (params.input) { ch_input = file(params.input, checkIfExists:true) } else { exit 1, 'Please provide an input samplesheet or path to Sarek results' }
if (params.mode.toLowerCase() == 'pcgr' && params.fasta) { ch_fasta = Channel.fromPath(params.fasta, checkIfExists:true) }

// Add Database as channel. String paths do not work on AWS.
// Only stage the genome db we are interested in.
ch_pcgr_dir = Channel.fromPath("${params.database}/data/${params.genome}")
ch_pcgr_dir.view()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK   } from '../subworkflows/local/input_check'
include { FORMAT_FILES  } from '../subworkflows/local/format_files'
include { MERGE_VCFS    } from '../subworkflows/local/merge_vcfs'

include { PCGR as RUN_PCGR } from '../modules/local/PCGR/Run/pcgr'
include { CPSR as RUN_CPSR } from '../modules/local/PCGR/Run/cpsr'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
//include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow PCGR {

    ch_versions = Channel.empty()

    // Read samplesheet/directory, create channel with meta.id, vcf, vcf_tbi, cna file
    INPUT_CHECK (
        ch_input
    )

    // TODO: Run CSPR and give outputs to pcgr. diverging these procs based on params.mode not sensible

    // Automatically add INFO and HEADER fields to somatic VCF files
    // Detect CNA tool used and format for PCGR
    // RUN PCGR
    if(params.mode.toLowerCase() == 'pcgr'){
        FORMAT_FILES(
            ch_fasta.collect(), INPUT_CHECK.out.ch_files
        )
        FORMAT_FILES.out.files
        MERGE_VCFS( FORMAT_FILES.out.files, ch_fasta.collect() )
        RUN_PCGR(
            MERGE_VCFS.out.pcgr_ready_vcf,
            ch_pcgr_dir
        )
    }

    if(params.mode.toLowerCase() == 'cpsr') RUN_CPSR( INPUT_CHECK.out.ch_files, ch_pcgr_dir )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
