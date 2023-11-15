/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// def valid_params = [
//     ena_metadata_fields : ['run_accession', 'experiment_accession', 'library_layout', 'fastq_ftp', 'fastq_md5']
// ]

// def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// // Validate input parameters
// WorkflowSra.initialise(params, log, valid_params)

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { MULTIQC } from '../modules/multiqc/main'

/*
========================================================================================
    BUILD WORKFLOW CHANNELS
========================================================================================
*/

// empty channel to populate with software versions
ch_versions = Channel.empty()

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow POSTPROCESSING {

    take:
    ch_multiqc // channel: path

    main:

    ch_multiqc.view()

    //
    // MODULE: combine QC into one report with multiQC
    //
    MULTIQC (
        ch_multiqc,
        params.multiqc_config,
        params.extra_multiqc_config,
        params.multiqc_logo
    )
    ch_versions = ch_versions.mix(MULTIQC.out.versions.ifEmpty(null))

    emit:
    // reads       = ch_trimmed_fastq      // channel: [ val(meta), fastq ]
    versions  = ch_versions           // channel: path

}

/*
========================================================================================
    THE END
========================================================================================
*/