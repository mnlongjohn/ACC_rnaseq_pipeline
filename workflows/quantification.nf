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

include { SALMON_INDEX  } from '../modules/salmon/index/main'
include { SALMON_QUANT  } from '../modules/salmon/quant/main'

/*
========================================================================================
    BUILD WORKFLOW CHANNELS
========================================================================================
*/

// empty channel to populate with software versions
ch_versions = Channel.empty()

// empty channel to populate with multiqc files
ch_multiqc = Channel.empty()

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow ALIGNMENT {

    take:
    ch_trimmed_fastq // channel: val( meta ), path([ fastq ])

    main:

    //
    // MODULE: Generate index file to quantify files against
    //
    SALMON_INDEX (
        params.fasta,
        params.transcript_fasta
    )
    .index
    .first()
    .set { ch_index }
    ch_versions = ch_versions.mix(SALMON_INDEX.out.versions.first().ifEmpty(null))

    //
    // MODULE: Qunatify fastq reads
    //
    SALMON_QUANT (
        ch_trimmed_fastq,
        ch_index,
        params.star_align_gtf,
        params.transcript_fasta,
        params.salmon_alignment_mode,
        params.salmon_lib_type
    )
    .results
    .set { ch_quant_results }
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first().ifEmpty(null))

    // add the star log file to the
    ch_multiqc = ch_multiqc.mix(SALMON_QUANT.out.json_info)

    emit:
    quant_results   = ch_quant_results       // channel: [ val(meta), fastq ]
    versions        = ch_versions           // channel: path
    multiqc         = ch_multiqc            // channel: [ val(meta), log ]

}
