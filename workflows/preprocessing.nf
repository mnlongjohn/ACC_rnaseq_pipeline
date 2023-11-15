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

include { CAT_FASTQ                             } from '../modules/command_line/cat_fastq/main'
include { FASTP                                 } from '../modules/fastp/main'
include { FASTQC as PRE_FASTQC                  } from '../modules/fastqc/main'
include { FASTQC as POST_FASTQC                 } from '../modules/fastqc/main'

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

workflow PREPROCESSING {

    take:
    fastq_files // channel: val( meta ), path([ fastq ])

    main:

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    // INPUT_CHECK (
    //     ch_input
    // )
    // .reads
    fastq_files
        .map {
            meta, fastq ->
                new_id = meta.id - ~/_L\d+/
                [ meta + [id: new_id], fastq ]
        }
        .groupTuple()
        .branch {
            meta, fastq ->
                single  : fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastq }
    // ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    //
    // MODULE: Analyze fastq QC before trimming
    //
    PRE_FASTQC (
        ch_cat_fastq
    )

    // add the pre fastqc zip file
    ch_multiqc = ch_multiqc.mix(PRE_FASTQC.out.zip)

    //
    // MODULE: Analyze fastq QC and trim reads
    //
    FASTP (
        ch_cat_fastq,
        params.adapter_fasta,
        params.save_trimmed_fail,
        params.save_merged
    )
    .reads
    .set { ch_trimmed_fastq }
    ch_versions = ch_versions.mix(FASTP.out.versions.first().ifEmpty(null))

    // add the fastp json file
    ch_multiqc = ch_multiqc.mix(FASTP.out.json)

    //
    // MODULE: Analyze fastq QC after trimming
    //
    POST_FASTQC (
        ch_trimmed_fastq
    )
    ch_versions = ch_versions.mix(POST_FASTQC.out.versions.first().ifEmpty(null))

    // add the post fastqc zip file
    ch_multiqc = ch_multiqc.mix(POST_FASTQC.out.zip)

    emit:
    reads       = ch_trimmed_fastq      // channel: [ val(meta), fastq ]
    versions    = ch_versions           // channel: path
    multiqc     = ch_multiqc            // channel: [ val(meta), log ]

}

/*
========================================================================================
    THE END
========================================================================================
*/