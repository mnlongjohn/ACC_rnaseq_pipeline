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

include { STAR_GENOMEGENERATE   } from '../modules/star/genomegenerate/main'
include { STAR_ALIGN            } from '../modules/star/align/main'

/*
========================================================================================
    BUILD WORKFLOW CHANNELS
========================================================================================
*/

// channel for the fasta files that we will index
Channel
    .fromPath(params.fasta, type: 'file', checkIfExists: true)
    .map { it -> tuple([id: it.simpleName], it) }
    .set{ ch_fasta }

// empty channel for the star genomegenerate gtf files required
channel
    .of( tuple([id: "blah"], params.star_index_gtf) )
    .first()
    .set { ch_star_index_gtf }

// empty channel for the star align gtf files required
channel
    .of( tuple([id: "blah"], params.star_align_gtf) )
    .first()
    .set { ch_star_align_gtf }

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
    // MODULE: Generate index file to align against
    //
    STAR_GENOMEGENERATE (
        ch_fasta,
        ch_star_index_gtf
    )
    .index
    .first()
    .set { ch_index }
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions.first().ifEmpty(null))

    //
    // MODULE: Analyze fastq QC and trim reads
    //
    STAR_ALIGN (
        ch_trimmed_fastq,
        ch_index,
        ch_star_align_gtf,
        params.star_ignore_sjdbgtf,
        params.seq_platform,
        params.seq_center
    )
    .bam_unsorted
    .set { ch_unsorted_bam }
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first().ifEmpty(null))

    // add the star log file to the
    ch_multiqc = ch_multiqc.mix(STAR_ALIGN.out.log_final)

    emit:
    bam       = ch_unsorted_bam       // channel: [ val(meta), fastq ]
    versions  = ch_versions           // channel: path
    multiqc   = ch_multiqc            // channel: [ val(meta), log ]

}
