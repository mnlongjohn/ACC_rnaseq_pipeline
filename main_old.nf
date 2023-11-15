#!/usr/bin/env nextflow
/*
========================================================================================
    Microbiome Analysis Pipeline
========================================================================================
    take in cleaned up fastq of short-read or long-read format and produce 
    the quantity of every refseq stored bacteria and/or virus along with an
    optional Krona plot. Eventually we will add phyolgency contruction as well 
    hopefully leveraging CAIR GPU clusters
========================================================================================
    Github : TOO BE ADDED
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// WorkflowMain.initialise(workflow, params, log)

// Read in fastq from --input file
Channel
    .fromFilePairs(params.input, size: params.single_end ? 1 : 2, checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find any .fastq files matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\n" }
    .map { it -> tuple([id: it[0], single_end: params.single_end], it[1]) }
    .set { fastq_files }

// add a channel to check for a fasta or index file and fail if none is found
// Read in fastq from --input file
// Channel
//     .fromFilePairs(params.input, size: params.single_end ? 1 : 2, checkIfExists: true)
//     .ifEmpty { exit 1, "Cannot find any .fastq files matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\n" }
//     .map { it -> tuple([id: it[0], single_end: params.single_end], it[1]) }
//     .set { fastq_files }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPROCESSING     } from './workflows/preprocessing'
include { ALIGNMENT         } from './workflows/alignment'
include { QUANTIFICATION    } from './workflows/quantification'
include { POSTPROCESSING    } from './workflows/postprocessing'

//
// WORKFLOW: Run main methmotif pipeline depending on the step provided
//
workflow     ACC_RNASEQ_PIPELINE {

    //
    // WORKFLOW: execute preprocessing steps
    //
    PREPROCESSING ( fastq_files )

    //
    // WORKFLOW: execute alignment steps
    //
    ALIGNMENT (PREPROCESSING.out.reads)

    //
    // WORKFLOW: execute quantification steps
    //
    QUANTIFICATION (PREPROCESSING.out.reads) 

    // create a channel containing the multiqc files
    PREPROCESSING.out.multiqc
        .mix(ALIGNMENT.out.multiqc)
        .map { it -> it[1] }
        .collect()
        .set{ ch_multiqc }

    //
    // WORKFLOW: execute postprocessing steps
    //
    POSTPROCESSING (ch_multiqc)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow {
    ACC_RNASEQ_PIPELINE ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/