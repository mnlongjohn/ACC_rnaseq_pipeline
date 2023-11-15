/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// def valid_params = [
//     aligners       : ['star_salmon'],
//     trimmers       : ['fastp'],
//     pseudoaligners : ['salmon'],
//     // rseqc_modules  : ['bam_stat', 'inner_distance', 'infer_experiment', 'junction_annotation', 'junction_saturation', 'read_distribution', 'read_duplication', 'tin']
// ]

// def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// // Validate input parameters
// WorkflowRnaseq.initialise(params, log, valid_params)

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input, params.multiqc_config,
    params.fasta, params.transcript_fasta, params.additional_fasta,
    params.gtf, params.gff, params.gene_bed,
    params.star_index, params.salmon_index
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
// ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
// ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
// ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

// // Header files for MultiQC
// ch_pca_header_multiqc        = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
// ch_clustering_header_multiqc = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
// ch_biotypes_header_multiqc   = file("$projectDir/assets/multiqc/biotypes_header.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
// include { BEDTOOLS_GENOMECOV                 } from '../modules/local/bedtools_genomecov'
// include { DESEQ2_QC as DESEQ2_QC_STAR_SALMON } from '../modules/local/deseq2_qc'
// include { DESEQ2_QC as DESEQ2_QC_SALMON      } from '../modules/local/deseq2_qc'
// include { DUPRADAR                           } from '../modules/local/dupradar'
// include { MULTIQC                            } from '../modules/local/multiqc'
// include { MULTIQC_CUSTOM_BIOTYPE             } from '../modules/local/multiqc_custom_biotype'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK    } from '../subworkflows/local/input_check'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
// include { ALIGN_STAR     } from '../subworkflows/local/align_star'
// include { QUANTIFY_SALMON as QUANTIFY_STAR_SALMON } from '../subworkflows/local/quantify_salmon'
// include { QUANTIFY_SALMON as QUANTIFY_SALMON      } from '../subworkflows/local/quantify_salmon'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
// include { SAMTOOLS_SORT               } from '../modules/nf-core/samtools/sort/main'
// include { SUBREAD_FEATURECOUNTS       } from '../modules/nf-core/subread/featurecounts/main'
// include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
// include { FASTQ_SUBSAMPLE_FQ_SALMON        } from '../subworkflows/nf-core/fastq_subsample_fq_salmon/main'
// include { FASTQ_FASTQC_UMITOOLS_FASTP      } from '../subworkflows/nf-core/fastq_fastqc_umitools_fastp/main'
// include { BAM_SORT_STATS_SAMTOOLS          } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'
// include { BAM_MARKDUPLICATES_PICARD        } from '../subworkflows/nf-core/bam_markduplicates_picard/main'
// include { BAM_RSEQC                        } from '../subworkflows/nf-core/bam_rseqc/main'
// include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME        } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main'
// include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main'
// include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD } from '../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig/main'
// include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE } from '../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// // Info required for completion email and summary
// def multiqc_report     = []
// def pass_mapped_reads  = [:]
// def pass_trimmed_reads = [:]
// def pass_strand_check  = [:]

workflow RNASEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
    PREPARE_GENOME (
        params.fasta,
        params.gtf,
        params.gff,
        params.additional_fasta,
        params.transcript_fasta,
        params.gene_bed,
        params.star_index,
        params.salmon_index,
        params.gencode,
        biotype,
        prepareToolIndices
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .reads
    .map {
        meta, fastq ->
            new_id = meta.id - ~/_T\d+/
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
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

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

    // // Branch FastQ channels if 'auto' specified to infer strandedness
    // ch_cat_fastq
    //     .branch {
    //         meta, fastq ->
    //             auto_strand : meta.strandedness == 'auto'
    //                 return [ meta, fastq ]
    //             known_strand: meta.strandedness != 'auto'
    //                 return [ meta, fastq ]
    //     }
    //     .set { ch_strand_fastq }
}