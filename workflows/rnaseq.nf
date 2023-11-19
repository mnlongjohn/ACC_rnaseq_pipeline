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

// Check alignment parameters
def prepareToolIndices  = []
if (!params.skip_alignment) { prepareToolIndices << params.aligner }
// if (!params.skip_pseudo_alignment && params.pseudo_aligner) { prepareToolIndices << params.pseudo_aligner }

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
include { ALIGN_STAR     } from '../subworkflows/local/align_star'
include { QUANTIFY_SALMON as QUANTIFY_STAR_SALMON } from '../subworkflows/local/quantify_salmon'
include { QUANTIFY_SALMON as QUANTIFY_SALMON      } from '../subworkflows/local/quantify_salmon'

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
include { FASTQ_SUBSAMPLE_FQ_SALMON        } from '../subworkflows/nf-core/fastq_subsample_fq_salmon/main'
include { FASTQ_FASTQC_UMITOOLS_FASTP      } from '../subworkflows/nf-core/fastq_fastqc_umitools_fastp/main'
include { BAM_SORT_STATS_SAMTOOLS          } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'
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

// Info required for completion email and summary
def multiqc_report     = []
def pass_mapped_reads  = [:]
def pass_trimmed_reads = [:]
def pass_strand_check  = [:]

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

    // Branch FastQ channels if 'auto' specified to infer strandedness
    ch_cat_fastq
        .branch {
            meta, fastq ->
                auto_strand : meta.strandedness == 'auto'
                    return [ meta, fastq ]
                known_strand: meta.strandedness != 'auto'
                    return [ meta, fastq ]
        }
        .set { ch_strand_fastq }

    //
    // SUBWORKFLOW: Sub-sample FastQ files and pseudo-align with Salmon to auto-infer strandedness
    //
    // Return empty channel if ch_strand_fastq.auto_strand is empty so salmon index isn't created
    PREPARE_GENOME.out.fasta
        .combine(ch_strand_fastq.auto_strand)
        .map { it.first() }
        .first()
        .set { ch_genome_fasta }

    FASTQ_SUBSAMPLE_FQ_SALMON (
        ch_strand_fastq.auto_strand,
        ch_genome_fasta,
        PREPARE_GENOME.out.transcript_fasta,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.salmon_index,
        !params.salmon_index && !('salmon' in prepareToolIndices)
    )
    ch_versions = ch_versions.mix(FASTQ_SUBSAMPLE_FQ_SALMON.out.versions)

    FASTQ_SUBSAMPLE_FQ_SALMON
        .out
        .json_info
        .join(ch_strand_fastq.auto_strand)
        .map { meta, json, reads ->
            return [ meta + [ strandedness: WorkflowRnaseq.getSalmonInferredStrandedness(json) ], reads ]
        }
        .mix(ch_strand_fastq.known_strand)
        .set { ch_strand_inferred_fastq }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with fastp
    //
    ch_filtered_reads      = Channel.empty()
    ch_fastqc_raw_multiqc  = Channel.empty()
    ch_fastqc_trim_multiqc = Channel.empty()
    ch_trim_log_multiqc    = Channel.empty()
    ch_trim_read_count     = Channel.empty()
    if (params.trimmer == 'fastp') {
        FASTQ_FASTQC_UMITOOLS_FASTP (
            ch_strand_inferred_fastq,
            params.skip_fastqc || params.skip_qc,
            params.with_umi,
            params.skip_umi_extract,
            params.umi_discard_read,
            params.skip_trimming,
            [],
            params.save_trimmed,
            params.save_trimmed,
            params.min_trimmed_reads
        )
        ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
        ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip
        ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip
        ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json
        ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_read_count
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)
    }

    //
    // Get list of samples that failed trimming threshold for MultiQC report
    //
    ch_trim_read_count
        .map {
            meta, num_reads ->
                pass_trimmed_reads[meta.id] = true
                if (num_reads <= params.min_trimmed_reads.toFloat()) {
                    pass_trimmed_reads[meta.id] = false
                    return [ "$meta.id\t$num_reads" ]
                }
        }
        .collect()
        .map {
            tsv_data ->
                def header = ["Sample", "Reads after trimming"]
                WorkflowRnaseq.multiqcTsvFromList(tsv_data, header)
        }
        .set { ch_fail_trimming_multiqc }

    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon
    //
    ch_genome_bam                 = Channel.empty()
    ch_genome_bam_index           = Channel.empty()
    ch_samtools_stats             = Channel.empty()
    ch_samtools_flagstat          = Channel.empty()
    ch_samtools_idxstats          = Channel.empty()
    ch_star_multiqc               = Channel.empty()
    ch_aligner_pca_multiqc        = Channel.empty()
    ch_aligner_clustering_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'star_salmon') {
        ALIGN_STAR (
            ch_filtered_reads,
            PREPARE_GENOME.out.star_index,
            PREPARE_GENOME.out.gtf,
            params.star_ignore_sjdbgtf,
            '',
            params.seq_center ?: '',
            PREPARE_GENOME.out.fasta.map { [ [:], it ] }
        )
        ch_genome_bam        = ALIGN_STAR.out.bam
        ch_genome_bam_index  = ALIGN_STAR.out.bai
        ch_transcriptome_bam = ALIGN_STAR.out.bam_transcript
        ch_samtools_stats    = ALIGN_STAR.out.stats
        ch_samtools_flagstat = ALIGN_STAR.out.flagstat
        ch_samtools_idxstats = ALIGN_STAR.out.idxstats
        ch_star_multiqc      = ALIGN_STAR.out.log_final
        if (params.bam_csi_index) {
            ch_genome_bam_index = ALIGN_STAR.out.csi
        }
        ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)

        //
        // SUBWORKFLOW: Count reads from BAM alignments using Salmon
        //
        QUANTIFY_STAR_SALMON (
            ch_transcriptome_bam,
            ch_dummy_file,
            PREPARE_GENOME.out.transcript_fasta,
            PREPARE_GENOME.out.gtf,
            true,
            params.salmon_quant_libtype ?: ''
        )
        ch_versions = ch_versions.mix(QUANTIFY_STAR_SALMON.out.versions)

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_STAR_SALMON (
                QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_aligner_pca_multiqc        = DESEQ2_QC_STAR_SALMON.out.pca_multiqc
            ch_aligner_clustering_multiqc = DESEQ2_QC_STAR_SALMON.out.dists_multiqc
            ch_versions = ch_versions.mix(DESEQ2_QC_STAR_SALMON.out.versions)
        }
    }
}