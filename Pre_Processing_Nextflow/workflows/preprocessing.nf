#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTX_TRIMMER                         } from '../modules/local/fastx_trimmer'
include { DOWNSAMPLING_BAM                      } from '../modules/local/downsampling_bam'
include { SAMTOOLS_VIEW_COUNT                   } from '../modules/local/count_reads'

// Quality control
include { MULTIQC                               } from '../modules/local/multiqc'
include { SATURATION                            } from '../modules/local/saturation'
include { CPG_COVERAGE                          } from '../modules/local/CpG_coverage'
include { CPG_ENRICHMENT                        } from '../modules/local/CpG_enrichment'
include { QC_SUMMARY                            } from '../modules/local/qc_summary'
include { FASTQC_SUMMARY                        } from '../modules/local/fastqc_summary' 
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// MODULE: Installed directly from nf-core/modules
//

include { FASTQC                                } from '../modules/nf-core/fastqc/main'
include { CUTADAPT                              } from '../modules/nf-core/cutadapt/main'
include { UMITOOLS_EXTRACT                      } from '../modules/nf-core/umitools/extract/main'
include { UMITOOLS_DEDUP                        } from '../modules/nf-core/umitools/dedup/main'
include { BOWTIE2_ALIGN                         } from '../modules/nf-core/bowtie2/align/main'
//include { BOWTIE2_BUILD               } from '../modules/nf-core/bowtie2/build/main'
include { SAMTOOLS_VIEW                         } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX                        } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT                         } from '../modules/nf-core/samtools/sort/main'
include { PICARD_MARKDUPLICATES                 } from '../modules/nf-core/picard/markduplicates/main'


// Quality control
include { SAMTOOLS_STATS                        } from '../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS                     } from '../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT                     } from '../modules/nf-core/samtools/flagstat/main'
include { DUMPSOFTWAREVERSIONS                  } from '../modules/nf-core/dumpsoftwareversions/main'
include { PICARD_COLLECTMULTIPLEMETRICS         } from '../modules/nf-core/picard/collectmultiplemetrics/main'

//
// SUBWORKFLOWS: Consisting of nf-core/modules, sometimes with slight adaptions
//
include { INPUT_CHECK                           } from '../subworkflows/input_check'
include { ALIGNEMENT_BOWTIE                     } from '../subworkflows/alignment_bowtie'
include { BAM_FILTER_SAMTOOLS                   } from '../subworkflows/bam_filter_samtools'
include { BAM_SORT_SAMTOOLS                     } from '../subworkflows/bam_sort_samtools'
include { BAM_STATS_SAMTOOLS                    } from '../subworkflows/bam_stats_samtools'
include { DEDUP_UMI_UMITOOLS                    } from '../subworkflows/dedup_umi_umitools'
include { MARK_DUPLICATES_PICARD                } from '../subworkflows/mark_duplicates_picard'
include { UMITOOLS_CUTADAPT                     } from '../subworkflows/umi_extract_trimming'
include { TRIMMING_FASTQC                       } from '../subworkflows/trimming'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN Preprocessing Workflow: rawdata (fastq) to ready-to-use bam-file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/   

// define channel for input (fastq files; seq raw data) and for human reference genome

reads_ch = Channel.fromPath("${params.meta_file}")                   
                .splitCsv(header:true, sep:',')
                .map {create_fastq_and_meta_channel(it) }
                .set { reads }   
   
//if ( params.alternative_dataset ) {
    def create_fastq_and_meta_channel(LinkedHashMap row) {
        // create meta map
        def meta = [:]
        meta.id             = row.SAMPLE_ID
        meta.sample         = row.SAMPLE_NAME
        meta.patient        = row.PATIENT_ID
        meta.sex            = row.SEX
        meta.phenotype      = row.PHENOTYPE
        meta.project        = row.PROJECT
        meta.library        = row.LIBRARY_TYPE
        meta.libraryshort   = row.LIBRARY_TYPE_SHORT
        meta.single_paired  = row.SEQUENCING_READ_TYPE
        meta.seq_platform   = row.INSTRUMENT_PLATFORM
        meta.file_id_R1     = row.file_accession_id_R1
        meta.file_id_R2     = row.file_accession_id_R2
        meta.single_end     = false

    // add path(s) of the fastq file(s) to the meta map 
        def fastq_meta = []
    
        if ((row.FASTQ_1) == false) {
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.FASTQ_1}"
        }
        if (meta.single_end) {
            fastq_meta = [ meta, [ (row.FASTQ_1="${params.rawfastq}/${row.FASTQ_1}") ] ]
        } else {
            if ((row.FASTQ_2) == false) {
                exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.FASTQ_2}"
            }
            fastq_meta = [ meta, [row.FASTQ_1="${params.rawfastq}/${row.file_accession_id_R1}/${row.FASTQ_1}", row.FASTQ_2="${params.rawfastq}/${row.file_accession_id_R2}/${row.FASTQ_2}" ] ]
        }
        return fastq_meta
    }
//} else {
        // Function to get list of [ meta, [ FASTQ_1, FASTQ_2, FASTQ_UMI, ... ] ]
//    def create_fastq_and_meta_channel(LinkedHashMap row) {
//            // create meta map
//           def meta = [:]
//            meta.id             = row.SAMPLE_ID
//            meta.sample         = row.SAMPLE_NAME
//            meta.patient        = row.PATIENT_ID
//            meta.source         = row.SOURCE
//            meta.tissue         = row.TISSUE_TYPE
//            meta.library        = row.LIBRARY_TYPE
//            meta.lane           = row.LANE_NO
//            meta.run            = row.RUN_ID
//            meta.index          = row.INDEX
//            meta.seq_center     = row.CENTER_NAME
//            meta.single_paired  = row.SEQUENCING_READ_TYPE
//            meta.seq_platform   = row.INSTRUMENT_PLATFORM
//            meta.instrument     = row.INSTRUMENT_MODEL
//            meta.single_end     = row.FASTQ_2.toBoolean()

        // add path(s) of the fastq file(s) to the meta map 
//        def fastq_meta = []
        
//        if(params.with_umi) {   
//            if ((row.FASTQ_1) == false) {
//                exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.FASTQ_1}"
//            }
//            if (meta.single_end) {
//                fastq_meta = [ meta, [ (row.FASTQ_1="${params.rawfastq}/${row.FASTQ_1}") ] ]
                //fastq_meta = [ meta, [ row.FASTQ_1="${params.rawfastq}/${row.FASTQ_1}", row.FASTQ_UMI="${params.rawfastq}/${row.FASTQ_UMI}" ] ]
//            } else {
//                if ((row.FASTQ_2) == false) {
//                    exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.FASTQ_2}"
//                }
//                fastq_meta = [ meta, [row.FASTQ_1="${params.rawfastq}/${row.FASTQ_1}", row.FASTQ_2="${params.rawfastq}/${row.FASTQ_2}", row.FASTQ_UMI="${params.rawfastq}/${row.FASTQ_UMI}"] ]
//            }
        
//        } else {
//            if ((row.FASTQ_1) == false) {
//                exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.FASTQ_1}"
//            }
//            if (meta.single_end) {
//                fastq_meta = [ meta, [ (row.FASTQ_1="${params.rawfastq}/${row.FASTQ_1}") ] ]
//            } else {
//                if ((row.FASTQ_2) == false) {
//                    exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.FASTQ_2}"
//                }
//                fastq_meta = [ meta, [row.FASTQ_1="${params.rawfastq}/${row.FASTQ_1}", row.FASTQ_2="${params.rawfastq}/${row.FASTQ_2}" ] ]
//            }
//        }
//        return fastq_meta   
//    }    
//}

//Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(', ')}"
}

// Prepare channels for fasta-File and for bowtie2-index files
if ( params.genome ) {
    genome_fasta = params.genomes[ params.genome ].fasta ?: false
    ch_fasta = Channel.value(file( "${genome_fasta}" ))
} else if ( params.fasta ) {
    ch_fasta = Channel.value(file( "${params.fasta}" ))
} else {
    ch_fasta = Channel.empty()
}

if ( params.genome ) {
    bwt2_index = params.genomes[ params.genome ].bowtie2 ?: false
    bwt2_indices = Channel.value(file( "${bwt2_index}/*" ))
} else if ( params.bwt2_index ) {
   bwt2_indices = Channel.empty()
   bwt2_indices = file("${params.bwt2_index}")
} else {
    bwt2_indices = Channel.empty()
}

//
/// Run Workflow
//

workflow PREPROCESSING {
    
    ch_versions = Channel.empty()

    // Step 0: Homogenize all reads to a length of 100bp (optional) 
    // MODULE: Run FASTX (fastx-toolkit)
    //
    
    ch_fastx_reads = Channel.empty()
    
    if(!params.skip_fastx) {
        FASTX_TRIMMER(reads)
        ch_fastx_reads = FASTX_TRIMMER.out.reads
        raw_reads = ch_fastx_reads
    }

    if(params.skip_fastx) {
        raw_reads = reads
    }
    
    // Step 1: (subworkflow) UMI extraction from I1 and addition to seq reads (optional); adapter and quality trimming
    // MODULE: Run UMI_EXTRACT, CUTADAPT
    //

    ch_reads_for_mapping = Channel.empty()
    ch_cutadapt_log = Channel.empty()

    if (params.with_umi && !params.skip_umi_extract) {
        UMITOOLS_CUTADAPT(raw_reads, params.skip_fastqc_raw, params.skip_fastqc_trim, params.skip_trimming)
        ch_reads_for_mapping = UMITOOLS_CUTADAPT.out.reads
    
        if (!params.skip_trimming) {
            ch_cutadapt_log = UMITOOLS_CUTADAPT.out.trim_log
        }

        ch_versions = ch_versions.mix(UMITOOLS_CUTADAPT.out.versions)
    }

    if (!params.with_umi && !params.skip_trimming) {   
            TRIMMING_FASTQC(raw_reads, params.skip_fastqc_raw, params.skip_fastqc_trim)
            ch_reads_for_mapping = TRIMMING_FASTQC.out.reads
        
            if (!params.skip_trimming) {
            ch_cutadapt_log = TRIMMING_FASTQC.out.trim_log
            }
            
            ch_versions = ch_versions.mix(TRIMMING_FASTQC.out.versions)
        }
    
    if ((params.skip_umi_extract || !params.with_umi ) & params.skip_trimming) {
        ch_reads_for_mapping = raw_reads 
    }    

    // Step 2: (subworkflow) alignment/mapping to human reference genome; GrCh37/hg19 , (sub-subworkflow) coordinate-sorting, indexing
    // MODULE: Run BOWTIE2_ALIGN 
    //
    ch_samtools_stats_bowtie2 = Channel.empty()
    ch_samtools_flagstats_bowtie2 = Channel.empty()
    ch_samtools_idxstats_bowtie2 = Channel.empty()
   
    ALIGNEMENT_BOWTIE(ch_reads_for_mapping, bwt2_indices, params.set_readgroups, params.save_unaligned, params.skip_samtools_stats_bowtie2, params.skip_fastqc_map)

    ch_bowtie2_log = ALIGNEMENT_BOWTIE.out.log_summary

    // quality control of aligned reads
    if (!params.skip_samtools_stats_bowtie2) {
        ch_samtools_stats_bowtie2 = ALIGNEMENT_BOWTIE.out.stats
        ch_samtools_flagstats_bowtie2 = ALIGNEMENT_BOWTIE.out.flagstat
        ch_samtools_idxstats_bowtie2 = ALIGNEMENT_BOWTIE.out.idxstats
    }
    
    ch_versions = ch_versions.mix(ALIGNEMENT_BOWTIE.out.versions)

    // Step 3: (subworkflow) quality filtering and indexing 
    // MODULE: Run SAMTOOLS_VIEW
    //

    ch_reads_for_umidedup = Channel.empty()
    ch_reads_for_markdup = Channel.empty()

    if (!params.skip_filtering) {
        BAM_FILTER_SAMTOOLS(ALIGNEMENT_BOWTIE.out.bam_bai, ch_fasta)
        
        ch_reads_for_umidedup = BAM_FILTER_SAMTOOLS.out.bam.join(BAM_FILTER_SAMTOOLS.out.bai, by: [0])
        ch_reads_for_markdup = BAM_FILTER_SAMTOOLS.out.bam
    }
    ch_versions = ch_versions.mix(BAM_FILTER_SAMTOOLS.out.versions)

    
    if (params.skip_filtering) {
        ch_reads_for_umidedup = ALIGNEMENT_BOWTIE.out.bam_bai
        ch_reads_for_markdup = ALIGNEMENT_BOWTIE.out.bam
    }

    // Step 4: DEDUPLICATION 
    // Deduplication with UMI (subworkflow DEDUP_UMI_UMITOOLS), indexing of sorted output
    //
   
    ch_samtools_stats_final = Channel.empty()
    ch_samtools_flagstats_final = Channel.empty()
    ch_samtools_idxstats_final = Channel.empty()

    ch_fastqc_final = Channel.empty()
    
    if (params.with_umi && params.skip_markduplicates) { 
        DEDUP_UMI_UMITOOLS(ch_reads_for_umidedup, params.umitools_dedup_stats, params.skip_samtools_stats_dedup, params.skip_fastqc_final) 
    
        if (!params.skip_samtools_stats_dedup) {
        ch_samtools_stats_final = DEDUP_UMI_UMITOOLS.out.stats
        ch_samtools_flagstats_final = DEDUP_UMI_UMITOOLS.out.flagstat
        ch_samtools_idxstats_final = DEDUP_UMI_UMITOOLS.out.idxstats
        }

    // quality control of final deduplicated BAM-files with FASTQC
        if (!params.skip_fastqc_final) {
            ch_fastqc_final = DEDUP_UMI_UMITOOLS.out.fastqc_final_zip
        }

        ch_versions = ch_versions.mix(DEDUP_UMI_UMITOOLS.out.versions)
    }

    // Deduplication with PICARD_MARKDUPLICATES (subworkflow MARK_DUPLICATES_PICARD), indexing of sorted output
    //
    ch_markduplicates_metrics = Channel.empty()
    
    if (!params.with_umi && !params.skip_markduplicates) {
        MARK_DUPLICATES_PICARD(ch_reads_for_markdup, params.skip_samtools_stats_dedup, params.skip_fastqc_final) 

        ch_markduplicates_metrics = MARK_DUPLICATES_PICARD.out.metrics
        
        if (!params.skip_samtools_stats_dedup) {
        ch_samtools_stats_final = MARK_DUPLICATES_PICARD.out.stats
        ch_samtools_flagstats_final = MARK_DUPLICATES_PICARD.out.flagstat
        ch_samtools_idxstats_final = MARK_DUPLICATES_PICARD.out.idxstats
        }

        // quality control of final deduplicated BAM-files with FASTQC
        if (!params.skip_fastqc_final) {
            ch_fastqc_final = MARK_DUPLICATES_PICARD.out.fastqc_final_zip
        }

        ch_versions = ch_versions.mix(MARK_DUPLICATES_PICARD.out.versions)
    }

    ///// Perform Downsampling to 10 Million paired reads (= 20 Mio total reads) by default
    
    if(!params.skip_downsampling) {
    
        if (params.with_umi && params.skip_markduplicates) { 
            DOWNSAMPLING_BAM(DEDUP_UMI_UMITOOLS.out.bam.join(DEDUP_UMI_UMITOOLS.out.bai, by: [0]) , params.downsamp)
            ch_downsamp_summary = DOWNSAMPLING_BAM2.out.summary_tsv
        }
     
        if (!params.with_umi && !params.skip_markduplicates) {
            DOWNSAMPLING_BAM(MARK_DUPLICATES_PICARD.out.bam.join(MARK_DUPLICATES_PICARD.out.bai, by: [0]) , params.downsamp)
            ch_downsamp_summary = DOWNSAMPLING_BAM2.out.summary_tsv
        }
    
        DOWNSAMPLING_BAM.out.summary_tsv.set{ch_downsamp_summary}
        ch_downsamp_summary.collectFile(name: 'readCount_downsamp.tsv', newLine:true, storeDir: "${params.preprocessing}/${params.runID}/downsamp/${params.downsamp}M/")
    
    } 
    
    /*
    ========================================================================================
        Quality Control for Pre-Processing
    ========================================================================================
    */

    // Run Picard CollectMultipleMetrics on mapped, coordinate-sorted BAM-files
    ch_picard_multiplemetrics = Channel.empty()

    if (!params.skip_picard_multiplemetrics) {
        PICARD_COLLECTMULTIPLEMETRICS(ALIGNEMENT_BOWTIE.out.bam , ch_fasta)
        ch_picard_multiplemetrics = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
    
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)
    }
   
    // Count reads
    ch_coverage_summary = Channel.empty()

    if (params.with_umi && params.skip_markduplicates) {
        SAMTOOLS_VIEW_COUNT(DEDUP_UMI_UMITOOLS.out.bam.join(DEDUP_UMI_UMITOOLS.out.bai, by: [0]))
        ch_coverage_summary = SAMTOOLS_VIEW_COUNT.out.summary_tsv
    }
    
    if (!params.with_umi && !params.skip_markduplicates) {
        SAMTOOLS_VIEW_COUNT(MARK_DUPLICATES_PICARD.out.bam.join(MARK_DUPLICATES_PICARD.out.bai, by: [0]))
        ch_coverage_summary = SAMTOOLS_VIEW_COUNT.out.summary_tsv
    }
     
    // Save single files for each sample into one common summary file
    SAMTOOLS_VIEW_COUNT.out.summary_tsv.collectFile(name: 'readCount_samtools.tsv', newLine:true, storeDir: "${params.preprocessing}/${params.runID}/QC/final_bam/")
                                        .set{readcount_summary}

    
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Run Methylation QC
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if(!params.skip_methylation_qc) {
    
    ///
    // Saturation analysis

    if (params.with_umi && params.skip_markduplicates) {
        SATURATION(DEDUP_UMI_UMITOOLS.out.bam.collect{it[1]}, params.single_end)
    }
    
    if (!params.with_umi && !params.skip_markduplicates) {
        SATURATION(MARK_DUPLICATES_PICARD.out.bam.collect{it[1]}, params.single_end)
    }

    ///
    // Evaluation of CpG Coverage

    if (params.with_umi && params.skip_markduplicates) {
        CPG_COVERAGE(DEDUP_UMI_UMITOOLS.out.bam.collect{it[1]}, params.single_end)
    }
    
    if (!params.with_umi && !params.skip_markduplicates) {
        CPG_COVERAGE(MARK_DUPLICATES_PICARD.out.bam.collect{it[1]}, params.single_end)
    }

    ///
    // Calculation of Methylation enrichment score

    if (params.with_umi && params.skip_markduplicates) {
        CPG_ENRICHMENT(DEDUP_UMI_UMITOOLS.out.bam.collect{it[1]}, params.single_end)
    }
    
    if (!params.with_umi && !params.skip_markduplicates) {
        CPG_ENRICHMENT(MARK_DUPLICATES_PICARD.out.bam.collect{it[1]}, params.single_end)
    }
    }
    
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Create Summaries
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    /// !! This part is currently not running properly !! 

    ///
    // Create Summary tsv-File 

    // Input directory for QC results
    // FASTQC_SUMMARY(ch_fastqc_final.collect{it[1]})
    
   // QC_SUMMARY(SATURATION.out.saturation_summary, CPG_COVERAGE.out.cpg_coverage_summary , CPG_ENRICHMENT.out.cpg_enrichment_summary , readcount_summary, FASTQC_SUMMARY.out.fastqc_summary, params.skip_methylation_qc)

    
    ///
    // collect software versions
    //

    //DUMPSOFTWAREVERSIONS (ch_versions.unique().collectFile(name: 'collated_versions.yml'))

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Run Multiqc for subsets of fastqc files
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    // Load up and check multiqc base config and (additional) custom configs
    def multiqc_report = []
    
    ch_multiqc_config        = Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

    if (!params.skip_multiqc_preprocess) {
        MULTIQC ( 
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            ch_fastqc_final.collect{it[1]}.ifEmpty([]),
            ch_picard_multiplemetrics.collect{it[1]}.ifEmpty([]),
            ch_cutadapt_log.collect{it[1]}.ifEmpty([]), 
            ch_bowtie2_log.collect{it[1]},
            ch_samtools_stats_bowtie2.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstats_bowtie2.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats_bowtie2.collect{it[1]}.ifEmpty([]),
            ch_samtools_stats_final.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstats_final.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats_final.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_metrics.collect{it[1]}.ifEmpty([])
        )
        
        multiqc_report = MULTIQC.out.report.toList()
    }

}