#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Define the default parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All of the default parameters are being set in `nextflow.config`

*/ 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Function which prints help message text
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/ 
def helpMessage() {
    log.info"""
MEDIP AND WGS PROCESSING  -  N F    v 1.0 
================================
Usage:
nextflow run generalprocessing_nextflow/main.nf <ARGUMENTS>
/// Required Arguments:
/ Input options:
    --rawfastq                    Folder containing paired-end FASTQ files ending with .fastq.gz;
                                containing either "_R1", "_R2" or "_I1" in the filename. (default: null)
or

    --meta_file                   Single file with the location of all input data. Must be formatted as a CSV with (at least) columns: SAMPLE_ID,R1,R2, I1 (default: $basedir/metadata_${runID}.csv"); 
                                additional columns in meta-file: SAMPLE_NAME, PATIENT_ID, SOURCE, TISSUE_TYPE, LIBRARY_TYPE, RUN_ID, LANE_NO, CENTER_NAME, INSTRUMENT_PLATFORM, INSTRUMENT_MODEL, INDEX, ILSE_NO, SEQUENCING_READ_TYPE

    --run_ID                      run_ID of sequencing run (either provided by sequencing center, or manually defined);(default: null)

/ Output options:
    --basedir                     MAIN directory for output files of entire analysis (default: null)
 
/ References:
    --reference                   Directory, storing reference genome to use for alignment, in FASTA format plus already created (Bowtie2) index (default: null)
    
/// Optional Arguments:
/ Output options:
    --preprocessing               MAIN directory for output files of preprocessing pipeline from raw-fastq files to ready-to-use bam-files (default: "$basedir/pre-processing")

/ References:
    --genome-fasta                Reference genome to use for alignment, in FASTA format (default: "$reference/genome.fa")

/ UMI handling and trimming:
    --with_UMI                    Libraries were sequenced with UMIs and UMI-sequence is available in I1.fastq.gz for each sample and has to be extracted plus added to *_{R1,R2}.fastq.gz (default: true)
    --skip_umi_extract            Libraries were not sequenced with UMIs/UMI sequence is not available in I1.fastq.gz (default: false)
  
    --umitools_bc_pattern         Pattern of UMI-Index, in case of NEBnext UDI-UMI adapters the UMI-INDEX contains 11 bases (default: 'NNNNNNNNNNN')
    --umitools_dedup_stats        (default: true)
    --skip_trimming               Skip trimming process and proceed with alignment process with untrimmed files (default: false)
    --save_trimmed                Retain trimmed reads and save as fastq.gz files (default: true)
    --nextseq_trim                Set nextseq_trim parameter for Cutadapt-Trimming to consider sequencing with 2-colour technology "NextSeq or NovaSeq" (default: 20)
    --adapter_forward             Sequences of forward sequencing adapter; default: Illumina adapter (default: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA')
    --adapter_reverse             Sequences of reverse sequencing adapter; default: Illumina adapter (default: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')

/ Alignement BOWTIE2:
    --save_unaligned              Retain and save unaligned/unmapped reads in addition to mapped reads after alignement process (default: false)
    --read_groups                 Assign readgroups to sequencing reads during alignement process (default: null)
    --bam_csi_index               Create .csi file with Samtools Index, instead of .bai file (default: false)
    
/ Samtools:
    --skip_samtools_stats         Skip statistics on SAM-Flags with samtools flagstat, idxstats, stats (default: false)
   
    """.stripIndent()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPROCESSING         } from './workflows/preprocessing'

//
// WORKFLOW: Run main nf-core/rnaseq analysis pipeline
//
workflow WGS_MEDIP_PIPELINE {
    
    PREPROCESSING ()
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
    WGS_MEDIP_PIPELINE ()
}  

/// Report workflow status
workflow.onComplete {
    println "Workflow run $workflow.runName completed at $workflow.complete with status " +
            "${ workflow.success ? 'success' : 'failure' }"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
