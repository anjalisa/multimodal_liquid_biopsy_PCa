#!/bin/bash
#
# Script Name: SeqData_PreProcess.sh <MAIN> <ID> <TYPE>
#
# Author: Florian Janke
# Last updated: 2022/03/15
#
# Description
#   (1) Creates directories to store pre-processed sequencing data and quality control data
#
# Run Information
#   (1) Script is called by XXX.sh
#
# --------------------------------------------------------------------------------------------------

module load R/3.6.2

MAIN=$1
ID=$2
TYPE=$3
PIPELINE="/omics/groups/OE0309/internal/janke/Pipeline/"


## Creates directories to store data ##
${PIPELINE}Scripts/Pre-Processing/Directory_PreProcess.sh ${MAIN} ${TYPE} ${ID}


## Obtains and prepares meta data sheet ##
Rscript ${PIPELINE}Scripts/Pre-Processing/MetaData_PreProcess.R ${MAIN} ${ID} ${TYPE}
META=${MAIN}tmp/${ID}/${ID}_meta.tsv
rm -rf ${MAIN}${ID}_meta.tsv

## Processing ##
# Adapter trimming using Trimmomatic #
${PIPELINE}Scripts/Pre-Processing/Trim_PreProcess.sh ${MAIN} ${META} ${ID}

# Alignment using bowtie2 #
${PIPELINE}Scripts/Pre-Processing/Bowtie2_PreProcess.sh ${MAIN} ${META} ${PIPELINE} ${ID}

# Filtering and deduplication; samtools and markduplicates #
${PIPELINE}Scripts/Pre-Processing/Filtering_PreProcess.sh ${MAIN} ${META} ${ID} ${TYPE}

# Clean-up temporary files #
rm -rf ${MAIN}tmp/${ID}/cutadapt
rm -rf ${MAIN}tmp/${ID}/bowtie2
rm -rf ${MAIN}tmp/${ID}/filtered


## Quality control ##
# Collects .fastqc files of unprocessed sequencing data #
${PIPELINE}Scripts/Pre-Processing/FastQCRAW_PreProcess.sh ${MAIN} ${META} ${TYPE}

# Generates .fastqc files of final .bam files #
${PIPELINE}Scripts/Pre-Processing/FastQC_PreProcess.sh ${MAIN} ${META} ${TYPE}

# Summarizes FastQC results #
Rscript ${PIPELINE}Scripts/Pre-Processing/FastQC_PreProcess.R ${MAIN} ${ID} ${TYPE}

# MeDIP-specific quality metrics #
if [ ${TYPE} = "MeDIP" ]; then

    # Estimates coverage saturation at the given sequencing depth #
    Rscript ${PIPELINE}Scripts/Pre-Processing/Saturation_PreProcess.R ${MAIN} ${META} ${ID} ${TYPE}

    # Calculates CpG enrichment scores per sample #
    Rscript ${PIPELINE}Scripts/Pre-Processing/CpG_enrichment_PreProcess.R ${MAIN} ${META} ${ID} ${TYPE}

    # Calculates CpG coverage scores per sample #
    Rscript ${PIPELINE}Scripts/Pre-Processing/CpG_coverage_Script.R ${MAIN} ${META} ${ID} ${TYPE}

fi

# Summarizes quality data #
Rscript ${PIPELINE}Scripts/Pre-Processing/QC_summary_PrepProcess.R ${MAIN} ${META} ${ID}
