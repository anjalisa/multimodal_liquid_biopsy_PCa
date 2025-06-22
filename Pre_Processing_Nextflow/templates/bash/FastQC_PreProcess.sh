#!/bin/bash
#
# Script Name: Bowtie2_PreProcess.sh <MAIN> <META> <PIPELINE> <ID>
#
# Author: Florian Janke
# Last updated: 2022/02/28
#
# Description
#   (1) Creates fastqc.zip files of filtered .bam files
#
# Run Information
#   (1) Script is called by SeqData_PreProcess.sh
#
# --------------------------------------------------------------------------------------------------

module load fastqc/0.11.5

MAIN=$1
META=$2
TYPE=$3

## FASTQC of final .bam files ##
for i in $(awk -F '\t' '{print $2}' ${META} | uniq)
do

    # Create file-specific directory if it does not exist yet #
    if [ ! -d "${MAIN}QC/fastqc/${i}" ]; then
        mkdir ${MAIN}QC/fastqc/${i}
        mkdir ${MAIN}QC/fastqc/${i}/${TYPE}
    fi

    # Check if fastqc was already run for this file #
    if [ -f "${MAIN}QC/fastqc/${i}/${TYPE}/${i}_filtered_fastqc.zip" ]; then
        continue
    fi

    # Run FASTQC #
    fastqc ${MAIN}bam/${TYPE}/${i}.bam -t 20 -o ${MAIN}QC/fastqc/${i}/${TYPE}

    # Remove .html file #
    rm ${MAIN}QC/fastqc/${i}/${TYPE}/${i}_fastqc.html

    # Rename .fastqc.zip file #
    mv ${MAIN}QC/fastqc/${i}/${TYPE}/${i}_fastqc.zip ${MAIN}QC/fastqc/${i}/${TYPE}/${i}_filtered_fastqc.zip

done