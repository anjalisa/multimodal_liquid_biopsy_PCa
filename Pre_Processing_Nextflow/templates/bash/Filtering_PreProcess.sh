#!/bin/bash
#
# Script Name: Filtering_PreProcess.sh <MAIN> <ID> <TYPE>
#
# Author: Florian Janke
# Last updated: 2022/02/28
#
# Description
#   (1) Filters .bam files (SAMTOOLS)
#   (2) Deduplicates .bam files (MARKDUPLICATES)
#   (3) Generates index file .bai
#   (4) Gets read coverage of filtered .bam files
#
# Run Information
#   (1) Script is called by SeqData_PreProcess.sh
#
# --------------------------------------------------------------------------------------------------

module load samtools/1.9
module load picard/2.25.1

MAIN=$1
META=$2
ID=$3
TYPE=$4


## .bam file filtering (SAMTOOLS), deduplication (MARKDUPLICATES) and .bai file generation ##
# Create file to store read coverage #
touch ${MAIN}tmp/${ID}/readCount_samtools.tsv

# Loop #
for i in $(awk -F '\t' '{print $2}' ${META} | uniq)
do

    # Check if filtering and deduplication was already run for this file (checks if final .bam file is already existing) #
    if [ -f "${MAIN}bam/${TYPE}/${i}.bam" ]; then
        continue
    fi

    # Filter .bam files using samtools view #
    samtools view -@ 10 -q 10 -F 12 -f 3 -b -o ${MAIN}tmp/${ID}/filtered/${i}.bam ${MAIN}tmp/${ID}/bowtie2/${i}.bam

    # Deduplicate filtered .bam files to create final .bams #
    java -jar /software/picard/2.25.1/bin/picard.jar MarkDuplicates -I ${MAIN}tmp/${ID}/filtered/${i}.bam -O ${MAIN}bam/${TYPE}/${i}.bam -M ${MAIN}QC/dedup_metrics/${i}_dedup_metrics.txt -REMOVE_DUPLICATES true

    # .bai file generation #
    samtools index ${MAIN}bam/${TYPE}/${i}.bam

    # Count reads #
    coverage=$(samtools view -c -@ 10 -q 20 -F 4 -F 2048 -F 256 -F 1024 ${MAIN}bam/${TYPE}/${i}.bam)
    echo -e "$i \t $coverage" >> ${MAIN}tmp/${ID}/readCount_samtools.tsv

done

