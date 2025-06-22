#!/bin/bash
#
# Script Name: Bowtie2_PreProcess.sh <MAIN> <META> <PIPELINE> <ID>
#
# Author: Florian Janke
# Last updated: 2022/02/28
#
# Description
#   (1) Aligns .fastq files to hg19
#   (2) Sorts .bam files
#   (3) Counts unfiltered .bam files
#
# Run Information
#   (1) Script is called by SeqData_PreProcess.sh
#
# --------------------------------------------------------------------------------------------------

module load bowtie2/2.3.5.1
module load samtools/1.9
module load sambamba/0.7.1

MAIN=$1
META=$2
PIPELINE=$3
ID=$4


# Alignment loop #
for i in `seq 1 2 $(cat ${META} | wc -l)`
do
    
    # Extract odcf IDs for R1 and R2 read from meta data file and store in variable #
    R1=$(awk -v var=$i -F '\t' 'FNR == var {print $1}' ${META})
    R2=$(awk -v var=$((i+1)) -F '\t' 'FNR == var {print $1}' ${META})

    # Extract user-defined IDs (one per sample) #
    user=$(awk -v var=$i -F '\t' 'FNR == var {print $2}' ${META})

    # Extract ILSE ID from meta_data file #
    ilse=$(awk -v var=$i -F '\t' 'FNR == var {print $3}' ${META})

    # Check if the current file was already processed #
    if [ -f "${MAIN}tmp/${ID}/bowtie2/${user}.bam" ]; then
        continue
    fi

    # Align and sort .bam file #
    bowtie2 -x ${PIPELINE}References/UCSC_hg19/hg19 -1 ${MAIN}tmp/${ID}/cutadapt/${R1} -2 ${MAIN}tmp/${ID}/cutadapt/${R2} --rg-id "${user}" --rg "SM:${ilse}" --rg "PL:ILLUMINA" --rg "LB:${R1::-12}" --phred33 --no-unal --minins 30 --maxins 700 --fr --threads 20 | samtools sort -@ 20 -m 2G -T ${MAIN}tmp/${ID}/bowtie2 -o - - > ${MAIN}tmp/${ID}/bowtie2/${user}.bam
    #bowtie2 -x ${PIPELINE}References/UCSC_hg19/hg19 -1 ${MAIN}tmp/${ID}/cutadapt/${R1} -2 ${MAIN}tmp/${ID}/cutadapt/${R2} --rg-id "${user}" --rg "SM:${ilse}" --rg "PL:ILLUMINA" --rg "LB:${R1::-12}" --phred33 --no-unal --minins 30 --maxins 700 --fr --threads 15 > ${MAIN}tmp/${ID}/bowtie2/${user}.bam
    #sambamba sort -t 30 -m 45G -o ${MAIN}tmp/${ID}/bowtie2/${user}.bam ${MAIN}tmp/${ID}/bowtie2/${user}.bam

done
