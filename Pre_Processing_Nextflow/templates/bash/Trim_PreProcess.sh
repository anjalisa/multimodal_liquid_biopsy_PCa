#!/bin/bash
#
# Script Name: Trim_PreProcess.sh <MAIN> <META> <ID>
#
# Author: Florian Janke
# Last updated: 2022/04/29
#
# Description
#   (1) Trims adapters from raw .fastq files
#
# Run Information
#   (1) Script is called by SeqData_PreProcess.sh
#
# --------------------------------------------------------------------------------------------------

module load trimmomatic/0.38

MAIN=$1
META=$2
ID=$3

### Adapter trimming loop ###
for i in `seq 1 2 $(cat ${META} | wc -l)`
do

    # Extract odcf IDs for R1 and R2 read from meta data file and store in variable #
    R1=$(awk -v var=$i -F '\t' 'FNR == var {print $1}' ${META})
    R2=$(awk -v var=$((i+1)) -F '\t' 'FNR == var {print $1}' ${META})

    # Extract user-defined IDs (one per sample) #
    user=$(awk -v var=$i -F '\t' 'FNR == var {print $2}' ${META})

    # Obtains path to run data; first tries midterm, then odcf server and finally a user provided path ('CUSTOM_PATH'; if present) #
    run=$(awk -v var=$i -F '\t' 'FNR == var {print $5}' ${META})

    if [ ! -d ${run} ]; then
        run=$(awk -v var=$i -F '\t' 'FNR == var {print $6}' ${META})
    fi

    if [ ! -d ${run} ]; then
        run=$(awk -v var=$i -F '\t' 'FNR == var {print $7}' ${META})
    fi

    # Check if the current file was already processed #
    if [ -f "${MAIN}tmp/${ID}/cutadapt/${R1}" ]; then
        continue
    fi

    # Create trimmomatic stat file to store statistics summary for each sample #
    #touch ${MAIN}QC/trimmomatic/${user}_stat.txt

    # Create empty log .txt files (log file is removed after the trimming) #
    #touch ${MAIN}tmp/${ID}/trimmomatic/${ID}_logfile.txt

    # Run Cutadapt #
    /home/jankef/.local/bin/cutadapt --trim-n --nextseq-trim=20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ${MAIN}tmp/${ID}/cutadapt/${R1} -p ${MAIN}tmp/${ID}/cutadapt/${R2} --overlap=1 --times=2 --error-rate=0.1 -j 10 --json=${MAIN}QC/cutadapt/${user}_stat.json ${run}${R1} ${run}${R2}

    #/software/trimmomatic/0.38/bin/trimmomatic.sh PE -threads 20 -phred33 -trimlog ${MAIN}tmp/${ID}/trimmomatic/${ID}_logfile.txt -summary ${MAIN}QC/trimmomatic/${user}_stat.txt ${run}${R1} ${run}${R2} ${MAIN}tmp/${ID}/trimmomatic/paired/${R1} ${MAIN}tmp/${ID}/trimmomatic/unpaired/${R1} ${MAIN}tmp/${ID}/trimmomatic/paired/${R2} ${MAIN}tmp/${ID}/trimmomatic/unpaired/${R2} ILLUMINACLIP:/software/trimmomatic/0.38/adapters/TruSeq3-PE-2.fa:2:30:10:3:keepBothReads LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:30

    # Remove log file #
    #rm ${MAIN}tmp/${ID}/trimmomatic/${ID}_logfile.txt

done