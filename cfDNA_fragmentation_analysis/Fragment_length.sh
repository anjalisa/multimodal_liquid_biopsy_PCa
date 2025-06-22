#!/bin/bash
#
# Script Name: Fragment_length.sh <MAIN> <PIPELINE>
#
# Author: Florian Janke (DKFZ)
# Last updated: 2022/03/03
#
# Description
#   (1) Creates directories to store fragment size information
#   (2) Obtains insert sizes of paired end sequencing data
#   (3) Summarizes insert sizes in frequency table
#
# Input
# final bam-files (mapped, deduplicated sequencing reads)  
# 
# --------------------------------------------------------------------------------------------------

module load samtools/1.9
module load R/4.0.0

MAIN="/path/to/bamfile/directory"
PIPELINE="/path/to/bioinformaticScripts/directory"


## Creates directories to store data ##
if [ ! -d ${MAIN}fragments ]; then
    mkdir ${MAIN}fragments
fi

if [ ! -d ${MAIN}fragments/size ]; then
    mkdir ${MAIN}fragments/size
fi

if [ ! -d ${MAIN}fragments/size/frequency_tables ]; then
    mkdir ${MAIN}fragments/size/frequency_tables
fi


## Defines <OUTPUT> directory ##
OUTPUT=${MAIN}fragments/size/frequency_tables/


## Obtains insert size of paired end sequencing data per sample ##
# Creates meta data sheet #
find ${MAIN}Sequencing_data/bam/ -maxdepth 1 -name '*bam' > ${OUTPUT}meta.tsv
META=${OUTPUT}meta.tsv

# Retrieves insert size #
for i in `seq 1 $(cat ${META} | wc -l)`
do

    # Gets path to sample and file name #
    path=$(awk -v var=$i -F '\t' 'FNR == var {print $1}' ${META})
    name=$(basename $path)
    name=${name::-4}

    # Checks if file was already processed #
    if [ -f ${OUTPUT}${name}.txt ]; then
        continue
    fi
    
    # Insert size #
    samtools view -f 66 -@ 10 ${MAIN}Sequencing_data/bam/${name}.bam | cut -f 9 > ${OUTPUT}${name}.txt

    # Prepares frequency table #
    Rscript ${PIPELINE}Scripts/Fragments/Fragment_length_FreqTable.R ${OUTPUT} ${name}

done

rm -rf ${OUTPUT}meta.tsv
