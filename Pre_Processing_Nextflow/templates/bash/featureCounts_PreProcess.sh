#!/bin/bash
#
# Script Name: featureCounts_PreProcess.sh <MAIN>
#
# Author: Florian Janke
# Last updated: 2022/05/18
#
# Description
#   (1) Counts features at genomic regions specified in provided .gtf file
#
# Run Information
#   (1) Script is called by SeqData_PreProcess.sh
#
# --------------------------------------------------------------------------------------------------

module load subread/1.5.3

MAIN=$1
PIPELINE="/omics/groups/OE0309/internal/janke/Pipeline/"

## Defines variables ##
INPUT=${MAIN}bam/MeDIP/
OUTPUT=${MAIN}featureCounts/


## Creates meta data file containing paths to all samples ##
find ${INPUT} -maxdepth 1 -name '*.bam' > ${OUTPUT}meta.tsv
META=${OUTPUT}meta.tsv


## Counts reads at specified .gtf file using featureCounts ##
for ((i=1; i <= $(cat ${META} | wc -l); i++))
do

    # Obtains path to and name of .bam file #
    path=$(awk -v var=$i -F '\t' 'FNR == var {print $1}' ${META})
    name=$(basename $path)
    name=${name::-4}

    # Check if featureCounts was already run for this file #
    if [ -f "${MAIN}featureCounts/${name}.txt" ]; then
        continue
    fi

    # Count at methylation blocks (Loyfer 2022) #
    featureCounts -a ${PIPELINE}References/Reference_methBlocks_hg19.gtf -o ${OUTPUT}${name}.txt ${path} -F GTF -g bin -t sequence_feature --largestOverlap -O -p -d 30 -D 700 -T 16

    # Remove .txt.summary file #
    rm ${OUTPUT}${name}.txt.summary

done


## Clean-up ##
rm -rf ${META}