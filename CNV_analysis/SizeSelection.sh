#!/bin/bash
#
# Script Name: SizeSelection.sh <INPUT> <OUTPUT> <LOW> <HIGH>
#
# Author: Florian Janke (Anja Riediger, adapted)
# Last updated: 2022/03/04 (adapted: 2023/01/26)
#
# Description
#   (1) Subsets .bam file by fragment length
#   (2) Upper and lower limits are set by <HIGH> and <LOW>
#
# Run Information
#   (1) 
#
# Input
# mapped, deduplicated, filtered Bam-Files
# --------------------------------------------------------------------------------------------------

module load samtools/1.9

BASEDIR="/path/to/base/directory/"

IN="/input/" # directory with bam-files for which size selection should be applied
OUT="/output/" # directory for output (size-selected) bam-files

LOW=$2 # lower size limit
HIGH1=$3 # upper size limit

FOLDERID=()  ## Folder names, in case multiple folders are processed

LOW=$1
HIGH1=$2
#HIGH1=(110 120 130 140 150 160 170)

##########

for folder in ${FOLDERID[@]}; do

    INPUT=("${BASEDIR}$folder${IN}")
    OUTPUT=("${BASEDIR}$folder${OUT}")

    ## Creates meta data sheet ##
    #find ${INPUT} -name '*.bam' > ${OUTPUT}meta_sizeselection.tsv
    #META=${OUTPUT}meta_sizeselection.tsv
    find ${INPUT} -name '*U*.bam' > ${OUTPUT}meta_sizeselection.tsv
    META=${OUTPUT}meta_sizeselection.tsv

    ################################# Perform size selection for LOW/HIGH1 #################################

    ## Creates directory to store size selected .bam files ##
    if [ ! -d ${OUTPUT}${LOW}_to_${HIGH1}bp ]; then
        mkdir ${OUTPUT}${LOW}_to_${HIGH1}bp
    fi

    OUTPUT1="${OUTPUT}${LOW}_to_${HIGH1}bp/"

    ## Subsets .bam files by fragment sizes within the limits set by <LOW> <HIGH> ##
    for i in `seq 1 $(cat ${META} | wc -l)`
    do

        # Saves sample path and name in variable #
        path=$(awk -v var=$i -F '\t' 'FNR == var {print $1}' ${META})
        name=$(basename $path)
        name=${name::-4}
        reducedname=${name::-6}

        # Checks if .bam file was already size-selected #
        if [ -f ${OUTPUT1}${reducedname}_${LOW}_${HIGH1}.bam ]; then
            continue
        fi

        # Subsets by fragment size #
        samtools view -h ${path} -@ 10 | \
            awk -v a="${LOW}" -v b="${HIGH1}" 'substr($0,1,1)=="@" || ($9>=a && $9<=b) || ($9<=-a && $9>=-b)' | \
            samtools view -b -@ 10 > ${OUTPUT1}${reducedname}_${LOW}_${HIGH1}.bam

        # Sort size-selected .bam file (coordinate-sorting)
        samtools sort ${OUTPUT1}${reducedname}_${LOW}_${HIGH1}.bam -o ${OUTPUT1}${reducedname}_${LOW}_${HIGH1}.bam -@ 10

        # Index sorted .bam file
        samtools index ${OUTPUT1}${reducedname}_${LOW}_${HIGH1}.bam -@ 10

    done

    ## Counts number of fragments per file within specified size range
    # Creates .tsv file to populate #
   touch $(dirname "${OUTPUT1}")/ReadCount_${LOW}_to_${HIGH1}bp.tsv

    # Populates .tsv file with coverage of size selected .bam files #
    for i in `seq 1 $(cat ${META} | wc -l)`
    do

        # Saves sample path and name in variable #
       path=$(awk -v var=$i -F '\t' 'FNR == var {print $1}' ${META})
        name=$(basename $path)
        name=${name::-4}
        reducedname=${name::-6}

        # Obtains coverage of size selected .bam #
        coverage=$(samtools view -c -@ 10 -q 10 -F 4 -F 2048 -F 256 -F 1024 ${OUTPUT1}${reducedname}_${LOW}_${HIGH1}.bam)

        # Save in .tsv files #
        echo -e "$reducedname \t $coverage" >> $(dirname "${OUTPUT1}")/ReadCount_${LOW}_to_${HIGH1}bp.tsv

    done


done