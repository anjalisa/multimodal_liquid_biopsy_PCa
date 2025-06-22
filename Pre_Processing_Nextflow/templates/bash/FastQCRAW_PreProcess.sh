#!/bin/bash
#
# Script Name: FastQCRAW_PreProcess.sh <MAIN> <ID> <TYPE>
#
# Author: Florian Janke
# Last updated: 2022/05/03
#
# Description
#   (1) Copies FASTQC files of raw sequencing data into <MAIN>/QC/fastqc/<TYPE>
#   (2) In case fastqc files do not exists they are created from scratch
#   (3) Renames files by user-defined ID and appends _unfiltered_fastqc.zip
#
# Run Information
#   (1) Script is called by SeqData_PreProcess.sh
#
# --------------------------------------------------------------------------------------------------

module load fastqc/0.11.5

MAIN=$1
META=$2
TYPE=$3


## LOOP ##
for ((i=1; i <= $(cat ${META} | wc -l); i++))
do

    # Obtains path to run data; first tries midterm and then odcf server #
    run=$(awk -v var=$i -F '\t' 'FNR == var {print $5}' ${META})
    path="CLUSTER"
    if [ ! -d ${run} ]; then
        run=$(awk -v var=$i -F '\t' 'FNR == var {print $6}' ${META})
        path="MIDTERM"
    fi

    if [ ! -d ${run} ]; then
        run=$(awk -v var=$i -F '\t' 'FNR == var {print $7}' ${META})
        path="CREATE_NEW"
    fi

    # Extract the i-th odcf and user-defined ID from the meta data .tsv file #
    odcf=$(awk -v var=$i -F '\t' 'FNR == var {print $1}' ${META})
    user=$(awk -v var=$i -F '\t' 'FNR == var {print $2}' ${META})

    # Get R1 or R2 read information #
    read=${odcf#*_}
    read=${read::-9}

    # Check if file was already copied #
    if [ -f "${MAIN}QC/fastqc/${user}/${TYPE}/${user}_${read}_unfiltered_fastqc.zip" ]; then
        continue
    fi

    # Creates directory to save files in #
    if [ ! -d "${MAIN}QC/fastqc/${user}" ]; then
        mkdir ${MAIN}QC/fastqc/${user}
    fi

    # Creates TYPE directory #
    if [ ! -d "${MAIN}QC/fastqc/${user}/${TYPE}" ]; then
        mkdir ${MAIN}QC/fastqc/${user}/${TYPE}
    fi

    # Locate fastqc.zip file of corresponding sample (get path) #
    if [ "${path}" = "CLUSTER" ]; then

        fastqc=$(dirname $(dirname ${run}))/view-by-pid
        fastqc=$(find $fastqc -name '*fastqc.zip' -a -name ${odcf::-9}*)

        # Copy fastqc.zip file to /QC/fastqc/ directory #
        cp -rf ${fastqc} ${MAIN}QC/fastqc/${user}/${TYPE}

        # Rename fastqc.zip file by the user-defined ID #
        mv ${MAIN}QC/fastqc/${user}/${TYPE}/${odcf::-9}_fastqc.zip ${MAIN}QC/fastqc/${user}/${TYPE}/${user}_${read}_unfiltered_fastqc.zip

    fi

    if [ "${path}" = "MIDTERM" ]; then
        fastqc=$(find $run -name '*fastqc.zip' -a -name ${odcf::-9}*)

        # Copy fastqc.zip file to /QC/fastqc/ directory #
        cp -rf ${fastqc} ${MAIN}QC/fastqc/${user}/${TYPE}

        # Rename fastqc.zip file by the user-defined ID #
        mv ${MAIN}QC/fastqc/${user}/${TYPE}/${odcf::-9}_fastqc.zip ${MAIN}QC/fastqc/${user}/${TYPE}/${user}_${read}_unfiltered_fastqc.zip

    fi

    # Create new fastqc.zip files in case they are not present on the cluster or midterm server #
    if [ "${path}" = "CREATE_NEW" ]; then

        # Run FASTQC #
        fastqc ${run}${odcf} -t 20 -o ${fastqc} ${MAIN}QC/fastqc/${user}/${TYPE}

        # Rename fastqc.zip file by the user-defined ID #
        mv ${MAIN}QC/fastqc/${user}/${TYPE}/${odcf::-9}_fastqc.zip ${MAIN}QC/fastqc/${user}/${TYPE}/${user}_${read}_unfiltered_fastqc.zip

        # Remove .html file #
        rm ${MAIN}QC/fastqc/${user}/${TYPE}/${odcf::-9}_fastqc.html

    fi


done