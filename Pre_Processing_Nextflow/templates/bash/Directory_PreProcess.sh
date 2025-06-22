#!/bin/bash
#
# Script Name: Directory_PreProcess.sh <MAIN> <TYPE> <ID>
#
# Author: Florian Janke
# Last updated: 2022/04/29
#
# Description
#   (1) Creates directories to store pre-processed sequencing data and quality control data
#
# Run Information
#   (1) Script is called by SeqData_PreProcess.sh
#
# --------------------------------------------------------------------------------------------------


MAIN=$1
TYPE=$2
ID=$3


## Creates directories to store ... ##
# ... quality control data #
if [ ! -d ${MAIN}QC ]; then
    mkdir ${MAIN}QC
    mkdir ${MAIN}QC/fastqc
    mkdir ${MAIN}QC/cutadapt
    mkdir ${MAIN}QC/dedup_metrics
    if [ ${TYPE} = "MeDIP" ]; then
        mkdir ${MAIN}QC/saturation
        mkdir ${MAIN}QC/saturation/${TYPE}
    fi
    mkdir ${MAIN}QC/summary
    mkdir ${MAIN}QC/summary/figures
fi


# ... final .bam files #
if [ ! -d ${MAIN}bam ]; then
    mkdir ${MAIN}bam
fi

if [ ! -d ${MAIN}bam/${TYPE} ]; then
    mkdir ${MAIN}bam/${TYPE}
fi

# ... meta data sheet for longterm storage #
if [ ! -d ${MAIN}meta_data ]; then
    mkdir ${MAIN}meta_data
fi

# ... temporary data #
if [ ! -d ${MAIN}tmp ]; then
    mkdir ${MAIN}tmp
fi

if [ ! -d ${MAIN}tmp/${ID} ]; then
    mkdir ${MAIN}tmp/${ID}

fi

if [ ! -d ${MAIN}tmp/${ID}/cutadapt ]; then
    mkdir ${MAIN}tmp/${ID}/cutadapt
fi

if [ ! -d ${MAIN}tmp/${ID}/bowtie2 ]; then
    mkdir ${MAIN}tmp/${ID}/bowtie2
fi

if [ ! -d ${MAIN}tmp/${ID}/filtered ]; then
    mkdir ${MAIN}tmp/${ID}/filtered
fi

# ... count matrix #
if [ ${TYPE} = "MeDIP" ]; then

    if [ ! -d ${MAIN}featureCounts ]; then
        mkdir ${MAIN}featureCounts
        mkdir ${MAIN}featureCounts/raw
    fi

fi



