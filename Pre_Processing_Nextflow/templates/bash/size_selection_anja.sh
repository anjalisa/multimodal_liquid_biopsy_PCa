#!/bin/bash

# Script for manually evaluating insert size and to perform size selection

# 1. evaluate insert size from mapped bam-files
# 2. evaluate insert size from final bamfile (only UCSC)

# 3. perform size selection of final bam files.  a) 20-150bp b) 20-175bp c) 180-800

# 4. Count total reads (flagstat) of downsampled files 

# 5. Perform downsampling of size selected bamfiles (to number of sample with smallest read count)

# 6. run MEDIP-QC with size selected files)

# 7. run hmmcopy with size selected, downsampled files
#########################
## Load needed modules
module load samtools/1.9

## Variables: Create directories and define important variables for the whole pipeline

MAIN="/omics/odcf/analysis/OE0562_projects/early_detection_prostate/Anja/cfMeDIP_run1_exp9/analysis_data_seq"                   # main directory for pipeline processing

SAMPLES="$(cat ${MAIN}/samples.txt)"                             # define variable which contains sample-ID of rowdata which has to be processed in the pipeline (variable is used for all loops)

MAP="${MAIN}/mapping_UCSC_hg19"                         # directory OUTPUT-files for mapping
FINAL_BAM="${MAIN}/final_bam_UCSC_hg19" 			# directory of final Bam-Files

mkdir ${MAIN}/size_selection			# create folder in which 

SIZE="${MAIN}/size_selection"			# define directory for output of insert size evaluation and size selection

mkdir ${SIZE}/20_150bp
mkdir ${SIZE}/20_175bp
mkdir ${SIZE}/180_800bp


### Evaluate insert size: for inital bam-file (after mapping and initial filtering): sorting has to be performed before evaluation; the final bam-file is already sorted

for i in ${SAMPLES}
do
samtools sort ${MAP}/${i}.filtered.bam -o ${MAP}/${i}.filtered_samtool_sorted.bam
samtools view -f66 ${MAP}/${i}.filtered_samtool_sorted.bam | cut -f9 | awk '{print sqrt($0^2)}' > ${SIZE}/${i}_initial_insertsizes.txt

samtools view -f66 ${FINAL_BAM}/${i}.final_UCSC_hg19.bam | cut -f9 | awk '{print sqrt($0^2)}' > ${SIZE}/${i}_final_insertsizes.txt
done
	
######### Size selection for final bam-file a) 20-150bp   plus evaluation of total reads (flagstat)
for i in ${SAMPLES}
do
samtools view -h ${FINAL_BAM}/${i}.final_UCSC_hg19.bam | awk 'substr($0,1,1)=="@" || ($9>= 20 && $9<=150) || ($9<=-20 && $9>=-150)' | samtools view -b > ${SIZE}/20_150bp/${i}_size_select_150.bam

echo "flagstat for size selected bamfile 20-150 ${i}" >> ${SIZE}/20_150bp/flagstat_size_select_150.txt 
samtools flagstat ${SIZE}/20_150bp/${i}_size_select_150.bam >> ${SIZE}/20_150bp/flagstat_size_select_150.txt
done 


########## Size selection for final bam-file b) 20-175bp   plus evaluation of total reads (flagstat)
for i in ${SAMPLES}
do
samtools view -h ${FINAL_BAM}/${i}.final_UCSC_hg19.bam | awk 'substr($0,1,1)=="@" || ($9>= 20 && $9<=175) || ($9<=-20 && $9>=-175)' | samtools view -b > ${SIZE}/20_175bp/${i}_size_select_175.bam

echo "flagstat for size selected bamfile 20-175 ${i}" >> ${SIZE}/20_175bp/flagstat_size_select_175.txt 
samtools flagstat ${SIZE}/20_175bp/${i}_size_select_175.bam >> ${SIZE}/20_175bp/flagstat_size_select_175.txt
done


######### Size selection for final bam-file  c) 180-800 .  plus evaluation of total reads (flagstat)
for i in ${SAMPLES}
do
samtools view -h ${FINAL_BAM}/${i}.final_UCSC_hg19.bam | awk 'substr($0,1,1)=="@" || ($9>= 180 && $9<=800) || ($9<=-180 && $9>=-800)' | samtools view -b > ${SIZE}/180_800bp/${i}_size_select_800.bam

echo "flagstat for size selected bamfile 180-800 ${i}" >> ${SIZE}/180_800bp/flagstat_size_select_800.txt 
samtools flagstat ${SIZE}/180_800bp/${i}_size_select_800.bam >> ${SIZE}/20_150bp/flagstat_size_select_800.txt
done

