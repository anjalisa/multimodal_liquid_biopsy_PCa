
#
# Author: Florian Janke ; Anja Riediger (Nextflow version)
# Last updated: 2022/02/28 ; 2022/26/08 (Nextflow version)
#
# Description
#   (1) Estimates whether sequence coverage is sufficiently covering the whole genome
#
# --------------------------------------------------------------------------------------------------

#!/usr/bin/env Rscript

# Install R packages, if not loaded via container
if (!require("BiocManager", quietly = TRUE)){
     install.packages("BiocManager", dependencies=TRUE, repos='http://cloud.r-project.org/')
}

if (!require("BSgenome.Hsapiens.UCSC.hg19")){
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
}

if (!require("MEDIPS")){
   BiocManager::install("MEDIPS")
}

if (!require("data.table")){
    install.packages("data.table", dependencies=TRUE, repos='http://cloud.r-project.org/')
}

# Load R packages
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)

# Define variables
args = commandArgs(trailingOnly=TRUE)

BAMFILES <- read.table(args[1], header=FALSE, col.names="Samples")       # Textfile with input files (= final BAM-files)
paired_single <- args[2]       # information of whether paired-end or single-end sequencing was performed
window_size <- args[3]

window_size <- as.numeric(window_size)

## Load meta data information ##
bamfiles <- as.character(BAMFILES$Samples)

# Check whether paired-end or single-end sequencing was performed
paired_end <- if(paired_single=='PAIRED') TRUE else FALSE

## Saturation analysis loop ##
for(bam in bamfiles){

  # Define bamname
  bamname <- gsub("*.final.bam", "", bam)
  
  # Run saturation analysis #
  sr <- MEDIPS.saturation(file = bam, BSgenome = "BSgenome.Hsapiens.UCSC.hg19", uniq = 0, window_size = window_size, paired = paired_end, chr.select = paste0("chr", c(1:22)))
  
  # Store saturation data per sample #
  x <- 2:length(bamfiles)
  sat.data <- data.frame(Sample_ID = bamname, saturation.r = sr$maxEstCor[2], read.cov = sr$maxEstCor[1])
  
  for(i in seq(along=x)){
    y <- i + 1 
    sat.data <- rbind(sat.data, data.frame(Sample_ID = bamname, saturation.r = sr$maxEstCor[2], read.cov = sr$maxEstCor[1])) 
  }
  
  #if(length(bamfiles) == 1){  
  # sat.data <- data.frame(Sample_ID = bamname, saturation.r = sr$maxEstCor[2], read.cov = sr$maxEstCor[1])  
  #}else{ 
  #  sat.data <- rbind(sat.data, data.frame(Sample_ID = bamname, saturation.r = sr$maxEstCor[2], read.cov = sr$maxEstCor[1])) 
  #}
  
  # Store saturation data per sample #
  sat <- data.frame(sr$estimation)
  colnames(sat) <- c("read.cov", "saturation.r")
  data.table::fwrite(sat, file = paste0("saturation_per_sample_", bamname, ".txt"))
  
  # Save data #
  write.table(sat.data, file = "saturation_summary.tsv", row.names = FALSE, sep = "\t")
}






