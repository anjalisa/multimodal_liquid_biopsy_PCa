#
# Author: Florian Janke ; Anja Riediger (Nextflow version)
# Last updated: 2022/02/28 ; 2022/07/08 (Nextflow version)
#
# Description
#   (1) Calculates CpG enrichment scores per sample
#
# --------------------------------------------------------------------------------------------------

#!/usr/bin/env Rscript

# Install R packages
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

if (!require("dplyr")){
    install.packages("dplyr", dependencies=TRUE, repos='http://cloud.r-project.org/')
}

## Load R packages
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(data.table)

## Define variables
args = commandArgs(trailingOnly=TRUE)

BAMFILES <- read.table(args[1], header=FALSE, col.names="Samples")       # Textfile with input files (= final BAM-files)
paired_single <- args[2]       # information of whether paired-end or single-end sequencing was performed

## Load meta data information ##
bamfiles <- as.character(BAMFILES$Samples)

# Check whether paired-end or single-end sequencing was performed
paired_end <- if(paired_single=='PAIRED') TRUE else FALSE

## CpG coverage analysis loop ##
for(bam in bamfiles){

  # Define bamname
  bamname <- gsub("*.final.bam", "", bam)
  
  # Run coverage analysis #
  cr <- MEDIPS.seqCoverage(file = bam, pattern = "CG", BSgenome = "BSgenome.Hsapiens.UCSC.hg19", chr.select = paste0("chr", c(1:22)), uniq = 0, paired = paired_end)
  
  # Summarize coverages per CpG ##
  cov.data <- data.frame(coverage = cr$cov.res, count = 1)
  
  # Calculate relative number of CpGs not covered by the sequencing data #
  cov <- 1-(nrow(cov.data[cov.data$coverage == 0, ]) / nrow(cov.data))
  
  # Calculate relative number of CpGs covered by more than 5 reads # 
  cov.data <- cov.data %>% dplyr::group_by(coverage) %>% dplyr::summarize(count = sum(count))
  cov.data.greater.5 <- sum(cov.data[c(which(cov.data$coverage == 6):nrow(cov.data)), ]$count) / sum(cov.data$count)
  
  # Calculate relative number of reads covering no CpG #
  cov.data.0 <- cr$numberReadsWO / (cr$numberReads)
  
  x <- 2:length(bamfiles)
  table <- data.frame(Sample_ID = bamname, CpGs.covered = cov, CpGs.greater.5 = cov.data.greater.5, reads.WO.CPG = cov.data.0) 
  
  for(i in seq(along=x)){
    y <- i + 1  
    table <- rbind(table, data.frame(Sample_ID = bamname, CpGs.covered = cov, CpGs.greater.5 = cov.data.greater.5, reads.WO.CPG = cov.data.0))  
  }
  #if(length(bamfiles) == 1){ 
  #table <- data.frame(Sample_ID = bamname, CpGs.covered = cov, CpGs.greater.5 = cov.data.greater.5, reads.WO.CPG = cov.data.0) 
  #} else {  
  #  table <- rbind(table, data.frame(Sample_ID = bamname, CpGs.covered = cov, CpGs.greater.5 = cov.data.greater.5, reads.WO.CPG = cov.data.0))  
 #}
  
  # Save data #
  write.table(table, file = "CpG_coverage_summary.tsv", row.names = FALSE, sep = "\t")
  
  # Create and save plots for each sample
  # Plot Coverage as pie; pie chart illustrates fraction of CpGs at coverage level (cov.level), histogram show overall coverage level
  jpeg(filename=paste0("Seq_coverage_pattern_pie_", bamname, ".jpg"), width = 650, height = 650)
  MEDIPS.plotSeqCoverage(seqCoverageObj=cr, type="pie", cov.level=c(0,1,2,3,4,5), main=paste0("Sequence_pattern_coverage_pie",":", bamname))
  dev.off()

}


