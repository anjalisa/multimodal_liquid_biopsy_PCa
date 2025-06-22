#
# Script Name: CpG_coverage_PreProcess.sh <MAIN> <META> <ID> <TYPE>
#
# Author: Florian Janke
# Last updated: 2022/02/28
#
# Description
#   (1) Calculates CpG enrichment scores per sample
#
# Run Information
#   (1) Script is called by SeqData_PreProcess.sh
#
# --------------------------------------------------------------------------------------------------

library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(dplyr)

args <- commandArgs()
MAIN <- args[6]
META <- args[7]
ID <- as.numeric(args[8])
TYPE <- args[9]

# Remove last backslash #
MAIN <- substr(MAIN, 1, nchar(MAIN)-1)


## Load meta data information ##
meta <- read.table(file.path(META), sep = "\t")[, c(2)]
meta <- as.character(meta[!duplicated(meta)])


## CpG coverage analysis loop ##
for(i in 1:length(meta)){
  
  # Run coverage analysis #
  cr <- MEDIPS.seqCoverage(file = file.path(MAIN, "bam", TYPE, paste0(meta[i], ".bam")), pattern = "CG", BSgenome = "BSgenome.Hsapiens.UCSC.hg19", chr.select = paste0("chr", c(1:22)), uniq = 0, paired = TRUE)
  
  # Summarize coverages per CpG ##
  cov.data <- data.frame(coverage = cr$cov.res, count = 1)
  
  # Calculate relative number of CpGs not covered by the sequencing data #
  cov <- 1-(nrow(cov.data[cov.data$coverage == 0, ]) / nrow(cov.data))
  
  # Calculate relative number of CpGs covered by more than 5 reads # 
  cov.data <- cov.data %>% dplyr::group_by(coverage) %>% dplyr::summarize(count = sum(count))
  cov.data.greater.5 <- sum(cov.data[c(which(cov.data$coverage == 6):nrow(cov.data)), ]$count) / sum(cov.data$count)
  
  # Calculate relative number of reads covering no CpG #
  cov.data.0 <- cr$numberReadsWO / (cr$numberReads)
  
  
  if(i == 1){
    
    table <- data.frame(Sample_ID = meta[i], CpGs.covered = cov, CpGs.greater.5 = cov.data.greater.5, reads.WO.CPG = cov.data.0)
    
  } else {
    
    table <- rbind(table, data.frame(Sample_ID = meta[i], CpGs.covered = cov, CpGs.greater.5 = cov.data.greater.5, reads.WO.CPG = cov.data.0))
    
  }
  
  # Save data #
  #write.table(table, file = file.path(MAIN, "tmp", ID, paste0("CpG_coverage_summary.tsv")), row.names = FALSE, sep = "\t")
  write.table(table, file = "/omics/odcf/analysis/OE0309_projects/nsclc_immunotherapy/jankef/Sequencing_data/QC/summary/CpG_coverage_summary.tsv", row.names = FALSE, sep = "\t")
  
}



