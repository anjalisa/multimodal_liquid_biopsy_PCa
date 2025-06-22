#
# Script Name: Saturation_PreProcess.sh <MAIN> <META> <ID> <TYPE>
#
# Author: Florian Janke
# Last updated: 2022/02/28
#
# Description
#   (1) Estimates whether sequence coverage is sufficiently covering the whole genome
#
# Run Information
#   (1) Script is called by SeqData_PreProcess.sh
#
# --------------------------------------------------------------------------------------------------

library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)

## Import path to main directory from shell script ##
args <- commandArgs()
MAIN <- args[6]
META <- args[7]
ID <- args[8]
TYPE <- args[9]

# Remove last backslash #
MAIN <- substr(MAIN, 1, nchar(MAIN)-1)


## Load meta data information ##
meta <- read.table(file.path(META), sep = "\t")[, c(2)]
meta <- as.character(meta[!duplicated(meta)])

## Saturation analysis loop ##
for(i in 1:length(meta)){
  
  # Run saturation analysis #
  sr <- MEDIPS.saturation(file = file.path(MAIN, "bam", TYPE, paste0(meta[i], ".bam")), BSgenome = "BSgenome.Hsapiens.UCSC.hg19", uniq = 0, window_size = 300, paired = TRUE, chr.select = paste0("chr", c(1:22)))
  
  # Store saturation data per sample #
  if(i == 1){
    
    sat.data <- data.frame(Sample_ID = meta[i], saturation.r = sr$maxEstCor[2], read.cov = sr$maxEstCor[1])
    
  } else {
    
    sat.data <- rbind(sat.data, data.frame(Sample_ID = meta[i], saturation.r = sr$maxEstCor[2], read.cov = sr$maxEstCor[1]))
    
  }
  
  # Store saturation data per sample #
  sat <- data.frame(sr$estimation)
  colnames(sat) <- c("read.cov", "saturation.r")
  fwrite(sat, file = file.path(MAIN, "QC/saturation", TYPE, paste0(meta[i], ".txt")))
  
  # Save data #
  write.table(sat.data, file = file.path(MAIN, "tmp", ID, paste0("saturation_summary.tsv")), row.names = FALSE, sep = "\t")
  #write.table(sat.data, file = "/omics/odcf/analysis/OE0309_projects/nsclc_immunotherapy/jankef/Sequencing_data/QC/summary/saturation_summary.tsv", row.names = FALSE, sep = "\t")
}






