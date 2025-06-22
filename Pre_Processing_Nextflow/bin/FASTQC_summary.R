#
# Script Name: FastQC_PreProcess.R <MAIN> <ID> <TYPE>
#
# Author: Florian Janke ; Anja Riediger (nextflow version)
# Last updated: 2022/03/02 ; 2022/08/26 (nextflow version)
#
# Description 
#   (1) Greps and summarizes relevant information from filtered and unfiltered .fastqc files
#
# --------------------------------------------------------------------------------------------------

#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
fastqc_zip <- read.table(args[1], header=FALSE, col.names="Samples") 

## Load meta data information ##
fastqc <- as.character(fastqc_zip$Samples)

### LOOP to extract and store relevant data (#total paired reads per sample, GC content and per sequence GC) ###
for(i in 1:length(fastqc)){
  
  ## List fastqc files of the same sample ##
  samp <- list.files(file.path(fastqc[i]))
  samp <- samp[order(samp)]
  
  for(j in 1:length(samp)){
    
    ## Create temporary directory to store unzipped fastqc file ##
    dir.create(file.path(${workdir}, "/QC/fastqc/temp"))
    
    ## Unzip fastqc file ##
    unzip(zipfile = samp[j], exdir = file.path(${workdir}, "/QC/fastqc/temp"))
    
    ## Load data and extract relevant information ##
    data <- read.table(file.path(${workdir}, "/QC/fastqc/temp/", "fastqc_data.txt"), sep = ";")
    
    # Split by >>END_MODULE #
    end <- which(data[, 1] == ">>END_MODULE")
    start <- c(1, end[1:(length(end)-1)] +1)
    
    # Extract 1st and 6th block into list #
    list <- list()
    for(n in c(1, 6)){
      
      tmp <- data.frame(V1 = data[c(start[n]:end[n]), ])
      tmp <- data.frame(tmp[c(1:(nrow(tmp)-1)), ])
      
      tmp <- data.frame(do.call("rbind", strsplit(as.character(tmp[,1]), "\t", fixed = TRUE)))
      colnames(tmp) <- as.character(tmp[1, ])
      tmp <- tmp[-1, ]
      
      list[colnames(tmp)[1]] <- list(tmp)
      
    }
    
    # Extract #total sequences and GC content of unfiltered data #
    #if(grepl("R1", samp[j])){
      
    #  R1 <- list[[1]]
    #  R1 <- R1[which(R1[, 1] %in% c("Total Sequences", "%GC")), ]
    #  R1[, 2] <- as.numeric(as.character(R1[, 2]))
      
    #}
    #if(grepl("R2", samp[j])){
      
    #  R2 <- list[[1]]
    #  R2 <- R2[which(R2[, 1] %in% c("Total Sequences", "%GC")), ]
    #  R2[, 2] <- as.numeric(as.character(R2[, 2]))
      
    #}
    
    # Extract #total sequences, GC content and per sequence GC of filtered data #
    #if(grepl("final", samp[j])){
      
      # Total sequences and GC content #
      filtered <- list[[1]]
      filtered <- filtered[which(filtered[, 1] %in% c("Total Sequences", "%GC")), ]
      filtered[, 2] <- as.numeric(as.character(filtered[, 2]))
      
      # Per sequence GC #
      #GC <- list[[2]]
      #colnames(GC) <- c("GC.content", "reads")
      #GC$GC.content <- c(0:100)
      #GC$reads <- as.numeric(as.character(GC$reads))
      #GC$Sample_ID <- meta[i]
    #}
    
    ## Remove temporary directory ##
    unlink(file.path(MAIN, "QC/fastqc/temp"), recursive = TRUE)
    
  }
  
  # Create summary table #
  if(i == 1){
    
    #table <- data.frame(Sample_ID = meta[i], Paired.reads.unfiltered = ((R1[1, 2] + R2[1, 2]) /2), Paired.reads.filtered = filtered[1, 2] /2, FPR = round((filtered[1, 2] /2) / ((R1[1, 2] + R2[1, 2]) /2), 3), GC.content = filtered[2, 2] /100)
    table <- data.frame(Sample_ID = fastqc[i], Paired.reads.filtered = filtered[1, 2] /2, GC.content = filtered[2, 2] /100)
    #GC.table <- GC
    
  } else {
    
    #table <- rbind(table, data.frame(Sample_ID = meta[i], Paired.reads.unfiltered = ((R1[1, 2] + R2[1, 2]) /2), Paired.reads.filtered = filtered[1, 2] /2, FPR = round((filtered[1, 2] /2) / ((R1[1, 2] + R2[1, 2]) /2), 3), GC.content = filtered[2, 2] /100))
    table <- rbind(table, data.frame(Sample_ID = fastqc[i], Paired.reads.filtered = filtered[1, 2] /2, GC.content = filtered[2, 2] /100))
    #GC.table <- rbind(GC.table, GC)
    
  }
  
  
}


### Export summary tables as .tsv files ###
write.table(table, file = "fastqc_summary.tsv", row.names = FALSE, sep = "\t")
#fwrite(GC.table, file = file.path(MAIN, "tmp", ID, paste0("perSequenceGC_summary.tsv")), row.names = FALSE, sep = "\t")



