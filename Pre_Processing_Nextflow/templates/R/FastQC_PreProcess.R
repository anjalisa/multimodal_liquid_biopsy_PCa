#
# Script Name: FastQC_PreProcess.R <MAIN> <ID> <TYPE>
#
# Author: Florian Janke
# Last updated: 2022/03/02
#
# Description 
#   (1) Greps and summarizes relevant information from filtered and unfiltered .fastqc files
#
# Run Information
#   (1) Script is called by SeqData_PreProcess_batch.sh
#
# --------------------------------------------------------------------------------------------------
library(ggplot2)
library(data.table)

args <- commandArgs()
MAIN <- args[6]
ID <- args[7]
TYPE <- args[8]

# Remove last backslash #
MAIN <- substr(MAIN, 1, nchar(MAIN)-1)


### Load meta data information ###
meta <- read.table(file.path(MAIN, "tmp", ID, paste0(ID, "_meta.tsv")), sep = "\t")[, c(2)]
meta <- as.character(meta[!duplicated(meta)])

### LOOP to extract and store relevant data (#total paired reads per sample, GC content and per sequence GC) ###
for(i in 1:length(meta)){
  
  ## List fastqc files of the same sample ##
  samp <- list.files(file.path(MAIN, "QC/fastqc", meta[i], TYPE))
  samp <- samp[order(samp)]
  
  for(j in 1:length(samp)){
    
    ## Create temporary directory to store unzipped fastqc file ##
    dir.create(file.path(MAIN, "QC/fastqc/temp"))
    
    ## Unzip fastqc file ##
    unzip(zipfile = file.path(MAIN, "QC/fastqc", meta[i], TYPE, samp[j]), exdir = file.path(MAIN, "QC/fastqc/temp"))
    
    ## Load data and extract relevant information ##
    data <- read.table(file.path(MAIN, "QC/fastqc/temp", list.files(file.path(MAIN, "QC/fastqc/temp")), "fastqc_data.txt"), sep = ";")
    
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
    if(grepl("R1", samp[j])){
      
      R1 <- list[[1]]
      R1 <- R1[which(R1[, 1] %in% c("Total Sequences", "%GC")), ]
      R1[, 2] <- as.numeric(as.character(R1[, 2]))
      
    }
    if(grepl("R2", samp[j])){
      
      R2 <- list[[1]]
      R2 <- R2[which(R2[, 1] %in% c("Total Sequences", "%GC")), ]
      R2[, 2] <- as.numeric(as.character(R2[, 2]))
      
    }
    
    # Extract #total sequences, GC content and per sequence GC of filtered data #
    if(grepl("_filtered", samp[j])){
      
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
    }
    
    ## Remove temporary directory ##
    unlink(file.path(MAIN, "QC/fastqc/temp"), recursive = TRUE)
    
  }
  
  # Create summary table #
  if(i == 1){
    
    table <- data.frame(Sample_ID = meta[i], Paired.reads.unfiltered = ((R1[1, 2] + R2[1, 2]) /2), Paired.reads.filtered = filtered[1, 2] /2, FPR = round((filtered[1, 2] /2) / ((R1[1, 2] + R2[1, 2]) /2), 3), GC.content = filtered[2, 2] /100)
    #GC.table <- GC
    
  } else {
    
    table <- rbind(table, data.frame(Sample_ID = meta[i], Paired.reads.unfiltered = ((R1[1, 2] + R2[1, 2]) /2), Paired.reads.filtered = filtered[1, 2] /2, FPR = round((filtered[1, 2] /2) / ((R1[1, 2] + R2[1, 2]) /2), 3), GC.content = filtered[2, 2] /100))
    #GC.table <- rbind(GC.table, GC)
    
  }
  
  
}


### Export summary tables as .tsv files ###
write.table(table, file = file.path(MAIN, "tmp", ID, paste0("fastqc_summary.tsv")), row.names = FALSE, sep = "\t")
#fwrite(GC.table, file = file.path(MAIN, "tmp", ID, paste0("perSequenceGC_summary.tsv")), row.names = FALSE, sep = "\t")



