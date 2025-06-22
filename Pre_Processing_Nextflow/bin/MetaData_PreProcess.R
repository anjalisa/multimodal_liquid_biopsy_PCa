#
# Script Name: MetaData_PreProcess.sh <MAIN> <ID> <TYPE>
#
# Author: Florian Janke
# Last updated: 2022/04/29
#
# Description 
#   (1) Searches for meta data sheet in <MAIN> or <MIDTERM>
#   (2) Copies meta data sheet to <MAIN>/meta_data/ for longterm storage
#   (3) Prepares meta data sheet for looping
#
# Run Information
#   (1) Script is called by SeqData_PreProcess.sh
#
# --------------------------------------------------------------------------------------------------

library(data.table)
library(filesstrings)

args <- commandArgs()
MAIN <- substr(args[6], 1, nchar(args[6])-1)
ID <- args[7]
TYPE <- args[8]


## Searches for meta data sheet in ... ##
# ... <MIDTERM> #
if(file.exists(file.path("/omics/gpcf/midterm", paste0(0, ID, "/data/", ID, "_meta.tsv")))) meta <- file.path("/omics/gpcf/midterm", paste0(0, ID, "/data/", ID, "_meta.tsv"))

# ... <MAIN> #
if(file.exists(file.path(MAIN, paste0(ID, "_meta.tsv")))) meta <- file.path(MAIN, paste0(ID, "_meta.tsv"))


## Prepares meta data for looping ##
# Loads meta data #
meta <- data.frame(fread(meta))

# Removes 'Undetermined' rows #
if(any(grepl("Undetermined", meta[, 1]))) meta <- meta[-grep("Undetermined", meta[, 1]), ]

# Converts 'SAMPLE_ID' to 'SAMPLE_NAME'; if necessary #
if(any(grepl("SAMPLE_ID", colnames(meta)))) colnames(meta)[which(colnames(meta) == "SAMPLE_ID")] <- "SAMPLE_NAME"

# Saves cleaned up meta file for longterm storage #
fwrite(meta, file = file.path(MAIN, "meta_data", paste0(ID, "_meta.tsv")), sep = "\t")

# Creates column with path to raw data on LSF cluster #
if(TYPE == "MeDIP") type <- "cfmedip_sequencing"
if(TYPE == "hMeSEAL") type <- "hme_seal_sequencing"
if(TYPE == "WGS") type <- "whole_genome_sequencing"
if(!"CUSTOM_PATH" %in% colnames(meta)){
  
  meta$CUSTOMER_TAGS <- gsub("ODCF_Project': '", "", str_sub(str_sub(meta$CUSTOMER_TAGS, 3, 50), -50, -3))
  meta$CLUSTER_PATH <- file.path("/omics/odcf/project", gsub("*_.*", "", meta$CUSTOMER_TAGS), tolower(gsub(paste0(gsub("*_.*", "", meta$CUSTOMER_TAGS), "_")[1], "", meta$CUSTOMER_TAGS)), "sequencing", type, "core", paste0("run", meta$RUN_ID, "/"))
  meta$MIDTERM_PATH <- file.path("/omics/gpcf/midterm", paste0("0", gsub("-.*", "", ID)), "data", meta$RUN_ID, str_sub(meta$FASTQ_FILE, -100, -13), "fastq/")
  
}
if("CUSTOM_PATH" %in% colnames(meta)){
  
  meta$CLUSTER_PATH <- "UNKNOWN"
  meta$MIDTERM_PATH <- "UNKNOWN"
  
}

# Checks if all required columns are present #
check <- c()
for(i in c("FASTQ_FILE", "SAMPLE_NAME", "ILSE_NO", "RUN_ID", "CLUSTER_PATH", "MIDTERM_PATH", "CUSTOM_PATH")){
  
  check <- c(check, any(grepl(i, colnames(meta))))
  
}

# Subsets and orders selected columns and saves file #
if(all(check)){
  
  meta <- meta[, colnames(meta) %in% c("FASTQ_FILE", "SAMPLE_NAME", "ILSE_NO", "RUN_ID", "CLUSTER_PATH", "MIDTERM_PATH", "CUSTOM_PATH")]
  meta <- meta[, c("FASTQ_FILE", "SAMPLE_NAME", "ILSE_NO", "RUN_ID", "CLUSTER_PATH", "MIDTERM_PATH", "CUSTOM_PATH")]
  
  ## Save meta file without header ##
  fwrite(meta, file = file.path(MAIN, "tmp", ID, paste0(ID, "_meta.tsv")), sep = "\t", col.names = FALSE)
  
}

