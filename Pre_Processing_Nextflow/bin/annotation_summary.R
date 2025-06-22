################################################
##     Feature annotation summary - HOMER     ##
################################################

###############################################
# Rscript ./Homer_Script.R <MAIN> <META> <ID> #
###############################################

#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(data.table)


#option_list <- list(make_option(c("-i", "--peak_files"), type="character", default=NULL, help="Comma-separated list of annotated files.", metavar="path"),
#                    make_option(c("-is", "--peak_files_shuffled"), type="character", default=NULL, help="Comma-separated list of annotated files.", metavar="string"),
#                    make_option(c("-s", "--sample_ids"), type="character", default=NULL, help="Comma-separated list of sample ids associated with peak files. Must be unique and in same order as peaks files input.", metavar="string"),
#                    make_option(c("-m", "--multiqc_txt"), type="character", default=NULL, help="Path to tab-delimited file, containing 22 columns from multiqc analysis ", metavar="string"),
#                    
#                    make_option(c("-o", "--outdir"), type="character", default='./', help="Output directory", metavar="path"),
#                    make_option(c("-p", "--outprefix"), type="character", default='macs2_peakqc', help="Output prefix", metavar="string"))

#opt_parser <- OptionParser(option_list=option_list)
#opt <- parse_args(opt_parser)

#peakfiles_unshuffled <- unlist(strsplit(opt$peak_files,","))
#peakfiles_shuffled <- unlist(strsplit(opt$peak_files_shuffled,","))
#SampleIDs <- unlist(strsplit(opt$sample_ids,","))
# reads <- data.frame(fread(opt$multiqc_txt))

args <- commandArgs(trailingOnly=TRUE)
peakfiles_unshuffled <- read.table(args[1], header=FALSE, col.names="Samples")
peakfiles_shuffled <- read.table(args[2], header=FALSE, col.names="Samples")
multiqc <- data.frame(fread(args[3]))

peakfiles_unshuffled <- as.character(peakfiles_unshuffled)
peakfiles_shuffled <- as.character(peakfiles_shuffled)

# Remove last backslash #
#MAIN <- substr(MAIN, 1, nchar(MAIN)-1)

## Load meta data information ##
#meta <- read.table(file.path(META), sep = "\t")[, c(2)]
#meta <- as.character(meta[!duplicated(meta)])

## Load number of paired reads per sample ##
#### multiqc_fastqc.txt files (first column (Sample), fifth column (Total Sequences))


## Get number of peaks per sample and number of peaks per feature ##
#for (sample in 1:length(SampleIDs)) {

for (idx in 1:length(peakfiles_unshuffled)) {
    sampleid = sub("*.bam", "", idx)
    #sampleid = SampleIDs[sample]

    data <- data.frame(fread(idx))

    features <- c("Promoter", "5' UTR", "1st Exon", "TES", "3' UTR", "Exon", "Intron", "Intergenic")
    
    # Adjust Annotation column to only include the features-of-interest #
    for(n in features){
      
      if(n == "Promoter") data[grep("TSS", data$Annotation), ]$Annotation <- n
      if(n == "5' UTR") data[grep("5' UTR", data$Annotation), ]$Annotation <- n
      if(n == "1st Exon") data[grep("exon 1 of", data$Annotation), ]$Annotation <- n
      if(n == "3' UTR") data[grep("3' UTR", data$Annotation), ]$Annotation <- n
      if(n == "TES") data[grep("TTS", data$Annotation), ]$Annotation <- n
      if(n == "Intron") data[grep("intron", data$Annotation), ]$Annotation <- n
      if(n == "Exon") data[grep("exon", data$Annotation), ]$Annotation <- n
      if(n == "Intergenic") data[grep("Intergenic", data$Annotation), ]$Annotation <- n
      
    }
     
     # Get number of peaks #
      input.reads <- reads[which(reads$Sample_ID == sampleid), which(colnames(reads) == "Total Sequences")]
      number.peaks <- data.frame(Sample_ID = sampleid, number.peaks = nrow(data), number.reads = input.reads, peaks.per.million = nrow(data) / input.reads *1000000)
      rm(input.reads)
    
    # Remove features not within features-of-interest #
    data <- data[data$Annotation %in% features, ]
    
    # Count peaks per feature #
    data$count <- 1
    data <- data %>% dplyr::group_by(Annotation) %>% dplyr::summarize(instances = sum(count, na.rm = TRUE))
    
    # Intermediate storage #
      data.samples <- data
      colnames(data.samples)[2] <- "Samples"
    
    rm(data)
}
    

for (idx in 1:length(peakfiles_shuffled)) {
  #sampleid = SampleIDs[idx]
  sampleid = sub("*.bam", "", idx)

    data <- data.frame(fread(idx))

    features <- c("Promoter", "5' UTR", "1st Exon", "TES", "3' UTR", "Exon", "Intron", "Intergenic")
    
    # Adjust Annotation column to only include the features-of-interest #
    for(n in features){
      
      if(n == "Promoter") data[grep("TSS", data$Annotation), ]$Annotation <- n
      if(n == "5' UTR") data[grep("5' UTR", data$Annotation), ]$Annotation <- n
      if(n == "1st Exon") data[grep("exon 1 of", data$Annotation), ]$Annotation <- n
      if(n == "3' UTR") data[grep("3' UTR", data$Annotation), ]$Annotation <- n
      if(n == "TES") data[grep("TTS", data$Annotation), ]$Annotation <- n
      if(n == "Intron") data[grep("intron", data$Annotation), ]$Annotation <- n
      if(n == "Exon") data[grep("exon", data$Annotation), ]$Annotation <- n
      if(n == "Intergenic") data[grep("Intergenic", data$Annotation), ]$Annotation <- n
      
    }
     
    # Remove features not within features-of-interest #
    data <- data[data$Annotation %in% features, ]
    
    # Count peaks per feature #
    data$count <- 1
    data <- data %>% dplyr::group_by(Annotation) %>% dplyr::summarize(instances = sum(count, na.rm = TRUE))
    
    # Intermediate storage #
      data.shuffled <- data
      colnames(data.shuffled)[2] <- "Random"
    
    rm(data)
}

  # Combine sample and shuffled data #
  data <- merge(data.samples, data.shuffled, by = "Annotation")
  rm(data.samples, data.shuffled)
  
  # Get log2 ratio #
  data$log2 <- log2(data$Samples / data$Random)
  
  # Annotate Sample ID #
  data$Sample_ID <- sampleid
  data <- data[, c(ncol(data), 1:(ncol(data)-1))]
  
  # Append samples #
  if(i == 1){
    
    table.annot <- data
    table.number <- number.peaks
    
  } else {
    
    table.annot <- rbind(table.annot, data)
    table.number <- rbind(table.number, number.peaks)
    
  }

#}

## Save feature annotation table ##
if (file.exists("Annotation_summary.tsv")) {
  
  annot.summary <- data.frame(fread(paste0("Annotation_summary.tsv")))
  annot.summary <- rbind(annot.summary, table.annot)
  
} else {
  
  annot.summary <- table.annot
  
}

fwrite(annot.summary, file = paste0("Annotation_summary.tsv"), sep = "\t")

## Save number of peak table ##
if (file.exists("Peak_number_summary.tsv")) {
  
  peak.summary <- data.frame(fread(paste0("Peak_number_summary.tsv")))
  peak.summary <- rbind(peak.summary, table.number)
  
} else {
  
  peak.summary <- table.number
  
}

fwrite(peak.summary, file = paste0("Peak_number_summary.tsv"), sep = "\t")
