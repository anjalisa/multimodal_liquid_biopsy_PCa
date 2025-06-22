#
# Author: Florian Janke ; Anja Riediger (Nextflow version)
# Last updated: 2022/02/28 ; 2022/07/08 (Nextflow version)
#
# Description
#   (1) Calculates CpG enrichment scores per sample
##
# --------------------------------------------------------------------------------------------------

#!/usr/bin/env Rscript

# Install R packages
#if (!require("BiocManager", quietly = TRUE)){
#     install.packages("BiocManager", dependencies=TRUE, repos='http://cloud.r-project.org/')
#}

#if (!require("BSgenome.Hsapiens.UCSC.hg19")){
#  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
#}

#if (!require("MEDIPS")){
#    BiocManager::install("MEDIPS")
#}

#if (!require("data.table")){
#    install.packages("data.table", dependencies=TRUE, repos='http://cloud.r-project.org/')
#}

#if (!require("gtools")){
#    install.packages("gtools", dependencies=TRUE, repos='http://cloud.r-project.org/')
#}


## Load R packages
library(BSgenome.Hsapiens.UCSC.hg19)
library(gtools)
library(data.table)
library(MEDIPS)

######################################################################
# CpG enrichment score - FUNCTION #
######################################################################
MEDIPS.CpGenrich.v2 <-function(file=NULL, BSgenome=NULL, extend=0, shift=0, uniq=1e-3, chr.select=NULL, paired=F){
  
  ## Proof correctness....
  if(is.null(BSgenome)){stop("Must specify a BSgenome library.")}
  
  ## Read region file		
  fileName=unlist(strsplit(file, "/"))[length(unlist(strsplit(file, "/")))]
  path=paste(unlist(strsplit(file, "/"))[1:(length(unlist(strsplit(file, "/"))))-1], collapse="/") 
  if(path==""){path=getwd()}		
  if(!fileName%in%dir(path)){stop(paste("File", fileName, " not found in", path, sep =" "))}	
  
  dataset=get(ls(paste("package:", BSgenome, sep = ""))[1])	
  #dataset=get(ls(paste0("package:", BSgenome)))

  if(!paired){GRange.Reads = getGRange(fileName, path, extend, shift, chr.select, dataset, uniq)}	else{GRange.Reads = getPairedGRange(fileName, path, extend, shift, chr.select, dataset, uniq)}
  
  ## Sort chromosomes
  if(length(unique(seqlevels(GRange.Reads)))>1){chromosomes=mixedsort(unique(seqlevels(GRange.Reads)))}
  if(length(unique(seqlevels(GRange.Reads)))==1){chromosomes=unique(seqlevels(GRange.Reads))}
  
  ## Get chromosome lengths for all chromosomes within data set.
  cat(paste("Loading chromosome lengths for ",BSgenome, "...\n", sep=""))		
  
  chr_lengths=as.numeric(seqlengths(dataset)[chromosomes])
  
  ranges(GRange.Reads) <- restrict(ranges(GRange.Reads),+1)
  
  ##Calculate CpG density for regions
  total=length(chromosomes)
  cat("Calculating CpG density for given regions...\n") 
  
  readsChars <- unlist(getSeq(dataset, GRange.Reads, as.character=TRUE))
  
  regions.CG = sum(vcountPattern("CG",readsChars))
  regions.C  = sum(vcountPattern("C",readsChars))
  regions.G  = sum(vcountPattern("G",readsChars))
  all.genomic= sum(width(readsChars))
  
  nReads <- length(readsChars)
  
  regions.relH=as.numeric(regions.CG)/as.numeric(all.genomic)*100
  regions.GoGe=(as.numeric(regions.CG)*as.numeric(all.genomic))/(as.numeric(regions.C)*as.numeric(regions.G))  
  
  CG <- DNAStringSet("CG")
  pdict0 <- PDict(CG)
  params <- new("BSParams", X = dataset, FUN = countPDict, simplify = TRUE, exclude = c("rand", "chrUn"))
  genome.CG=sum(bsapply(params, pdict = pdict0))			
  params <- new("BSParams", X = dataset, FUN = alphabetFrequency, exclude = c("rand", "chrUn"), simplify=TRUE)
  alphabet=bsapply(params)
  genome.l=sum(as.numeric(alphabet))
  genome.C=as.numeric(sum(alphabet[2,]))
  genome.G=as.numeric(sum(alphabet[3,]))
  genome.relH=genome.CG/genome.l*100
  genome.GoGe=(genome.CG*genome.l)/(genome.C*genome.G);
  
  ##Calculate CpG density for reference genome
  
  enrichment.score.relH=regions.relH/genome.relH	
  enrichment.score.GoGe=regions.GoGe/genome.GoGe	
  
  gc()
  return(list(genome=BSgenome, regions.CG=regions.CG, regions.C=regions.C, regions.G=regions.G, regions.relH=regions.relH, regions.GoGe=regions.GoGe, genome.C=genome.C, genome.G=genome.G, genome.CG=genome.CG, genome.relH=genome.relH, genome.GoGe=genome.GoGe, enrichment.score.relH=enrichment.score.relH, enrichment.score.GoGe=enrichment.score.GoGe))
  
}
############################################################################################################################################

## Define variables
args = commandArgs(trailingOnly=TRUE)

BAMFILES <- read.table(args[1], header=FALSE, col.names="Samples")       # Textfile with input files (= final BAM-files)
paired_single <- args[2]       # information of whether paired-end or single-end sequencing was performed

## Load meta data information ##
bamfiles <- as.character(BAMFILES$Samples)

# Check whether paired-end or single-end sequencing was performed
paired_end <- if(paired_single=='PAIRED') TRUE else FALSE

## CpG enrichment score calculation ##
# Chromosomes to include #
chr.select <- paste0("chr", c(1:22))

# Genome build to be used #
BSgenome <- "BSgenome.Hsapiens.UCSC.hg19"

# Loop # #### adapted ####
for(i in 1:length(bamfiles)){
  
  # Define bamname
  bamname <- gsub("*.final.bam", "", bamfiles[i])
  
  # Run command for Methylation enrichment score calculation 
  file <- MEDIPS.CpGenrich.v2(bamfiles[i], BSgenome, chr.select = chr.select, paired = paired_end)
  
  if(i == 1){
    data <- data.frame(Sample_ID = bamname, genome=file$genome, regions.CG=file$regions.CG, regions.C=file$regions.C, regions.G=file$regions.G, regions.relH=file$regions.relH, regions.GoGe=file$regions.GoGe, genome.C=file$genome.C, genome.G=file$genome.G, genome.CG=file$genome.CG, genome.relH=file$genome.relH, genome.GoGe=file$genome.GoGe, enrichment.score.relH=file$enrichment.score.relH, enrichment.score.GoGe=file$enrichment.score.GoGe)
  } else {
    data <- rbind(data, data.frame(Sample_ID = bamname, genome=file$genome, regions.CG=file$regions.CG, regions.C=file$regions.C, regions.G=file$regions.G, regions.relH=file$regions.relH, regions.GoGe=file$regions.GoGe, genome.C=file$genome.C, genome.G=file$genome.G, genome.CG=file$genome.CG, genome.relH=file$genome.relH, genome.GoGe=file$genome.GoGe, enrichment.score.relH=file$enrichment.score.relH, enrichment.score.GoGe=file$enrichment.score.GoGe))
  }

 message("... Sample ", bamfiles[i], " of ", length(bamfiles), " ...")
  
  ## Save data ##
  if(i==length(bamfiles)) {
    write.table(data, file = "CpG_enrich_summary.tsv", row.names = FALSE, sep = "\t")
  }
}