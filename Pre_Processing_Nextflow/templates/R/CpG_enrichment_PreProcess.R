#
# Script Name: CpG_enrichment_PreProcess.sh <MAIN> <META> <ID> <TYPE>
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

library(BSgenome.Hsapiens.UCSC.hg19)
library(gtools)
library(data.table)
library(MEDIPS)

args <- commandArgs()
MAIN <- args[6]
META <- args[7]
ID <- as.numeric(args[8])
TYPE <- args[9]

###################################
# CpG enrichment score - FUNCTION #
###################################
MEDIPS.CpGenrich.v2 <-function(file=NULL, BSgenome=NULL, extend=0, shift=0, uniq=1e-3, chr.select=NULL, paired=F){
  
  ## Proof correctness....
  if(is.null(BSgenome)){stop("Must specify a BSgenome library.")}
  
  ## Read region file		
  fileName=unlist(strsplit(file, "/"))[length(unlist(strsplit(file, "/")))]
  path=paste(unlist(strsplit(file, "/"))[1:(length(unlist(strsplit(file, "/"))))-1], collapse="/") 
  if(path==""){path=getwd()}		
  if(!fileName%in%dir(path)){stop(paste("File", fileName, " not found in", path, sep =" "))}	
  
  dataset = get(ls(paste("package:", BSgenome, sep = "")))	
  
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


# Remove last backslash #
MAIN <- substr(MAIN, 1, nchar(MAIN)-1)

## Load ILSE meta data information ##
files <- read.table(file.path(META), sep = "\t")[, c(2)]
files <- as.character(files[!duplicated(files)])


## CpG enrichment score calculation ##
# Chromosomes to include #
chr.select <- paste0("chr", c(1:22))

# Genome build to be used #
BSgenome <- "BSgenome.Hsapiens.UCSC.hg19"
dataset <- BSgenome.Hsapiens.UCSC.hg19

# Loop #
for(i in 1:length(files)){
  
  name <- paste0(files[i], ".bam")
  
  file <- MEDIPS.CpGenrich.v2(file.path(MAIN, "bam", TYPE, name), BSgenome, chr.select = chr.select, paired = TRUE)
  
  if(i == 1){
    
    data <- data.frame(Sample_ID = gsub(".bam", "", name), enrichment.score = file$enrichment.score.relH)
    
  } else {
    
    data <- rbind(data, data.frame(Sample_ID = gsub(".bam", "", name), enrichment.score = file$enrichment.score.relH))
    
  }
  
  message("... Sample ", i, " of ", length(files), " ...")
  
  ## Save data ##
  write.table(data, file = file.path(MAIN, "tmp", ID, paste0("CpG_enrich_summary.tsv")), row.names = FALSE, sep = "\t")
  #write.table(data, file = "/omics/odcf/analysis/OE0309_projects/nsclc_immunotherapy/jankef/Sequencing_data/QC/summary/CpG_enrich_summary.tsv", row.names = FALSE, sep = "\t")
  
}









