#!/usr/bin/env Rscript

#######################################################################################
# Create empirical blacklist from control data for further calculation of t-MAD score
########################################################################################
#### https://github.com/sdchandra/tMAD

library(CNAclinic)
library(QDNAseq)

args = commandArgs(trailingOnly=TRUE)

downsampled_controls = read.table(args[1], header=FALSE, col.names="Samples")           # "/path/to/sampled/control_bamfiles"
total_reads_needed = args[2]       # Number of total downsampled reads
binSizes = c(50,500,1000)              # Reads binned into x Kbp windows, example "c(30)" means 30 Kbp windows : binSizes <- c(30) # Reads binned into 30 Kbp windows

###########################################
###########################################
binSizes <- as.numeric(binSizes)

total_reads_label <- paste0(total_reads_needed, "M")
bamfiles <- as.vector(downsampled_controls$Samples)

for(b in 1:length(binSizes)){
    
    binSize <- binSizes[b]

    for(t in 1:length(total_reads_label)){

        userMadeBins <- QDNAseq::getBinAnnotations(binSize=binSize,
                                                   genome="hg19")
        
        
        readCounts <- QDNAseq::binReadCounts(bins=userMadeBins,
                                             bamfiles=bamfiles,
                                             cache=TRUE,
                                             pairedEnds=TRUE)
        
        ctrl <- readCounts
        
        readCounts <- QDNAseq::applyFilters(readCounts, residual=FALSE, 
                                            blacklist=FALSE,
                                            mappability=FALSE, 
                                            bases=FALSE,
                                            chromosomes = c("X", "Y", "MT", "M"))
        
        userMadeBins$residual <- QDNAseq::iterateResiduals(readCounts)
        
        chromosomes = c("X", "Y", "MT", "M")
        
        # Create a residual filter from cfDNA controls
        condition <- rep(TRUE, times=nrow(readCounts))
        condition <- !(Biobase::fData(readCounts)$chromosome %in% chromosomes)
        condition <- condition & !is.na(Biobase::fData(readCounts)$gc)
        
        residuals <- userMadeBins$residual
        cutoff <- TRUE * matrixStats::madDiff(residuals, na.rm=TRUE)
        residualsMissing <- aggregate(residuals,
                            by=list(chromosome=Biobase::fData(readCounts)$chromosome),
                            function(x) all(is.na(x)))
        chromosomesWithResidualsMissing <-
                residualsMissing$chromosome[residualsMissing$x]
        chromosomesToInclude <-
                setdiff(chromosomesWithResidualsMissing, chromosomes)
        if (length(chromosomesToInclude) > 0) {
            message("Note: Residual filter missing for chromosomes: ",
            paste(chromosomesToInclude, collapse=", "))
                residuals[Biobase::fData(readCounts)$chromosome %in% chromosomesToInclude] <- 0
        }

        # If FALSE, filter genomic bin from analysis, if TRUE keep bin
        condition <- condition & !is.na(residuals)

        saveRDS(condition, 
                file=paste0("control_blacklist_", total_reads_label[t], 
                            "_", binSize, "Kbp",".rds"))   
        
    }
}