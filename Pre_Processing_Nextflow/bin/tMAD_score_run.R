#!/usr/bin/env Rscript

###############################################
# Calculate t-MAD score
###############################################
#### https://github.com/sdchandra/tMAD

library(CNAclinic)

args = commandArgs(trailingOnly=TRUE)

###########################################
# Arguments to change
###########################################

# 1) Path to all your downsampled (sampl_10Million_***) BAM files

downsampled_tumors  = read.table(args[1], header=FALSE, col.names="Samples")               # "/path/to/sampled/bamfiles"

# 2) What is the total read count expected
 
total_reads_needed = args[2]

# 3) Which genomic bin size should be used as resolution

selected_binSize = args[3]
selected_binSize = as.numeric(selected_binSize)

# 4) Path to the control that should be used to normalise the samples 
# This could be a single control, or a control made from a panel of controls
# It needs to have the same total number of reads as the test samples

downsampled_control = args[4]                                          # "path/to/sampled/control/samp_10Million_XXX.bam"
downsampled_control = file.path(downsampled_control)
# 5) The path to blacklist created in the previous step

path_to_blacklist = args[5]                                           # "path/to/saved/blacklist"
path_to_blacklist = file.path(path_to_blacklist)
# 6) User defined name for the control used to normalise

control_name = args[6]                                                 # "CONTROL_XXX"
control_name = as.character(control_name)

# 7) Name of downsampled tumor file, without .bam
#tumor_name = args[7]

###########################################
# Check the code below for necessary changes
###########################################

total_reads_label <- paste0(total_reads_needed, "M")

controlSample <- c(control_name, "none")

segType <- c("CH")

for(t in 1:length(total_reads_label)){
    
    if(length(total_reads_label) != length(selected_binSize))
        stop("Check variables")
    
    binSize <- selected_binSize[t]

    # cfDNA_blacklist <- readRDS(paste0(path_to_blacklist, "control_blacklist_", total_reads_label[t], "_", binSize, "Kbp"))
    cfDNA_blacklist <- readRDS(paste0(path_to_blacklist))
    
    # bamfiles <- Sys.glob(paste0(downsample_dir, "/samp_", total_reads_label[t], "_*.bam"))
    bamfiles <- as.vector(downsampled_tumors$Samples)
    
    bamnames <- unlist(lapply(strsplit(bamfiles, ".bam"), function(x){ x[[length(x)]]}))        # changed by anja

    bamfiles <- c(downsampled_control, bamfiles)
    
    bamnames <-  c(control_name, bamnames)
    
    
    for(j in 1:length(controlSample)){
        
        # Doing only median normalization for these samples        
        control <- controlSample[j]
        
        processedData <- NULL
        
        if(control == "none"){
            
            outfile_suffix <- paste0(total_reads_label[t], "_",
                                     binSize, "_control_none")
            
            processedData <- processForSegmentation(
                                bamfiles=bamfiles,
                                binSize=binSize,
                                chromosomesFilter=c("X", "Y", "MT", "M"),
                                cache=FALSE,
                                isPaired=TRUE,
                                saveCountData=FALSE)


            # Run segmentation
            if(segType == "CH"){
                CNAData <- runSegmentation(processedData, genome="hg19",
                                           segmentType=c("CBS", "HMM"),
                                           summaryMethod="mean")
                
                saveRDS(CNAData, file=paste0("CNAData_", outfile_suffix, ".rds"))
            }
        }else if(control == control_name){
            
            if(length(control_name) != 1)
                stop("Check control specified")
            
            which_control <- which(bamnames == control_name)
            refSamples <- rep(control_name, length(bamnames))
            
            # drop the control sample after normalizing test samples 
            # as we will have log2R=0 if it is normalized by its own bin counts
            refSamples[which_control] <- 'drop'
            
            outfile_suffix <- paste0(total_reads_label[t], "_",
                                     binSize, "_control_",
                                     control)
            
            processedData <- processForSegmentation(cache=FALSE,
                                bamfiles=bamfiles,
                                bamnames=bamnames,
                                refSamples=refSamples,
                                binSize=binSize,
                                isPaired=TRUE,
                                skipMedianNormalization=TRUE,
                                saveCountData=FALSE)
            
            

            # Run segmentation
            if(segType == "CH"){
                CNAData <- runSegmentation(processedData, genome="hg19",
                                           segmentType=c("CBS", "HMM"),
                                           summaryMethod="mean")
                
                saveRDS(CNAData, file=paste0("CNAData_", outfile_suffix, ".rds"))
            }
        }
        
        #sampleNames <- unlist(lapply(sampleNames(CNAData),function(x){strsplit(x, ".vs.")[[1]][1]}))
        #sampleNames(CNAData) <- sampleNames
        
        stopifnot(length(cfDNA_blacklist) == length(usebin(CNAData)))
        
        segData <-  segSummary(CNAData)
        
        segData <- segData[cfDNA_blacklist & usebin(CNAData), ]
        
        segData[abs(segData) >=5 ] <- NA
        
        tMAD <- apply(segData, 2,  function(x){ 
            abs(mad(x = x, center=0, na.rm = TRUE))
        })
        
        outData <- data.frame(sampleNames=sampleNames(CNAData),
                              tMAD,
                              stringsAsFactors=FALSE)
        
        outData$sampleNames <- gsub(paste0("*.vs..", control_name), "", outData$sampleNames)
        colnames(outData)[1] <- c("Sample_ID")

        write.csv(outData, file=paste0("tMAD_", outfile_suffix, ".csv"), 
                  quote=F, row.names=F)
        
    }
}