#!/usr/bin/env Rscript
####################
#     ichorCNA     #
####################
# file:   ichorCNA.R
# author: Gavin Ha, Ph.D.
#         Fred Hutchinson Cancer Research Center
# contact: <gha@fredhutch.org>
# # website: https://GavinHaLab.org
#
# author: Justin Rhoades, Broad Institute
#
# ichorCNA website: https://github.com/GavinHaLab/ichorCNA
# date:   January 6, 2020
#
# description: Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.
# This script is the main script to run the HMM.

# Adapted by: Anja Riediger

library(optparse)

option_list <- list(
  make_option(c("--WIG"), type = "character", help = "Path to tumor WIG file. Required."),
  make_option(c("--NORMWIG"), type = "character", default=NULL, help = "Path to normal WIG file. Default: [%default]"),
  make_option(c("--gcWig"), type = "character", help = "Path to GC-content WIG file; Required"),
  make_option(c("--mapWig"), type = "character", default=NULL, help = "Path to mappability score WIG file. Default: [%default]"),
  make_option(c("--repTimeWig"), type = "character", default=NULL, help ="Path to replication timing WIG file. Default: [%default]"),
  make_option(c("--normalPanel"), type="character", default=NULL, help="Median corrected depth from panel of normals. Default: [%default]"),
  make_option(c("--sex"), type = "character", default = NULL, help = "User specified gender: male or female [Default: %default]"),
  make_option(c("--exons.bed"), type = "character", default=NULL, help = "Path to bed file containing exon regions. Default: [%default]"),
  make_option(c("--id"), type = "character", default="test", help = "Patient ID. Default: [%default]"),
  make_option(c("--centromere"), type="character", default=NULL, help = "File containing Centromere locations; if not provided then will use hg19 version from ichorCNA package. Default: [%default]"),
  make_option(c("--minMapScore"), type = "numeric", default=0.9, help="Include bins with a minimum mappability score of this value. Default: [%default]."),
  make_option(c("--rmCentromereFlankLength"), type="numeric", default=1e5, help="Length of region flanking centromere to remove. Default: [%default]"),
  make_option(c("--normal"), type="character", default="0.5", help = "Initial normal contamination; can be more than one value if additional normal initializations are desired. Default: [%default]"),
  make_option(c("--normal.init"), type="character", default="c(0.5, 0.5)", help = "Specific initialization of normal contamination for multiple samples. Default: [%default]"),
  make_option(c("--scStates"), type="character", default="NULL", help = "Subclonal states to consider. Default: [%default]"),
  make_option(c("--scPenalty"), type="numeric", default=0.1, help = "Penalty for subclonal state transitions. 0.1 penalizes subclonal states by ~10%. Default: [%default]"),
  make_option(c("--normal2IgnoreSC"), type="numeric", default=1.0, help="Ignore subclonal analysis when normal proportion is >= this value. Default: [%default]"),
  make_option(c("--coverage"), type="numeric", default=NULL, help = "PICARD sequencing coverage. Default: [%default]"),
  make_option(c("--likModel"), type="character", default="t", help="Likelihood model to use. \"t\" or \"gaussian\". Use \"gaussian\" for faster runtimes. Default: [%default]"),
  make_option(c("--lambda"), type="character", default="NULL", help="Initial Student's t precision; must contain 4 values (e.g. c(1500,1500,1500,1500)); if not provided then will automatically use based on variance of data. Default: [%default]"),
  make_option(c("--lambdaScaleHyperParam"), type="numeric", default=3, help="Hyperparameter (scale) for Gamma prior on Student's-t precision. Default: [%default]"),
  #	make_option(c("--kappa"), type="character", default=50, help="Initial state distribution"),
  make_option(c("--ploidy"), type="character", default="2", help = "Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired. Default: [%default]"),
  make_option(c("--maxCN"), type="numeric", default=7, help = "Total clonal CN states. Default: [%default]"),
  make_option(c("--estimateNormal"), type="logical", default=TRUE, help = "Estimate normal. Default: [%default]"),
  make_option(c("--estimateScPrevalence"), type="logical", default=TRUE, help = "Estimate subclonal prevalence. Default: [%default]"),
  make_option(c("--estimatePloidy"), type="logical", default=TRUE, help = "Estimate tumour ploidy. Default: [%default]"),
  make_option(c("--maxFracCNASubclone"), type="numeric", default=0.7, help="Exclude solutions with fraction of subclonal events greater than this value. Default: [%default]"),
  make_option(c("--maxFracGenomeSubclone"), type="numeric", default=0.5, help="Exclude solutions with subclonal genome fraction greater than this value. Default: [%default]"),
  make_option(c("--minSegmentBins"), type="numeric", default=50, help="Minimum number of bins for largest segment threshold required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction."),
  make_option(c("--altFracThreshold"), type="numeric", default=0.05, help="Minimum proportion of bins altered required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction. Default: [%default]"),
  make_option(c("--chrNormalize"), type="character", default="c(1:22)", help = "Specify chromosomes to normalize GC/mappability biases. Default: [%default]"),
  make_option(c("--chrTrain"), type="character", default="c(1:22)", help = "Specify chromosomes to estimate params. Default: [%default]"),
  make_option(c("--chrs"), type="character", default="c(1:22,\"X\")", help = "Specify chromosomes to analyze. Default: [%default]"),
  make_option(c("--genomeBuild"), type="character", default="hg19", help="Geome build. Default: [%default]"),
  make_option(c("--genomeStyle"), type = "character", default = "NCBI", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
  make_option(c("--normalizeMaleX"), type="logical", default=TRUE, help = "If male, then normalize chrX by median. Default: [%default]"),
  make_option(c("--minTumFracToCorrect"), type="numeric", default=0.05, help = "Tumor-fraction correction of bin and segment-level CNA if sample has minimum estimated tumor fraction. [Default: %default]"),
  make_option(c("--fracReadsInChrYForMale"), type="numeric", default=0.001, help = "Threshold for fraction of reads in chrY to assign as male. Default: [%default]"),
  make_option(c("--includeHOMD"), type="logical", default=FALSE, help="If FALSE, then exclude HOMD state. Useful when using large bins (e.g. 1Mb). Default: [%default]"),
  make_option(c("--txnE"), type="numeric", default=0.9999999, help = "Self-transition probability. Increase to decrease number of segments. Default: [%default]"),
  make_option(c("--txnStrength"), type="numeric", default=1e7, help = "Transition pseudo-counts. Exponent should be the same as the number of decimal places of --txnE. Default: [%default]"),
  make_option(c("--multSampleTxnStrength"), type="numeric", default=1, help="Strength of same state transition between multiple samples. Default: [%default]"),
  make_option(c("--plotFileType"), type="character", default="pdf", help = "File format for output plots. Default: [%default]"),
	make_option(c("--plotYLim"), type="character", default="c(-2,2)", help = "ylim to use for chromosome plots. Default: [%default]"),
  make_option(c("--outDir"), type="character", default="./", help = "Output Directory. Default: [%default]"),
  make_option(c("--libdir"), type = "character", default=NULL, help = "Script library path. Usually exclude this argument unless custom modifications have been made to the ichorCNA R package code and the user would like to source those R files. Default: [%default]"),
  make_option(c("--cores"), type="numeric", default = 1, help = "Number of cores to use for EM. Default: [%default]")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=0, stringsAsFactors=F)

#library(ichorCNA)
library(HMMcopy)
library(GenomicRanges)
library(GenomeInfoDb)
library(data.table)
library(foreach)
library(doMC)
library(dplyr)
library(ggplot2)
#library(BSgenome.Hsapiens.UCSC.hg19)
#library(BSgenome.Hsapiens.UCSC.hg38)
options(stringsAsFactors=FALSE)
options(bitmapType='cairo')

patientID <- opt$id
tumour_file <- opt$WIG
normal_file <- opt$NORMWIG
gcWig <- opt$gcWig
mapWig <- opt$mapWig
repTimeWig <- opt$repTimeWig
normal_panel <- opt$normalPanel
sex <- opt$sex
exons.bed <- opt$exons.bed  # "0" if none specified
centromere <- opt$centromere
minMapScore <- opt$minMapScore
flankLength <- opt$rmCentromereFlankLength
normal <- eval(parse(text = opt$normal))
normal.init <- eval(parse(text = opt$normal.init))
scStates <- eval(parse(text = opt$scStates))
subclone.penalty <- opt$scPenalty
likModel <- opt$likModel
lambda <- eval(parse(text = opt$lambda))
lambdaScaleHyperParam <- opt$lambdaScaleHyperParam
estimateNormal <- opt$estimateNormal
estimatePloidy <- opt$estimatePloidy
estimateScPrevalence <- opt$estimateScPrevalence
maxFracCNASubclone <- opt$maxFracCNASubclone
maxFracGenomeSubclone <- opt$maxFracGenomeSubclone
minSegmentBins <- opt$minSegmentBins
altFracThreshold <- opt$altFracThreshold
ploidy <- eval(parse(text = opt$ploidy))
coverage <- opt$coverage
maxCN <- opt$maxCN
txnE <- opt$txnE
txnStrength <- opt$txnStrength
multSampleTxnStrength <- opt$multSampleTxnStrength
normalizeMaleX <- as.logical(opt$normalizeMaleX)
includeHOMD <- as.logical(opt$includeHOMD)
fracReadsInChrYForMale <- opt$fracReadsInChrYForMale
chrXMedianForMale <- -0.1
minTumFracToCorrect <- opt$minTumFracToCorrect     # adapted form original script
normal2IgnoreSC <- opt$normal2IgnoreSC
outDir <- opt$outDir
libdir <- opt$libdir
plotFileType <- opt$plotFileType
plotYLim <- eval(parse(text=opt$plotYLim))
outImage <- paste0(outDir,"/", patientID,".RData")
genomeBuild <- opt$genomeBuild
genomeStyle <- opt$genomeStyle
chrs <- as.character(eval(parse(text = opt$chrs)))
chrTrain <- as.character(eval(parse(text=opt$chrTrain))); 
chrNormalize <- as.character(eval(parse(text=opt$chrNormalize))); 
seqlevelsStyle(chrs) <- genomeStyle
seqlevelsStyle(chrNormalize) <- genomeStyle
seqlevelsStyle(chrTrain) <- genomeStyle
cores <- opt$cores 


## load ichorCNA library or source R scripts
if (!is.null(libdir) && libdir != "None"){
	source(paste0(libdir,"/R/utils.R"))
  source(paste0(libdir,"/R/segmentation.R"))
	source(paste0(libdir,"/R/EM.R"))
  source(paste0(libdir,"/R/output.R"))
	source(paste0(libdir,"/R/plotting.R"))
} else {
    library(ichorCNA)
}

## load seqinfo 

####################################
##### FUNCTION GET SEQINFO ######
####################################
#getSeqInfo <- function(genomeBuild = "hg19", genomeStyle = "NCBI"){
#	bsg <- paste0("BSgenome.Hsapiens.UCSC.", genomeBuild)
#	if (!require(bsg, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)) {
#		seqinfo <- Seqinfo(genome=genomeBuild)
#	} else {
#		seqinfo <- seqinfo(get(bsg))
#	}
#	seqlevelsStyle(seqinfo) <- genomeStyle
#	seqinfo <- keepSeqlevels(seqinfo, value = chrs)
#	#seqinfo <- cbind(seqnames = seqnames(seqinfo), as.data.frame(seqinfo))
#	return(seqinfo)	
#}

########
seqinfo <- getSeqInfo(genomeBuild, genomeStyle, chrs)

if (substr(tumour_file,nchar(tumour_file)-2,nchar(tumour_file)) == "wig") {
  wigFiles <- data.frame(cbind(patientID, tumour_file))
} else {
  wigFiles <- read.delim(tumour_file, header = FALSE, as.is = TRUE)
}

## FILTER BY EXONS IF PROVIDED ##
## add gc and map to GRanges object ##
if (is.null(exons.bed) || exons.bed == "None" || exons.bed == "NULL"){
  targetedSequences <- NULL
}else{
  targetedSequences <- fread(exons.bed) 
}

## load PoN
if (is.null(normal_panel) || normal_panel == "None" || normal_panel == "NULL"){
	normal_panel <- NULL
}

if (is.null(centromere) || centromere == "None" || centromere == "NULL"){ # no centromere file provided
	centromere <- system.file("extdata", "GRCh37.p13_centromere_UCSC-gapTable.txt", 
			package = "ichorCNA")
}
centromere <- read.delim(centromere,header=T,stringsAsFactors=F,sep="\t")
save.image(outImage)

## LOAD GC/MAP/REPTIME WIG FILES ###
message("Reading GC and mappability files")
gc <- wigToGRanges(gcWig)
if (is.null(gc)){
    stop("GC wig file not provided but is required")
}
map <- wigToGRanges(mapWig)
if (is.null(map)){
  message("No mappability wig file input, excluding from correction")
}
repTime <- wigToGRanges(repTimeWig)
if (is.null(repTime)){
  message("No replication timing wig file input, excluding from correction")
}else{
  if (mean(repTime$value, na.rm = TRUE) > 1){
    repTime$value <- repTime$value / 100 ## values in [0,1] - for LNCaP_repTime_10kb_hg38.txt
  }
}


## LOAD IN WIG FILES ##
numSamples <- nrow(wigFiles)
S <- numSamples
tumour_copy <- list()
counts <- list()

for (i in 1:numSamples) {
  id <- wigFiles[i,1]
  ## create output directories for each sample ##
  dir.create(paste0(outDir, "/", id, "/"), recursive = TRUE)
  ### LOAD TUMOUR AND NORMAL FILES ###
  message("Loading tumour file:", wigFiles[i,1])
  tumour_reads <- wigToGRanges(wigFiles[i,2])
  
  message("Correcting Tumour")
  
  counts[[id]] <- loadReadCountsFromWig(tumour_reads, chrs = chrs, gc = gc, map = map, repTime = repTime,
                                       centromere = centromere, flankLength = flankLength, 
                                       targetedSequences = targetedSequences, chrXMedianForMale = chrXMedianForMale,
                                       genomeStyle = genomeStyle, fracReadsInChrYForMale = fracReadsInChrYForMale,
                                       chrNormalize = chrNormalize, mapScoreThres = minMapScore)
  gender <- counts[[id]]$gender
  
  if ((!is.null(sex) && sex != "None" && sex != "NULL") && gender$gender != sex ){ #compare with user-defined sex
    message("Estimated gender (", gender$gender, ") doesn't match to user-defined gender (", sex, "). Use ", sex, " instead.")
    gender$gender <- sex
  }

  ## load in normal file if provided 
  if (!is.null(normal_file) && normal_file != "None" && normal_file != "NULL"){
	message("Loading normal file:", normal_file)
	normal_reads <- wigToGRanges(normal_file)
	message("Correcting Normal")
	counts.normal <- loadReadCountsFromWig(normal_reads, chrs=chrs, gc=gc, map=map, repTime = repTime,
  			centromere=centromere, flankLength = flankLength, targetedSequences=targetedSequences,
  			genomeStyle = genomeStyle, chrNormalize = chrNormalize, mapScoreThres = minMapScore)
  	normal_copy <- counts.normal$counts #as(counts$counts, "GRanges")
  	counts[[id]]$counts$cor.gc.normal <- counts.normal$counts$cor.gc
  	counts[[id]]$counts$cor.map.normal <- counts.normal$counts$cor.map
  	counts[[id]]$counts$cor.rep.normal <- counts.normal$counts$cor.rep
  	gender.normal <- counts[[id]]$gender
  }else{
	  normal_copy <- NULL
  }
  ### DETERMINE GENDER ###
  ## if normal file not given, use chrY, else use chrX
  message("Determining gender...", appendLF = FALSE)
  gender.mismatch <- FALSE
  if (!is.null(normal_copy)){
	if (gender$gender != gender.normal$gender){ #use tumour # use normal if given
	# check if normal is same gender as tumour
	  gender.mismatch <- TRUE
	}
  }
  message("Gender ", gender$gender)

  ## NORMALIZE GENOME-WIDE BY MATCHED NORMAL OR NORMAL PANEL (MEDIAN) ##
  tumour_copy[[id]] <- counts[[id]]$counts
  
  tumour_copy[[id]] <- normalizeByPanelOrMatchedNormal(tumour_copy[[id]], chrs = chrs, 
      normal_panel = normal_panel, normal_copy = normal_copy, 
      gender = gender$gender, normalizeMaleX = normalizeMaleX)
	
	### OUTPUT FILE ###
	### PUTTING TOGETHER THE COLUMNS IN THE OUTPUT ###
	outMat <- as.data.frame(tumour_copy[[id]])
	#outMat <- outMat[,c(1,2,3,12)]
	outMat <- outMat[,c("seqnames","start","end","copy")]
	colnames(outMat) <- c("chr","start","end","log2_TNratio_corrected")
	outFile <- paste0(outDir,"/",id,".correctedDepth.txt")
	message(paste("Outputting to:", outFile))
	write.table(outMat, file=outFile, row.names= FALSE, col.names = TRUE, quote = FALSE, sep="\t")

} ## end of for each sample

chrInd <- as.character(seqnames(tumour_copy[[1]])) %in% chrTrain
## get positions that are valid
valid <- tumour_copy[[1]]$valid & !is.na(tumour_copy[[1]]$copy)
if (length(tumour_copy) >= 2) {
  for (i in 2:length(tumour_copy)){ 
    valid <- valid & tumour_copy[[i]]$valid & !is.na(tumour_copy[[i]]$copy)
  } 
}
save.image(outImage)

  ### RUN HMM ###
registerDoMC()
options(cores = cores)
message("ichorCNA: Using ", getDoParWorkers(), " cores.")

## store the results for different normal and ploidy solutions ##
ptmTotalSolutions <- proc.time() # start total timer
results <- list()
numCombinations <- (length(normal) * length(ploidy)) ^ S
loglik <- as.data.frame(matrix(NA, nrow = numCombinations, ncol = 7, 
                 dimnames = list(c(), c("n_0", "phi_0", "n_est", "phi_est", 
                 												"Frac_genome_subclonal", "Frac_CNA_subclonal", "loglik"))))
fracGenomeSub <- as.data.frame(matrix(NA, nrow = numCombinations, ncol = S))
fracAltSub <- as.data.frame(matrix(NA, nrow = numCombinations, ncol = S))
# prepare normal/ploidy restarts for different solutions
normal.restarts <- expand.grid(rep(list(normal), S))
if (!is.null(normal.init) && numSamples > 1){
  normal.restarts <- unique(rbind(normal.restarts, normal.init))
}
#ploidy.restarts <- expand.grid(rep(list(ploidy), S))
counter <- 1
compNames <- rep(NA, nrow(loglik))
mainName <- vector('list', S) #rep(NA, nrow(normal.restarts) * nrow(ploidy.restarts))
if (numSamples > 1){
  maxCN <- min(maxCN, 4)
  message("Using maxCN=4")
}
#### restart for purity and ploidy values ####
for (i in 1:length(ploidy)){
  p <- rep(ploidy[i], S)
  for (j in 1:nrow(normal.restarts)){
    n <- as.numeric(normal.restarts[j, ])
    ## skip restarts where normal=0.95 and ploidy not diploid (2)
    if (sum(n == 0.95 & p != 2) > 0) {
      next
    }
    message("Running EM for normal=", paste(n, collapse=","), ", ploidy=", paste(p, collapse=","))
    # if (!grepl("gauss", likModel, ignore.case = TRUE)){
    #   likModel <- "Gaussian"
    #   message("Switching to Gaussian likelihood model.")
    # }
    
    logR <- as.data.frame(lapply(tumour_copy, function(x) { x$copy })) # NEED TO EXCLUDE CHR X #
    param <- getDefaultParameters(logR[valid & chrInd, , drop=F], n_0 = n, maxCN = maxCN, 
      includeHOMD = includeHOMD, 
      ct.sc=scStates, normal2IgnoreSC = normal2IgnoreSC, ploidy_0 = floor(p), 
      e=txnE, e.subclone = subclone.penalty, e.sameState = multSampleTxnStrength, 
      strength=txnStrength, likModel = likModel)

      ############################################
      ######## CUSTOM PARAMETER SETTINGS #########
      
      ##############################################################################
      # 0.1x cfDNA # ### used for plasma and urinary cfDNA samples
      ##############################################################################
      
      #if (is.null(lambda)){
      #  logR.var <- 1 / ((apply(logR, 2, sd, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
      #  param$lambda <- rep(logR.var, length(param$ct))
      #  param$lambda[param$ct %in% c(2)] <- logR.var 
      #  param$lambda[param$ct %in% c(1,3)] <- logR.var 
      #  param$lambda[param$ct >= 4] <- logR.var / 5
      #  param$lambda[param$ct == max(param$ct)] <- logR.var / 15
      #  param$lambda[param$ct.sc.status] <- logR.var / 10
      #}else{
      #  param$lambda[param$ct %in% c(2)] <- lambda[2]
      #  param$lambda[param$ct %in% c(1)] <- lambda[1]
      #  param$lambda[param$ct %in% c(3)] <- lambda[3]
      #  param$lambda[param$ct >= 4] <- lambda[4]
      #  param$lambda[param$ct == max(param$ct)] <- lambda[2] / 15
      #  param$lambda[param$ct.sc.status] <- lambda[2] / 10
      #}
      #param$alphaLambda <- rep(lambdaScaleHyperParam, length(param$ct))  
      
      ##############################################################################
      # 1x bulk tumors # ### used for PCa tissue and matched buffy coat samples
      ##############################################################################

      #param$lambda[param$ct %in% c(2)] <- 2000
      #param$lambda[param$ct %in% c(1)] <- 1750
      #param$lambda[param$ct %in% c(3)] <- 1750
      #param$lambda[param$ct >= 4] <- 1500
      #param$lambda[param$ct == max(param$ct)] <- 1000 / 25
      #param$lambda[param$ct.sc.status] <- 1000 / 75
      #param$alphaLambda[param$ct.sc.status] <- 4
      #param$alphaLambda[param$ct %in% c(1,3)] <- 5
      #param$alphaLambda[param$ct %in% c(2)] <- 5
      #param$alphaLambda[param$ct == max(param$ct)] <- 4
      
      #############################################
      ################ RUN HMM ####################
      #############################################
      hmmResults.cor <- HMMsegment(tumour_copy, valid, dataType = "copy", 
                                 param = param, chrTrain = chrTrain, maxiter = 20,
                                 estimateNormal = estimateNormal, estimatePloidy = estimatePloidy, 
                                 estimatePrecision = TRUE, estimateVar = TRUE, 
                                 estimateSubclone = estimateScPrevalence, estimateTransition = TRUE,
                                 estimateInitDist = TRUE, likChangeConvergence = 1e-4, verbose = TRUE)

    for (s in 1:numSamples){
  		iter <- hmmResults.cor$results$iter
  		id <- names(hmmResults.cor$cna)[s]
  
      ## correct integer copy number based on estimated purity and ploidy
   		correctedResults <- correctIntegerCN(cn = hmmResults.cor$cna[[s]],
     				segs = hmmResults.cor$results$segs[[s]],
      			purity = 1 - hmmResults.cor$results$n[s, iter], ploidy = hmmResults.cor$results$phi[s, iter],
      			cellPrev = 1 - hmmResults.cor$results$sp[s, iter],
      			maxCNtoCorrect.autosomes = maxCN, maxCNtoCorrect.X = maxCN, minPurityToCorrect = minTumFracToCorrect,
      			gender = gender$gender, chrs = chrs, correctHOMD = includeHOMD)
  		hmmResults.cor$results$segs[[s]] <- correctedResults$segs
  	  hmmResults.cor$cna[[s]] <- correctedResults$cn
  		## convert full diploid solution (of chrs to train) to have 1.0 normal or 0.0 purity
  		## check if there is an altered segment that has at least a minimum # of bins
  		segsS <- hmmResults.cor$results$segs[[s]]
  		segsS <- segsS[segsS$chr %in% chrTrain, ]
  		segAltInd <- which(segsS$event != "NEUT")
  		maxBinLength = -Inf
  		if (length(segAltInd) > 0){
  			maxInd <- which.max(segsS$end[segAltInd] - segsS$start[segAltInd] + 1)
  			maxSegRD <- GRanges(seqnames=segsS$chr[segAltInd[maxInd]], 
  								ranges=IRanges(start=segsS$start[segAltInd[maxInd]], end=segsS$end[segAltInd[maxInd]]))
  			hits <- findOverlaps(query=maxSegRD, subject=tumour_copy[[s]][valid, ])
  			maxBinLength <- length(subjectHits(hits))
  		}
		## check if there are proportion of total bins altered 
		# if segment size smaller than minSegmentBins, but altFrac > altFracThreshold, then still estimate TF
  		cnaS <- hmmResults.cor$cna[[s]]
  		altInd <- cnaS[cnaS$chr %in% chrTrain, "event"] == "NEUT"
  		altFrac <- sum(!altInd, na.rm=TRUE) / length(altInd)
  		if ((maxBinLength <= minSegmentBins) & (altFrac <= altFracThreshold)){
  			hmmResults.cor$results$n[s, iter] <- 1.0
  		}
        
        ## plot solution ##
        mainName[[s]] <- c(mainName[[s]], paste0("Solution: ", counter, ", Sample: ", id, ", n: ", n[s], ", p: ", p[s], 
        ", log likelihood: ", signif(hmmResults.cor$results$loglik[hmmResults.cor$results$iter], digits = 4)))
      ## plot individual samples 
      if (numSamples == 1){ # if only single-sample analysis
        outPlotFile <- paste0(outDir, "/", id, "/", id, "_genomeWide_", "n", n[s], "-p", p[s])      
        plotGWSolution(hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType=plotFileType, 
              logR.column = "logR", call.column = "Corrected_Call",
        			 plotYLim=plotYLim, estimateScPrevalence=estimateScPrevalence, seqinfo=seqinfo, main=mainName[[s]][counter])
      }
    }
    iter <- hmmResults.cor$results$iter
    results[[counter]] <- hmmResults.cor
    loglik[counter, "loglik"] <- signif(hmmResults.cor$results$loglik[iter], digits = 4)
    subClonalBinCount <- unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$subclone.status) }))
    fracGenomeSub[counter, ] <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ nrow(x) }))
    fracAltSub[counter, ] <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$copy.number != 2) }))
    fracAltSubVect <- lapply(fracAltSub[counter, ], function(x){if (is.na(x)){0}else{x}})
    loglik[counter, "Frac_genome_subclonal"] <- paste0(signif(fracGenomeSub[counter, ], digits=2), collapse=",")
    loglik[counter, "Frac_CNA_subclonal"] <- paste0(signif(as.numeric(fracAltSubVect), digits=2), collapse=",")
    loglik[counter, "n_0"] <- paste0(n, collapse = ",")
    loglik[counter, "phi_0"] <- paste0(p, collapse = ",")
    loglik[counter, "n_est"] <- paste(signif(hmmResults.cor$results$n[, iter], digits = 2), collapse = ",")
    loglik[counter, "phi_est"] <- paste(signif(hmmResults.cor$results$phi[, iter], digits = 4), collapse = ",")
    
    counter <- counter + 1
  }
}
  ## remove solutions witn a likelihood of NA (i.e. no solution)
  loglik <- loglik[!is.na(loglik$loglik), ]

  ## get total time for all solutions ##
  elapsedTimeSolutions <- proc.time() - ptmTotalSolutions
  message("Total ULP-WGS HMM Runtime: ", format(elapsedTimeSolutions[3] / 60, digits = 2), " min.")
  
  ### SAVE R IMAGE ###
  #save.image(outImage)
  #save(tumour_copy, results, loglik, file=paste0(outDir,"/",id,".RData"))
  
  ### SELECT SOLUTION WITH LARGEST LIKELIHOOD ###

 # if (estimateScPrevalence){ ## sort but excluding solutions with too large % subclonal 
 #   fracInd <- which(loglik[, "Frac_CNA_subclonal"] <= maxFracCNASubclone & 
 #                      loglik[, "Frac_genome_subclonal"] <= maxFracGenomeSubclone)
 #   if (length(fracInd) > 0){ ## if there is a solution satisfying % subclonal
 #     ind <- fracInd[order(loglik[fracInd, "loglik"], decreasing=TRUE)]
 #   }else{ # otherwise just take largest likelihood
 #     ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
 #   }
 # }else{#sort by likelihood only
 #   ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
 # }

  ### SORT SOLUTIONS BY DECREASING LIKELIHOOD ###
fracAltSub[is.nan(fracAltSub$V1), ] <- 0
if (estimateScPrevalence){ ## sort but excluding solutions with too large % subclonal 
  if (numSamples > 1){
      fracInd <- which(rowSums(fracAltSub <= maxFracCNASubclone) == S &
                      rowSums(fracGenomeSub <= maxFracGenomeSubclone) == S)
  }else{
      fracInd <- which(fracAltSub <= maxFracCNASubclone &
                      fracGenomeSub <= maxFracGenomeSubclone)
  }
	if (length(fracInd) > 0){ ## if there is a solution satisfying % subclonal
    # sort by highest loglik and then by lowest fracCNAsubclonal (if ties)
		ind <- fracInd[order(-loglik[fracInd, "loglik"], loglik[fracInd, "Frac_CNA_subclonal"])]
    # sort by highest loglik for solutions that don't pass fracCNAsubclonal threshold
    fracInd.exclude <- setdiff(1:nrow(loglik), ind)
    ind.excludeSC <- fracInd.exclude[order(-loglik[fracInd.exclude, "loglik"])]
    ind <- c(ind, ind.excludeSC)
	}else{ # otherwise just take largest likelihood
		ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
	}
}else{#sort by likelihood only
  ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
}
  
  ## print combined PDF of all solutions
  #new loop by order of solutions (ind)
# individual samples 
for (s in 1:numSamples){
  id <- names(hmmResults.cor$cna)[s]
  outPlotFile <- paste0(outDir, "/", id, "/", id, "_genomeWide_all_sols")
  for (i in 1:length(ind)) {
    hmmResults.cor <- results[[ind[i]]]
    turnDevOff <- FALSE
    turnDevOn <- FALSE
    if (i == 1){
    	turnDevOn <- TRUE
    }
    if (i == length(ind)){
    	turnDevOff <- TRUE
    }
    plotGWSolution(hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType="pdf", 
                       logR.column = "logR", call.column = "Corrected_Call",
                       plotYLim=plotYLim, estimateScPrevalence=estimateScPrevalence, 
                       seqinfo = seqinfo,
                       turnDevOn = turnDevOn, turnDevOff = turnDevOff, main=mainName[[s]][ind[i]])
  }
}

# multisample, combined plots
if (numSamples > 1){
  outPlotFile <- paste0(outDir, "/", patientID, "_all_sols_multiSample")
  if (plotFileType == "png"){ 
    outPlotFile <- paste0(outPlotFile, ".png")
    png(outPlotFile,width=20,height=numSamples*6,units="in",res=300)
   }else{
    outPlotFile <- paste0(outPlotFile, ".pdf")
    pdf(outPlotFile,width=20,height=numSamples*6)
  } 
  for (i in 1:length(ind)) {
    hmmResults.cor <- results[[ind[i]]]
    par(mfrow=c(numSamples, 1))
    for (s in 1:numSamples){
      plotGWSolution(hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType="pdf", 
                       logR.column = "logR", call.column = "Corrected_Call",
                       plotYLim=plotYLim, estimateScPrevalence=estimateScPrevalence, 
                       seqinfo = seqinfo, spacing = 4, 
                       cex.text = 1.25, cex=0.5,
                       turnDevOn = FALSE, turnDevOff = FALSE, main=mainName[[s]][ind[i]])
    }
  }
  dev.off()
}

  ## output text files
  hmmResults.cor <- results[[ind[1]]]
  hmmResults.cor$results$loglik <- as.data.frame(loglik)
  hmmResults.cor$results$gender <- gender$gender
  hmmResults.cor$results$chrYCov <- gender$chrYCovRatio
  hmmResults.cor$results$chrXMedian <- gender$chrXMedian
  hmmResults.cor$results$coverage <- coverage
  
  outputHMM(cna = hmmResults.cor$cna, segs = hmmResults.cor$results$segs, 
            results = hmmResults.cor$results, patientID = patientID, outDir=outDir)
  outFile <- paste0(outDir, "/", patientID, ".params.txt")
  outputParametersToFile(hmmResults.cor, file = outFile)
  
  ## plot solutions for all samples 
  plotSolutions(hmmResults.cor, tumour_copy, chrs, outDir, counts, numSamples=numSamples,
              logR.column = "logR", call.column = "event", likModel = likModel,
              plotFileType=plotFileType, plotYLim=plotYLim, seqinfo = seqinfo,
              estimateScPrevalence=estimateScPrevalence, maxCN=maxCN)



#################################### PLOTTING ####################################

####
if (substr(tumour_file,nchar(tumour_file)-2,nchar(tumour_file)) == "wig") {
  wigFiles <- data.frame(cbind(patientID, tumour_file))
} else {
  wigFiles <- read.delim(tumour_file, header=F, as.is=T)
}

numSamples <- nrow(wigFiles)

for (i in 1:numSamples) {
  id <- wigFiles[i,1]
  
  # Input path
  #path.in <- file.path(outDir, paste0(id))

    # Extract TFx and ploidy #
    params <- data.frame(fread(paste0(outDir,"/", id, ".params.txt"), sep = ";"))
    params <- params[c(6:7), ]
    params <- data.frame(do.call("rbind", strsplit(as.character(params), "\t", fixed = TRUE)))
    params[, 1] <- gsub(":", "", params[, 1])
    params[, 2] <- as.numeric(as.character(params[, 2]))
    
    # Prepare for plotting #
    cna <- data.frame(fread(paste0(outDir,"/", id, ".cna.seg")))
    txt <- data.frame(fread(paste0(outDir,"/", id, ".seg.txt")))
    
    txt <- txt[, c(2:4, 6)]
    colnames(txt) <- c("chr", "start", "end", "seg")      
    cna <- cna[, c(1:4, 6, 9)]
    cna$bin <- seq(1:nrow(cna))
    colnames(cna) <- c("chr", "start", "end", "copy_number", "logR", "annot", "bin")
    #cna <- cna[!cna$chr == "X", ]        ## adapted from original script (# added)
    cna$copy_number = as.character(cna$copy_number)
    cna <- data.frame(cna)
    txt <- merge(cna, txt, by = c("chr", "end"))
    txt <- txt[, c(7, 9)]
    txt <- txt[order(txt$bin), ]
    y <- data.frame(bin = 1, seg = txt[1, 2])
    
    
    txt$start <- txt$bin
    txt$end <- txt$bin
    y <- data.frame(start = txt$start)
    x <- data.frame(start = -1)
    y <- rbind(x, y)
    y <- data.frame(start = y[-c(nrow(y)), ])
    txt$start <- y$start
    txt$start <- txt$start + 1
    txt$bin <- NULL
    
    # Creating data necessary for visualization #
    chr <- cna %>% dplyr::group_by(chr) %>% dplyr::summarize(count = n()) #make sure plyr package is not loaded!!!
    chr$chr <-  sub("chr*", "", chr$chr)    ## adapted from original script
    chr$chr <- factor(chr$chr, levels=c(seq(1:22), "X", "Y")) # adapted from original script
    chr <- chr[order(chr$chr), ]
    chr$counts_half <- chr$count / 2
    chr$cumsum <- cumsum(chr$count)
    chr$tick <- chr$cumsum -chr$count + chr$counts_half
    
    cna$chr <- sub("chr*", "", cna$chr)    ## adapted from original script
    
    # Defining colors #
    colors <- c("3" = "#8a0000", "4" = "#ff3737", "5" = "#ff3737", "2" = "#0000ff", "1" = "#006300")
    
    # Plotting #
    base_size <- 14
    TFx <- format(round(params[1, 2], digits=2), nsmall = 2)
    
    min(cna$logR, na.rm = TRUE)
    
    #!!! adapted from original script: c(seq(1:22),"X", "Y")
    ggplot() +
      geom_point(data = cna, aes(x = bin, y = logR, color = copy_number), size = 0.75) +
      scale_color_manual(values = colors) +
      geom_vline(xintercept = c(0, chr$cumsum), linetype = "dashed", alpha = 0.5, size = 0.4) +
      geom_segment(aes(x = 0, xend = nrow(cna[cna$chr %in% c(seq(1:22), "X", "Y"), ]), y = 0, yend = 0), alpha = 0.4) +
      scale_x_continuous(name = "Chromosomes", limits = c(0, nrow(cna[cna$chr %in% c(seq(1:22), "X", "Y"), ])), breaks = chr[!chr$chr %in% c(15,17,19,21,13, 20, 11), ]$tick, labels = chr[!chr$chr %in% c(15,17,19,21,13, 20, 11), ]$chr) +
      scale_y_continuous(limits = c(-2, 2), name = "Copy number (log2 ratio)") +
      labs(title = id, subtitle = paste0("TFx: ", TFx)) +
      geom_segment(aes(y = txt$seg, yend = txt$seg, x = txt$start, xend = txt$end), size = 0.8) +
      theme(
        axis.line =         element_blank(),
        axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, vjust = 0.5, angle = 0),
        axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, hjust = 1),
        axis.title.x =      element_text(size = base_size * 0.8, lineheight = 0.9, hjust = 0.5),
        axis.title.y =      element_text(size = base_size * 0.9, angle = 90, vjust = 2),
        axis.ticks.x =      element_blank(), 
        axis.ticks.length.x = unit(0.3, "lines"),
        axis.ticks.length.y = unit(0.2, "lines"),
        panel.background =  element_rect(fill = "white", colour = NA), 
        panel.border =      element_rect(fill = NA, colour="grey20"), 
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank(),
        strip.background =  element_rect(fill = "grey80", colour = "grey50"), 
        strip.text.x =      element_text(size = base_size * 0.8),
        strip.text.y =      element_text(size = base_size * 0.8, angle = -90),
        plot.background =   element_rect(colour = NA),
        plot.title =        element_text(size = base_size * 1.2),
        legend.position = "none"
      )
    
    # Saving plot #
    ggsave(
      file.path(outDir, paste0(id, ".svg")),
      plot = ggplot2::last_plot(),
      device = "png",
      scale = 1.25,
      width = 20,
      height = 8,
      units = "cm",
      dpi = 300,
      bg = "transparent",
      limitsize = TRUE)
    
    ggsave(
      file.path(outDir, paste0(id, ".png")),
      plot = ggplot2::last_plot(),
      device = "png",
      scale = 1.25,
      width = 20,
      height = 8,
      units = "cm",
      dpi = 200,
      bg = "transparent",
      limitsize = TRUE)

# Summary table #
    params <- data.frame(Sample_ID = id, TFx = params[1, 2], Ploidy = params[2, 2])
    
    if(i == 1){
      
      table <- params
      
    } else {
      
      table <- rbind(table, params)
      
    }
   
}

## Saving summary table ##
if(file.exists(file.path(outDir, paste0("ichorCNA_summary.tsv")))){
  
  summary <- data.frame(fread(file.path(outDir, paste0("ichorCNA_summary.tsv"))))
  summary <- rbind(summary, table)
  
} else {
  
  summary <- table
  
}

fwrite(summary, file = file.path(outDir, paste0("ichorCNA_summary.tsv")), sep = "\t")
