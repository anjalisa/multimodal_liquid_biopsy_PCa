#
# Script name: tMAD_scores.R
#
# Author: Florian Janke
# Last updated: 20220725
#
# Description:
#   (1) Takes a vector of paths to .bam files and vector of paths to files used as panel of normals
#   (2) Downsamples all files to specified read depth
#   (3) Creates blacklist for t-MAD calculation
#   (4) Calculates t-MAD scores
#
# Arguments to provide:
#   (1) <bam.files>         Paths to .bam files to be processed
#   (2) <references>        Paths to .bam files to be used as panel of normals (has to be part of <bam.files>)
#   (3) <bam.blacklist>     Paths to .bam files to be used for blacklist creation (has to be part of <bam.files>)
#   (4) <downsample>        Number of reads to downsample to (paired-reads are given here!)
#   (5) <bin.size>          Bin size used for genome segregation (in kb)
#   (6) <path.out>          Path to output t-MAD scores and downsampled .bam files
#
# -------------------------------------------------------------------------------------------------------


tMAD_scores <- function(bam.files, references, bam.blacklist, downsample, blacklist = TRUE, bin.size, path.out) {

  library(QDNAseq.hg19)
  library(CNAclinic)
  library(stringr)
  library(data.table)

  # Create downsampling directory (if necessary)
  if(!dir.exists(file.path(path.out, "downsampling"))) dir.create(file.path(path.out, "downsampling"))
  if(!dir.exists(file.path(path.out, "downsampling", paste0(downsample/1000000, "M")))) dir.create(file.path(path.out, "downsampling", paste0(downsample/1000000, "M")))

  # Downsampling
  for(i in 1:length(bam.files)){

    # Check if file was already downsampled
    if(file.exists(file.path(path.out, "downsampling", paste0(downsample/1000000, "M"), basename(bam.files[i])))) next

    # Get coverage
    command <- paste0("/software/samtools/1.9/bin/samtools view -c -@ 10 -q 10 -F 4 -F 2048 -F 256 -F 1024 ", bam.files[i])
    coverage <- as.numeric(system(command, intern=TRUE, wait=TRUE))

    # Calculate downsampling factor
    factor <- 8 +(downsample *2 /coverage)

    # Downsample
    command <- paste0("/software/samtools/1.9/bin/samtools view -s ", factor, " -q 10 -F 4 -F 2048 -F 1024 -b -@ 10 ", bam.files[i], " > ", file.path(path.out, "downsampling", paste0(downsample/1000000, "M"), basename(bam.files[i])))
    system(command, intern=FALSE, wait=TRUE)

    # Sort downsampled .bam file
    command <- paste0("/software/samtools/1.9/bin/samtools sort ", file.path(path.out, "downsampling", paste0(downsample/1000000, "M"), basename(bam.files[i])), " -o ", file.path(path.out, "downsampling", paste0(downsample/1000000, "M"), basename(bam.files[i])), " -@ 10")
    system(command, intern=FALSE, wait=TRUE)

    # Index .bam file
    command <- paste0("/software/samtools/1.9/bin/samtools index ", file.path(path.out, "downsampling", paste0(downsample/1000000, "M"), basename(bam.files[i])), " -@ 10")
    system(command, intern=FALSE, wait=TRUE)

    # Message
    message("Downsampling ", i, " of ", length(bam.files), " ... !")

  }

  # Create control reference file (i.e., merge, sort and downsample)
  if(!file.exists(file.path(path.out, "downsampling", paste0(downsample/1000000, "M"), "Reference.bam"))){

    # Merge
    command <- paste0("/software/samtools/1.9/bin/samtools merge ", file.path(path.out, "downsampling", paste0(downsample/1000000, "M"), "Reference_tmp.bam"), " ", paste0(file.path(path.out, "downsampling", paste0(downsample/1000000, "M"), basename(references)), collapse = " "), " -@ 10")
    system(command, intern=FALSE, wait=TRUE)

    # Sort
    command <- paste0("/software/samtools/1.9/bin/samtools sort -n ", file.path(path.out, "downsampling", paste0(downsample/1000000, "M"), "Reference_tmp.bam"), " -o ", file.path(path.out, "downsampling", paste0(downsample/1000000, "M"), "Reference_tmp.bam -@ 10"))
    system(command, intern=FALSE, wait=TRUE)

    # Downsample (downsampling if off by a factor of 2 after merging)
    command <- paste0("/software/samtools/1.9/bin/samtools view -c -@ 10 -q 10 -F 4 -F 2048 -F 256 -F 1024 ", file.path(path.out, "downsampling", paste0(downsample/1000000, "M"), "Reference_tmp.bam"))
    coverage <- as.numeric(system(command, intern=TRUE, wait=TRUE))
    factor <- 8 +(downsample /coverage)

    command <- paste0("/software/samtools/1.9/bin/samtools view -s ", factor, " -q 10 -F 4 -F 2048 -F 1024 -b -@ 10 ", file.path(path.out, "downsampling", paste0(downsample/1000000, "M"), "Reference_tmp.bam"), " > ", file.path(path.out, "downsampling", paste0(downsample/1000000, "M"), "Reference.bam"))
    system(command, intern=FALSE, wait=TRUE)

    # Remove temporary reference file
    file.remove(file.path(path.out, "downsampling", paste0(downsample/1000000, "M"), "Reference_tmp.bam"))

  }


  # Create blacklist
  userMadeBins <- QDNAseq::getBinAnnotations(binSize = bin.size, genome="hg19")
  readCounts <- QDNAseq::binReadCounts(bins = userMadeBins, bamfiles = file.path(path.out, "downsampling", paste0(downsample/1000000, "M"), basename(bam.blacklist)), cache = TRUE, pairedEnds = TRUE)
  readCounts <- QDNAseq::applyFilters(readCounts, residual = FALSE, blacklist = FALSE, mappability = FALSE, bases = FALSE, chromosomes = c("X", "Y", "MT", "M"))
  userMadeBins$residual <- QDNAseq::iterateResiduals(readCounts)
  chromosomes = c("X", "Y", "MT", "M")
  condition <- rep(TRUE, times=nrow(readCounts))
  condition <- !(Biobase::fData(readCounts)$chromosome %in% chromosomes)
  condition <- condition & !is.na(Biobase::fData(readCounts)$gc)
  residuals <- userMadeBins$residual
  cutoff <- TRUE * matrixStats::madDiff(residuals, na.rm=TRUE)
  residualsMissing <- aggregate(residuals, by=list(chromosome=Biobase::fData(readCounts)$chromosome), function(x) all(is.na(x)))
  chromosomesWithResidualsMissing <- residualsMissing$chromosome[residualsMissing$x]
  chromosomesToInclude <- setdiff(chromosomesWithResidualsMissing, chromosomes)
  if (length(chromosomesToInclude) > 0) {
    message("Note: Residual filter missing for chromosomes: ", paste(chromosomesToInclude, collapse=", "))
    residuals[Biobase::fData(readCounts)$chromosome %in% chromosomesToInclude] <- 0
  }
  blacklist <- condition & !is.na(residuals)
  rm(userMadeBins, residualsMissing, residuals, readCounts, chromosomes, chromosomesToInclude, chromosomesWithResidualsMissing, condition)

  # Calculate t-MAD scores
  segType <- c("CH")
  bamnames <- c("Reference", gsub(".bam", "", basename(bam.files)))
  control <- control_name <- "Reference"

  if(length(control_name) != 1) stop("Check control specified")

  which_control <- which(bamnames == control_name)
  refSamples <- rep(control_name, length(bamnames))

  # drop the control sample after normalizing test samples 
  # as we will have log2R=0 if it is normalized by its own bin counts
  refSamples[which_control] <- 'drop'


#outfile_suffix <- paste0(total_reads_label, "_",
  #                         binSize, "kb")

  processedData <- processForSegmentation(cache=FALSE,
                                          bamfiles = c(file.path(path.out, "downsampling", paste0(downsample/1000000, "M"), "Reference.bam"), file.path(path.out, "downsampling", paste0(downsample/1000000, "M"), basename(bam.files))),
                                          bamnames = bamnames,
                                          refSamples = refSamples,
                                          binSize = bin.size,
                                          isPaired = TRUE,
                                          skipMedianNormalization = TRUE,
                                          saveCountData = FALSE)



  # Run segmentation
  if(segType == "CH"){
    CNAData <- runSegmentation(processedData, genome="hg19",
                               segmentType = c("CBS", "HMM"),
                               summaryMethod = "mean")

  }

  sampleNames <- unlist(lapply(sampleNames(CNAData), function(x){strsplit(x, ".vs.")[[1]][1]}))
  stopifnot(length(blacklist) == length(usebin(CNAData)))
  segData <-  segSummary(CNAData)
  segData <- segData[blacklist & usebin(CNAData), ]
  segData[abs(segData) >=5 ] <- NA

  tMAD <- apply(segData, 2,  function(x){
    abs(mad(x = x, center=0, na.rm = TRUE))
  })

  outData <- data.frame(sampleNames=sampleNames(CNAData), tMAD, stringsAsFactors=FALSE)
  outData$sampleNames <- gsub(".bam.*", "", outData$sampleNames)
  colnames(outData)[1] <- c("Sample_ID")
  outData$Sample_ID <- gsub("\\.", "-", gsub(".vs..*", "", outData$Sample_ID))

  # Save t-MAD scores
  fwrite(outData, file = file.path(path.out, "t-MAD_scores.txt"))

  return(outData)


}
  
