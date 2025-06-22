#!/usr/bin/env Rscript

library(HMMcopy)
library(GenomicRanges)
library(optparse)
library(ichorCNA)
library(ggplot2)

options(stringsAsFactors=FALSE, scipen=0)
options(bitmapType='cairo')

option_list <- list(
	make_option(c("--gcWig"), type = "character", help = "GC Wig file for reference genome"),
	make_option(c("--mapWig"), type = "character", default=NULL, help = "Mappabiliy Wig file for reference genome"),
	make_option(c("--repTimeWig"), type = "character", default=NULL, help = "Path to replication timing WIG file.")
	make_option(c("-f", "--filelist"), type = "character", help = "List of of wig files."),
	make_option(c("-o", "--outfile"), type = "character", help = "Output file."),
	make_option(c("-c", "--centromere"), type="character", help = "File containing Centromere locations"),
	make_option(c("--rmCentromereFlankLength"), type="numeric", default=1e5, help="Length of region flanking centromere to remove. Default: [%default]"),
	make_option(c("--chrs"), type="character", default="c(1:22,\"X\")", help = "Specify chromosomes to analyze."),
	make_option(c("--genomeBuild"), type="character", default="hg19", help="Genome build. Default: [%default]"),
    make_option(c("--genomeStyle"), type = "character", default = "NCBI", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
	make_option(c("--chrNormalize"), type="character", default="c(1:22)", help = "Specify chromosomes to normalize GC/mappability biases"),
	make_option(c("--minMapScore"), type = "numeric", default=0.0, help="Include bins with a minimum mappability score of this value. Default: [%default]."),
  	make_option(c("--maleChrXLogRThres"), type="numeric", default=-0.80, help = "ChrX Log ratio threshold to confirm as male gender."),
	make_option(c("--fracReadsInChrYForMale"), type="numeric", default=0.001, help = "Threshold for fraction of reads in chrY to assign as male. Default: [%default]"),
	make_option(c("-e", "--exons.bed"), type = "character", default=NULL, help = "Path to bed file containing exon regions."),
	make_option(c("--method"), type = "character", default="median", help="Median or Mean.")
	make_option(c("--ylim"), type = "character", default="c(-2,2)", help="Y-limits for plotting of mean/median log ratios")
	make_option(c("--plotChrPanels"), type = "logical", default=FALSE, help="Plot PoN values")
)

parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)

#id <- opt$id
gcWig <- opt$gcWig
mapWig <- opt$mapWig
repTimeWig <- opt$repTimeWig
filelist <- opt$filelist
exons.bed <- opt$exons.bed  # "0" if none specified
centromere <- opt$centromere
method <- opt$method
outfile <- opt$outfile
genomeStyle <- opt$genomeStyle
genomeBuild <- opt$genomeBuild
ylim <- eval(parse(text = opt$ylim))
maleChrXLogRThres <- opt$maleChrXLogRThres
chrs <- as.character(eval(parse(text = opt$chrs)))
chrNormalize <- as.character(eval(parse(text=opt$chrNormalize))); 
flankLength <- opt$rmCentromereFlankLength
minMapScore <- opt$minMapScore
fracReadsInChrYForMale <- opt$fracReadsInChrYForMale					
plotChrPanels <- opt$plotChrPanels					
					
#####

## load functions from MESA package
my_path <- ${projectDir}/assets/ichorCNA/Rscripts_v050  # set your path
source_files <- list.files(my_path, "*.R$")  # locate all .R files
map(paste0(my_path, source_files), source)  # source all your R scripts!


createPanelOfNormals(gcWig, mapWig, repTimeWig = repTimeWig, filelist, outfile, centromere, flankLength = flankLength,
					chrs = chrs, genomeStyle = genomeStyle, genomeBuild = genomeBuild, 
					chrNormalize = chrNormalize, minMapScore = minMapScore, maleChrXLogRThres = maleChrXLogRThres, 
					fracReadsInChrYForMale = fracReadsInChrYForMale, exons.bed = exons.bed, method = method, ylim = ylim, 
					plotChrPanels = plotChrPanels) 
