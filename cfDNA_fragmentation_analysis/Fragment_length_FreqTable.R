#
# Script Name: Fragment_length_FreqTable.R <OUTPUT> <NAME>
#
# Author: Florian Janke (DKFZ Heidelberg)
# Last updated: 2022/03/03
#
# Description
#   Greps insert sizes extracted from bam-files (with samtools) and prepares a frequency table
#
# Input
# text file which stores insert sizes, previously created with bash-script "Fragment_length.sh"
# --------------------------------------------------------------------------------------------------

library(data.table)

args <- commandArgs()
OUTPUT <- args[6]
NAME <- args[7]

## Loading insert sizes ##
data <- data.frame(fread(file.path(OUTPUT, paste0(NAME, ".txt"))))


## Preparing frequency table ##
data[, 1] <- abs(data[, 1])
data <- data.frame(table(data))
colnames(data) <- c("length", "freq")
data$length <- as.numeric(as.character(data$length))


# Saving frequency table #
fwrite(data, file = file.path(OUTPUT, paste0(NAME, ".txt")))