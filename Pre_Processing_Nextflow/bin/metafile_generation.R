#######
library(dplyr)                        # Load dplyr package
library(readr)                        # Load readr package
library(tidyverse)                    # load tidyverse package

#lese metafile="/Users/riediger/Documents/Lab_work_DKFZ/Versuche/exp15_cfMeDIP-Seq_run4/sequencing/27217/27217_meta.tsv
metadata_27215<- read.table(file="/Users/riediger/Documents/Lab_work_DKFZ/Versuche/exp15_cfMeDIP-Seq_run4/sequencing/27215/27215_meta.tsv", header=T,sep='\t')

metadata_27216<- read.table(file="/Users/riediger/Documents/Lab_work_DKFZ/Versuche/exp15_cfMeDIP-Seq_run4/sequencing/27216/27216_meta.tsv", header=T,sep='\t')

metadata_27217<- read.table(file="/Users/riediger/Documents/Lab_work_DKFZ/Versuche/exp15_cfMeDIP-Seq_run4/sequencing/27217/27217_meta.tsv", header=T,sep='\t')

metadata_27218<- read.table(file="/Users/riediger/Documents/Lab_work_DKFZ/Versuche/exp15_cfMeDIP-Seq_run4/sequencing/27218/27218_meta.tsv", header=T,sep='\t')

## kombiniere alle Metadaten-Files
metadata_all <- bind_rows(metadata_27215,metadata_27216, metadata_27217, metadata_27218)

# selektiere Columns
metadata_all_reduced <- select(metadata_all, FASTQ_FILE, READ, LANE_NO, SAMPLE_NAME, ILSE_NO, PATIENT_ID, TISSUE_TYPE, RUN_ID, LANE_NO, CENTER_NAME, INSTRUMENT_PLATFORM, INSTRUMENT_MODEL, INDEX, ILSE_NO, SEQUENCING_READ_TYPE)
metadata_all_reduced <- filter(metadata_all_reduced, !grepl("Undetermined", FASTQ_FILE))

### füge die Spalten "SOURCE" und "SAMPLE_TYPE" hinzu 
metadata_all_select_filt_newcolumn <- add_column(metadata_all_reduced, .after="SAMPLE_NAME", SOURCE=NA, LIBRARY_TYPE=NA, LIBRARY_TYPE_SHORT=NA)

# schreibe "U1" in Spalte "SOURCE", wenn bei SAMPLE_NAME ein "US1 oder US2" vorhanden ist; und "U3", wenn "US3 oder US4" vorhanden ist
metadata_all_select_filt_newcolumn[grep("US1|US2", metadata_all_select_filt_newcolumn$SAMPLE_NAME),]$SOURCE <-"U1"
metadata_all_select_filt_newcolumn[grep("US3|US4", metadata_all_select_filt_newcolumn$SAMPLE_NAME),]$SOURCE <-"U3"

# schreibe "P1" in Spalte "SOURCE", wenn bei SAMPLE_NAME ein "P1 oder P2" vorhanden ist; und "P3", wenn "P3 oder P4" vorhanden sind
metadata_all_select_filt_newcolumn[grep("P1|P2", metadata_all_select_filt_newcolumn$SAMPLE_NAME),]$SOURCE <-"P1"
metadata_all_select_filt_newcolumn[grep("P3|P4", metadata_all_select_filt_newcolumn$SAMPLE_NAME),]$SOURCE <-"P3"

# schreibe "B1" in Spalte "SOURCE", wenn bei SAMPLE_NAME ein "B1 oder B2" vorhanden ist und "B3" wenn "B3 oder B4" vorhanden ist
metadata_all_select_filt_newcolumn[grep("B1|B2", metadata_all_select_filt_newcolumn$SAMPLE_NAME),]$SOURCE <-"B1"
metadata_all_select_filt_newcolumn[grep("B3|B4", metadata_all_select_filt_newcolumn$SAMPLE_NAME),]$SOURCE <-"B3"


# schreibe "MEDIP" in Spalte "LIBRARY_TYPE", wenn bei SAMPLE_NAME ein "IP" vorhanden ist und "M" in LIBRARY_TYPE_SHORT
metadata_all_select_filt_newcolumn[grep("IP", metadata_all_select_filt_newcolumn$SAMPLE_NAME),]$LIBRARY_TYPE <-"MEDIP"
metadata_all_select_filt_newcolumn[grep("IP", metadata_all_select_filt_newcolumn$SAMPLE_NAME),]$LIBRARY_TYPE_SHORT <-"M"

# schreibe "WGS" in Spalte "LIBRARY_TYPE", wenn bei SAMPLE_NAME ein "IC" vorhanden ist und "W" in LIBRARY_TYPE_SHORT
metadata_all_select_filt_newcolumn[grep("IC", metadata_all_select_filt_newcolumn$SAMPLE_NAME),]$LIBRARY_TYPE <-"WGS"
metadata_all_select_filt_newcolumn[grep("IC", metadata_all_select_filt_newcolumn$SAMPLE_NAME),]$LIBRARY_TYPE_SHORT <-"W"

### create new column with name "SAMPLE_ID" for the processing pipeline
metadata_all_final <- add_column(metadata_all_select_filt_newcolumn, .before="FASTQ_FILE", SAMPLE_ID=NA)
metadata_all_final$SAMPLE_ID <- paste(metadata_all_final$TISSUE_TYPE, metadata_all_final$PATIENT_ID, metadata_all_final$SOURCE, metadata_all_final$LIBRARY_TYPE_SHORT, sep="-")
metadata_all_final <- metadata_all_final %>% select(SAMPLE_ID, FASTQ_FILE, READ, SAMPLE_NAME, SOURCE, PATIENT_ID, LIBRARY_TYPE, LIBRARY_TYPE_SHORT, TISSUE_TYPE, INDEX, ILSE_NO, RUN_ID, LANE_NO, INSTRUMENT_MODEL, INSTRUMENT_PLATFORM, CENTER_NAME, SEQUENCING_READ_TYPE)

write_csv(metadata_all_final, file="/Users/riediger/Documents/Lab_work_DKFZ/Versuche/exp14_cfMeDIP-Seq_run3/metadata_all.csv")

### Konvertiere die Liste in wide-Format
metadata_all_wider <- pivot_wider(metadata_all_final, names_from = READ, values_from = FASTQ_FILE) 
metadata_all_wider <- rename(metadata_all_wider, c("FASTQ_1" = "1", "FASTQ_2" = "2", "FASTQ_UMI" = "I1"))
metadata_all_wider <- metadata_all_wider %>% select(SAMPLE_ID, FASTQ_1, FASTQ_2, FASTQ_UMI, SAMPLE_NAME, SOURCE, PATIENT_ID, LIBRARY_TYPE, LIBRARY_TYPE_SHORT, TISSUE_TYPE, INDEX, ILSE_NO, RUN_ID, LANE_NO, INSTRUMENT_MODEL, INSTRUMENT_PLATFORM, CENTER_NAME, SEQUENCING_READ_TYPE)

write_csv(metadata_all_wider, file="/Users/riediger/Documents/Lab_work_DKFZ/Versuche/exp14_cfMeDIP-Seq_run3/metadata_all_wider.csv")

# create a short version of metadata_wider
metadata_all_wider_short <- select(metadata_all_wider, SAMPLE_ID, FASTQ_1, FASTQ_2, FASTQ_UMI)
write_csv(metadata_all_wider_short, file="/Users/riediger/Documents/Lab_work_DKFZ/Versuche/exp14_cfMeDIP-Seq_run3/metadata.csv")
