#
# Script Name: cfDNA_fragmentation_features.R
#
# Author: Anja Riediger (DKFZ Heidelberg)
# Last updated: 2025/06/17
#
# Description
#   Analysis of cfDNA fragmentation features based on relative frequency distributions of insert sizes 
# 
# Input
# insert size frequency tables: previoulsy created with R-script "Fragment_length_FreqTable.R" 
# ----------------------------------------------------------------------------------------------------------

### load R packages
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(stringr)
library(dplyr)                        
library(readr)                        
library(tidyverse)                    
library(gtools)
library(readxl)
library(rstatix)
library(matrixStats)
library(DescTools)
library(PMCMRplus)
library(liver)
library(zoo)

##### RELEVANT FUNCTIONS
calculate_peaks <- function(x) {
  peaks <-
    x %>%
    dplyr::filter(length %in% .data$length[which((rollapply(.data$freq, width = 5, FUN = function(x) x[3] == max(x), align = "center", fill = NA)))])
  return(peaks)
}

calculate_valleys <- function(x) {
  valleys <-
    x %>%
    dplyr::filter(length %in% .data$length[which((rollapply(.data$freq, width = 5, FUN = function(x) x[3] == min(x), align = "center", fill = NA)))])
  return(valleys)
}

##
calculate_interpeak_dist <- function(x) {
  dist <-
    x %>%
    mutate(interpeak_dist = .data$length - lag(.data$length))
  return(dist)
}

##
calculate_intervalley_dist <- function(x) {
  dist <-
    x %>%
    mutate(intervalley_dist = .data$length - lag(.data$length))
  return(dist)
}

##############
# define directory for output/results
outdir_results <- "/path/to/output/directory"

# define intput; load insert size frequency tables previously created with R script "Fragment_length_FreqTable.R"
path <- "/path/to/frequency_tables"
files <- list.files(path)

# Combine samples into one data.frame
for(i in 1:length(files)){
  
  # Load data
  data <- data.frame(fread(file.path(path, files[i])))
  
  # Add sample ID information
  data$Sample_ID <- gsub(".final.txt","", files[i])
  
  cat("Processing: " , unique(data$Sample_ID))
  
  # Get relative frequency
  data$freq <- data$freq / sum(data$freq, na.rm = TRUE)
  
  # add cumulative frequency
  table_cum_freq <- data %>%                             # Apply group_by & mutate functions
    group_by(Sample_ID) %>%
    dplyr::mutate(cum_freq = cumsum(freq))
  
  table_cum_freq <- as.data.frame(table_cum_freq)
  
  ######## Calculate general parameter: 
  # max frequency, fragment length at maximum frequency (=modal size), cumulative frequency at maximum relative frequency
  # mean fragment length
  # relative and cumulative frequency at 167bp and at 334bp
  # fragment length at cum-frequency = 0.5 >> = median fragment length
  
  max_frequency <- max(table_cum_freq$freq)
  max_frequency_cumfreq <- mean(table_cum_freq[table_cum_freq$freq == max(table_cum_freq$freq), ]$cum_freq)
  mode_fragment_length <- mean(table_cum_freq[table_cum_freq$freq == max(table_cum_freq$freq),]$length)
  
  mean_fragment_length <- sum(table_cum_freq$freq* table_cum_freq$length)
  
  rel_freq_167bp <- table_cum_freq[table_cum_freq$length == 167, ]$freq
      #cum_freq_167bp <- table_cum_freq[table_cum_freq$length == 167, ]$cum_freq
  rel_freq_334bp <- table_cum_freq[table_cum_freq$length == 334, ]$freq
      #cum_freq_334bp <- table_cum_freq[table_cum_freq$length == 334, ]$cum_freq
  
  ### define the fragment length, when cum freq is 0.5 (or at least the closest to 0.5, since it does not always reach exactly 0.5)
  cumfreq05_fragment_length <- table_cum_freq$length[min(which(table_cum_freq$cum_freq >= 0.5))]
  
  ##########################################################################
  ## Calculate proportions
  #########################################################################
  # frequency table ranges from 30 - 700 bp
  
  P30_60 <- round(sum(data[data$length >= 30 & data$length <= 60,  ]$freq), 5)
  P30_100 <- round(sum(data[data$length <= 100, ]$freq), 5)
  P30_150 <- round(sum(data[data$length <= 150, ]$freq), 5)
  P30_180 <- round(sum(data[data$length <= 180, ]$freq), 5)
  P90_150 <- round(sum(data[data$length >= 90 & data$length <= 150, ]$freq), 5)
  P160_180 <- round(sum(data[data$length >= 160 & data$length <= 180, ]$freq), 5)
  P163_169 <- round(sum(data[data$length >= 163 & data$length <= 169, ]$freq), 5)
  P180_220 <- round(sum(data[data$length >= 180 & data$length <= 220, ]$freq), 5)
  P150_300 <- round(sum(data[data$length >= 150 & data$length <= 300, ]$freq), 5)
  P250_320 <- round(sum(data[data$length >= 250 & data$length <= 320, ]$freq), 5)
  P250_420 <- round(sum(data[data$length >= 250 & data$length <= 420, ]$freq), 5)
  P324_344 <- round(sum(data[data$length >= 324 & data$length <= 344, ]$freq), 5)
  P420_700 <- round(sum(data[data$length >= 420 & data$length <= 700, ]$freq), 5)
  
  ##################################################################################
  # Perform analysis with different ratios of proportions
  ##################################################################################
  ## create table with the following ratios: 
  
  ratio_P30_150_P160_180 = round((sum(data[data$length <= 150, ]$freq) / sum(data[data$length >= 160 & data$length <= 180, ]$freq)), 5)
  ratio_P30_150_P163_169 = round((sum(data[data$length <= 150, ]$freq) / sum(data[data$length >= 163 & data$length <= 169, ]$freq)), 5)
  ratio_P30_100_P160_180 = round((sum(data[data$length <= 100, ]$freq) / sum(data[data$length >= 160 & data$length <= 180, ]$freq)), 5)
  ratio_P30_100_P163_169 = round((sum(data[data$length <= 100, ]$freq) / sum(data[data$length >= 163 & data$length <= 169, ]$freq)), 5)
  ratio_P30_150_P150_300 = round((sum(data[data$length <= 150, ]$freq) / sum(data[data$length >= 150 & data$length <= 300, ]$freq)), 5)
  ratio_P90_150_P163_169 = round((sum(data[data$length >= 90 & data$length <= 150, ]$freq) / sum(data[data$length >= 163 & data$length <= 169, ]$freq)), 5)
  ratio_P160_180_P250_420 = round((sum(data[data$length >= 160 & data$length <= 180, ]$freq) / sum(data[data$length >= 250 & data$length <= 420, ]$freq)), 5)
  
  #################################################################################
  ## 10bp-Oscillation
  #################################################################################
  ### ### ### ### 10bp oscillation in 30-150bp ### ### ### ### ### ### 
  
  # filter for fragment lengths in the range 30-150bp
  data_30_150 <- table_cum_freq[table_cum_freq$length <=150 ,]
  
  # define peaks/valleys by selecting the positions y such that y was the largest/lowest value in the interval [y − 2, y + 2]
  local_maxima_30_150 <- data_30_150 %>%
    calculate_peaks
  
  local_minima_30_150 <- data_30_150 %>%
    calculate_valleys
  
  # Initialize an dataframe to store the filtered rows; add the first row of the current dataframe to this new df, because the first row can´t be compared to one row above 
  # local_maxima_filtered_df <- data.frame(length = numeric(0), freq = numeric(0), Sample_ID , cum_freq = numeric(0))
  local_maxima_30_150_filtered <- local_maxima_30_150[1,]
  local_minima_30_150_filtered <- local_minima_30_150[1,]
  
  # Loop through each row starting from the second row
  for (n in 2:nrow(local_maxima_30_150)) {
    # Compare the current row's 'length' with the last row's 'length' from xx_filtered_df, and see whether the difference is >6
    if (local_maxima_30_150[n,]$length - local_maxima_30_150_filtered[nrow(local_maxima_30_150_filtered),]$length >= 6) {
      # If the current 'length' is at least >8 than the last row in xx_filtered_df, add the current row to the filtered dataframe
      local_maxima_30_150_filtered <- rbind(local_maxima_30_150_filtered, local_maxima_30_150[n, ])
    } else {
      next
    }
  }
  
  # Loop through each row starting from the second row
  for (n in 2:nrow(local_minima_30_150)) {
    # Compare the current row's 'length' with the last row's 'length' from xx_filtered_df, and see whether the difference is >6
    if (local_minima_30_150[n,]$length - local_minima_30_150_filtered[nrow(local_minima_30_150_filtered),]$length >= 6) {
      # If the current 'length' is at least >8 than the last row in xx_filtered_df, add the current row to the filtered dataframe
      local_minima_30_150_filtered <- rbind(local_minima_30_150_filtered, local_minima_30_150[n, ])
    } else {
      next
    }
  }
  
  # calculate interpeak-/intervalley distance
  local_maxima_30_150_filtered <- local_maxima_30_150_filtered %>%
    calculate_interpeak_dist 
  
  local_minima_30_150_filtered <- local_minima_30_150_filtered %>%
    calculate_intervalley_dist
  
  # first row always has interpeak distance = NA (because it can´t be compared to one row above) >> add 0 instead of NA
  local_maxima_30_150_filtered[1,]$interpeak_dist <- 0
  local_minima_30_150_filtered[1,]$intervalley_dist <- 0
  
  # count number of local maxima and local minima
  number_local_maxima_30_150 <- nrow(local_maxima_30_150_filtered)
  number_local_minima_30_150 <- nrow(local_minima_30_150_filtered)
  
  # Calculate the sum of "frequency" values at the filtered local maxima positions
  sum_local_maxima_30_150 <- sum(local_maxima_30_150_filtered$freq)
  # cum_sum_local_maxima_30_150 <- sum(local_maxima_30_150_filtered$cum_freq)
  
  # Calculate the sum of "frequency" values at the filtered local minima positions
  sum_local_minima_30_150 <- sum(local_minima_30_150_filtered$freq)
  # cum_sum_local_minima_30_150 <- sum(local_minima_30_150_filtered$cum_freq)
  
  # Calculate the mean interpeak distance
  mean_interpeak_distance_30_150 <- mean(local_maxima_30_150_filtered$interpeak_dist)
  mean_intervalley_distance_30_150 <- mean(local_minima_30_150_filtered$intervalley_dist)
  
  # Calculate the cum freq/density between 30bp and 150bp
  cum_freq_30_150 <- data_30_150[data_30_150$length ==150,]$cum_freq - data_30_150[data_30_150$length ==30,]$cum_freq
  
  # Subtract the sum of the minima from the sum of the height of the maxima ( = oscillation factor). The larger this difference, the more distinct are the peaks.
  oscillation_factor_freq_30_150 <- sum_local_maxima_30_150 - sum_local_minima_30_150
  # oscillation_factor_cum_freq_30_150 <- cum_sum_local_maxima_30_150 - cum_sum_local_minima_30_150
  
  # Create two dataframes with locations for local maxima and local minima for each sample
  Sample_name <-gsub(".final.txt","", files[i])
  
  data_local_maxima_30_150 <- data.frame(t(select(local_maxima_30_150_filtered, length)))
  rownames(data_local_maxima_30_150) <- Sample_name
  colnames(data_local_maxima_30_150) <- paste0("peak", seq_along(data_local_maxima_30_150))
  
  data_local_minima_30_150 <- data.frame(t(select(local_minima_30_150_filtered, length)))
  rownames(data_local_minima_30_150) <- Sample_name
  colnames(data_local_minima_30_150) <- paste0("valley", seq_along(data_local_minima_30_150))
  
  if(i == 1) {
    data_local_maxima_table_30_150 <- data_local_maxima_30_150
    data_local_minima_table_30_150 <- data_local_minima_30_150
  }
  if(i != 1) {
    data_local_maxima_table_30_150 <- dplyr::bind_rows(data_local_maxima_table_30_150, data_local_maxima_30_150)
    data_local_minima_table_30_150 <- dplyr::bind_rows(data_local_minima_table_30_150, data_local_minima_30_150)
  }
  #######################################################################
  ### ### ### ### 10bp oscillation in 150 - 300bp ### ### ### ### ### ### 
  
  # filter for fragment lengths in the range 30-150bp
  data_150_300 <- table_cum_freq[table_cum_freq$length >=150 & table_cum_freq$length <= 300,]
  
  # define peaks/valleys by selecting the positions y such that y was the largest/lowest value in the interval [y − 2, y + 2]
  local_maxima_150_300 <- data_150_300 %>%
    calculate_peaks
  
  local_minima_150_300 <- data_150_300 %>%
    calculate_valleys
  
  # Initialize an dataframe to store the filtered rows; add the first row of the current dataframe to this new df, because the first row can´t be compared to one row above 
  # local_maxima_filtered_df <- data.frame(length = numeric(0), freq = numeric(0), Sample_ID , cum_freq = numeric(0))
  local_maxima_150_300_filtered <- local_maxima_150_300[1,]
  local_minima_150_300_filtered <- local_minima_150_300[1,]
  
  if(nrow(local_maxima_150_300) > 2){
    # Loop through each row starting from the second row
    for (n in 2:nrow(local_maxima_150_300)) {
      # Compare the current row's 'length' with the last row's 'length' from xx_filtered_df, and see whether the difference is >6
      if (local_maxima_150_300[n,]$length - local_maxima_150_300_filtered[nrow(local_maxima_150_300_filtered),]$length >= 6) {
        # If the current 'length' is at least >8 than the last row in xx_filtered_df, add the current row to the filtered dataframe
        local_maxima_150_300_filtered <- rbind(local_maxima_150_300_filtered, local_maxima_150_300[n, ])
      } else {
        next
      }
    }
    
    # calculate interpeak distance
    local_maxima_150_300_filtered <- local_maxima_150_300_filtered %>%
      calculate_interpeak_dist 
    
    # first row always has interpeak distance = NA (because it can´t be compared to one row above) >> add 0 instead of NA
    local_maxima_150_300_filtered[1,]$interpeak_dist <- 0
    
    # count number of local maxima and local minima
    number_local_maxima_150_300 <- nrow(local_maxima_150_300_filtered)
    
    # Calculate the sum of "frequency" values at the filtered local maxima positions
    sum_local_maxima_150_300 <- sum(local_maxima_150_300_filtered$freq)
    # cum_sum_local_maxima_150_300 <- sum(local_maxima_150_300_filtered$cum_freq)
    
    # Calculate the mean interpeak distance
    mean_interpeak_distance_150_300 <- mean(local_maxima_150_300_filtered$interpeak_dist)
    
  } else {
    
    # count number of local maxima and local minima
    number_local_maxima_150_300 <- nrow(local_maxima_150_300_filtered)
    
    # Calculate the sum of "frequency" values at the filtered local maxima positions
    sum_local_maxima_150_300 <- NA
    
    # Calculate the mean interpeak distance
    mean_interpeak_distance_150_300 <- NA
    
  }
  
  if(nrow(local_minima_150_300) > 2){
    # Loop through each row starting from the second row
    for (n in 2:nrow(local_minima_150_300)) {
      # Compare the current row's 'length' with the last row's 'length' from xx_filtered_df, and see whether the difference is >6
      if (local_minima_150_300[n,]$length - local_minima_150_300_filtered[nrow(local_minima_150_300_filtered),]$length >= 6) {
        # If the current 'length' is at least >8 than the last row in xx_filtered_df, add the current row to the filtered dataframe
        local_minima_150_300_filtered <- rbind(local_minima_150_300_filtered, local_minima_150_300[n, ])
      } else {
        next
      }
    }
    
    # calculate intervalley distance  
    local_minima_150_300_filtered <- local_minima_150_300_filtered %>%
      calculate_intervalley_dist
    
    # first row always has interpeak distance = NA (because it can´t be compared to one row above) >> add 0 instead of NA
    local_minima_150_300_filtered[1,]$intervalley_dist <- 0
    
    # count number of local maxima and local minima
    number_local_minima_150_300 <- nrow(local_minima_150_300_filtered)
    
    # Calculate the sum of "frequency" values at the filtered local minima positions
    sum_local_minima_150_300 <- sum(local_minima_150_300_filtered$freq)
    # cum_sum_local_minima_150_300 <- sum(local_minima_150_300_filtered$cum_freq)
    
    # Calculate the mean interpeak distance
    mean_intervalley_distance_150_300 <- mean(local_minima_150_300_filtered$intervalley_dist)
    
  } else {
    
    # count number of local maxima and local minima
    number_local_minima_150_300 <- nrow(local_minima_150_300_filtered)
    
    # Calculate the sum of "frequency" values at the filtered local minima positions
    sum_local_minima_150_300 <- NA
    
    # Calculate the mean interpeak distance
    mean_intervalley_distance_150_300 <- NA
    
  }
  
  
  # Calculate the cum freq/density between 30bp and 150bp
  cum_freq_150_300 <- data_150_300[data_150_300$length == 300,]$cum_freq - data_150_300[data_150_300$length ==150,]$cum_freq
  
  # Subtract the sum of the minima from the sum of the height of the maxima ( = oscillation factor). The larger this difference, the more distinct are the peaks.
  oscillation_factor_freq_150_300 <- sum_local_maxima_150_300 - sum_local_minima_150_300
  # oscillation_factor_cum_freq_150_300 <- cum_sum_local_maxima_150_300 - cum_sum_local_minima_150_300
  
  # Create two dataframes with locations for local maxima and local minima for each sample
  Sample_name <-gsub(".final.txt","", files[i])
  
  data_local_maxima_150_300 <- data.frame(t(select(local_maxima_150_300_filtered, length)))
  rownames(data_local_maxima_150_300) <- Sample_name
  colnames(data_local_maxima_150_300) <- paste0("peak", seq_along(data_local_maxima_150_300))
  
  data_local_minima_150_300 <- data.frame(t(select(local_minima_150_300_filtered, length)))
  rownames(data_local_minima_150_300) <- Sample_name
  colnames(data_local_minima_150_300) <- paste0("valley", seq_along(data_local_minima_150_300))
  
  if(i == 1) {
    data_local_maxima_table_150_300 <- data_local_maxima_150_300
    data_local_minima_table_150_300 <- data_local_minima_150_300
  }
  if(i != 1) {
    data_local_maxima_table_150_300 <- dplyr::bind_rows(data_local_maxima_table_150_300, data_local_maxima_150_300)
    data_local_minima_table_150_300 <- dplyr::bind_rows(data_local_minima_table_150_300, data_local_minima_150_300)
  }
  
  ############################## Combine all features ##############################
  if (i == 1){
    
    features_final <- data.frame(Sample_ID = basename(gsub(".final.txt","", files[i])),
                           max_frequency = max_frequency, 
                           max_frequency_cumfreq =  max_frequency_cumfreq,
                           mode_fragment_length = mode_fragment_length, 
                           mean_fragment_length = mean_fragment_length, 
                           rel_freq_167bp = rel_freq_167bp, 
                           rel_freq_334bp = rel_freq_334bp, 
                           cumfreq05_fragment_length = cumfreq05_fragment_length, 
                           P30_60 = P30_60, 
                           P30_100 = P30_100, 
                           P30_150 = P30_150, 
                           P30_180 = P30_180,
                           P90_150 = P90_150,
                           P160_180 = P160_180, 
                           P163_169 = P163_169,
                           P180_220 = P180_220, 
                           P150_300 = P150_300, 
                           P324_344 = P324_344,
                           P250_320 = P250_320, 
                           P250_420 = P250_420, 
                           P420_700 = P420_700, 
                           ratio_P30_150_P160_180 = ratio_P30_150_P160_180, 
                           ratio_P30_150_P163_169 = ratio_P30_150_P163_169, 
                           ratio_P30_100_P160_180 = ratio_P30_100_P160_180,
                           ratio_P30_100_P163_169 = ratio_P30_100_P163_169, 
                           ratio_P30_150_P150_300 = ratio_P30_150_P150_300, 
                           ratio_P90_150_P163_169 = ratio_P90_150_P163_169, 
                           ratio_P160_180_P250_420 = ratio_P160_180_P250_420,
                           number_local_maxima_30_150 = number_local_maxima_30_150 , 
                           number_local_minima_30_150 = number_local_minima_30_150 , 
                           oscillation_factor_freq_30_150 = oscillation_factor_freq_30_150,  
                           mean_interpeak_distance_30_150 = mean_interpeak_distance_30_150, 
                           mean_intervalley_distance_30_150 = mean_intervalley_distance_30_150, 
                           number_local_maxima_150_300 = number_local_maxima_150_300 , 
                           number_local_minima_150_300 = number_local_minima_150_300 , 
                           oscillation_factor_freq_150_300 = oscillation_factor_freq_150_300, 
                           mean_interpeak_distance_150_300 = mean_interpeak_distance_150_300, 
                           mean_intervalley_distance_150_300 = mean_intervalley_distance_150_300 
                           )
  }
  
  if (i != 1){
    
    features_final <- rbind(features_final, data.frame(Sample_ID = basename(gsub(".final.txt","", files[i])), 
                                           max_frequency = max_frequency, 
                                           max_frequency_cumfreq =  max_frequency_cumfreq,
                                           mode_fragment_length = mode_fragment_length, 
                                           mean_fragment_length = mean_fragment_length, 
                                           rel_freq_167bp = rel_freq_167bp, 
                                           rel_freq_334bp = rel_freq_334bp, 
                                           cumfreq05_fragment_length = cumfreq05_fragment_length, 
                                           P30_60 = P30_60, 
                                           P30_100 = P30_100, 
                                           P30_150 = P30_150, 
                                           P30_180 = P30_180,
                                           P90_150 = P90_150,
                                           P160_180 = P160_180, 
                                           P163_169 = P163_169,
                                           P180_220 = P180_220, 
                                           P150_300 = P150_300, 
                                           P324_344 = P324_344,
                                           P250_320 = P250_320, 
                                           P250_420 = P250_420, 
                                           P420_700 = P420_700, 
                                           ratio_P30_150_P160_180 = ratio_P30_150_P160_180, 
                                           ratio_P30_150_P163_169 = ratio_P30_150_P163_169, 
                                           ratio_P30_100_P160_180 = ratio_P30_100_P160_180,
                                           ratio_P30_100_P163_169 = ratio_P30_100_P163_169, 
                                           ratio_P30_150_P150_300 = ratio_P30_150_P150_300, 
                                           ratio_P90_150_P163_169 = ratio_P90_150_P163_169, 
                                           ratio_P160_180_P250_420 = ratio_P160_180_P250_420,
                                           number_local_maxima_30_150 = number_local_maxima_30_150 , 
                                           number_local_minima_30_150 = number_local_minima_30_150 , 
                                           oscillation_factor_freq_30_150 = oscillation_factor_freq_30_150,  
                                           mean_interpeak_distance_30_150 = mean_interpeak_distance_30_150, 
                                           mean_intervalley_distance_30_150 = mean_intervalley_distance_30_150, 
                                           number_local_maxima_150_300 = number_local_maxima_150_300 , 
                                           number_local_minima_150_300 = number_local_minima_150_300 , 
                                           oscillation_factor_freq_150_300 = oscillation_factor_freq_150_300, 
                                           mean_interpeak_distance_150_300 = mean_interpeak_distance_150_300, 
                                           mean_intervalley_distance_150_300 = mean_intervalley_distance_150_300 
                                          ))
    
    
  }
  
}

### Annotate dataframe for cfDNA fragmentation features with sample metadata

### load metadata (metadata should include a column "tumorstatus" which defines the tumorstatus of the respective sample (e.g. tumor, control) and a column "group" for the three cohorts: localized_PCa, advanced_PCa, controls)
meta_data_all <- read.csv(file="xxx", header=TRUE, sep=',')

## combine/join "table" and "meta_data_all"
features_final <- full_join(features, meta_data_all, by="Sample_ID")

############## Statistical Overview ################
# calculate median and range (min, max) for each fragmentation feature and each tumorstatus (tumor, control)
summary_stats_tumorstatus <- select(features_final, c(max_frequency:mean_intervalley_distance_150_300, tumorstatus)) %>%
  group_by(tumorstatus) %>%
  summarise(across(max_frequency:mean_intervalley_distance_150_300, list( min = ~min(.x, na.rm = TRUE), 
                                                                          max = ~max(.x, na.rm = TRUE), 
                                                                          mean = ~mean(.x, na.rm = TRUE),
                                                                          median = ~median(.x, na.rm = TRUE))))

############## Statistical testing ################
# List of unique features
features_extracted <- names(select(features_final, c(max_frequency:mean_intervalley_distance_150_300)))

####
# Step 1: Perform Wilcoxon test for each feature between pairs of groups (tumor vs. control)
wilcox_results_tumorstatus <- list()

for (feature in features_extracted) {
  wilcox_results_tumorstatus[[feature]] <- list(
    tumor_vs_control_pvalue = wilcox.test(features_final[[feature]][features_final$tumorstatus == "tumor"], 
                                          features_final[[feature]][features_final$tumorstatus == "control"], paired = FALSE)$p.value
                                          )
}

# Convert the list to a data frame
wilcox_results_tumorstatus_df <- data.frame(do.call(rbind.data.frame, wilcox_results_tumorstatus))

# Ajdust for multiple testing with Benjamini Hochberg Method
wilcox_results_tumorstatus_df$p_adjust_BH <- p.adjust(wilcox_results_tumorstatus_df$tumor_vs_control_pvalue, method = "BH")

####
# Step 2: Perform Kruskal-Wallis test for each feature across several groups (localized cancer, advanced cancer, controls)

kruskal_results <- numeric(length(features_extracted))
names(kruskal_results) <- features_extracted

for (feature in features_extracted) {
  kruskal_results[feature] <- kruskal.test(features_final[[feature]] ~ features_final$group)$p.value
}
kruskal_results_df <- as.data.frame(kruskal_results)

####
# Step 3: Perform PostHoc testing (Dunn`s test with Benjamini Hochberg for adjustment of multiple testing`)
dunn_result <- data.frame()

for (feature in features_extracted) {
  dunn <- features_final %>% dunn_test(as.formula(paste(feature, "~ group")), p.adjust.method = "BH")
  
  dunn_result <- rbind(dunn_result, dunn)
}


#######################################################################################################################################
## Plotting
#######################################################################################################################################
## (A) Different proportions of cfDNA fragment length ranges (in relation to all fragments with 30–700 bp length) in all tumor and control samples. 
# Box plots: center lines indicate the median, and boxes illustrate the interquartile range with Tukey whiskers.

# transform to longer-format
proportions_longer <- features_final %>% pivot_longer(cols = 9:21, names_to = "fragment_length_range", values_to = "proportion" )

#### Ploting with ggplot ####
levels_proportion = c( "P30_60", "P30_100", "P30_150" , "P30_180", "P90_150", "P160_180", "P163_169","P180_220", "P150_300","P324_344", "P250_320","P250_420", "P420_700" ) # define levels
base_size = 22 # define factor for sizes in theme()

ggplot(data = proportions_longer, aes(x = factor(fragment_length_range, levels = levels_proportion), y = proportion, fill=tumorstatus)) +
  geom_boxplot() + 
  #geom_point(alpha = 0.2, position = position_jitterdodge(), color = "grey30", size = 2.0) +
  scale_fill_manual(values = c("tumor" = "coral2", "control" = "grey80")) + 
  labs(title= "Proportions of different cfDNA fragment length ranges",
       x = "fragment length range", y = "proportion of total") +
  scale_x_discrete(guide = guide_axis(angle = 90)) + guides(fill=guide_legend(title="tumorstage")) +
  theme(
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = base_size * 1.2 , lineheight = 0.9, vjust = 0.5, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = base_size * 1.1, lineheight = 0.9, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = base_size * 1.2, angle = 90, vjust = 0.8),
    axis.ticks.length = unit(0.3, "lines"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.border = element_rect(fill = NA, colour = "transparent", linewidth = 0.75),
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(color = "black", linewidth = 0.75, fill = "#9ec5c8"),
    strip.text.x = element_text(size = base_size * 1, angle = 0, vjust = 1),
    plot.background = element_rect(colour = NA, fill = "transparent"),
    plot.title = element_text(size = base_size * 0.8), 
    plot.subtitle = element_text(size = base_size * 0.85),
    legend.position = "right",
    legend.key = element_rect(colour = "transparent", fill = "transparent"),
    legend.title = element_text(size = base_size *1.1),
    legend.key.size = unit(0.55, 'cm'),
    legend.text = element_text(size=base_size * 1.0)
  )

#### 
## (B) Different ratios of proportions of cfDNA fragment length ranges (in relation to all fragments with 30–700 bp length) in all tumor and control samples. 
# Box plots: center lines indicate the median, and boxes illustrate the interquartile range with Tukey whiskers.
 
# transform to longer-format
ratios_longer <- features_final %>% pivot_longer(cols = 22:28, names_to = "ratio_between_fragment_length_ranges", values_to = "ratio" )

#### Ploting with ggplot ####
levels_ratios = c("ratio_P30_100_P160_180" ,"ratio_P30_100_P163_169" , "ratio_P30_150_P160_180","ratio_P30_150_P163_169", "ratio_P30_150_P150_300","ratio_P90_150_P163_169", "ratio_P160_180_P250_420") # define levels
base_size = 22 # define factor for sizes in theme()

ggplot(data = ratios_longer, aes(x = factor(ratio_between_fragment_length_ranges, levels = levels_ratios), y = ratio, fill=tumorstatus)) +
  geom_boxplot() + 
  #geom_point(alpha = 0.2, position = position_jitterdodge(), color = "grey30", size = 2.0) +
  scale_fill_manual(values = c("tumor" = "coral2", "control" = "grey80")) + 
  labs(title= "Ratios of different proportions of cfDNA fragment length ranges", 
       x = "ratio_proportion1_proportion2 [proportion of fragment length ranges to all (30-700bp)]", y = "ratio (proportion 1 / proportion 2) ") +
  scale_x_discrete(guide = guide_axis(angle = 60)) + guides(fill=guide_legend(title="tumorstage")) + 
  theme(
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = base_size * 1.0 , lineheight = 0.9, vjust = 0.5, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = base_size * 1.0, lineheight = 0.9, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = base_size * 1.2, angle = 90, vjust = 0.8),
    axis.ticks.length = unit(0.3, "lines"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.border = element_rect(fill = NA, colour = "transparent", linewidth = 0.75),
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(color = "black", linewidth = 0.75, fill = "#9ec5c8"),
    strip.text.x = element_text(size = base_size * 1, angle = 0, vjust = 1),
    plot.background = element_rect(colour = NA, fill = "transparent"),
    plot.title = element_text(size = base_size * 0.7), 
    plot.subtitle = element_text(size = base_size * 0.85),
    legend.position = "right",
    legend.key = element_rect(colour = "transparent", fill = "transparent"),
    legend.title = element_text(size = base_size *1.1),
    legend.key.size = unit(0.55, 'cm'),
    legend.text = element_text(size=base_size * 1.0)
  )
  

#######################################################################################################################################
## Assessment of 10bp-oscillation patterns
#######################################################################################################################################

### calculate positions of local maxima and minima (30-150bp fragment length range)
data_local_maxima_table_30_150  # dataframe including local maxima in 30-150bp fragment length range
data_local_minima_table_30_150  # dataframe including local minima in 30-150bp fragment length range

round(colMeans(data_local_maxima_table_30_150, na.rm = TRUE, dims = 1),0) # positions of local maxima
round(colMeans(data_local_minima_table_30_150, na.rm = TRUE, dims = 1),0) # positions of local minima

### calculate positions of local maxima and minima (150-300bp fragment length range
data_local_maxima_table_150_300  # dataframe including local maxima in 30-150bp fragment length range
data_local_minima_table_150_300  # dataframe including local minima in 30-150bp fragment length range

round(colMeans(data_local_maxima_table_150_300, na.rm = TRUE, dims = 1),0) # positions of local maxima
round(colMeans(data_local_minima_table_150_300, na.rm = TRUE, dims = 1),0) # positions of local minima


########## Plotting: 10 bp-oscillation score ##########
# Box plots: center lines indicate the median, and boxes illustrate the interquartile range with Tukey whiskers.
 
base_size = 22 # define factor for sizes in theme()

ggplot(data = features_final, aes(x = group, y = oscillation_factor_freq_30_150, fill=tumorstatus)) +
  geom_boxplot(outlier.shape = NA) + geom_point(alpha = 0.7, position = position_jitterdodge(), color = "grey25", size = 2.5) +
  scale_fill_manual(values = c("tumor"  = "coral2", "control"  = "grey80")) + 
  labs(title= "10bp-oscillation score (cfDNA) in tumor and control samples (30-150bp fragment length range)", y = "oscillation score") +
  theme(
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = base_size * 1.45 , lineheight = 0.9, vjust = 0.5, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = base_size * 1.1, lineheight = 0.9, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = base_size * 1.2, angle = 90, vjust = 0.5),
    axis.ticks.length = unit(0.3, "lines"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.border = element_rect(fill = NA, colour = "transparent", linewidth = 0.75),
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(color = "black", linewidth = 0.75, fill = "#9ec5c8"),
    strip.text.x = element_text(size = base_size * 1, angle = 0, vjust = 1),
    plot.background = element_rect(colour = NA, fill = "transparent"),
    plot.title = element_text(size = base_size * 0.65), 
    plot.subtitle = element_text(size = base_size * 0.85),
    legend.position = "none",
    legend.key = element_rect(colour = "transparent", fill = "transparent"),
    legend.title = element_text(size = 20),
    legend.key.size = unit(0.55, 'cm'),
    legend.text = element_text(size=18)
  )