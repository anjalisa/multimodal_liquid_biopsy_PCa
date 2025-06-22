#
# Script Name: cfDNA_fragment_distribution.R
#
# Author: Anja Riediger (DKFZ Heidelberg)
# Last updated: 2025/06/17
#
# Description
#   Assessment of cfDNA fragmentation based on relative and cummulative frequency distributions of insert sizes 
#   Comparison of cumulative frequency distributions between tumor and control samples with Kolmogorov-Smirnov Testing
# 
# Input
# insert size frequency tables: previoulsy created with R-script "Fragment_length_FreqTable.R" 
# ----------------------------------------------------------------------------------------------------------

# load R packages
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(stringr)
library(dplyr)                        
library(readr)                        
library(tidyverse)                   
library(plyranges)
library(readxl)
library(rstatix)
library(matrixStats) 

#####################
# define directory for output/results
outdir_results <- "/path/to/output/directory"

# define intput; load insert size frequency tables previously created with R script "Fragment_length_FreqTable.R"
path <- "/path/to/frequency_tables"
files <- list.files(path)

#####
# Load insert size frequency tables and combine data from all samples into one dataframe
for(i in 1:length(files)){
  
  # Load data
  data <- data.frame(fread(file.path(path, files[i])))
  
  # Change sample ID information
  data$Sample_ID <- gsub(".final.txt","", files[i])
  
  # Get relative frequency
  data$freq <- data$freq / sum(data$freq, na.rm = TRUE)
  
  # add cumulative frequency
  table_cum_freq <- data %>%                             # Apply group_by & mutate functions
    group_by(Sample_ID) %>%
    dplyr::mutate(cum_freq = cumsum(freq))
  
  table_cum_freq <- as.data.frame(table_cum_freq)
  
  # Append dataframe with fragment length, frequency
  if(i == 1) table <- table_cum_freq
  if(i != 1) table <- rbind(table, table_cum_freq)
  
}

### Annotate dataframe with sample metadata
# load metadata (metadata should include a column "group" which defines the tumorstatus of the respective sample (e.g. localized prostate cancer (PCa), advanced PCa, control))
meta_data_all <- read.csv(file="xxx", header=TRUE, sep=',')

# combine/join "table" and "meta_data_all"
table_annotated <- full_join(table, meta_data_all, by = "Sample_ID")

#### Extract the median relative frequency distribution of cfDNA fragments in tumor and control samples, respectively
## Calculate medians for each combination of "length" and "group" and transform data table to wider-format   
table_annotated_wider_summary_group <- table_annotated  %>%
  group_by(length, group) %>%
  summarise(median_freq = median(freq, na.rm = TRUE)) %>%
  pivot_wider(names_from = group, values_from = median_freq)


##############
###### Plotting: median relative frequency distributions of tumor and control samples #############

base_size = 22 # define factor for sizes in theme()

# relative Frequency / genomewide cfDNA fragmentation plots: Create a ggplot curve with median values
ggplot(table_annotated, aes(x = length, y = freq, group = group, color = factor(group, levels = c("control", "localized_PCa", "advanced_PCa")))) +
  stat_summary(fun = "median", geom = "line", linewidth = 1.5, alpha = 0.6) + scale_color_manual(values = c(control = "darkgreen", localized_PCa = "blue3", advanced_PCa = "darkred")) +
  labs(x = "fragment length [bp]", y = "relative frequency", title = "CfDNA fragmentation profiles for controls, localized PCa and advanced PCa patients", subtitle= "median relative frequency (referred to 30-700bp) for each group") + 
  guides(color=guide_legend(title="tumorstage")) + coord_cartesian(xlim = c(30,700)) +
    geom_vline(xintercept = c(167, 167*2, 167*3), color = "#525252", linetype = "dashed")  +
    scale_x_continuous(limits = c(30,700), breaks = c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700)) + 
    # scale_x_continuous(limits = c(30,250), breaks = c(50, 100, 150, 200, 250)) +                                  # plot cfDNA fragment distribution in 30-250 bp length range
    # scale_x_continuous(limits = c(250,700), breaks = c(250,300, 350, 400, 450, 500, 550, 600, 650, 700)) +        # plot cfDNA fragment distribution in 250-700 bp length range
    # scale_x_continuous(limits = c(30,150), breaks = c(30,40,50,60,70,80,90,100,110,120,130,140,150)) +            # plot cfDNA fragment distribution in 30-150 bp length range
  theme(
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = base_size * 1.2 , lineheight = 0.9, vjust = 0.5, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = base_size * 1.2, lineheight = 0.9, hjust = 1),
    axis.title.x = element_text(size = base_size * 1.4, angle = 0, vjust = 0.8),
    axis.title.y = element_text(size = base_size * 1.4, angle = 90, vjust = 0.8),
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
    plot.subtitle = element_text(size = base_size * 0.75),
    legend.position = "top",  legend.margin=margin(), legend.box="vertical", legend.text.align = 0, 
    legend.justification = c(0.05, 0.5),
    legend.box.just = "left",  # Align items inside the legend box to the left
    legend.spacing.x = unit(0.4, 'cm'),  # Adjust horizontal spacing
    legend.box.spacing = unit(0.2, 'cm'),  # Adjust space between legend box and plot
    legend.key = element_rect(colour = "transparent", fill = "transparent"),
    legend.title = element_text(size = base_size *1.2),
    legend.key.size = unit(0.55, 'cm'),
    legend.text = element_text(size=base_size * 1.1)
  )

  

### comparison between plasma and urinary cfDNA fragmentation
# relative Frequency / genomewide cfDNA fragmentation plots: Create a ggplot curve with median values
base_size = 22 # define factor for sizes in theme()

ggplot() +
  stat_summary(data =table_annotated[table_annotated$source == "plasma",], aes(x = length, y = freq, group = group, color = factor(group, levels = c("control", "localized_PCa", "advanced_PCa")), linetype = "plasma"), 
               fun = "median", geom = "line", linewidth = 1.2, alpha = 0.6) + 
  stat_summary(data =table_annotated[table_annotated$source == "urine",], aes(x = length, y = freq, group = group, color = factor(group, levels = c("control", "localized_PCa", "advanced_PCa")), linetype = "urine"), 
               fun = "median", geom = "line", linewidth = 1.2, alpha = 0.6) + scale_color_manual(values = c(control = "darkgreen", localized_PCa = "blue3", advanced_PCa = "darkred")) +
  labs(x = "fragment length [bp]", y = "relative frequency", title = "Fragmentation profiles of plasma and urinary cfDNA for controls, localized PCa and advanced PCa patients", subtitle= "median relative frequency (referred to 30-700bp) for each group") + 
  guides(color=guide_legend(title="tumorstage")) + coord_cartesian(xlim = c(30,700)) + scale_linetype_manual(name = "liquid biopsy source", values = c("plasma" = "solid", "urine" = "twodash")) +
  geom_vline(xintercept = c(167, 167*2, 167*3), color = "#525252", linetype = "dashed")  +
  scale_x_continuous(limits = c(30,700), breaks = c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700)) + 
  # scale_x_continuous(limits = c(30,250), breaks = c(50, 100, 150, 200, 250)) +                                  # plot cfDNA fragment distribution in 30-250 bp length range
  # scale_x_continuous(limits = c(250,700), breaks = c(250,300, 350, 400, 450, 500, 550, 600, 650, 700)) +        # plot cfDNA fragment distribution in 250-700 bp length range
  # scale_x_continuous(limits = c(30,150), breaks = c(30,40,50,60,70,80,90,100,110,120,130,140,150)) +            # plot cfDNA fragment distribution in 30-150 bp length range
  theme(
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = base_size * 1.2 , lineheight = 0.9, vjust = 0.5, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = base_size * 1.2, lineheight = 0.9, hjust = 1),
    axis.title.x = element_text(size = base_size * 1.4, angle = 0, vjust = 0.8),
    axis.title.y = element_text(size = base_size * 1.4, angle = 90, vjust = 0.8),
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
    plot.subtitle = element_text(size = base_size * 0.75),
    legend.position = "top",  legend.margin=margin(), legend.box="vertical", legend.text.align = 0, 
    legend.justification = c(0.05, 0.5),
    legend.box.just = "left",  # Align items inside the legend box to the left
    legend.spacing.x = unit(0.4, 'cm'),  # Adjust horizontal spacing
    legend.box.spacing = unit(0.2, 'cm'),  # Adjust space between legend box and plot
    legend.key = element_rect(colour = "transparent", fill = "transparent"),
    legend.title = element_text(size = base_size *1.2),
    legend.key.size = unit(0.55, 'cm'),
    legend.text = element_text(size=base_size * 1.1)
  )

##################################################################
### Perfom Kolmogorov Smirnov Testing with cumulative Frequencies
##################################################################
## annotated data tables were already created above

#### Test for differences in the distributions, using Kolmogorov-Smirnov test
## Calculate median for each combination of "length" and "group" (control, localized PCa, advanced PCa) and transform data table to wider_format   
table_annotated_wider_cum_freq_median_group <- table_annotated  %>%
  group_by(length, group) %>%
  summarise(median_cum_freq = median(cum_freq, na.rm = TRUE)) %>%
  pivot_wider(names_from = group, values_from = median_cum_freq)

### Calculate median for each combination of "length" and "group" (control, tumor) and transform data table to wider_format   
table_annotated_wider_cum_freq_median_tumorstatus <- table_annotated  %>%
  group_by(length, tumorstatus) %>%
  summarise(median_cum_freq = median(cum_freq, na.rm = TRUE)) %>%
  pivot_wider(names_from = tumorstatus, values_from = median_cum_freq)


# Perform Kolmogorov-Smirnov testing to compare cumulative frequency distributions of cfDNA fragmentation between tumor and control samples (for different fragment length ranges)
ks_test_result_30_700bp_control_tumor_cum_freq <- ks.test(table_annotated_wider_cum_freq_median_tumorstatus$control, table_annotated_wider_cum_freq_median_tumorstatus$tumor)
ks_test_result_30_150bp_control_tumor_cum_freq <- ks.test(table_annotated_wider_cum_freq_median_tumorstatus[table_annotated_wider_cum_freq_median_tumorstatus$length <=150,]$control, table_annotated_wider_cum_freq_median_tumorstatus[table_annotated_wider_cum_freq_median_tumorstatus$length <=150,]$tumor)
ks_test_result_30_250bp_control_tumor_cum_freq <- ks.test(table_annotated_wider_cum_freq_median_tumorstatus[table_annotated_wider_cum_freq_median_tumorstatus$length <=250,]$control, table_annotated_wider_cum_freq_median_tumorstatus[table_annotated_wider_cum_freq_median_tumorstatus$length <=250,]$tumor)
ks_test_result_250_700bp_control_tumor_cum_freq <- ks.test(table_annotated_wider_cum_freq_median_tumorstatus[table_annotated_wider_cum_freq_median_tumorstatus$length > 250,]$control, table_annotated_wider_cum_freq_median_tumorstatus[table_annotated_wider_cum_freq_median_tumorstatus$length >250,]$tumor)

##############
###### Plotting: median cumulative frequency distributions of tumor and control samples #############
base_size = 22    # define factor for sizes in theme()

# cumulative Frequency / genomewide Fragmentation plots: Create a ggplot curve with median values
ggplot(table_annotated[table_annotated$LIBRARY_TYPE =="WGS" & table_annotated$pre_post == "pre" & table_annotated$source == "plasma",], aes(x = length, y = cum_freq, group = group_extended, color = factor(group_extended, levels = c("control", "localized_PCa", "advanced_PCa")))) +
  stat_summary(fun = "median", geom = "line", linewidth = 1.5, alpha = 0.6) + scale_color_manual(values = c(control = "darkgreen", localized_PCa = "blue3", advanced_PCa = "darkred")) +
  labs(x = "fragment length [bp]", y = "cumulative frequency", title = "CfDNA fragmentation profiles for controls, localized PCa and advanced PCa patients", subtitle= "median cumulative frequency (referred to 30-700bp) for each group") + 
  guides(color=guide_legend(title="tumorstage")) + coord_cartesian(xlim = c(30,700)) +
  geom_vline(xintercept = c(167, 167*2, 167*3), color = "#525252", linetype = "dashed")  +
  scale_x_continuous(limits = c(30,700), breaks = c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700)) + 
  # scale_x_continuous(limits = c(30,250), breaks = c(50, 100, 150, 200, 250)) +                                  # plot cfDNA fragment distribution in 30-250 bp length range
  # scale_x_continuous(limits = c(250,700), breaks = c(250,300, 350, 400, 450, 500, 550, 600, 650, 700)) +        # plot cfDNA fragment distribution in 250-700 bp length range
  # scale_x_continuous(limits = c(30,150), breaks = c(30,40,50,60,70,80,90,100,110,120,130,140,150)) +            # plot cfDNA fragment distribution in 30-150 bp length range
  theme(
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = base_size * 1.2 , lineheight = 0.9, vjust = 0.5, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = base_size * 1.2, lineheight = 0.9, hjust = 1),
    axis.title.x = element_text(size = base_size * 1.4, angle = 0, vjust = 0.8),
    axis.title.y = element_text(size = base_size * 1.4, angle = 90, vjust = 0.8),
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
    plot.subtitle = element_text(size = base_size * 0.75),
    legend.position = "top",  legend.margin=margin(), legend.box="vertical", legend.text.align = 0, 
    legend.justification = c(0.05, 0.5),
    legend.box.just = "left",  # Align items inside the legend box to the left
    legend.spacing.x = unit(0.4, 'cm'),  # Adjust horizontal spacing
    legend.box.spacing = unit(0.2, 'cm'),  # Adjust space between legend box and plot
    legend.key = element_rect(colour = "transparent", fill = "transparent"),
    legend.title = element_text(size = base_size *1.2),
    legend.key.size = unit(0.55, 'cm'),
    legend.text = element_text(size=base_size * 1.1)
  )



### comparison between plasma and urinary cfDNA fragmentation
# cumulative frequency / genomewide cfDNA fragmentation plots: Create a ggplot curve with median values
base_size = 22 # define factor for sizes in theme()

ggplot() +
  stat_summary(data =table_annotated[table_annotated$source == "plasma",], aes(x = length, y = cum_freq, group = group, color = factor(group, levels = c("control", "localized_PCa", "advanced_PCa")), linetype = "plasma"), 
               fun = "median", geom = "line", linewidth = 1.2, alpha = 0.6) + 
  stat_summary(data =table_annotated[table_annotated$source == "urine",], aes(x = length, y = cum_freq, group = group, color = factor(group, levels = c("control", "localized_PCa", "advanced_PCa")), linetype = "urine"), 
               fun = "median", geom = "line", linewidth = 1.2, alpha = 0.6) + scale_color_manual(values = c(control = "darkgreen", localized_PCa = "blue3", advanced_PCa = "darkred")) +
  labs(x = "fragment length [bp]", y = "cumulative frequency", title = "Fragmentation profiles of plasma and urinary cfDNA for controls, localized PCa and advanced PCa patients", subtitle= "median cumulative frequency (referred to 30-700bp) for each group") + 
  guides(color=guide_legend(title="tumorstage")) + coord_cartesian(xlim = c(30,700)) + scale_linetype_manual(name = "liquid biopsy source", values = c("plasma" = "solid", "urine" = "twodash")) +
  geom_vline(xintercept = c(167, 167*2, 167*3), color = "#525252", linetype = "dashed")  +
  scale_x_continuous(limits = c(30,700), breaks = c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700)) + 
  # scale_x_continuous(limits = c(30,250), breaks = c(50, 100, 150, 200, 250)) +                                  # plot cfDNA fragment distribution in 30-250 bp length range
  # scale_x_continuous(limits = c(250,700), breaks = c(250,300, 350, 400, 450, 500, 550, 600, 650, 700)) +        # plot cfDNA fragment distribution in 250-700 bp length range
  # scale_x_continuous(limits = c(30,150), breaks = c(30,40,50,60,70,80,90,100,110,120,130,140,150)) +            # plot cfDNA fragment distribution in 30-150 bp length range
   theme(
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = base_size * 1.2 , lineheight = 0.9, vjust = 0.5, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = base_size * 1.2, lineheight = 0.9, hjust = 1),
    axis.title.x = element_text(size = base_size * 1.4, angle = 0, vjust = 0.8),
    axis.title.y = element_text(size = base_size * 1.4, angle = 90, vjust = 0.8),
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
    plot.subtitle = element_text(size = base_size * 0.75),
    legend.position = "top",  legend.margin=margin(), legend.box="vertical", legend.text.align = 0, 
    legend.justification = c(0.05, 0.5),
    legend.box.just = "left",  # Align items inside the legend box to the left
    legend.spacing.x = unit(0.4, 'cm'),  # Adjust horizontal spacing
    legend.box.spacing = unit(0.2, 'cm'),  # Adjust space between legend box and plot
    legend.key = element_rect(colour = "transparent", fill = "transparent"),
    legend.title = element_text(size = base_size *1.2),
    legend.key.size = unit(0.55, 'cm'),
    legend.text = element_text(size=base_size * 1.1)
  )