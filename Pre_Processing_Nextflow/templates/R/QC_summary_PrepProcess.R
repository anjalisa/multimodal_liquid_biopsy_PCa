#
# Script Name: QC_summary_PreProcess.sh <MAIN> <META> <ID> <TYPE>
#
# Author: Florian Janke
# Last updated: 2022/03/02
#
# Description
#   (1) Summarizes qualtiy metrics
#
# Run Information
#   (1) Script is called by SeqData_PreProcess.sh
#
# --------------------------------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(reshape2)
library(matrixStats)

args <- commandArgs()
MAIN <- args[6]
META <- args[7]
ID <- args[8]
TYPE <- args[9]

# Remove last backslash #
MAIN <- substr(MAIN, 1, nchar(MAIN)-1)


## Combine quality metrics ##
# Get list of files to combine #
files <- list.files(file.path(MAIN, "tmp", ID))
files <- files[grep(".tsv", files)]
files <- files[-grep("meta", files)]

# Combine #
for(i in 1:length(files)){
  
  tmp <- data.frame(fread(file.path(MAIN, "tmp", ID, files[i])))
  if(files[i] == "readCount_samtools.tsv") colnames(tmp) <- c("Sample_ID", "unpaired.reads")
  
  if(i == 1){
    
    table <- tmp
    
  } else {
    
    table <- merge(table, tmp, by = "Sample_ID")
    
  }
}

# Save #
if(file.exists(file.path(MAIN, "QC", "summary", "QC_summary.tsv"))){
  
  summary <- data.frame(fread(file.path(MAIN, "QC/summary/QC_summary.tsv")))
  summary <- rbind(summary, table)
  
  
} else {
  
  summary <- table
  
}
fwrite(summary, file = file.path(MAIN, "QC", "summary", "QC_summary.tsv"), sep = "\t")


## Summary table displaying median and range of each quality metric ##
summary.2 <- summary[, -which(colnames(summary) %in% c("unpaired.reads"))]
rownames(summary.2) <- summary.2$Sample_ID
summary.2$Sample_ID <- NULL
summary.2[, c(1:2)] <- summary.2[, c(1:2)] /1000000
summary.2 <- data.frame(t(summary.2))

summary.2 <- data.frame(metric = rownames(summary.2), median = rowMedians(as.matrix(summary.2), na.rm = TRUE),
                        min = rowMins(as.matrix(summary.2), na.rm = TRUE), max = rowMaxs(as.matrix(summary.2), na.rm = TRUE))

# Save summary table #
fwrite(summary.2, file = file.path(MAIN, "QC", "summary", "Metric_summary.tsv"), sep = "\t")


## Generate figures ##
base_size <- 14

# Coverage #
cov <- summary[, which(colnames(summary) %in% c("Sample_ID", "Paired.reads.unfiltered", "Paired.reads.filtered"))]
cov[, c(2:3)] = cov[, c(2:3)] / 1000000
cov$Sample_ID <- factor(cov$Sample_ID, levels = cov[order(cov$Paired.reads.unfiltered), ]$Sample_ID)
cov$Paired.reads.unfiltered <- cov$Paired.reads.unfiltered - cov$Paired.reads.filtered
colnames(cov)[c(2:3)] <- c("Total", "Filtered")
cov <- reshape2::melt(cov)

ggplot(data = cov, aes(fill = variable, y = value, x = Sample_ID)) + 
  geom_bar(position="stack", stat="identity") +
  scale_y_continuous(name = "Paired-reads (x1e06)") +
  scale_x_discrete(name = "") +
  scale_fill_manual(values = c("gray90", "#89a7b2")) +
  coord_flip() +
  theme(
    axis.line =         element_line(colour = "grey20"),
    axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, vjust = 0.5, angle = 0, hjust = 0.5),
    axis.text.y =       element_text(size = base_size * 0.6, lineheight = 0.9, hjust = 1),
    axis.title.x =      element_text(size = base_size * 0.8, lineheight = 0.9, hjust = 0.5),
    axis.title.y =      element_text(size = base_size * 0.9, angle = 90, vjust = 2),
    axis.ticks.length = unit(0.3, "lines"),
    panel.background =  element_rect(fill = "transparent", colour = NA), 
    panel.border =      element_blank(), 
    panel.grid.major =  element_blank(),
    panel.grid.minor =  element_blank(),
    strip.background =  element_blank(), 
    strip.text.x =      element_blank(),
    strip.text.y =      element_blank(),
    strip.text =        element_blank(),
    plot.background =   element_rect(colour = NA, fill = "transparent"),
    plot.title =        element_text(size = base_size * 1),
    legend.position = "right",
    legend.title = element_blank(),
    legend.key = element_blank()
  )


ggsave(
  file.path(MAIN, "QC/summary/figures/Coverage_overview.svg"),
  plot = ggplot2::last_plot(),
  device = "png",
  scale = 1.25,
  width = 14,
  height = 12,
  units = "cm",
  dpi = 300,
  bg = "transparent",
  limitsize = TRUE)

ggsave(
  file.path(MAIN, "QC/summary/figures/Coverage_overview.png"),
  plot = ggplot2::last_plot(),
  device = "png",
  scale = 1.25,
  width = 14,
  height = 12,
  units = "cm",
  dpi = 200,
  bg = "transparent",
  limitsize = TRUE)


if(TYPE == "MeDIP"){
  
  # CpG enrichment #
  enrich.Score <- summary[order(summary$enrichment.score), ]
  enrich.Score$SampleNo <- c(1:nrow(enrich.Score))
  ggplot(data = enrich.Score, aes(x = SampleNo, y = enrichment.score)) +
    geom_point(size = 2) +
    scale_y_continuous(limits = c(0.8, 3), breaks = seq(1, 3, 0.5), name = "CpG enrichment score") +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    theme(
      axis.line =         element_line(colour = "grey20"),
      axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, vjust = 0.5, angle = 0, hjust = 0.5),
      axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, hjust = 1),
      axis.title.x =      element_text(size = base_size * 0.8, lineheight = 0.9, hjust = 0.5),
      axis.title.y =      element_text(size = base_size * 0.9, angle = 90, vjust = 2),
      axis.ticks.length = unit(0.3, "lines"),
      panel.background =  element_rect(fill = "transparent", colour = NA), 
      panel.border =      element_blank(), 
      panel.grid.major =  element_blank(),
      panel.grid.minor =  element_blank(),
      strip.background =  element_blank(), 
      strip.text.x =      element_blank(),
      strip.text.y =      element_blank(),
      strip.text =        element_blank(),
      plot.background =   element_rect(colour = NA, fill = "transparent"),
      plot.title =        element_text(size = base_size * 1),
      legend.position = "bottom"
    )
  
  ggsave(
    file.path(MAIN, "QC/summary/figures/CpG_enrichment.svg"),
    plot = ggplot2::last_plot(),
    device = "png",
    scale = 1.25,
    width = 12,
    height = 10,
    units = "cm",
    dpi = 300,
    bg = "transparent",
    limitsize = TRUE)
  
  ggsave(
    file.path(MAIN, "QC/summary/figures/CpG_enrichment.png"),
    plot = ggplot2::last_plot(),
    device = "png",
    scale = 1.25,
    width = 12,
    height = 10,
    units = "cm",
    dpi = 200,
    bg = "transparent",
    limitsize = TRUE)
  
  
  ## Saturation analysis ##
  summary$read.cov <- summary$read.cov /1000000
  ggplot(data = summary, aes(x = read.cov, y = saturation.r)) +
    geom_point(size = 2) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), name = "Estimated correlation") +
    scale_x_continuous(name = "Number of paired reads (x1e06)") +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    theme(
      axis.line =         element_line(colour = "grey20"),
      axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, vjust = 0.5, angle = 0, hjust = 0.5),
      axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, hjust = 1),
      axis.title.x =      element_text(size = base_size * 0.8, lineheight = 0.9, hjust = 0.5),
      axis.title.y =      element_text(size = base_size * 0.9, angle = 90, vjust = 2),
      axis.ticks.length = unit(0.3, "lines"),
      panel.background =  element_rect(fill = "transparent", colour = NA), 
      panel.border =      element_blank(), 
      panel.grid.major =  element_blank(),
      panel.grid.minor =  element_blank(),
      strip.background =  element_blank(), 
      strip.text.x =      element_blank(),
      strip.text.y =      element_blank(),
      strip.text =        element_blank(),
      plot.background =   element_rect(colour = NA, fill = "transparent"),
      plot.title =        element_text(size = base_size * 1),
      legend.position = "bottom"
    )
  
  
  ggsave(
    file.path(MAIN, "QC/summary/figures/Saturation_overview.svg"),
    plot = ggplot2::last_plot(),
    device = "png",
    scale = 1.25,
    width = 12,
    height = 10,
    units = "cm",
    dpi = 300,
    bg = "transparent",
    limitsize = TRUE)
  
  ggsave(
    file.path(MAIN, "QC/summary/figures/Saturation_overview.png"),
    plot = ggplot2::last_plot(),
    device = "png",
    scale = 1.25,
    width = 12,
    height = 10,
    units = "cm",
    dpi = 200,
    bg = "transparent",
    limitsize = TRUE)
  
}

