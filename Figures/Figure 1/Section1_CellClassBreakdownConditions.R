
# close all, clear all
graphics.off()  
rm(list = ls())

# Packages
library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(scales)
library(dplyr)
library(circlize)
library(cowplot)
library(patchwork)
library(reticulate)
library(stringr)
library(NICHES)

## MEMO: For section 1 of paper; getting exact percentages (that are not shown directly on the stacked bar plots) for lineage in PD and BAL

# Set working directory
setwd("~/Desktop/Datasets/Starting Populations")

# Load data
load("BAL.integrated_03-29-2025.Robj")
load("PD.integrated_03-29-2025.Robj")

# Quickly check structure
str(BAL.integrated@meta.data)
str(PD.integrated@meta.data)
table(BAL.integrated$CellClass_Merged)
table(PD.integrated$CellClass_Merged)

# Get the full set of cell classes from PD
all_classes <- names(table(PD.integrated$CellClass_Merged))

# Calculate proportions for BAL
bal_counts <- table(BAL.integrated$CellClass_Merged)
bal_percent <- round(100 * prop.table(bal_counts), 1)

# Calculate proportions for PD
pd_counts <- table(PD.integrated$CellClass_Merged)
pd_percent <- round(100 * prop.table(pd_counts), 1)

# Align to all_classes (fill missing with NA)
bal_percent_full <- bal_percent[all_classes]
pd_percent_full <- pd_percent[all_classes]

# Replace NAs in BAL (if some classes don't exist) with 0
bal_percent_full[is.na(bal_percent_full)] <- 0

# Combine into one summary data frame
lineage_summary <- data.frame(
  CellClass = all_classes,
  BAL_Percent = as.numeric(bal_percent_full),
  PD_Percent = as.numeric(pd_percent_full))

# View result
print(lineage_summary)

### But if we want stats (mean, standard dev)
# Define function to get % of each class per replicate
get_replicate_percentages <- function(seurat_obj, class_col = "CellClass_Merged", replicate_col = "Orig_ID") {
  meta <- seurat_obj@meta.data
  meta %>%
    group_by(.data[[replicate_col]], .data[[class_col]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(.data[[replicate_col]]) %>%
    mutate(percent = 100 * n / sum(n)) %>%
    ungroup()
}
pd_perc_reps <- get_replicate_percentages(PD.integrated)
bal_perc_reps <- get_replicate_percentages(BAL.integrated)

# Summarize across replicates for PD
pd_summary <- pd_perc_reps %>%
  group_by(CellClass_Merged) %>%
  summarise(
    mean_percent = round(mean(percent), 1),
    sd = round(sd(percent), 1),
    n = n(),
    .groups = "drop")
# Same for BAL
bal_summary <- bal_perc_reps %>%
  group_by(CellClass_Merged) %>%
  summarise(
    mean_percent = round(mean(percent), 1),
    sd = round(sd(percent), 1),
    n = n(),
    .groups = "drop")

# Take a look at summaries
pd_summary
bal_summary
