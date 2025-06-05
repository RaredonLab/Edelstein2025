
## MEMO: Get percent breakdown by cell type in global engineered object

# close all, clear all
graphics.off()  
rm(list = ls())

# Load packages
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

# Set wd to be true and load datasets
setwd("~/Desktop/Datasets/Datasets")
load("eng.subset.integrated_HK.Robj") # Organoids merged

# Set idents to node aligned meta-data slot
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType.NodeAligned
# Extract metadata
meta <- eng.subset.integrated_HK@meta.data

# Calculate total cells per replicate
rep_total <- meta %>%
  group_by(Condition, Orig_ID) %>%
  summarise(TotalCells = n(), .groups = "drop")

#  Count cells per CellType.NodeAligned per replicate
rep_counts <- meta %>%
  group_by(Condition, Orig_ID, CellType.NodeAligned) %>%
  summarise(CellTypeCount = n(), .groups = "drop")

# Join and calculate percent per replicate
rep_percent <- left_join(rep_counts, rep_total, by = c("Condition", "Orig_ID")) %>%
  mutate(Percent = 100 * CellTypeCount / TotalCells)

# Summarize across replicates to get mean ± SD
summary_table <- rep_percent %>%
  group_by(Condition, CellType.NodeAligned) %>%
  summarise(
    MeanPercent = mean(Percent),
    SD = sd(Percent),
    N = n(),
    .groups = "drop")

# View the final summary table
print(summary_table)
View(summary_table)
