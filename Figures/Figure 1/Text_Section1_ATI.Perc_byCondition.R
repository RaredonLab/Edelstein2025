
# close all, clear all
graphics.off()  
rm(list = ls())

library(Seurat)
library(ggplot2)
library(SeuratObject)
library(dplyr)
library(tidyr)

## Memo: getting ATI percentages across conditions for our first section (to provide explicit data)

# Set wd
setwd("~/Desktop/Datasets")

# Load data
load("eng.subset.integrated_HK.Robj")

# Look at structure
eng.subset.integrated_HK
str(eng.subset.integrated_HK@meta.data)
table(eng.subset.integrated_HK$CellType.combined.Integrated)
table(eng.subset.integrated_HK$Condition)

# Extract metadata from our obj
meta <- eng.subset.integrated_HK@meta.data

# Calculate % ATI_Like per condition
ati_percent_df <- meta %>%
  group_by(Condition) %>%
  summarise(
    total_cells = n(),
    ati_cells = sum(CellType.combined.Integrated == "ATI_Like"),
    ati_percent = round(100 * ati_cells / total_cells, 1))

# View results
print(ati_percent_df)

# What about for condition level objects?
load("PD_3D.integrated_HK.Robj")
load("Mixed_3D.integrated_HK.Robj")
load("BAL_3D.integrated_HK.Robj")

# Check structure
str(PD_3D.integrated_HK@meta.data)
str(Mixed_3D.integrated_HK@meta.data)
str(BAL_3D.integrated_HK@meta.data)

table(PD_3D.integrated_HK$CellType)
table(Mixed_3D.integrated_HK$CellType)
table(BAL_3D.integrated_HK$CellType)

# Total cell counts (can also be calculated via nrow(object@meta.data))
PD_total <- sum(table(PD_3D.integrated_HK$CellType))         # = 8724
Mixed_total <- sum(table(Mixed_3D.integrated_HK$CellType))   # = 6344
BAL_total <- sum(table(BAL_3D.integrated_HK$CellType))       # = 9676

# Function to calculate % ATI_Like directly from Seurat object
get_ati_percent <- function(seurat_obj, celltype_col = "CellType", ati_label = "ATI_Like") {
  meta <- seurat_obj@meta.data
  total_cells <- nrow(meta)
  ati_cells <- sum(meta[[celltype_col]] == ati_label)
  ati_percent <- round(100 * ati_cells / total_cells, 1)
  return(ati_percent)
}

# Run above function for each condition-level object
pd_ati <- get_ati_percent(PD_3D.integrated_HK)
mixed_ati <- get_ati_percent(Mixed_3D.integrated_HK)
bal_ati <- get_ati_percent(BAL_3D.integrated_HK)

# Combine results
ATI_cross.condition = 
  data.frame(Condition = c("PD_3D", "Mixed_3D", "BAL_3D"), ATI_Like_Percent = c(pd_ati, mixed_ati, bal_ati))
print(ATI_cross.condition)
