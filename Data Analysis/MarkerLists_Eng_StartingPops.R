
# close all, clear all
graphics.off()  
rm(list = ls())

# Load packages
library(NICHES)
library(Seurat)
library(dplyr)
library(patchwork)
library(writexl)

## MEMO: Generating marker lists by cell type for supplemental data

##Set working directory to be true and load data
## First load engineered object
setwd("~/Desktop/Datasets/Engineered Global Objects")
load("eng.subset.integrated_NodeAligned.Robj")

# Now load the starting populations
setwd("~/Desktop/Datasets/Starting Populations")
load("BAL.integrated_03-29-2025.Robj")
load("PD.integrated_03-29-2025.Robj")

# Check structure of all objects
str(eng.subset.integrated_HK@meta.data)
str(BAL.integrated@meta.data)
str(PD.integrated@meta.data)

# Check cell type annotations for each
table(eng.subset.integrated_HK$CellType.NodeAligned)
table(BAL.integrated$CellType)
table(PD.integrated$CellType)

## Create marker list for engineered object (organoids; merged, all conditions)
# Set idents to the CellType.NodeAligned 
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType.NodeAligned
eng.subset.celltype.mark = FindAllMarkers(eng.subset.integrated_HK, min.pct = 0.1,logfc.threshold = 0.1)
eng.subset.celltype.mark$ratio = eng.subset.celltype.mark$pct.1/eng.subset.celltype.mark$pct.2
eng.subset.celltype.mark$power = eng.subset.celltype.mark$ratio*eng.subset.celltype.mark$avg_log2FC
# Set wd to where we want this spreadsheet saved
setwd("~/Desktop/Edelstein2025/Biorxiv Submission/Marker Lists")
write_xlsx(eng.subset.celltype.mark, path = "eng.subset.celltype.mark.xlsx")
# Can save as csv, too (if you want))
write.csv(eng.subset.celltype.mark, file = "eng.subset.celltype.mark.csv", row.names = TRUE)

## Create marker list for starting BAL
# Checking CellType to confirm
table(BAL.integrated$CellType)
# Set idents to CellType
Idents(BAL.integrated) = BAL.integrated$CellType
BAL.celltype.mark = FindAllMarkers(BAL.integrated, min.pct = 0.1,logfc.threshold = 0.1)
BAL.celltype.mark$ratio = BAL.celltype.mark$pct.1/BAL.celltype.mark$pct.2
BAL.celltype.mark$power = BAL.celltype.mark$ratio*BAL.celltype.mark$avg_log2FC
# Set wd to where we want this spreadsheet saved
setwd("~/Desktop/Edelstein2025/Biorxiv Submission/Marker Lists")
write_xlsx(BAL.celltype.mark, path = "BAL.celltype.mark.xlsx")
# Can save as csv, too (if you want)
write.csv(BAL.celltype.mark, file = "BAL.celltype.mark.csv", row.names = TRUE)

## Create marker list for starting Pulmonary Dissociation (PD)
# Checking CellType to confirm
table(PD.integrated$CellType)
# Set idents to CellType
Idents(PD.integrated) = PD.integrated$CellType
PD.celltype.mark = FindAllMarkers(PD.integrated, min.pct = 0.1,logfc.threshold = 0.1)
PD.celltype.mark$ratio = PD.celltype.mark$pct.1/PD.celltype.mark$pct.2
PD.celltype.mark$power = PD.celltype.mark$ratio*PD.celltype.mark$avg_log2FC
# Set wd to where we want this spreadsheet saved
setwd("~/Desktop/Edelstein2025/Biorxiv Submission/Marker Lists")
write_xlsx(PD.celltype.mark, path = "PD.celltype.mark.xlsx")
# Can save as csv, too (if you want)
write.csv(PD.celltype.mark, file = "PD.celltype.mark.csv", row.names = TRUE)



