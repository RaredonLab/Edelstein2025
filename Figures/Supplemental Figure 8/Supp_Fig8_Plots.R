
# close all, clear all
graphics.off()  
rm(list = ls())

library(Seurat)
library(ggplot2)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(cowplot)
library(grid)
library(forcats)
library(patchwork)
library(viridis)  # If not installed: install.packages("viridis")

## MEMO: Supplemental Figure 7: Engineered Object (Global) Cluster Justifications

setwd("~/Desktop/Datasets/Engineered Global Objects")
load("eng.subset.integrated_NodeAligned.Robj")
table(eng.subset.integrated_HK$CellType.NodeAligned)

Cond.cols = c("PD_3D" = "#4B2C77",
              "Mixed_3D" = "#F05C60",
              "BAL_3D" = "#02AEBC")
igf1_split.condition_vln = VlnPlot(eng.subset.integrated_HK, features = c("Igf1"), split.by = "Condition",cols = Cond.cols)
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 8")
ggsave("igf1_split.condition_vln.png", plot = igf1_split.condition_vln, width = 11, height = 5.5, dpi = 600)
igf1_split.condition_vln

Il1b_split.condition_vln = VlnPlot(eng.subset.integrated_HK, features = c("Il1b"), split.by = "Condition",cols = Cond.cols)
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 8")
ggsave("Il1b_split.condition_vln.png", plot = Il1b_split.condition_vln, width = 11, height = 5.5, dpi = 600)
