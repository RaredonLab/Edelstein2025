
# close all, clear all
graphics.off()  
rm(list = ls())
# Load packages
require(Seurat)
require(dplyr)

## UPDATED FOR NODE ALIGNMENT - 05.21.2025

# Load ComplexHeatmap_3 (file is ComplexHeatmap_3.continuous) function
# Also load ComplexHeatmap_2 (for two vars, one being pseudotime continuous)

# Define object as whatever object we are working with
# Load local functions
# source("/Users/sophieedelstein/Desktop/SE Single Cell/Organoids_BASC_Mono/CustomHeatmapOnly3.R")
# Or you can source it locally

# Set working directory to be true to get object
setwd("~/Desktop/Datasets/Engineered Global Objects")
load("eng.subset.integrated_NodeAligned.Robj")

str(eng.subset.integrated_HK@meta.data)
table(eng.subset.integrated_HK$CellType.NodeAligned)

# Define object as whatever object we are working with
object = eng.subset.integrated_HK
# Inspect data to get a handle on the metadata slot names
names(eng.subset.integrated_HK@meta.data)

# Look at the ordering of different metadata slots
table(eng.subset.integrated_HK$CellClass.NodeAligned)
table(eng.subset.integrated_HK$CellType.NodeAligned)
table(eng.subset.integrated_HK$Condition)

# Re-order the metadata for plotting
eng.subset.integrated_HK$CellClass.NodeAligned <- factor(eng.subset.integrated_HK$CellClass.NodeAligned,
                                                                 levels = c('Epithelium')) # Modify the order here if needed
eng.subset.integrated_HK$CellType.NodeAligned <- factor(eng.subset.integrated_HK$CellType.NodeAligned,
                                                              levels = c('Secretory','RAS_Like','Stressed_Progenitor')) # Modify the levels here
eng.subset.integrated_HK$Condition <- factor(eng.subset.integrated_HK$Condition,
                                             levels = c('BAL_3D','Mixed_3D','PD_3D')) # Modify the levels here
# Look at the effect of re-ordering
table(eng.subset.integrated_HK$CellType.NodeAligned, eng.subset.integrated_HK$CellClass.NodeAligned, eng.subset.integrated_HK$Condition)

# Color palette - pulling these from ColorPalette_OrganoidProj.R
col.pal <- list()
col.pal$CellClass.NodeAligned <- c('#D982C6') #SE Custom
names(col.pal$CellClass.NodeAligned) <- c('Epithelium')
col.pal$CellType.NodeAligned <- c('#2E8B57','#595db0','#c71585','#DB7093','#4682B4','#b22222')
names(col.pal$CellType.NodeAligned) <- c('Secretory','RAS_Like','Stressed_Progenitor')
col.pal$Condition <- c('#00b0be','#ff585e','#4a2377') #SE Custom
names(col.pal$Condition) <- c('BAL_3D','PD_3D','Mixed_3D')

# Scale the object
# downsampled <- ScaleData(downsampled, features = rownames(downsampled))
# Create a marker list
Idents(eng.subset.integrated_HK) <- eng.subset.integrated_HK$CellType
mark <- FindAllMarkers(eng.subset.integrated_HK,
                       only.pos = FALSE,
                       min.pct = 0.1,
                       logfc.threshold = 0.1)
mark$ratio <- mark$pct.1/mark$pct.2
mark$power <- mark$ratio * mark$avg_log2FC
# Look at marker list
View(mark)

# Always remember to scale data or code will prompt a bug
eng.subset.integrated_HK <- ScaleData(eng.subset.integrated_HK,features = rownames(eng.subset.integrated_HK))

# Define the cell types you want to include in the heatmap
selected_celltypes <- c("Stressed_Progenitor", "RAS_Like", "Secretory")  # example
# Subset the Seurat object
subset_obj <- subset(eng.subset.integrated_HK, subset = CellType.NodeAligned %in% selected_celltypes)
# Define matching color palette
subset_palette <- col.pal$CellType[names(col.pal$CellType.NodeAligned) %in% selected_celltypes]

GOI <- c('Sox9','Tmem255b','Rasd2','Cox4i2','Stc1','Espn','Apln','Icam5','Ptges',
         'Scgb3a1','Agr2','Ndst3','Muc5b','Bpifb1','Bpifa5','Scgb3a2','Scgb1a1',
         'Sftpc','Lamp3','Napsa')
# Define row labels (these are what we want to show labels for)
row.labels <- c('Sox9','Tmem255b','Rasd2','Cox4i2','Stc1','Espn','Apln','Icam5','Ptges',
                'Scgb3a1','Agr2','Ndst3','Muc5b','Bpifb1','Bpifa5','Scgb3a2','Scgb1a1',
                'Sftpc','Lamp3','Napsa')
# Now plot
png(file = 'plot10.png', width = 11, height = 4.5, units = 'in', res = 600)
ComplexHeatMap_2.inferno(
  object = subset_obj,
  data.type = 'RNA',
  primary = 'CellType.NodeAligned',
  secondary = 'Condition',
  primary.cols = col.pal$CellType.NodeAligned,
  secondary.cols = col.pal$Condition,
  features = GOI,
  labels = c('Cell Type', 'Condition'),
  selected.row.anotations = row.labels,
  selected.label.size = 8,
  use.scale.data = TRUE,
  range.frac = 0.5,
  row.dendrogram = FALSE)
dev.off()


############################## USING ENTIRE INTEGRATED OBJECT ######################################
# Set working directory to be true to get object
setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets/")
load("eng.subset.integrated_HK.Robj")

# Define object as whatever object we are working with
object = eng.subset.integrated_HK
# Inspect data to get a handle on the metadata slot names
names(eng.subset.integrated_HK@meta.data)

# Look at the ordering of different metadata slots
table(eng.subset.integrated_HK$CellClass.NodeAligned)
table(eng.subset.integrated_HK$CellType.NodeAligned)
table(eng.subset.integrated_HK$Condition)

# Re-order the metadata for plotting
eng.subset.integrated_HK$CellClass.NodeAligned <- factor(eng.subset.integrated_HK$CellClass.NodeAligned,
                                                                 levels = c('Epithelium')) # Modify the order here if needed
eng.subset.integrated_HK$CellType.NodeAligned <- factor(eng.subset.integrated_HK$CellType.NodeAligned,
                                                                 levels = c('Secretory','RAS_Like','Stressed_Progenitor')) # Modify the levels here
eng.subset.integrated_HK$Condition <- factor(eng.subset.integrated_HK$Condition,
                                             levels = c('BAL_3D','Mixed_3D','PD_3D')) # Modify the levels here
# Look at the effect of re-ordering
table(eng.subset.integrated_HK$CellClass.NodeAligned, eng.subset.integrated_HK$CellType.NodeAligned, eng.subset.integrated_HK$Condition)

# Color palette - pulling these from ColorPalette_OrganoidProj.R
col.pal <- list()
col.pal$CellClass.NodeAligned <- c('#D982C6') #SE Custom
names(col.pal$CellClass.NodeAligned) <- c('Epithelium')
col.pal$CellType.NodeAligned <- c('#2E8B57','#595db0','#c71585')
names(col.pal$CellType.NodeAligned) <- c('Secretory','RAS_Like','Stressed_Progenitor')
col.pal$Condition <- c('#00b0be','#ff585e','#4a2377') #SE Custom
names(col.pal$Condition) <- c('BAL_3D','PD_3D','Mixed_3D')

# Scale the object
# downsampled <- ScaleData(downsampled, features = rownames(downsampled))
# Create a marker list
Idents(eng.subset.integrated_HK) <- eng.subset.integrated_HK$CellClass.combined.Integrated
mark <- FindAllMarkers(eng.subset.integrated_HK,
                       only.pos = FALSE,
                       min.pct = 0.1,
                       logfc.threshold = 0.1)
mark$ratio <- mark$pct.1/mark$pct.2
mark$power <- mark$ratio * mark$avg_log2FC
# Look at marker list
View(mark)

# Always remember to scale data or code will prompt a bug
eng.subset.integrated_HK <- ScaleData(eng.subset.integrated_HK,features = rownames(eng.subset.integrated_HK))

# Define the cell types you want to include in the heatmap
# selected_celltypes <- c("Stressed_Progenitor", "RAS_Like", "Secretory")  # example
# Subset the Seurat object
# subset_obj <- subset(eng.subset.integrated_HK, subset = CellType %in% selected_celltypes)
# Define matching color palette
# subset_palette <- col.pal$CellType[names(col.pal$CellType) %in% selected_celltypes]

GOI <- c('Sox9','Tmem255b','Rasd2','Cox4i2','Stc1','Espn','Apln','Icam5','Ptges',
         'Scgb3a1','Agr2','Ndst3','Muc5b','Bpifb1','Bpifa5','Scgb3a2','Scgb1a1',
         'Sftpc','Lamp3','Napsa')
# Define row labels (these are what we want to show labels for)
row.labels <- c('Sox9','Tmem255b','Rasd2','Cox4i2','Stc1','Espn','Apln','Icam5','Ptges',
                'Scgb3a1','Agr2','Ndst3','Muc5b','Bpifb1','Bpifa5','Scgb3a2','Scgb1a1',
                'Sftpc','Lamp3','Napsa')
# Now plot
setwd("~/Desktop")
png(file = 'plot10.png', width = 10, height = 4, units = 'in', res = 600)
ComplexHeatMap_2.inferno(object = eng.subset.integrated_HK,
                         data.type = 'RNA',
                         primary = 'CellType.NodeAligned',
                         secondary = 'Condition',
                         primary.cols = col.pal$CellType.NodeAligned,
                         secondary.cols = col.pal$Condition,
                         features = GOI,
                         labels = c('Cell Type', 'Condition'),
                         selected.row.anotations = row.labels,
                         selected.label.size = 8,
                         use.scale.data = TRUE,
                         range.frac = 0.5,
                         row.dendrogram = FALSE)
dev.off()
