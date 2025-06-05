
## MEMO: Making heatmap for select ligands and receptors (Figure 6)

# close all, clear all
graphics.off()  
rm(list = ls())

# Load required packages
require(Seurat)
require(dplyr)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)

# Load ComplexHeatmap_3 (file is ComplexHeatmap_3.continuous) function
# Also load ComplexHeatmap_2 (for two vars, one being pseudotime continuous)

# Define object as whatever object we are working with
# Load local functions
# source("/Users/sophieedelstein/Desktop/SE Single Cell/Organoids_BASC_Mono/CustomHeatmapOnly3.R")
# Or you can source it locally

# Set working directory and load object for CTC-based analysis
setwd("~/Desktop/Manuscript Figures/Figure 6")
load("regen.subset_CTC_bySample.Robj")

# Inspect structure and relevant metadata
regen.circuit.object
str(regen.circuit.object@meta.data)
table(regen.circuit.object$CellType.regen.spec)

# Assign object to working variable
object = regen.circuit.object
# Inspect data to get a handle on the metadata slot names
names(regen.circuit.object@meta.data)

# Look at the ordering of different metadata slots
table(regen.circuit.object$CellType.regen.spec)
table(regen.circuit.object$CellClass.sub.regen)
table(regen.circuit.object$Condition)

# Ensure consistent ordering of metadata variables for plotting
regen.circuit.object$CellClass.sub.regen <- factor(regen.circuit.object$CellClass.sub.regen,
                                                                 levels = c('Epithelium', "Immune","Mesenchyme")) # Modify the order here if needed
regen.circuit.object$CellType.regen.spec <- factor(regen.circuit.object$CellType.regen.spec,
                                                              levels = c('Hillock_Like','Polarized_Mac','Rspo3+_Mes','Pdgfrb+_Pericyte')) # Modify the levels here
regen.circuit.object$Condition <- factor(regen.circuit.object$Condition,
                                             levels = c('BAL_3D','Mixed_3D','PD_3D')) # Modify the levels here
# Look at the effect of re-ordering
table(regen.circuit.object$CellClass.sub.regen, regen.circuit.object$CellType.regen.spec)

# Color palette - pulling these from ColorPalette_OrganoidProj.R
col.pal <- list()
col.pal$CellClass.sub.regen <- c('#D982C6','#87B37A','#F4A261') #SE Custom
names(col.pal$CellClass.sub.regen) <- c('Epithelium', 'Immune','Mesenchyme')
col.pal$CellType.regen.spec <- c('#ff0000','#ff9e80','#89a5a5','#F4AFB4')
names(col.pal$CellType.regen.spec) <- c('Hillock_Like','Polarized_Mac','Rspo3+_Mes','Pdgfrb+_Pericyte')
col.pal$Condition <- c('#00b0be','#ff585e','#4a2377') #SE Custom
names(col.pal$Condition) <- c('BAL_3D','PD_3D','Mixed_3D')

# Create a marker list for ligands 
Idents(regen.circuit.object) <- regen.circuit.object$CellType.regen.spec
mark_ligands <- FindAllMarkers(regen.circuit.object,features = lig.list,
                       only.pos = FALSE,
                       min.pct = 0.1,
                       logfc.threshold = 0.1)
mark_ligands$ratio <- mark_ligands$pct.1/mark$pct.2
mark_ligands$power <- mark_ligands$ratio * mark$avg_log2FC
# Look at marker list
View(mark_ligands)

# Create a marker list for receptors
Idents(regen.circuit.object) <- regen.circuit.object$CellType.regen.spec
mark_receptors <- FindAllMarkers(regen.circuit.object,features = rec.list,
                               only.pos = FALSE,
                               min.pct = 0.1,
                               logfc.threshold = 0.1)
mark_receptors$ratio <- mark_receptors$pct.1/mark$pct.2
mark_receptors$power <- mark_receptors$ratio * mark$avg_log2FC
# Look at marker list
View(mark_receptors)

# Always remember to scale data or code will prompt a bug
regen.circuit.object <- ScaleData(regen.circuit.object,features = rownames(regen.circuit.object))
# Start w ligands for both, then receptors
GOI <- c('C1qb', 'Il1b','Alox5ap','Icam2','Ccl3','Cxcl2','Il18','Igf1', # mac ligands
         'Rspo3','Vcam1','Col4a3','Lama2','Dcn','Cxcl12','Cxcl13','Fgf10','Ncam1','Cxcl6', # mes ligands
          'Bdnf','Col8a1','Col4a1', # lig pericytes
         'Cd4','Tlr7','Csf1r','Cd40','Itgax','Il2rg','Tlr2','Tlr6','Csf3r','Sdc3','Mertk','Asgr2', # mac receptors
         'Pdgfra','Cdh2','Fzd1','Tgfbr3','Ror2','Smo','Plaur', # mes lig
         'Pdgfrb','Aqp1')
# Define row labels (these are what we want to show labels for)
row.labels <- c('C1qb', 'Il1b','Alox5ap','Icam2','Ccl3','Cxcl2', # mac ligands
                'Rspo3','Vcam1','Col4a3','Lama2','Dcn','Cxcl12','Cxcl13','Fgf10','Cxcl6', # mes ligands
                'Bdnf','Col4a1','Apoe', # lig pericytes
                'Cd4','Tlr7','Csf1r','Cd40','Itgax','Il2rg','Tlr2','Csf3r', # mac receptors
                'Pdgfra','Cdh2','Fzd1','Tgfbr3','Smo', # mes lig
                'Pdgfrb','Aqp1')
# Now plot
png(file = 'Fig6_Lig_Receptor_Heatmap.png', width = 7.5, height = 9, units = 'in', res = 600)
ComplexHeatMap_3.inferno(
  object = regen.circuit.object,
  data.type = 'RNA',
  primary = 'CellClass.sub.regen',
  secondary = 'CellType.regen.spec',
  tertiary = "Condition",
  primary.cols = col.pal$CellClass.sub.regen,
  secondary.cols = col.pal$CellType.regen.spec,
  tertiary.cols = col.pal$Condition,
  features = GOI,
  labels = c('Cell Class', 'Cell Type', 'Condition'),
  selected.row.anotations = row.labels,
  selected.label.size = 8,
  use.scale.data = TRUE,
  range.frac = 0.5,
  row.dendrogram = FALSE)
dev.off()
