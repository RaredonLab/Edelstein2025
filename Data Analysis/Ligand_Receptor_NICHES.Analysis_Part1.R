
# close all, clear all
graphics.off()  
rm(list = ls())

# Load packages
library(NICHES)
library(Seurat)
library(dplyr)
library(patchwork)
library(writexl)

# Set WD to be true and load NICHES object
setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets/NICHES_Engineering")
load("eng.CTC.final.Robj")
# Fix WD and load the transcriptomic object
setwd("~/Desktop/Datasets/Engineered Global Objects")
load("eng.subset.integrated_NodeAligned.Robj")

# Object
eng.subset.integrated_HK

# Inspect metadata of both objects
str(eng.subset.integrated_HK@meta.data)
str(eng.CTC.final@meta.data)

# Visualize cell class annotations
a = DimPlot(eng.subset.integrated_HK, group.by = "CellClass.integrated", reduction = "umap.rpca")
b = DimPlot(eng.subset.integrated_HK, group.by = "CellClass.combined.Integrated", reduction = "umap.rpca")
a | b

# Load ground truth ligand-receptor database from NICHES (rat-specific)
ligands <- unique(NICHES::ncomms8866_rat$Ligand.ApprovedSymbol)
receptors <- unique(NICHES::ncomms8866_rat$Receptor.ApprovedSymbol)

# Identify which ligands/receptors are present in the dataset
lig.list <- ligands[ligands %in% rownames(eng.subset.integrated_HK)]
rec.list <- receptors[receptors %in% rownames(eng.subset.integrated_HK)]

# Look at both, notice the length has not changed
lig.list
rec.list

# Set identities based on combined integrated cell class
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellClass.combined.Integrated
# Find ligand markers across cell classes
glbl.mark.class.lig = FindAllMarkers(eng.subset.integrated_HK, features = lig.list, 
                                     only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)
# Add custom metrics for visualization or ranking
glbl.mark.class.lig$ratio = glbl.mark.class.lig$pct.1/glbl.mark.class.lig$pct.2
glbl.mark.class.lig$power = glbl.mark.class.lig$ratio*glbl.mark.class.lig$avg_log2FC
# View
View(glbl.mark.class.lig)
# Set wd and then save outputs
setwd("~/Desktop/Edelstein2025/Manuscript Figures/Basal vs. Hillock LR Analysis for MSBR")
write_xlsx(glbl.mark.class.lig, path = "glbl.mark.class.lig.xlsx")
save(glbl.mark.class.lig, file = "glbl.mark.class.lig.Robj")

# Generate marker list with idents set to cell class - for RECEPTORS
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellClass.combined.Integrated
glbl.mark.class.recep = FindAllMarkers(eng.subset.integrated_HK, features = rec.list, 
                                     only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)
# Add custom metrics for visualization or ranking
glbl.mark.class.recep$ratio = glbl.mark.class.recep$pct.1/glbl.mark.class.recep$pct.2
glbl.mark.class.recep$power = glbl.mark.class.recep$ratio*glbl.mark.class.recep$avg_log2FC
View(glbl.mark.class.recep)
# Set wd and then save outputs
setwd("~/Desktop/Edelstein2025/Manuscript Figures/Basal vs. Hillock LR Analysis for MSBR")
write_xlsx(glbl.mark.class.recep, path = "glbl.mark.class.recep.xlsx")
save(glbl.mark.class.recep, file = "glbl.mark.class.recep.Robj")

## Load the proximal subset 
setwd("~/Desktop/Datasets/Subset Population Objects/Proximal Subset")
load("eng.subset.prox.int.Robj")

# Explore metadata and cell types
str(eng.subset.prox.int@meta.data)
DimPlot(eng.subset.prox.int, reduction = "umap.rpca", group.by = "CellType_Prox_Subset")

# Set idents for subset cell types
Idents(eng.subset.prox.int) = eng.subset.prox.int$CellType_Prox_Subset
# Ligand marker analysis (proximal subset)
prox.markers.lig = FindAllMarkers(eng.subset.prox.int, features = lig.list,
                                       only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)
prox.markers.lig$ratio = prox.markers.lig$pct.1/prox.markers.lig$pct.2
prox.markers.lig$power = prox.markers.lig$ratio*prox.markers.lig$avg_log2FC
# View
View(prox.markers.lig)
# Then save lists
setwd("~/Desktop/Edelstein2025/Manuscript Figures/Basal vs. Hillock LR Analysis for MSBR")
write_xlsx(prox.markers.lig, path = "prox.markers.lig.xlsx")
save(prox.markers.lig, file = "prox.markers.lig.Robj")

## Do the same, but not for receptors
Idents(eng.subset.prox.int) = eng.subset.prox.int$CellType_Prox_Subset
prox.markers.recep = FindAllMarkers(eng.subset.prox.int, features = rec.list,
                                  only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)
prox.markers.recep$ratio = prox.markers.recep$pct.1/prox.markers.recep$pct.2
prox.markers.recep$power = prox.markers.recep$ratio*prox.markers.recep$avg_log2FC
View(prox.markers.recep)
# Save
setwd("~/Desktop/Edelstein2025/Manuscript Figures/Basal vs. Hillock LR Analysis for MSBR")
write_xlsx(prox.markers.recep, path = "prox.markers.recep.xlsx")
save(prox.markers.recep, file = "prox.markers.recep.Robj")

### Exploratory FeaturePlots: Ephrin family and Erbb3 ### 
# Cell to cell contact molecules (ephrins, locally attractive/repulsive)
FeaturePlot(eng.subset.prox.int, features = c("Ephb3"), reduction = "umap.rpca")
FeaturePlot(eng.subset.prox.int, features = c("Erbb3"), reduction = "umap.rpca") # receptor
FeaturePlot(eng.subset.integrated_HK, features = "Erbb3", reduction = "umap.rpca", min.cutoff = 1.5, order = T)

# Go look at the Fantom5 database
View(NICHES::ncomms8866_rat)
# Make a list of the mechanisms
mech.list = NICHES::ncomms8866_rat
View(mech.list)

# Filter
mech.list = mech.list[mech.list$Ligand.ApprovedSymbol %in% lig.list &
                      mech.list$Receptor.ApprovedSymbol %in% rec.list, ]
# look at filtered list
View(mech.list)

# Extract ligands binding Erbb3
loi <- mech.list[mech.list$Receptor.ApprovedSymbol == 'Erbb3',]$Ligand.ApprovedSymbol
# View
loi
# Plot each Erbb3-binding ligand
FeaturePlot(eng.subset.integrated_HK, features = "Areg", reduction = "umap.rpca",order=T)
FeaturePlot(eng.subset.integrated_HK, features = "Btc", reduction = "umap.rpca",order=T)
FeaturePlot(eng.subset.integrated_HK, features = "Cdh1", reduction = "umap.rpca",order=T)
FeaturePlot(eng.subset.integrated_HK, features = "Egf", reduction = "umap.rpca",order=T)
FeaturePlot(eng.subset.integrated_HK, features = "Ereg", reduction = "umap.rpca",order=T)
FeaturePlot(eng.subset.integrated_HK, features = "Nrg1", reduction = "umap.rpca",order=T)

# IL-1 signaling analysis
FeaturePlot(eng.subset.prox.int, features = c("Il1a"), reduction = "umap.rpca")
# Visualize Il1a and Il1r1 expression patterns
FeaturePlot(eng.subset.integrated_HK, features = "Il1a", reduction = "umap.rpca",order=T)
FeaturePlot(eng.subset.integrated_HK, features = "Il1r1", reduction = "umap.rpca",order=T,min.cutoff = 1,pt.size = 1)
# Split by condition
FeaturePlot(eng.subset.integrated_HK, features = "Il1a", reduction = "umap.rpca",order=T,split.by = 'Condition')
FeaturePlot(eng.subset.integrated_HK, features = "Il1r1", reduction = "umap.rpca",order=T,min.cutoff = 1,pt.size = 1)
VlnPlot(eng.subset.integrated_HK,'Il1a',group.by = 'Condition',adjust = 5)

# Identify known receptors for Il1a
roi <- mech.list[mech.list$Ligand.ApprovedSymbol == 'Il1a',]$Receptor.ApprovedSymbol
roi
# Visualize other Il1a-associated receptors
FeaturePlot(eng.subset.integrated_HK, features = "Il1r1", reduction = "umap.rpca",order=T)
FeaturePlot(eng.subset.integrated_HK, features = "Il1r2", reduction = "umap.rpca",order=T)
FeaturePlot(eng.subset.integrated_HK, features = "Il1rap", reduction = "umap.rpca",order=T)
