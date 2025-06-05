
# close all, clear all
graphics.off()  
rm(list = ls())

# Load libraries
library(NICHES)
library(Seurat)
library(ggplot2)
library(dplyr)

# load global object (transriptomic)
setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets")
load("eng.subset.integrated_HK.Robj")

# Load NICHES object
setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets/NICHES_Engineering")
load("eng.CTC.final.Robj")

setwd("/Users/see27/Desktop/Single Cell/Final Objects")
load("eng.subset.integrated_HK.Robj")

# object
eng.subset.integrated_HK

# Look at meta-data
str(eng.subset.integrated_HK@meta.data)
table(eng.subset.integrated_HK$CellType)
eng.subset.integrated_HK$CellType <- as.character(eng.subset.integrated_HK$CellType)
eng.subset.integrated_HK$CellType[eng.subset.integrated_HK$CellType == "M1_Mac_Like"] <- "Pro_Inflamm_Mac"
eng.subset.integrated_HK$CellType[eng.subset.integrated_HK$CellType == "M2_Mac_Like"] <- "Anti_Inflamm_Mac"
eng.subset.integrated_HK$CellType <- factor(eng.subset.integrated_HK$CellType)  # convert back to factor if needed
# Confirm the update
table(eng.subset.integrated_HK$CellType)
# Resave
setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets/")
save(eng.subset.integrated_HK, file = "eng.subset.integrated_HK.Robj")

### Subset out mesenchyme
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellClass.combined.Integrated
mes.sub = subset(eng.subset.integrated_HK, ident = "Mesenchyme")
# Now re-scale, etc
mes.sub <- NormalizeData(mes.sub) # Don't have to do but doing it anyway
# Identify highly variable features
mes.sub <- FindVariableFeatures(mes.sub)
# Scale the data
mes.sub <- ScaleData(mes.sub)

# Perform PCA
mes.sub <- RunPCA(mes.sub, features = VariableFeatures(object = mes.sub))
# Visualize PCA with an Elbow Plot
ElbowPlot(mes.sub, ndims = 200)
PCHeatmap(mes.sub,cells=200,balanced=T,dims=1:9) 
PCHeatmap(mes.sub,cells=200,balanced=T,dims=10:18) 
PCHeatmap(mes.sub,cells=200,balanced=T,dims=19:27)

# Run UMAP and clustering - SE picking specific PCs based on what we know about structure
mes.sub <- RunUMAP(mes.sub, dims = c(1:3, 7,9))
mes.sub <- FindNeighbors(mes.sub, dims = c(1:3, 7,9))
mes.sub <- FindClusters(mes.sub, resolution = 0.2)
DimPlot(mes.sub)
FeaturePlot(mes.sub, features = c("Pdgfrb","Col1a1","Acta2","Gucy1b1","Rspo3"))
DimPlot(mes.sub, group.by = "CellType.combined.Integrated", reduction = "umap")

# Grab your ground truth ligand and receptor lists
ligands <- unique(NICHES::ncomms8866_rat$Ligand.ApprovedSymbol)
receptors <- unique(NICHES::ncomms8866_rat$Receptor.ApprovedSymbol)

# See where ligand and receptor ground truth intersects with our data
lig.list <- ligands[ligands %in% rownames(eng.subset.integrated_HK)]
rec.list <- receptors[receptors %in% rownames(eng.subset.integrated_HK)]
# Look at both, notice the length has not changed
lig.list
rec.list

# Generate marker list with idents set to cell class - for ligands
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellClass.combined.Integrated
glbl.mark.class.lig = FindAllMarkers(eng.subset.integrated_HK, features = lig.list, 
                                     only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)
glbl.mark.class.lig$ratio = glbl.mark.class.lig$pct.1/glbl.mark.class.lig$pct.2
glbl.mark.class.lig$power = glbl.mark.class.lig$ratio*glbl.mark.class.lig$avg_log2FC
View(glbl.mark.class.lig)
setwd("~/Desktop/Manuscript Figures/Basal vs. Hillock LR Analysis for MSBR")
save(glbl.mark.class.lig, file = "glbl.mark.class.lig.Robj")

# Generate marker list with idents set to cell class - for receptors
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellClass.combined.Integrated
glbl.mark.class.recep = FindAllMarkers(eng.subset.integrated_HK, features = rec.list, 
                                       only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)
glbl.mark.class.recep$ratio = glbl.mark.class.recep$pct.1/glbl.mark.class.recep$pct.2
glbl.mark.class.recep$power = glbl.mark.class.recep$ratio*glbl.mark.class.recep$avg_log2FC
View(glbl.mark.class.recep)
setwd("~/Desktop/Manuscript Figures/Basal vs. Hillock LR Analysis for MSBR")
save(glbl.mark.class.recep, file = "glbl.mark.class.recep.Robj")

#### Load the proximal subset 
setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets/Subset Population Objects/Proximal Subset")
load("eng.subset.prox.int.Robj")
# look at the object
str(eng.subset.prox.int@meta.data)
DimPlot(eng.subset.prox.int, reduction = "umap.rpca", group.by = "CellType_Prox_Subset")

Idents(eng.subset.prox.int) = eng.subset.prox.int$CellType_Prox_Subset
prox.markers.lig = FindAllMarkers(eng.subset.prox.int, features = lig.list,
                                  only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)
prox.markers.lig$ratio = prox.markers.lig$pct.1/prox.markers.lig$pct.2
prox.markers.lig$power = prox.markers.lig$ratio*prox.markers.lig$avg_log2FC
View(prox.markers.lig)
setwd("~/Desktop/Manuscript Figures/Basal vs. Hillock LR Analysis for MSBR")
save(prox.markers.lig, file = "prox.markers.lig.Robj")

# Same thing for receptors
Idents(eng.subset.prox.int) = eng.subset.prox.int$CellType_Prox_Subset
prox.markers.recep = FindAllMarkers(eng.subset.prox.int, features = rec.list,
                                    only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)
prox.markers.recep$ratio = prox.markers.recep$pct.1/prox.markers.recep$pct.2
prox.markers.recep$power = prox.markers.recep$ratio*prox.markers.recep$avg_log2FC
View(prox.markers.recep)
setwd("~/Desktop/Manuscript Figures/Basal vs. Hillock LR Analysis for MSBR")
save(prox.markers.recep, file = "prox.markers.recep.Robj")

Idents(eng.subset.prox.int) = eng.subset.prox.int$CellType_Prox_Subset
prox.markers.lig_2 = FindMarkers(eng.subset.prox.int, features = lig.list,
                                    ident.1 = c("Hillock_Luminal","Hillock_Basal"),
                                 ident.2 = NULL,
                                  only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)
prox.markers.lig_2$ratio = prox.markers.lig_2$pct.1/prox.markers.lig_2$pct.2
prox.markers.lig_2$power = prox.markers.lig_2$ratio*prox.markers.lig_2$avg_log2FC
View(prox.markers.lig_2)
setwd("~/Desktop/Manuscript Figures/Basal vs. Hillock LR Analysis for MSBR")
save(prox.markers.lig_2, file = "prox.markers.lig_2.Robj")

# Same thing for receptors
Idents(eng.subset.prox.int) = eng.subset.prox.int$CellType_Prox_Subset
prox.markers.recep_2 = FindMarkers(eng.subset.prox.int, features = rec.list, 
                                   ident.1 = c("Hillock_Luminal", "Hillock_Basal"),
                                   ident.2 = NULL,
                                    only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)
prox.markers.recep_2$ratio = prox.markers.recep_2$pct.1/prox.markers.recep_2$pct.2
prox.markers.recep_2$power = prox.markers.recep_2$ratio*prox.markers.recep_2$avg_log2FC
View(prox.markers.recep_2)
setwd("~/Desktop/Manuscript Figures/Basal vs. Hillock LR Analysis for MSBR")
save(prox.markers.recep_2, file = "prox.markers.recep_2.Robj")

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

### Luminal Hillock Specific ligand: Il1a
FeaturePlot(eng.subset.prox.int, features = c("Il1a"), reduction = "umap.rpca", order = T, label = T, pt.size = 0.75)
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType.combined.Integrated
FeaturePlot(eng.subset.integrated_HK, features = "Il1a", reduction = "umap.rpca", min.cutoff = 1.0, order = T, pt.size = 0.5, label = T)

# What receptors can bind Il1a?
roi <- mech.list[mech.list$Ligand.ApprovedSymbol == 'Il1a',]$Receptor.ApprovedSymbol
roi
# Look at in global object
a = FeaturePlot(eng.subset.integrated_HK, features = "Il1r1", reduction = "umap.rpca",order=T, min.cutoff = 1.0, pt.size = 0.5, label = T)
b = FeaturePlot(eng.subset.integrated_HK, features = "Il1r2", reduction = "umap.rpca",order=T, pt.size = 0.25, max.cutoff = 1, label = T)
c = FeaturePlot(eng.subset.integrated_HK, features = "Il1rap", reduction = "umap.rpca",order=T, min.cutoff = 1, label = T, pt.size = 0.25)
a | b | c
# Not telling us anything about hillock to immune though 

### Luminal Hillock Specific ligand: Serpina1
FeaturePlot(eng.subset.prox.int, features = c("Serpina1"), reduction = "umap.rpca", order = T, label = T, pt.size = 0.75)
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType.combined.Integrated
FeaturePlot(eng.subset.integrated_HK, features = "Serpina1", reduction = "umap.rpca", order = T, pt.size = 0.5, label = T, min.cutoff = 0.5)

# What receptors can bind Serpina1?
roi <- mech.list[mech.list$Ligand.ApprovedSymbol == 'Serpina1',]$Receptor.ApprovedSymbol
roi
# Look at in global object
FeaturePlot(eng.subset.integrated_HK, features = "Lrp1", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T, min.cutoff = 1)

### Luminal Hillock Specific ligand: Slpi
DimPlot(eng.subset.prox.int, split.by = "Orig_ID", ncol = 3, reduction = "umap.rpca")
FeaturePlot(eng.subset.prox.int, features = c("Slpi"), reduction = "umap.rpca", order = T, label = T, pt.size = 0.75)
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType.combined.Integrated
FeaturePlot(eng.subset.integrated_HK, features = "Slpi", reduction = "umap.rpca", order = T, pt.size = 0.5, label = T, min.cutoff = 0.5)

# What receptors can bind Slpi?
roi <- mech.list[mech.list$Ligand.ApprovedSymbol == 'Slpi',]$Receptor.ApprovedSymbol
roi
# Only one that can bind Slpi is Cd4
FeaturePlot(eng.subset.integrated_HK, features = "Cd4", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)

# What receptors can bind Slpi?
loi <- mech.list[mech.list$Receptor.ApprovedSymbol == 'Cd4',]$Ligand.ApprovedSymbol
loi
# Only one that can bind Slpi is Cd4
FeaturePlot(eng.subset.integrated_HK, features = "Cd4", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)

FeaturePlot(eng.subset.integrated_HK, features = c("Il1a"), reduction = "umap.rpca")

### Luminal Hillock Specific ligand: Plat
DimPlot(eng.subset.prox.int, split.by = "Orig_ID", ncol = 3, reduction = "umap.rpca")
FeaturePlot(eng.subset.prox.int, features = c("Plat"), reduction = "umap.rpca", order = T, label = T, pt.size = 0.75)
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType.combined.Integrated
FeaturePlot(eng.subset.integrated_HK, features = "Plat", reduction = "umap.rpca", order = T, pt.size = 0.5, label = T, min.cutoff = 0.5)

# What receptors can bind Plat?
roi <- mech.list[mech.list$Ligand.ApprovedSymbol == 'Plat',]$Receptor.ApprovedSymbol
roi
# Only one that can bind Slpi is Cd4
FeaturePlot(eng.subset.integrated_HK, features = "Itgb2", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)
FeaturePlot(eng.subset.integrated_HK, features = "Lrp1", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)

### Luminal Hillock Specific ligand: Ephb3
DimPlot(eng.subset.prox.int, split.by = "Orig_ID", ncol = 3, reduction = "umap.rpca")
FeaturePlot(eng.subset.prox.int, features = c("Plau"), reduction = "umap.rpca", order = T, label = T, pt.size = 0.75)
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType.combined.Integrated
FeaturePlot(eng.subset.integrated_HK, features = "Tyro3", reduction = "umap.rpca", order = T, pt.size = 0.5, label = T, min.cutoff = 0.5)

# What receptors can bind Defb1?
roi <- mech.list[mech.list$Ligand.ApprovedSymbol == 'Cd24',]$Receptor.ApprovedSymbol
roi
# Only one that can bind Slpi is Cd4
FeaturePlot(eng.subset.integrated_HK, features = "Selp", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)
FeaturePlot(eng.subset.integrated_HK, features = "Efnb2", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)
FeaturePlot(eng.subset.integrated_HK, features = "Efnb3", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)

### Luminal Hillock Specific ligand: Defb1
FeaturePlot(eng.subset.prox.int, features = c("Defb1"), reduction = "umap.rpca", order = T, label = T, pt.size = 0.75)
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType.combined.Integrated
FeaturePlot(eng.subset.integrated_HK, features = "Defb1", reduction = "umap.rpca", order = T, pt.size = 0.5, label = T, min.cutoff = 0.5)

# What receptors can bind Defb1?
roi <- mech.list[mech.list$Ligand.ApprovedSymbol == 'Defb1',]$Receptor.ApprovedSymbol
roi
# Only one that can bind Defb1 is Ccr6 but nothing is hearing it because Ccr6 not present
FeaturePlot(eng.subset.integrated_HK, features = "Ccr6", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)
# NOPE

### Hillock Specific receptor: Tspan1
FeaturePlot(eng.subset.prox.int, features = c("Tspan1"), reduction = "umap.rpca", order = T, label = T, pt.size = 0.75)
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType.combined.Integrated
FeaturePlot(eng.subset.integrated_HK, features = "Tspan1", reduction = "umap.rpca", order = T, pt.size = 0.5, label = T, min.cutoff = 0.5)

# What receptors can bind Defb1?
loi <- mech.list[mech.list$Receptor.ApprovedSymbol == 'Tspan1',]$Ligand.ApprovedSymbol
loi
# Only one that can bind Defb1 is Ccr6 but nothing is hearing it because Ccr6 not present
FeaturePlot(eng.subset.integrated_HK, features = "Mdk", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)
# NOPE

### Basal Hillock Specific Receptor? Nrg3
FeaturePlot(eng.subset.prox.int, features = c("Nrg4"), reduction = "umap.rpca", order = T, label = T, pt.size = 0.75)
VlnPlot(eng.subset.prox.int, features = c("Nrg4"), adjust = 3)
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType.combined.Integrated
FeaturePlot(eng.subset.integrated_HK, features = "Nrg4", reduction = "umap.rpca", order = T, pt.size = 0.5, label = T, min.cutoff = 0.5)

# What receptors can bind Nrg4?
roi <- mech.list[mech.list$Ligand.ApprovedSymbol == 'Nrg4',]$Receptor.ApprovedSymbol
roi
# Look at in global object
a = FeaturePlot(eng.subset.integrated_HK, features = "Egfr", reduction = "umap.rpca",order=T, min.cutoff = 1.0, pt.size = 0.5, label = T)
b = FeaturePlot(eng.subset.integrated_HK, features = "Erbb2", reduction = "umap.rpca",order=T, pt.size = 0.25, max.cutoff = 1, label = T)
c = FeaturePlot(eng.subset.integrated_HK, features = "Erbb4", reduction = "umap.rpca",order=T, min.cutoff = 1, label = T, pt.size = 0.25)
a | b | c


############## We can also bring in the immune subset
setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets/Subset Population Objects")
load("immune_subset.clean.Robj")

# look at the object
str(immune_subset.clean@meta.data)
DimPlot(immune_subset.clean, reduction = "umap", group.by = "CellType.sub")
FeaturePlot(immune_subset.clean, features = c("Cd4","Pparg","Mrc1","Prodh2"), reduction = "umap")

# Get ligands specific to the immune polarized populations
Idents(immune_subset.clean) = immune_subset.clean$CellType.sub
imm.markers.lig = FindAllMarkers(immune_subset.clean, features = lig.list,
                                 only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)
imm.markers.lig$ratio = imm.markers.lig$pct.1/imm.markers.lig$pct.2
imm.markers.lig$power = imm.markers.lig$ratio*imm.markers.lig$avg_log2FC
View(imm.markers.lig)
setwd("~/Desktop/Manuscript Figures/Basal vs. Hillock LR Analysis for MSBR")
save(imm.markers.lig, file = "imm.markers.lig.Robj")

# Same thing for receptors
Idents(immune_subset.clean) = immune_subset.clean$CellType.sub
imm.markers.recep = FindAllMarkers(immune_subset.clean, features = rec.list,
                                   only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)
imm.markers.recep$ratio = imm.markers.recep$pct.1/imm.markers.recep$pct.2
imm.markers.recep$power = imm.markers.recep$ratio*imm.markers.recep$avg_log2FC
View(imm.markers.recep)
setwd("~/Desktop/Manuscript Figures/Basal vs. Hillock LR Analysis for MSBR")
save(imm.markers.recep, file = "imm.markers.recep.Robj")

# Look at ligands first (Pro Inflamm)
FeaturePlot(immune_subset.clean, features = c("Vcan"), reduction = "umap") # way more specific
FeaturePlot(eng.subset.integrated_HK, features = "Vcan", reduction = "umap.rpca", order = T)
FeaturePlot(immune_subset.clean, features = c("Cxcl3"), reduction = "umap")
FeaturePlot(eng.subset.integrated_HK, features = "Cxcl3", reduction = "umap.rpca", order = T)
FeaturePlot(immune_subset.clean, features = c("Osm"), reduction = "umap")
FeaturePlot(eng.subset.integrated_HK, features = "Osm", reduction = "umap.rpca", order = T)

# Find receptors that bind Vcan
loi <- mech.list[mech.list$Ligand.ApprovedSymbol == 'Vcan',]$Receptor.ApprovedSymbol
loi

FeaturePlot(eng.subset.integrated_HK, features = "Cd44", reduction = "umap.rpca",order=T) # interesting b/c everything but luminal hillock
FeaturePlot(eng.subset.integrated_HK, features = "Egfr", reduction = "umap.rpca",order=T)
FeaturePlot(eng.subset.integrated_HK, features = "Itga4", reduction = "umap.rpca",order=T)
FeaturePlot(eng.subset.integrated_HK, features = "Itgb1", reduction = "umap.rpca",order=T)
FeaturePlot(eng.subset.integrated_HK, features = "Tlr1", reduction = "umap.rpca",order=T)
FeaturePlot(eng.subset.integrated_HK, features = "Tlr2", reduction = "umap.rpca",order=T)
FeaturePlot(eng.subset.integrated_HK, features = "Sell", reduction = "umap.rpca",order=T)
FeaturePlot(eng.subset.integrated_HK, features = "Selp", reduction = "umap.rpca",order=T)

# Find receptors that bind Osm (this has come up in some of my other analyses)
loi <- mech.list[mech.list$Ligand.ApprovedSymbol == 'Osm',]$Receptor.ApprovedSymbol
loi

FeaturePlot(eng.subset.integrated_HK, features = "Il6st", reduction = "umap.rpca",order=T, min.cutoff = 0.5) # Higher in the hillock territory
FeaturePlot(eng.subset.integrated_HK, features = "Lifr", reduction = "umap.rpca",order=T) # Specific to distal epi
FeaturePlot(eng.subset.integrated_HK, features = "Osmr", reduction = "umap.rpca",order=T)

################ Mesenchyme subset
# look at the object
str(mes.sub@meta.data)
DimPlot(mes.sub, reduction = "umap", group.by = "CellType.combined.Integrated")

# Get ligands specific to the immune polarized populations
Idents(mes.sub) = mes.sub$CellType.combined.Integrated
mes.sub.lig = FindMarkers(mes.sub, features = lig.list, ident.1 = "Rspo3+_Mes",
                                 only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)
mes.sub.lig$ratio = mes.sub.lig$pct.1/mes.sub.lig$pct.2
mes.sub.lig$power = mes.sub.lig$ratio*mes.sub.lig$avg_log2FC
View(mes.sub.lig)
setwd("~/Desktop/Manuscript Figures/Basal vs. Hillock LR Analysis for MSBR")
save(mes.sub.lig, file = "mes.sub.lig.Robj")

# Same thing for receptors
Idents(immune_subset.clean) = immune_subset.clean$CellType.sub
imm.markers.recep = FindAllMarkers(immune_subset.clean, features = rec.list,
                                   only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)
imm.markers.recep$ratio = imm.markers.recep$pct.1/imm.markers.recep$pct.2
imm.markers.recep$power = imm.markers.recep$ratio*imm.markers.recep$avg_log2FC
View(imm.markers.recep)
setwd("~/Desktop/Manuscript Figures/Basal vs. Hillock LR Analysis for MSBR")
save(imm.markers.recep, file = "imm.markers.recep.Robj")
# ABOVE DID NOT WORK BECAUSE WE ONLY HAVE ONE POPULATION
# Use only ligands expressed in your subset
lig.expr_mes <- GetAssayData(mes.sub, slot = "data")[lig.list, , drop = FALSE]
# Calculate percent expressed
pct.exp <- Matrix::rowMeans(lig.expr_mes > 0)
# Calculate average expression
avg.exp <- Matrix::rowMeans(lig.expr_mes)
# Combine into a data frame
lig.summary_mes <- data.frame(
  gene = rownames(lig.expr_mes),
  pct.exp = pct.exp,
  avg.exp = avg.exp)
# Add a simple "power" score if desired (e.g., pct.exp * avg.exp)
lig.summary_mes$power <- lig.summary_mes$pct.exp * lig.summary_mes$avg.exp
View(lig.summary_mes)

## Receptors
# Use only receptors expressed in your subset
rec.expr_mes <- GetAssayData(mes.sub, slot = "data")[rec.list, , drop = FALSE]
# Calculate percent expressed
pct.exp.rec <- Matrix::rowMeans(rec.expr_mes > 0)
# Calculate average expression
avg.exp.rec <- Matrix::rowMeans(rec.expr_mes)
# Combine into a data frame
rec.summary_mes <- data.frame(
  gene = rownames(rec.expr_mes),
  pct.exp = pct.exp.rec,
  avg.exp = avg.exp.rec)
# Add a simple "power" score if desired (e.g., pct.exp * avg.exp)
rec.summary_mes$power <- rec.summary_mes$pct.exp * rec.summary_mes$avg.exp
View(rec.summary_mes)

### Mesenchyme ligand: Fn1
DimPlot(mes.sub, split.by = "Orig_ID", ncol = 3, reduction = "umap.rpca")
FeaturePlot(mes.sub, features = c("Fn1"), reduction = "umap", order = T, label = T, pt.size = 0.75)
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType.combined.Integrated
FeaturePlot(eng.subset.integrated_HK, features = "Fn1", reduction = "umap.rpca", order = T, pt.size = 0.5, label = T, min.cutoff = 2.0)

# What receptors can bind Fn1?
roi <- mech.list[mech.list$Ligand.ApprovedSymbol == 'Fn1',]$Receptor.ApprovedSymbol
roi
# Bind Fn1: Itgb1, Cd44, Itgav, Robo4, Tshr, Il17rc
FeaturePlot(eng.subset.integrated_HK, features = "Il17rc", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)
FeaturePlot(eng.subset.integrated_HK, features = "Itgav", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T, min.cutoff = 1.0)
FeaturePlot(eng.subset.integrated_HK, features = "Flt4", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)
FeaturePlot(eng.subset.integrated_HK, features = "Cd44", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)
FeaturePlot(eng.subset.integrated_HK, features = "Plaur", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)
FeaturePlot(eng.subset.integrated_HK, features = "Itgb3", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)
FeaturePlot(eng.subset.integrated_HK, features = "Itgb1", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)
# Consider Fn1-Itgav, Fn1-Itgb3

### Mesenchyme ligand: Vim
DimPlot(mes.sub, split.by = "Orig_ID", ncol = 3, reduction = "umap.rpca")
FeaturePlot(mes.sub, features = c("Vim"), reduction = "umap", order = T, label = T, pt.size = 0.75)
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType.combined.Integrated
FeaturePlot(eng.subset.integrated_HK, features = "Vim", reduction = "umap.rpca", order = T, pt.size = 0.5, label = T, min.cutoff = 2.0)
# Super high for immune/mes

# What receptors can bind Vim?
roi <- mech.list[mech.list$Ligand.ApprovedSymbol == 'Vim',]$Receptor.ApprovedSymbol
roi
# Cd44 is only receptor to bind Vim
FeaturePlot(eng.subset.integrated_HK, features = "Cd44", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)
# Not really helpful

### Mesenchyme ligand: Rspo3
DimPlot(mes.sub, split.by = "Orig_ID", ncol = 3, reduction = "umap.rpca")
FeaturePlot(mes.sub, features = c("Rspo3"), reduction = "umap", order = T, label = T, pt.size = 0.75)
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType.combined.Integrated
FeaturePlot(eng.subset.integrated_HK, features = "Rspo3", reduction = "umap.rpca", order = T, pt.size = 0.5, label = T)

# What receptors can bind Vim?
roi <- mech.list[mech.list$Ligand.ApprovedSymbol == 'Rspo3',]$Receptor.ApprovedSymbol
roi
# Receptors to bind: "Fzd8" "Lgr4" "Lgr5" "Lgr6" "Lrp6" "Sdc4"
FeaturePlot(eng.subset.integrated_HK, features = "Fzd8", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)
FeaturePlot(eng.subset.integrated_HK, features = "Lgr4", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)
FeaturePlot(eng.subset.integrated_HK, features = "Lgr5", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)
FeaturePlot(eng.subset.integrated_HK, features = "Lgr6", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)
FeaturePlot(eng.subset.integrated_HK, features = "Lrp6", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)
FeaturePlot(eng.subset.integrated_HK, features = "Sdc4", reduction = "umap.rpca",order=T, pt.size = 0.5, label = T)

Idents(regen.circuit.object) = regen.circuit.object$CellType.regen.spec
# Ligands
VlnPlot(regen.circuit.object, features = c("Defb1","Plat","Adm","Cd24","Cgn","Serpina1","Slpi"))
# Receptor
VlnPlot(regen.circuit.object, features = c("Tspan1","Ephb3","Tyro3","F3","Cldn4","Erbb3"))

library(Seurat)
library(ggplot2)
library(patchwork)

genes <- c("Tspan1","Ephb3","Tyro3","F3","Cldn4","Erbb3")

# Create plots with y-axis text removed
vln_plots <- lapply(seq_along(genes), function(i) {
  p <- VlnPlot(
    regen.circuit.object,
    features = genes[i],
    group.by = "CellType.regen.spec",
    pt.size = 0.1,
    cols = regen.cols
  ) +
    ggtitle(genes[i]) +
    ylab(NULL) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      axis.title.y = element_blank()
    )
  
  # Remove x-axis from all but bottom plot
  if (i != length(genes)) {
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
      )
  }
  
  return(p)
})

# Stack and add a single shared y-axis label
final.rec = wrap_plots(vln_plots, ncol = 1) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(),
      axis.title.y = element_text(size = 14)
    )
  ) &
  theme(axis.title.y = element_blank()) &
  patchwork::plot_layout() &
  ylab("Log-normalized expression")  # Single y-axis label for the whole stack

setwd("~/Desktop/Manuscript Figures/Figure 6")
ggsave("rec.up.hillock.png", plot = final.rec, width = 5, height = 10, dpi = 600)

genes <- c("Areg","Ereg","Erbb3","Ephb3","Tgfa")

# Create plots with y-axis text removed
vln_plots <- lapply(seq_along(genes), function(i) {
  p <- VlnPlot(
    regen.circuit.object,
    features = genes[i],
    group.by = "CellType.regen.spec",
    pt.size = 0.1,
    cols = regen.cols
  ) +
    ggtitle(genes[i]) +
    ylab(NULL) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      axis.title.y = element_blank()
    )
  
  # Remove x-axis from all but bottom plot
  if (i != length(genes)) {
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
      )
  }
  
  return(p)
})

# Stack and add a single shared y-axis label
final.lig = wrap_plots(vln_plots, ncol = 1) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(),
      axis.title.y = element_text(size = 14)
    )
  ) &
  theme(axis.title.y = element_blank()) &
  patchwork::plot_layout() &
  ylab("Log-normalized expression") 
# Single y-axis label for the whole stack
final.lig
setwd("~/Desktop/Manuscript Figures/Figure 6")
ggsave("lig.up.hillock.png", plot = final.lig, width = 5, height = 10, dpi = 600)

lig.rec.vln = final.lig | final.rec
setwd("~/Desktop/Manuscript Figures/Figure 6")
ggsave("lig.rec.vln.hillock.png", plot = lig.rec.vln, width = 6, height = 9, dpi = 600)

