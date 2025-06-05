
# close all, clear all
graphics.off()  
rm(list = ls())

options(future.globals.maxSize = 16000 * 1024^2)

# Set seed
set.seed(2)

## MEMO: NICHES Subset Analysis

# Packages
library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(scales)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(cowplot)
library(patchwork)
library(SeuratWrappers)
library(reticulate)
library(stringr)
library(NICHES)

# Color palettes at start for later visualization
cols.cellclass = c(
  "Epithelium" = "#D982C6",
  "Immune" = "#87B37A",
  "Mesenchyme" = "#F4A261",
  "Endothelium" = "#2A9D8F")
Condition.cols = c(
  "PD" = "#29339B",
  "BAL" = "#e9c716",
  "PD_3D" = "#4a2377",
  "Mixed_3D" = "#ff585e",
  "BAL_3D" = "#00b0be")

# Load data
setwd("/Volumes/Home-1/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets")
load("BAL_3D.integrated_HK.Robj") # Load BAL_3D
load("PD_3D.integrated_HK.Robj") # Load PD_3D
load("Mixed_3D.integrated_HK.Robj") # Load Mixed_3D

setwd("/Users/see27/Desktop/Single Cell/Final Objects")
load("eng.subset.integrated_HK.Robj") # Load merged engineered data just in case

# Object
eng.subset.integrated_HK

# Look at meta-data
str(eng.subset.integrated_HK@meta.data)
table(eng.subset.integrated_HK$CellType)

# Fix cell type annotations 
eng.subset.integrated_HK$CellType <- as.character(eng.subset.integrated_HK$CellType)
eng.subset.integrated_HK$CellType[eng.subset.integrated_HK$CellType == "M1_Mac_Like"] <- "Pro_Inflamm_Mac"
eng.subset.integrated_HK$CellType[eng.subset.integrated_HK$CellType == "M2_Mac_Like"] <- "Anti_Inflamm_Mac"
eng.subset.integrated_HK$CellType <- factor(eng.subset.integrated_HK$CellType)  # convert back to factor if needed
# Confirm the update
table(eng.subset.integrated_HK$CellType)
# Resave
setwd("/Volumes/Home-1/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets")
save(eng.subset.integrated_HK, file = "eng.subset.integrated_HK.Robj")
setwd("/Users/see27/Desktop/Single Cell/Final Objects")
save(eng.subset.integrated_HK, file = "eng.subset.integrated_HK.Robj")

# Load the NICHES object from Nuoya (this is for our engineered combined - in case we need to reference it)
setwd("/Volumes/Home-1/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets/NICHES_Engineering")
load("eng.CTC.final.Robj")

# Plot example
DimPlot(eng.CTC.final,group.by = 'clusters',shuffle=T, raster = F, cols = full.col.pal,label = T, reduction = 'umap') # CellToCell regular embedding and clustering
DimPlot(eng.CTC.final,group.by = 'integrated.clusters',shuffle=T, raster = F, cols = full.col.pal,label = T, reduction = 'umap.rpca') # Embedding integrated by condition
DimPlot(eng.CTC.final,group.by = 'integrated.clusters.2',shuffle=T, raster = F, cols = full.col.pal,label = T, reduction = 'umap.rpca.2') # Embedding integrated by sample

str(eng.subset.integrated_HK@meta.data)
# Subset out just the cells we want from eng.CTC?
str(eng.subset.integrated_HK@meta.data)
table(eng.CTC.final$CellClass.combined.Integrated.Joint)
table(eng.subset.integrated_HK$CellType)

# Pull out: Hillock_Luminal, Hillock_Basal, Rspo3+_Mes, Pdgfrb+_Pericyte, Anti_Inflamm_Mac, Pro_Inflamm_Mac (BUT FROM THE ORIGINAL SINGLE CELL OBJECT NOT NICHES)
selected_types <- c("Hillock_Luminal", "Hillock_Basal", "Rspo3+_Mes",
                    "Pdgfrb+_Pericyte", "Anti_Inflamm_Mac", "Pro_Inflamm_Mac")
regen.circuit.object <- subset(eng.subset.integrated_HK, 
                               subset = CellType %in% selected_types)

# Now re-scale, etc
regen.circuit.object <- NormalizeData(regen.circuit.object) # Don't have to do but doing it anyway
# Identify highly variable features
regen.circuit.object <- FindVariableFeatures(regen.circuit.object)
# Scale the data
regen.circuit.object <- ScaleData(regen.circuit.object)

# Perform PCA
regen.circuit.object <- RunPCA(regen.circuit.object, features = VariableFeatures(object = regen.circuit.object))
# Visualize PCA with an Elbow Plot
ElbowPlot(regen.circuit.object, ndims = 200)
PCHeatmap(regen.circuit.object,cells=200,balanced=T,dims=1:9) # Don't like 8
PCHeatmap(regen.circuit.object,cells=200,balanced=T,dims=10:18) # Things start to look iffy around PC 13
PCHeatmap(regen.circuit.object,cells=200,balanced=T,dims=19:27)

# Run UMAP and clustering - SE picking specific PCs based on what we know about structure
regen.circuit.object <- RunUMAP(regen.circuit.object, dims = c(1:5, 7, 9, 10:12,14))
regen.circuit.object <- FindNeighbors(regen.circuit.object, dims = c(1:5, 7, 9, 10:12,14))
regen.circuit.object <- FindClusters(regen.circuit.object, resolution = 0.4)

Idents(regen.circuit.object) = regen.circuit.object$seurat_clusters
aa = DimPlot(regen.circuit.object, label = T) + ggtitle("Clusters")
bb = DimPlot(regen.circuit.object, group.by = "CellType", label = T) # These are old annotations
cc = DimPlot(regen.circuit.object, group.by = "Condition", label = F, cols = Condition.cols) # Condition; cool, cool, this checks out
dd = DimPlot(regen.circuit.object, group.by = "CellClass.sub.regen", label = T, cols = cols.cellclass) # These are old annotations
# Taking a look at all plots together
(aa | bb) / (cc | dd)

FeaturePlot(regen.circuit.object, features = c("Epcam","Ptprc","Col1a1"), order = T)
cell.epi = WhichCells(regen.circuit.object, idents = c(0,3,5,6,9))
cell.immu = WhichCells(regen.circuit.object, idents = c(2,4,7,1))
cell.mes = WhichCells(regen.circuit.object, idents = c(8,10))
regen.circuit.object = SetIdent(regen.circuit.object, cells = cell.epi, value = 'Epithelium')
regen.circuit.object = SetIdent(regen.circuit.object, cells = cell.mes, value = 'Mesenchyme')
regen.circuit.object = SetIdent(regen.circuit.object, cells = cell.immu, value = 'Immune')
# Stash the cell class
regen.circuit.object$CellClass.sub.regen = Idents(regen.circuit.object)
table(regen.circuit.object$CellClass.sub.regen)
# confirming that everything is labeled
sum(is.na(regen.circuit.object$CellClass.sub.regen))
DimPlot(regen.circuit.object, label = T)

Idents(regen.circuit.object) = regen.circuit.object$seurat_clusters
# Cell Type annos - just adding another meta-data layer
regen.circuit.object$CellType.regen.spec <- NA
c0 <- WhichCells(regen.circuit.object, idents = "0") # Hillock
c1 <- WhichCells(regen.circuit.object, idents = "1") # Polarized_Mac
c2 <- WhichCells(regen.circuit.object, idents = "2") # Polarized_Mac
c3 <- WhichCells(regen.circuit.object, idents = "3") # Hillock
c4 <- WhichCells(regen.circuit.object, idents = "4") # Polarized_Mac
c5 <- WhichCells(regen.circuit.object, idents = "5") # Hillock
c6 <- WhichCells(regen.circuit.object, idents = "6") # Hillock
c7 <- WhichCells(regen.circuit.object, idents = "7") # Polarized_Mac
c8 <- WhichCells(regen.circuit.object, idents = "8") #Rspo_3+ Mes
c9 <- WhichCells(regen.circuit.object, idents = "9") # Hillock
c10 <- WhichCells(regen.circuit.object, idents = "10") # Pericytes

regen.circuit.object$CellType.regen.spec[c0] <- 'Hillock_Like'
regen.circuit.object$CellType.regen.spec[c1] <- 'Polarized_Mac'
regen.circuit.object$CellType.regen.spec[c2] <- 'Polarized_Mac'
regen.circuit.object$CellType.regen.spec[c3] <- 'Hillock_Like'
regen.circuit.object$CellType.regen.spec[c4] <- 'Polarized_Mac'
regen.circuit.object$CellType.regen.spec[c5] <- 'Hillock_Like'
regen.circuit.object$CellType.regen.spec[c6] <- 'Hillock_Like'
regen.circuit.object$CellType.regen.spec[c7] <- 'Polarized_Mac'
regen.circuit.object$CellType.regen.spec[c8] <- 'Rspo_3+_Mes' 
regen.circuit.object$CellType.regen.spec[c9] <- 'Hillock_Like'
regen.circuit.object$CellType.regen.spec[c10] <- 'Pdgfrb+_Pericyte'

table(regen.circuit.object$CellType.regen.spec)
DimPlot(regen.circuit.object, group.by = 'CellType.regen.spec', label = T)

setwd("~/Desktop/NICHES Work Manuscript/Regen Circuit NICHES Subset Analysis")
setwd("~/Desktop/Manuscript Figures/Figure 6")
save(regen.circuit.object, file = "regen.circuit.object.Robj")

####### Now, let's run NICHES ##########
# Check data structure
table(regen.circuit.object$Orig_ID)
table(regen.circuit.object$Condition)
table(regen.circuit.object$CellType)
table(regen.circuit.object$CellClass.sub.regen)
table(regen.circuit.object$CellType.regen.spec)

# Split by conditions first
BAL_3D.regen <- subset(regen.circuit.object, subset = Condition == 'BAL_3D')
Mixed_3D.regen <- subset(regen.circuit.object, subset = Condition == 'Mixed_3D')
PD_3D.regen <- subset(regen.circuit.object, subset = Condition == 'PD_3D')

### NICHES: based on conditions, no imputation

# Condition: BAL_3D
BAL_3D.list <- SplitObject(BAL_3D.regen,split.by='Orig_ID') # splitting by replicate here
NICHES.list <- list()

for(i in 1:length(BAL_3D.list)){
  print(i)
  CellType_Dist <- table(BAL_3D.list[[i]]$CellType.regen.spec)
  to.delete <- names(CellType_Dist[CellType_Dist<=1])
  BAL_3D.list[[i]] <- subset(BAL_3D.list[[i]],subset = CellType.regen.spec %in% to.delete,invert=T)
  
  NICHES.list[[i]] <- RunNICHES(BAL_3D.list[[i]],
                                LR.database = "fantom5",
                                species = "rat",
                                assay = "RNA",
                                cell_types = "CellType.regen.spec", # RUNNING WITH OUR SUBSET CELL TYPES
                                meta.data.to.map = names(BAL_3D.list[[i]]@meta.data),
                                SystemToCell = T,
                                CellToCell = T,
                                CellToSystem = T)
}

## Save CellToCell signals
temp.list <- list()
for(i in 1:length(NICHES.list)){
  temp.list[[i]] <- NICHES.list[[i]]$CellToCell # Isolate CellToCell Signaling
  gc()
}
BAL_3D.CTC.bySample <- merge(temp.list[[1]],temp.list[2:length(temp.list)])
BAL_3D.CTC.bySample <- JoinLayers(BAL_3D.CTC.bySample)

##### Condition: PD_3D
PD_3D.list <- SplitObject(PD_3D.regen,split.by='Orig_ID')
NICHES.list <- list()

for(i in 1:length(PD_3D.list)){
  print(i)
  CellType_Dist <- table(PD_3D.list[[i]]$CellType.regen.spec)
  to.delete <- names(CellType_Dist[CellType_Dist<=1])
  PD_3D.list[[i]] <- subset(PD_3D.list[[i]],subset = CellType.regen.spec %in% to.delete,invert=T)
  
  NICHES.list[[i]] <- RunNICHES(PD_3D.list[[i]],
                                LR.database = "fantom5",
                                species = "rat",
                                assay = "RNA",
                                cell_types = "CellType.regen.spec",
                                meta.data.to.map = names(PD_3D.list[[i]]@meta.data),
                                SystemToCell = T,
                                CellToCell = T,
                                CellToSystem = T)
}


## Save CellToCell signals
temp.list <- list()
for(i in 1:length(NICHES.list)){
  temp.list[[i]] <- NICHES.list[[i]]$CellToCell # Isolate CellToCell Signaling
  gc()
}
PD_3D.CTC.bySample <- merge(temp.list[[1]],temp.list[2:length(temp.list)])
PD_3D.CTC.bySample <- JoinLayers(PD_3D.CTC.bySample)

#### Condition: Mixed_3D
Mixed_3D.list <- SplitObject(Mixed_3D.regen,split.by='Orig_ID')

NICHES.list <- list()

for(i in 1:length(Mixed_3D.list)){
  print(i)
  CellType_Dist <- table(Mixed_3D.list[[i]]$CellType.regen.spec)
  to.delete <- names(CellType_Dist[CellType_Dist<=1])
  Mixed_3D.list[[i]] <- subset(Mixed_3D.list[[i]],subset = CellType.regen.spec %in% to.delete,invert=T)
  
  NICHES.list[[i]] <- RunNICHES(Mixed_3D.list[[i]],
                                LR.database = "fantom5",
                                species = "rat",
                                assay = "RNA",
                                cell_types = "CellType.regen.spec",
                                meta.data.to.map = names(Mixed_3D.list[[i]]@meta.data),
                                SystemToCell = T,
                                CellToCell = T,
                                CellToSystem = T)
}

## Save CellToCell signals
temp.list <- list()
for(i in 1:length(NICHES.list)){
  temp.list[[i]] <- NICHES.list[[i]]$CellToCell # Isolate CellToCell Signaling
  gc()
}
Mixed_3D.CTC.bySample <- merge(temp.list[[1]],temp.list[2:length(temp.list)])
Mixed_3D.CTC.bySample <- JoinLayers(Mixed_3D.CTC.bySample)

# Merge conditions together
eng.CTC.bySample_regen <- merge(BAL_3D.CTC.bySample, c(PD_3D.CTC.bySample,Mixed_3D.CTC.bySample))
eng.CTC.bySample_regen <- JoinLayers(eng.CTC.bySample_regen)

### Filter out low-information signals
VlnPlot(eng.CTC.bySample_regen,c('nFeature_CellToCell'),group.by = 'Orig_ID.Sending',raster=F,pt.size = 0.1) + ggtitle('Before Filtration')+ylab('nFeature_CellToCell')
eng.CTC.bySample_regen <- subset(eng.CTC.bySample_regen,nFeature_CellToCell>40)
# Check again after filtering
VlnPlot(eng.CTC.bySample_regen,c('nFeature_CellToCell'),group.by = 'Orig_ID.Sending',raster=F,pt.size = 0.1)+ggtitle('After Filtration')+ylab('nFeature_CellToCell')

# Make copy just in case
regen.subset_CTC_bySample <- eng.CTC.bySample_regen

# Embedding and clustering
regen.subset_CTC_bySample <- ScaleData(regen.subset_CTC_bySample)
regen.subset_CTC_bySample <- FindVariableFeatures(regen.subset_CTC_bySample)

regen.subset_CTC_bySample <- RunPCA(regen.subset_CTC_bySample,npcs = 100)
setwd("~/Desktop/NICHES Work Manuscript/Regen Circuit NICHES Subset Analysis")
pdf(file='Regen_Subset.CTC.bySample.PCs.pdf',width=10,height=8)
ElbowPlot(regen.subset_CTC_bySample,ndims = 100)
PCHeatmap(regen.subset_CTC_bySample,cells=200,balanced=T,dims=1:9) # Look good through PC 9
PCHeatmap(regen.subset_CTC_bySample,cells=200,balanced=T,dims=10:18)
PCHeatmap(regen.subset_CTC_bySample,cells=200,balanced=T,dims=19:27)
PCHeatmap(regen.subset_CTC_bySample,cells=200,balanced=T,dims=28:36)
PCHeatmap(regen.subset_CTC_bySample,cells=200,balanced=T,dims=37:45)
PCHeatmap(regen.subset_CTC_bySample,cells=200,balanced=T,dims=46:54)
PCHeatmap(regen.subset_CTC_bySample,cells=200,balanced=T,dims=55:63)
PCHeatmap(regen.subset_CTC_bySample,cells=200,balanced=T,dims=64:72)
PCHeatmap(regen.subset_CTC_bySample,cells=200,balanced=T,dims=73:81)
PCHeatmap(regen.subset_CTC_bySample,cells=200,balanced=T,dims=82:90)
PCHeatmap(regen.subset_CTC_bySample,cells=200,balanced=T,dims=91:99)
dev.off()

# Embedding, finding nearest neighbors, and finding clusters (handpicking PCs)
regen.subset_CTC_bySample <- RunUMAP(regen.subset_CTC_bySample,dims = c(1:10,13:14,17:18,21,23))
regen.subset_CTC_bySample <- FindNeighbors(regen.subset_CTC_bySample,dims = c(1:10,13:14,17:18,21,23))
regen.subset_CTC_bySample <- FindClusters(regen.subset_CTC_bySample,resolution=1.3) # Deliberately choosing this high resolution so that the circular cluster in the middle can be 
DimPlot(regen.subset_CTC_bySample,label=T,raster = F,cols=full.col.pal)

# Plot specifics for condition sending, cell types (sending and receiving)
DimPlot(regen.subset_CTC_bySample,label=T,raster = F, cols=full.col.pal,split.by = 'Condition.Sending',ncol=3)
DimPlot(regen.subset_CTC_bySample,group.by = 'CellType.regen.spec.Sending',raster=F,shuffle = T,cols=celltype.col.pal)
DimPlot(regen.subset_CTC_bySample,group.by = 'CellType.regen.spec.Receiving',raster=F,shuffle = T,cols=celltype.col.pal)

# Want to see three plots side-by-side
cols.cellclass = c(
  "Epithelium" = "#D982C6",
  "Immune" = "#87B37A",
  "Mesenchyme" = "#F4A261")
bb <- DimPlot(
  regen.subset_CTC_bySample,
  group.by = 'CellClass.sub.regen.Sending',
  raster = FALSE,
  shuffle = TRUE,
  cols = cols.cellclass) +
  ggtitle("Sending Lineage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()
cc <- DimPlot(
  regen.subset_CTC_bySample,
  group.by = 'CellClass.sub.regen.Receiving',
  raster = FALSE,
  shuffle = TRUE,
  cols = cols.cellclass) +
  ggtitle("Receiving Lineage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()
send.rec.subset = bb | cc
print(send.rec.subset)

# Save regular clusters, in column "clusters"
regen.subset_CTC_bySample$clusters <- regen.subset_CTC_bySample$seurat_clusters
# Not going to integrate right now because we want to preserve differences; save object
setwd("~/Desktop/Regen Circuit - for MSBR")
save(regen.subset_CTC_byCondition, file = 'regen.circuit_CTC.Robj')
save(regen.circuit.object, file = "regen.circuit.object.Robj")
# Now we have both our subset object and our NICHES data for the subset object



