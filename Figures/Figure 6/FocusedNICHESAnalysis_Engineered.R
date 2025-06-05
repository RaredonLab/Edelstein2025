
options(future.globals.maxSize = 16000 * 1024^2)

# Running NICHES on condition objects

# Set seed
set.seed(2)

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

# Load data
setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets/")
load("BAL_3D.integrated_HK.Robj")
load("PD_3D.integrated_HK.Robj")
load("Mixed_3D.integrated_HK.Robj")
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

# Load the NICHES object for Nuoya
setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets/NICHES_Engineering")
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

# Visualize
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

Idents(regen.circuit.object) = regen.circuit.object$seurat_clusters
aa = DimPlot(regen.circuit.object, label = T) + ggtitle("Clusters")
bb = DimPlot(regen.circuit.object, group.by = "CellType", label = T) # These are old annotations
cc = DimPlot(regen.circuit.object, group.by = "Condition", label = F, cols = Condition.cols) # Condition; cool, cool, this checks out
dd = DimPlot(regen.circuit.object, group.by = "CellClass.sub.regen", label = T, cols = cols.cellclass) # These are old annotations
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

setwd("~/Desktop/Regen Circuit NICHES Subset Analysis")
setwd("~/Desktop/Manuscript Figures/Figure 6")
save(regen.circuit.object, file = "regen.circuit.object.Robj")
load("regen.circuit.object.Robj")

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
                                cell_types = "CellType.regen.spec",
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
setwd("~/Desktop/Regen Circuit NICHES Subset Analysis") # Set wd to save this
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
aa = DimPlot(regen.subset_CTC_bySample,label=T,raster = F,cols=full.col.pal)+ggtitle("NICHES Subset Clusters") + theme(plot.title = element_text(hjust = 0.5))
bb = DimPlot(regen.subset_CTC_bySample,group.by = 'CellClass.sub.regen.Sending',raster=F,shuffle = T,cols=cols.cellclass) + ggtitle("Sending Lineage") + theme(plot.title = element_text(hjust = 0.5))
cc = DimPlot(regen.subset_CTC_bySample,group.by = 'CellClass.sub.regen.Receiving',raster=F,shuffle = T,cols=cols.cellclass) + ggtitle("Receiving Lineage") + theme(plot.title = element_text(hjust = 0.5))
bb | cc
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

# Combine plots side-by-side
send.rec.subset = bb | cc
print(send.rec.subset)
setwd("~/Desktop/Manuscript Figures/Figure 6")
ggsave("send.rec.subset.png", plot = send.rec.subset, width = 12, height = 6, dpi = 600)

dd = DimPlot(regen.subset_CTC_bySample,group.by = 'CellType.regen.spec.Sending',raster=F,shuffle = T,cols=celltype.col.pal)
ee = DimPlot(regen.subset_CTC_bySample,group.by = 'CellType.regen.spec.Receiving',raster=F,shuffle = T,cols=celltype.col.pal)
regen.circuit_NICHES_2 = aa | dd | ee
print(regen.circuit_NICHES_2)

# Save regular clusters, in column "clusters"
regen.subset_CTC_bySample$clusters <- regen.subset_CTC_bySample$seurat_clusters
# Not going to integrate right now because we want to preserve differences; save object
setwd("~/Desktop/Manuscript Figures/Figure 6")
save(regen.subset_CTC_bySample, file = 'regen.circuit_CTC.Robj')

table(regen.subset_CTC_bySample$CellType.regen.spec.Sending)
# Run marker test
Idents(regen.subset_CTC_bySample) = regen.subset_CTC_bySample$CellType.regen.spec.Joint
mark_cell.sending = FindAllMarkers(regen.subset_CTC_bySample, min.pct = 0.1,logfc.threshold = 0.1)
mark_cell.sending$ratio = mark_cell.sending$pct.1/mark_cell.sending$pct.2
mark_cell.sending$power = mark_cell.sending$ratio*mark_cell.sending$avg_log2FC
View(mark_cell.sending)
# Look at top 10 genes per cluster
top10_mark = mark %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(regen.subset_CTC_bySample, features = top10_mark$gene) 
View(top10_mark)

# Okay bring back our side-by-side plots for reference
regen.circuit_NICHES
# Make one for orig object
a = DimPlot(eng.CTC.final,group.by = 'clusters',shuffle=T, raster = F, cols = full.col.pal,label = T, reduction = 'umap') # CellToCell regular embedding and clustering
b = DimPlot(eng.CTC.final,group.by = 'CellType.Sending',raster=F,shuffle = T,cols=celltype.col.pal,reduction='umap')
c = DimPlot(eng.CTC.final,group.by = 'CellType.Receiving',raster=F,shuffle = T,cols=celltype.col.pal,reduction='umap')
glbl.combined = a | b | c
print(glbl.combined)

# Let's first look at Luminal Hillock -> Pro_Inflamm_Mac (aka cluster 14)
FeaturePlot(regen.subset_CTC_bySample, features = c("Il1rn—Il1r2"), order = T)
FeaturePlot(regen.subset_CTC_bySample, features = c("Wnt4—Fzd6"), label = T)

# Let's first look at Rspo3+_Mes -> Luminal Hillock (aka cluster 12)
FeaturePlot(regen.subset_CTC_bySample, features = c("Fn1—Itgb6"))
FeaturePlot(regen.subset_CTC_bySample, features = c("Wnt4—Fzd6"), label = T)

Idents(regen.subset_CTC_bySample) = regen.subset_CTC_bySample$CellType.Sending
VlnPlot(regen.subset_CTC_bySample, features = c("Fn1—Itgb6")) # specific here; does it hold true in broader object?
VlnPlot(regen.subset_CTC_bySample, features = c("Col1a1—Ddr1")) # specific here; does it hold true in broader object?
Idents(regen.subset_CTC_bySample) = regen.subset_CTC_bySample$CellType.Receiving
VlnPlot(regen.subset_CTC_bySample, features = c("Col1a1—Ddr1")) # specific here; does it hold true in broader object?

# In global object
FeaturePlot(eng.CTC.final, features = c("Fn1—Itgb6"))
Idents(eng.CTC.final) = eng.CTC.final$CellType.combined.Integrated.Sending
VlnPlot(eng.CTC.final, features = c("Fn1—Itgb6")) # Not holding in global
VlnPlot(eng.CTC.final, features = c("Col1a1—Ddr1")) # Not holding in global

# Not a whole lot here (above) so what about the Macrophages to Hillock (we have two combos here)
# Path 1: Pro_Inflamm_Mac --> Hillock_Luminal (Cluster 9)
# Path 2: Anti_Inflamm_Mac --> Hillock_Basal (Cluster 13)

# Identify conserved and non-conserved?
Idents(regen.subset_CTC_bySample) = regen.subset_CTC_bySample$clusters
markers_lum.pro_vs_bas.anti <- FindMarkers(regen.subset_CTC_bySample, 
                              ident.1 = "9", 
                              ident.2 = "13", 
                              logfc.threshold = 0.1, 
                              min.pct = 0.1, 
                              only.pos = FALSE)
View(markers_lum.pro_vs_bas.anti)
markers_lum.pro_vs_bas.anti$ratio = markers_lum.pro_vs_bas.anti$pct.1/mark$pct.2
markers_lum.pro_vs_bas.anti$power = markers_lum.pro_vs_bas.anti$ratio*mark$avg_log2FC
# Genes up in pro-inflamm/luminal (avg_log2FC > 0)
# Genes up in anti-inflamm/basal (avg_log2FC < 0)

# Find markers individually
markers_lum.pro <- FindMarkers(regen.subset_CTC_bySample, ident.1 = "9", logfc.threshold = 0.1)
View(markers_lum.pro)
markers_lum.pro$ratio = markers_lum.pro$pct.1/mark$pct.2
markers_lum.pro$power = markers_lum.pro$ratio*mark$avg_log2FC

markers_bas.anti <- FindMarkers(regen.subset_CTC_bySample, ident.1 = "13", logfc.threshold = 0.1)
View(markers_bas.anti)
markers_bas.anti$ratio = markers_bas.anti$pct.1/mark$pct.2
markers_bas.anti$power = markers_bas.anti$ratio*mark$avg_log2FC

# Identify conserved?
common_pairs <- intersect(rownames(markers_lum.pro), rownames(markers_bas.anti))
# Refine
# Define a cutoff
sig_common <- common_pairs[
  markers_lum.pro[common_pairs, "p_val_adj"] < 0.05 & 
    markers_bas.anti[common_pairs, "p_val_adj"] < 0.05]

# Divergent pairs
unique_to_lum.pro <- setdiff(rownames(markers_lum.pro), rownames(markers_bas.anti))
unique_to_bas.anti <- setdiff(rownames(markers_bas.anti), rownames(markers_lum.pro))
unique_to_lum.pro
unique_to_bas.anti

Idents(regen.subset_CTC_bySample) = regen.subset_CTC_bySample$clusters

## Il1b—Adrb2 (Macrophage --> Luminal Hillock) # LIKE THIS MECHANISM
FeaturePlot(regen.subset_CTC_bySample, features = c("Il1b—Adrb2"), order = T, label = T)
# Check ligand and receptor in sc object
FeaturePlot(regen.circuit.object, features = c("Il1b","Adrb2"), order = T, label = T)
# Checking in the object itself
Idents(regen.circuit.object) = regen.circuit.object$CellType
VlnPlot(regen.circuit.object, features = c("Il1b","Adrb2"))
# In eng global transcriptomic object
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType
VlnPlot(eng.subset.integrated_HK, features = c("Il1b","Adrb2"))

Idents(eng.CTC.final) = eng.CTC.final$clusters
FeaturePlot(eng.CTC.final, features = c("Il1b—Adrb2"), order = T, label = T, reduction = "umap")
Idents(eng.CTC.final) = eng.CTC.final$CellType.Sending
FeaturePlot(eng.CTC.final, features = c("Wnt4—Fzd6"), order = T, reduction = "umap")
VlnPlot(eng.CTC.final, features = c("Il1b","Adrb2"))

## Tgfb1—Itgb6 - not specific enough
FeaturePlot(regen.subset_CTC_bySample, features = c("Tgfb1—Itgb6"), order = T, label = T)
# Check ligand and receptor in sc object
FeaturePlot(regen.circuit.object, features = c("Tgfb1","Itgb6"), order = T, label = T)
# Checking in the object itself
Idents(regen.circuit.object) = regen.circuit.object$CellType
VlnPlot(regen.circuit.object, features = c("Tgfb1","Itgb6"))
# In eng global transcriptomic object
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType
VlnPlot(eng.subset.integrated_HK, features = c("Tgfb1","Itgb6"))
Idents(eng.CTC.final) = eng.CTC.final$CellType.Sending
FeaturePlot(eng.CTC.final, features = c("Tgfb1—Itgb6"), order = T, reduction = "umap")

## Osm—Il6st (Macrophage --> Luminal Hillock)
FeaturePlot(regen.subset_CTC_bySample, features = c("Osm—Il6st"), order = T, label = T)
# Check ligand and receptor in sc object
FeaturePlot(regen.circuit.object, features = c("Osm","Il6st"), order = T, label = T)
# Checking in the object itself
Idents(regen.circuit.object) = regen.circuit.object$CellType
VlnPlot(regen.circuit.object, features = c("Osm","Il6st"))
# In eng global transcriptomic object
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType
VlnPlot(eng.subset.integrated_HK, features = c("Osm","Il6st"))
Idents(eng.CTC.final) = eng.CTC.final$CellType.Sending
FeaturePlot(eng.CTC.final, features = c("Osm—Il6st"), order = T, reduction = "umap")

## Rtn4—Gjb2 (Macrophage --> Luminal Hillock)
FeaturePlot(regen.subset_CTC_bySample, features = c("Rtn4—Gjb2"), order = T, label = T)
# Check ligand and receptor in sc object
FeaturePlot(regen.circuit.object, features = c("Rtn4","Gjb2"), order = T, label = T)
# Checking in the object itself
Idents(regen.circuit.object) = regen.circuit.object$CellType
VlnPlot(regen.circuit.object, features = c("Rtn4","Gjb2"))
# In eng global transcriptomic object
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType
VlnPlot(eng.subset.integrated_HK, features = c("Rtn4","Gjb2"))
Idents(eng.CTC.final) = eng.CTC.final$CellType.Sending
FeaturePlot(eng.CTC.final, features = c("Rtn4—Gjb2"), order = T, reduction = "umap")

## Nrg1—Erbb2 (Macrophage --> Luminal Hillock)
FeaturePlot(regen.subset_CTC_bySample, features = c("Nrg1—Erbb2"), order = T, label = T)
# Check ligand and receptor in sc object
FeaturePlot(regen.circuit.object, features = c("Nrg1","Erbb2"), order = T, label = T)
# Checking in the object itself
Idents(regen.circuit.object) = regen.circuit.object$CellType
VlnPlot(regen.circuit.object, features = c("Nrg1","Erbb2"))
# In eng global transcriptomic object
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType
VlnPlot(eng.subset.integrated_HK, features = c("Nrg1","Erbb2"))
Idents(eng.CTC.final) = eng.CTC.final$CellType.Sending
FeaturePlot(eng.CTC.final, features = c("Nrg1—Erbb2"), order = T, reduction = "umap")


############ Now looking at Anti_Inflamm_Mac --> Basal_Hillock ##############

## Lpl—Sdc1 (Macrophage --> Basal Hillock)
FeaturePlot(regen.subset_CTC_bySample, features = c("Lpl—Sdc1"), order = T, label = T)
# Check ligand and receptor in sc object
FeaturePlot(regen.circuit.object, features = c("Lpl","Sdc1"), order = T, label = T)
# Checking in the object itself
Idents(regen.circuit.object) = regen.circuit.object$CellType
VlnPlot(regen.circuit.object, features = c("Lpl","Sdc1"))
# In eng global transcriptomic object
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType
VlnPlot(eng.subset.integrated_HK, features = c("Lpl","Sdc1"))
Idents(eng.CTC.final) = eng.CTC.final$CellType.Sending
FeaturePlot(eng.CTC.final, features = c("Lpl—Sdc1"), order = T, reduction = "umap")

## Psen1—Notch3 (Macrophage --> Basal Hillock)
FeaturePlot(regen.subset_CTC_bySample, features = c("Psen1—Notch3"), order = T, label = T)
# Check ligand and receptor in sc object
FeaturePlot(regen.circuit.object, features = c("Psen1","Notch3"), order = T, label = T)
# Checking in the object itself
Idents(regen.circuit.object) = regen.circuit.object$CellType
VlnPlot(regen.circuit.object, features = c("Psen1","Notch3"))
# In eng global transcriptomic object
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType
VlnPlot(eng.subset.integrated_HK, features = c("Psen1","Notch3"))
Idents(eng.CTC.final) = eng.CTC.final$CellType.Sending
FeaturePlot(eng.CTC.final, features = c("Psen1—Notch3"), order = T, reduction = "umap")

## Psap—Celsr1 - not specific enough; anti-inflamm to pretty much all immune
FeaturePlot(regen.subset_CTC_bySample, features = c("Psap—Celsr1"), order = T, label = T)
# Check ligand and receptor in sc object
FeaturePlot(regen.circuit.object, features = c("Psap","Celsr1"), order = T, label = T)
# Checking in the object itself
Idents(regen.circuit.object) = regen.circuit.object$CellType
VlnPlot(regen.circuit.object, features = c("Psap","Celsr1"))
# In eng global transcriptomic object
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType
VlnPlot(eng.subset.integrated_HK, features = c("Psap","Celsr1"))
Idents(eng.CTC.final) = eng.CTC.final$CellType.Sending
FeaturePlot(eng.CTC.final, features = c("Psap—Celsr1"), order = T, reduction = "umap")

regen.cols = c("Hillock_Like" = "#ff0000", "Pdgfrb+_Pericyte" = "#F4AFB4",
               "Polarized_Mac" = "#ff9e80", "Rspo3+_Mes" = "#89a5a5")
  
VlnPlot(regen.circuit.object, features = c("Serpina1","Lrp1"), group.by = "CellType.regen.spec", pt.size = 0, adjust = 1, cols = regen.cols)
VlnPlot(regen.circuit.object, features = c("Plat","Itgb2"), group.by = "CellType.regen.spec", pt.size = 0, adjust = 1, cols = regen.cols)
VlnPlot(regen.circuit.object, features = c("Slpi","Cd4"), group.by = "CellType.regen.spec", pt.size = 0, adjust = 1, cols = regen.cols)
VlnPlot(regen.circuit.object, features = c("Il1a","Il1r2"), group.by = "CellType.regen.spec", pt.size = 0.1, adjust = 3, cols = regen.cols)
VlnPlot(regen.circuit.object, features = c("Cd24","Selp"), group.by = "CellType.regen.spec", pt.size = 0, adjust = 3, cols = regen.cols)

VlnPlot(regen.circuit.object, features = c("Serpina1"), group.by = "CellType.regen.spec", pt.size = 0, adjust = 1, cols = regen.cols)


(VlnPlot(
  regen.circuit.object,
  features = c("Il1a"),
  group.by = "CellType.regen.spec",
  pt.size = 0.1,
  adjust = 0.5,   # or even 3
  cols = regen.cols)
  + ylim(0, 1))

ggplot(expr_df, aes(x = CellType, y = Expression, color = CellType)) +
  geom_jitter(width = 0.2, size = 0.4, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", width = 0.6) +
  facet_wrap(~ Gene, scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = regen.cols)



