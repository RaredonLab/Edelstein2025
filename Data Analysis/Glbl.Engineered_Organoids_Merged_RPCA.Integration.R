

# Memo: Merging and integrating HK cleaned engineered samples
# Condition: All engineered cleaned by HK
# Date: 03-12-2025

# close all, clear all
graphics.off()  
rm(list = ls())

# Packages
library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(scales)
library(dplyr)
library(circlize)
library(cowplot)
library(patchwork)
library(reticulate)
library(stringr)
library(NICHES)

# Set global options
options(future.globals.maxSize = 10000 * 1024^2)

# Set seed
set.seed(2)

setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets")
load("BAL_3D.integrated_HK.Robj")
load("Mixed_3D.integrated_HK.Robj")
load("PD_3D.integrated_HK.Robj")

# Look at all objects
BAL_3D.integrated_HK # 9676 cells
Mixed_3D.integrated_HK # 6344 cells
PD_3D.integrated_HK # 8724 cells

# Set wd 
setwd("~/Desktop/Single Cell/Organoid Project/HK Cleaned Engineered Subset/Integration")
# Merge the conditions
eng.subset_HK.cleaned <- merge(BAL_3D.integrated_HK, y = c(Mixed_3D.integrated_HK, PD_3D.integrated_HK))

# Join layers
eng.subset_HK.cleaned = JoinLayers(eng.subset_HK.cleaned)

# Scale (now that we have merged)
eng.subset_HK.cleaned = ScaleData(eng.subset_HK.cleaned)
# Normalizing the data
eng.subset_HK.cleaned = NormalizeData(eng.subset_HK.cleaned)
# Finding variable features
eng.subset_HK.cleaned = FindVariableFeatures(eng.subset_HK.cleaned)

# Setting wd to local folder to share PCA outputs
setwd("~/Desktop/Single Cell/Organoid Project/HK Cleaned Engineered Subset/Integration")
# Principle component analysis
eng.subset_HK.cleaned = RunPCA(eng.subset_HK.cleaned, npcs = 100)
pdf(file='eng.subset_HK.cleaned_PCs_v1.pdf',width=10,height=8)
# Visualizing the PCAs
ElbowPlot(eng.subset_HK.cleaned,ndims = 100)
PCHeatmap(eng.subset_HK.cleaned,cells=200,balanced=T,dims=1:9)
PCHeatmap(eng.subset_HK.cleaned,cells=200,balanced=T,dims=10:18)
PCHeatmap(eng.subset_HK.cleaned,cells=200,balanced=T,dims=19:27)
PCHeatmap(eng.subset_HK.cleaned,cells=200,balanced=T,dims=28:36)
PCHeatmap(eng.subset_HK.cleaned,cells=200,balanced=T,dims=37:45)
dev.off()

# UMAP
eng.subset_HK.cleaned = RunUMAP(eng.subset_HK.cleaned, dims = 1:36)
DimPlot(eng.subset_HK.cleaned, label = TRUE)
# Cluster data (creating nearest neighbor graph, not clustering)
eng.subset_HK.cleaned = FindNeighbors(eng.subset_HK.cleaned, dims = 1:36)
# Defining clusters
eng.subset_HK.cleaned = FindClusters(eng.subset_HK.cleaned, res = 0.32)

# Merged plots
mergedDim = DimPlot(eng.subset_HK.cleaned, label = TRUE) + ggtitle("Engineered Subset Merged (HK)") +
  theme(plot.title = element_text(hjust = 0.5))
mergedDim_Split = DimPlot(eng.subset_HK.cleaned, label = TRUE, split.by = "Orig_ID", ncol = 3)
merged_eng.HK = mergedDim + mergedDim_Split
merged_eng.HK

# Lineage
FeaturePlot(eng.subset_HK.cleaned, features = c("Epcam", "Cdh5", "Col1a1", "Ptprc"), label = TRUE)

# Look at object
eng.subset_HK.cleaned

# Look at object meta-data
colnames(eng.subset_HK.cleaned@meta.data)
table(eng.subset_HK.cleaned@meta.data$Condition)
table(eng.subset_HK.cleaned@meta.data$Orig_ID)
table(eng.subset_HK.cleaned@meta.data$CellClass)
table(eng.subset_HK.cleaned@meta.data$Phase)
table(eng.subset_HK.cleaned@meta.data$CellType)

# Get a rough idea of current embedding
DimPlot(eng.subset_HK.cleaned,group.by = 'CellClass.integrated')
DimPlot(eng.subset_HK.cleaned,group.by = 'CellType', split.by = 'Orig_ID', ncol = 3)

# Split object first
# tmp[["RNA"]] <- split(tmp[["RNA"]], f = tmp$Orig_ID, layers = "scale.data")
eng.subset_HK.cleaned[["RNA"]] <- split(eng.subset_HK.cleaned[["RNA"]], f = eng.subset_HK.cleaned$Orig_ID, 
                                        layers = c("data", "counts", "scale.data")) # need to split all layers for IntegrateLayers to work

# Embed data, via individual 'layers'
eng.subset_HK.cleaned <- NormalizeData(eng.subset_HK.cleaned)
eng.subset_HK.cleaned <- FindVariableFeatures(eng.subset_HK.cleaned)
eng.subset_HK.cleaned <- ScaleData(eng.subset_HK.cleaned)
eng.subset_HK.cleaned <- RunPCA(eng.subset_HK.cleaned)

# Now we can run integration
eng.subset_HK.cleaned.integrated <- IntegrateLayers(object = eng.subset_HK.cleaned, method = RPCAIntegration,
                                  orig.reduction ="pca", new.reduction = "integrated.rpca", verbose = F)
eng.subset_HK.cleaned.integrated[["RNA"]] <- JoinLayers(eng.subset_HK.cleaned.integrated[["RNA"]])

# RPCA is going to limit batch effect here
# 

# Choose number of PCs
ElbowPlot(eng.subset_HK.cleaned.integrated,ndims=100)+ggtitle('Engineered Subset ElbowPlot')
PCHeatmap(eng.subset_HK.cleaned.integrated,cells=200,balanced=T,dims=1:9)
PCHeatmap(eng.subset_HK.cleaned.integrated,cells=200,balanced=T,dims=10:18)
PCHeatmap(eng.subset_HK.cleaned.integrated,cells=200,balanced=T,dims=19:27)
PCHeatmap(eng.subset_HK.cleaned.integrated,cells=200,balanced=T,dims=28:36)
PCHeatmap(eng.subset_HK.cleaned.integrated,cells=200,balanced=T,dims=37:45) # number of PCs = 25 is appropriate!
PCHeatmap(eng.subset_HK.cleaned.integrated,cells=200,balanced=T,dims=45:50)

eng.subset_HK.cleaned.integrated <- RunUMAP(eng.subset_HK.cleaned.integrated, reduction='integrated.rpca', 
                          dims=1:25, reduction.name = 'umap.rpca')
DimPlot(eng.subset_HK.cleaned.integrated,reduction = 'umap.rpca')

# Cluster
eng.subset_HK.cleaned.integrated <- FindNeighbors(eng.subset_HK.cleaned.integrated, reduction='integrated.rpca', dims= 1:25)
eng.subset_HK.cleaned.integrated <- FindClusters(eng.subset_HK.cleaned.integrated, resolution = 0.4, cluster.name ='rpca_clusters')
d2 = DimPlot(eng.subset_HK.cleaned.integrated, reduction = "umap.rpca") + ggtitle("Engineered Subset Integrated (For Pub)") +
  theme(plot.title = element_text(hjust = 0.5))

# final inspection
integrated_split = DimPlot(eng.subset_HK.cleaned.integrated, reduction = "umap.rpca",
                           group.by = c("CellType", "rpca_clusters"),split.by = 'Condition') 
d2 + integrated_split
ggsave("integrated_split.png", plot = integrated_split, width = 20, height = 12, dpi = 600)

p1 = DimPlot(eng.subset_HK.cleaned.integrated, reduction = "umap.rpca", label = TRUE) + 
  ggtitle("Integrated Engineered Subset") +
  theme(plot.title = element_text(hjust = 0.5))
p1


p2 = DimPlot(eng.subset_HK.cleaned.integrated, reduction = "umap.rpca", label = T, group.by = "CellClass.integrated") + 
  ggtitle("Integrated Engineered Subset by CellClass") +
  theme(plot.title = element_text(hjust = 0.5))
p2
total = p1 + p2
total
ggsave("umap.global_cellclass.png", plot = total, width = 20, height = 10, dpi = 600)

# Classing by Canonical Markers
fp_integrated = FeaturePlot(eng.subset_HK.cleaned.integrated, features = c("Epcam","Ptprc","Cdh5","Col1a1"), reduction = "umap.rpca", label = T)
fp_integrated
ggsave("fp_integrated.png", plot = fp_integrated, width = 10, height = 8, dpi = 600)

# Stash integrated clusters
eng.subset_HK.cleaned.integrated$stash.integrated = Idents(eng.subset_HK.cleaned.integrated)

# Label integrated clusters by cell class (lineage)
cell.epi = WhichCells(eng.subset_HK.cleaned.integrated, idents = c(0,1,2,4,5,6,7,8,9,10,11))
cell.imm = WhichCells(eng.subset_HK.cleaned.integrated, idents = c(3,13))
cell.mes = WhichCells(eng.subset_HK.cleaned.integrated, idents = c(12))
# Adding lineage annotations
eng.subset_HK.cleaned.integrated = SetIdent(eng.subset_HK.cleaned.integrated, cells = cell.epi, value = 'Epithelium')
eng.subset_HK.cleaned.integrated = SetIdent(eng.subset_HK.cleaned.integrated, cells = cell.imm, value = 'Immune')
eng.subset_HK.cleaned.integrated = SetIdent(eng.subset_HK.cleaned.integrated, cells = cell.mes, value = 'Mesenchyme')

# Stash the cell class
eng.subset_HK.cleaned.integrated$CellClass.combined.Integrated = Idents(eng.subset_HK.cleaned.integrated)
table(eng.subset_HK.cleaned.integrated$CellClass.combined.Integrated)
# confirming that everything is labeled
sum(is.na(eng.subset_HK.cleaned.integrated$CellClass.combined.Integrated))

# Restores the cell identities that were previously stashed or stored elsewhere.
Idents(eng.subset_HK.cleaned.integrated) = eng.subset_HK.cleaned.integrated$stash.integrated
table(Idents(eng.subset_HK.cleaned.integrated))

# Now, getting a sense of where things are before we start trying to define cell types by cluster
FeaturePlot(eng.subset_HK.cleaned.integrated, features = c("Krt13","Krt5","Sox9","Sftpc","Napsa","Lamp3","Klf4","Scgb3a2"), reduction = "umap.rpca", order = T, label = T)
FeaturePlot(eng.subset_HK.cleaned.integrated, features = c("Scgb1a1","Sftpc","Scgb3a2"), reduction = "umap.rpca", label = T, order = T)

# Developing marker list
mark1 = FindAllMarkers(eng.subset_HK.cleaned.integrated, min.pct = 0.1,logfc.threshold = 0.1)
mark1$ratio = mark1$pct.1/mark1$pct.2
mark1$power = mark1$ratio*mark1$avg_log2FC
# Top markers
top10_mark1 = mark1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(eng.subset_HK.cleaned.integrated, features = top10_mark1$gene) 
# Viewing marker list and top 10
View(mark1)
View(top10_mark1)

# Set idents to rpca_clusters
Idents(eng.subset_HK.cleaned.integrated) = eng.subset_HK.cleaned.integrated$rpca_clusters

FeaturePlot(eng.subset_HK.cleaned.integrated, features = c("Scgb1a1","Scgb3a2","Sftpc"),reduction = "umap.rpca", label = T)
FeaturePlot(eng.subset_HK.cleaned.integrated, features = c("Top2a","Mki67","Ube2c"),reduction = "umap.rpca", label = T)
FeaturePlot(eng.subset_HK.cleaned.integrated, features = c("Sox9","Cox4i2","Ptges"),reduction = "umap.rpca", label = T)

# Create empty meta-data slot for annotations
eng.subset_HK.cleaned.integrated$CellType.combined.Integrated = NA
# Apply annotations to integrated obj
c0 <- WhichCells(eng.subset_HK.cleaned.integrated, idents = "0")
c1 <- WhichCells(eng.subset_HK.cleaned.integrated, idents = "1")
c2 <- WhichCells(eng.subset_HK.cleaned.integrated, idents = "2")
c3 <- WhichCells(eng.subset_HK.cleaned.integrated, idents = "3")
c4 <- WhichCells(eng.subset_HK.cleaned.integrated, idents = "4")
c5 <- WhichCells(eng.subset_HK.cleaned.integrated, idents = "5")
c6 <- WhichCells(eng.subset_HK.cleaned.integrated, idents = "6")
c7 <- WhichCells(eng.subset_HK.cleaned.integrated, idents = "7")
c8 <- WhichCells(eng.subset_HK.cleaned.integrated, idents = "8")
c9 <- WhichCells(eng.subset_HK.cleaned.integrated, idents = "9")
c10 <- WhichCells(eng.subset_HK.cleaned.integrated, idents = "10")
c11 <- WhichCells(eng.subset_HK.cleaned.integrated, idents = "11")
c12 <- WhichCells(eng.subset_HK.cleaned.integrated, idents = "12")
c13 <- WhichCells(eng.subset_HK.cleaned.integrated, idents = "13")

eng.subset_HK.cleaned.integrated$CellType.combined.Integrated[c0] <- 'ATII_Like'
eng.subset_HK.cleaned.integrated$CellType.combined.Integrated[c1] <- 'Cycling_Distal_Epi'
eng.subset_HK.cleaned.integrated$CellType.combined.Integrated[c2] <- 'Stressed_Progenitor'
eng.subset_HK.cleaned.integrated$CellType.combined.Integrated[c3] <- 'Polarized_Mac'
eng.subset_HK.cleaned.integrated$CellType.combined.Integrated[c4] <- 'Hillock_Like'
eng.subset_HK.cleaned.integrated$CellType.combined.Integrated[c5] <- 'Basal_Like'
eng.subset_HK.cleaned.integrated$CellType.combined.Integrated[c6] <- 'ATI_Like'
eng.subset_HK.cleaned.integrated$CellType.combined.Integrated[c7] <- 'RAS_Like'
eng.subset_HK.cleaned.integrated$CellType.combined.Integrated[c8] <- 'Secretory' 
eng.subset_HK.cleaned.integrated$CellType.combined.Integrated[c9] <- 'Cycling_Proximal_Epi'
eng.subset_HK.cleaned.integrated$CellType.combined.Integrated[c10] <- 'Hillock_Like'
eng.subset_HK.cleaned.integrated$CellType.combined.Integrated[c11] <- 'Ciliated' 
eng.subset_HK.cleaned.integrated$CellType.combined.Integrated[c12] <- 'Rspo3+_Mes'
eng.subset_HK.cleaned.integrated$CellType.combined.Integrated[c13] <- 'Cycling_Immune'

# Check the distribution to confirm all annotations were applied
table(eng.subset_HK.cleaned.integrated$CellType.combined.Integrated)
View(eng.subset_HK.cleaned.integrated@meta.data)
sum(is.na(eng.subset_HK.cleaned.integrated$CellType.combined.Integrated))

# Visualize with annotations
DimPlot(eng.subset_HK.cleaned.integrated, reduction = "umap.rpca", group.by = "CellType.combined.Integrated", label = T)

# Removing unnecessary replicate data, but rename first for simplicity
eng.subset_HK.cleaned.integrated
eng.subset_HK.cleaned.integrated[["RNA"]]$scale.data.BAL_3D_1 <- NULL
eng.subset_HK.cleaned.integrated[["RNA"]]$scale.data.BAL_3D_2 <- NULL
eng.subset_HK.cleaned.integrated[["RNA"]]$scale.data.BAL_3D_3 <- NULL
eng.subset_HK.cleaned.integrated[["RNA"]]$scale.data.PD_3D_1 <- NULL
eng.subset_HK.cleaned.integrated[["RNA"]]$scale.data.PD_3D_2 <- NULL
eng.subset_HK.cleaned.integrated[["RNA"]]$scale.data.PD_3D_3 <- NULL
eng.subset_HK.cleaned.integrated[["RNA"]]$scale.data.Mixed_3D_1 <- NULL
eng.subset_HK.cleaned.integrated[["RNA"]]$scale.data.Mixed_3D_2 <- NULL
eng.subset_HK.cleaned.integrated[["RNA"]]$scale.data.Mixed_3D_3 <- NULL
# Look at object; clean again, much better.
eng.subset_HK.cleaned.integrated

# Rename object before saving
eng.subset.integrated_HK <- eng.subset_HK.cleaned.integrated
# Saving object
save(eng.subset.integrated_HK, file='eng.subset.integrated_HK.Robj')


