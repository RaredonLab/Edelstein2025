
# Subset of Proximal Cluster from Global Object
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

# Set wd to be true and load data
setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets")
load("eng.subset.integrated_HK.Robj") # Merged engineered subset (global)

# Set colors for class
CellClass.cols <- c(
  "Epithelium" = "#D982C6",
  "Immune" = "#87B37A",
  "Mesenchyme" = "#F4A261",
  "Endothelium" = "#2A9D8F",
  "General" = "gray70")

# Define CellType color palette
CellType.cols <- c(
  "ATI_Like" = "#b22222",
   "ATII_Like" = "#4682B4",
  "Secretory" = "#2E8B57",
  "Ciliated" = "#9B489B",
  "Stressed_Progenitor" = "#c71585",
  "Hillock_Like" = "#ff0000","Hillock_Luminal" = "#ff0000",
  "Hillock_Basal" = "#2800c7",
  "Basal_Like" = "#c6e308", "Basal" = "#c6e308",
  "B" = "#41571b",
  "Mac_Alv" = "#00FBFF",
  "Mac_Inter" = "#FF481B",
  "gCaps" = "#A0522D",
  "Tuft" = "#15DDB5",
  "Mesothelium" = "#4A314D",
  "Neutrophils" = "#fcfd1d",
  "pDCs" = "#E9165D",
  "Monocytes" = "#3cd500",
  "NK" = "#F51BEA",
  "T" = "#95B8D1",
  "Arterial" = "#7E5109",
  "Venous" = "#660708",
  "RAS_Like" = "#595db0",
  "Cycling_Epithelium" = "#DB7093",
  "M1_Mac_Like" = "#9497fd",
  "M2_Mac_Like" = "#ff9e80","Polarized_Mac" = "#ff9e80",
  "Rspo3+_Mes" = "#89a5a5","Actc1_Mural" = "#1179fa",
  "Cycling_Distal_Epi" = "#93d0ff",
  "Cycling_Proximal_Epi" = "#c4c2ff",
  "Pdgfrb+_Pericyte" = "#F4AFB4",
  "Cell_Cycle" = "#e0bf1b","Cycling_Immune" = "#e0bf1b",
  "Fzd7+_Stressed" = "#2F004F",
  "LEPs" = "#E2AEDD") # Add LEPs explicitly if needed

# Check structure of data
str(eng.subset.integrated_HK@meta.data)
table(eng.subset.integrated_HK$CellType.combined.Integrated)
# Get a sense of current embedding
DimPlot(eng.subset.integrated_HK, reduction = "umap.rpca", label = T, group.by = "CellType.combined.Integrated", cols = CellType.cols)

# Specify legend order
legend_order <- c(
  "ATI_Like", "ATII_Like", "Basal_Like", "Secretory", "Ciliated", "Stressed_Progenitor",
  "Hillock_Like","RAS_Like","Cycling_Proximal_Epi","Cycling_Distal_Epi",
  "Polarized_Mac", "Cycling_Immune", "Rspo3+_Mes")

# Create a global UMAP
glbl.umap <- DimPlot(
  eng.subset.integrated_HK,
  reduction = "umap.rpca",
  label = TRUE,
  group.by = "CellType.combined.Integrated",
  cols = CellType.cols,
  repel = TRUE) +
  theme(
    plot.title = element_blank(),         # removes the title
    legend.position = "",           # moves legend to the bottom
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank())
glbl.umap
setwd("~/Desktop/Manuscript Figures/Figure 3")
ggsave("glbl.umap.png", plot = glbl.umap, width = 6, height = 6, dpi = 600)

# Want to pull out 4,5,9,10
# Confirm the cluster IDs are treated as characters or factors
levels(eng.subset.integrated_HK@meta.data$rpca_clusters)
# Want to pull out 4,5,9,10
clusters_to_keep <- c("4","5","9","10")
# Subset the object
eng.subset.integ.prox_sub <- subset(eng.subset.integrated_HK, 
                                 subset = rpca_clusters %in% clusters_to_keep)

# Scale (now that we have subset)
eng.subset.integ.prox_sub = ScaleData(eng.subset.integ.prox_sub)
# Normalizing the data
eng.subset.integ.prox_sub = NormalizeData(eng.subset.integ.prox_sub)
# Finding variable features
eng.subset.integ.prox_sub = FindVariableFeatures(eng.subset.integ.prox_sub)

# Principle component analysis
eng.subset.integ.prox_sub = RunPCA(eng.subset.integ.prox_sub, npcs = 100)
ElbowPlot(eng.subset.integ.prox_sub,ndims = 100)
PCHeatmap(eng.subset.integ.prox_sub,cells=200,balanced=T,dims=1:9)
PCHeatmap(eng.subset.integ.prox_sub,cells=200,balanced=T,dims=10:18)
PCHeatmap(eng.subset.integ.prox_sub,cells=200,balanced=T,dims=19:27) # 21 looks good?
PCHeatmap(eng.subset.integ.prox_sub,cells=200,balanced=T,dims=28:36)
PCHeatmap(eng.subset.integ.prox_sub,cells=200,balanced=T,dims=37:45)

# UMAP
eng.subset.integ.prox_sub = RunUMAP(eng.subset.integ.prox_sub, dims = c(1:18))
DimPlot(eng.subset.integ.prox_sub, label = TRUE)
# Cluster data (creating nearest neighbor graph, not clustering) 
eng.subset.integ.prox_sub = FindNeighbors(eng.subset.integ.prox_sub, dims = c(1:18))
# Defining clusters
eng.subset.integ.prox_sub = FindClusters(eng.subset.integ.prox_sub, res = 0.5)
DimPlot(eng.subset.integ.prox_sub, label = TRUE)

# For cleanliness of data, getting rid of 8,10,11
# redo subset of object
eng.subset.integ.prox_sub_1 = subset(eng.subset.integ.prox_sub, idents = c("8","10","11"), invert = TRUE)
# Scale (now that we have subset further)
eng.subset.integ.prox_sub_1 = ScaleData(eng.subset.integ.prox_sub_1)
# Normalizing the data
eng.subset.integ.prox_sub_1 = NormalizeData(eng.subset.integ.prox_sub_1)
# Finding variable features
eng.subset.integ.prox_sub_1 = FindVariableFeatures(eng.subset.integ.prox_sub_1)

# Principle component analysis
eng.subset.integ.prox_sub_1 = RunPCA(eng.subset.integ.prox_sub_1, npcs = 100)
ElbowPlot(eng.subset.integ.prox_sub_1,ndims = 100)
PCHeatmap(eng.subset.integ.prox_sub_1,cells=200,balanced=T,dims=1:9)
PCHeatmap(eng.subset.integ.prox_sub_1,cells=200,balanced=T,dims=10:18)
PCHeatmap(eng.subset.integ.prox_sub_1,cells=200,balanced=T,dims=19:27) # 22 looks good?
PCHeatmap(eng.subset.integ.prox_sub_1,cells=200,balanced=T,dims=28:36)
PCHeatmap(eng.subset.integ.prox_sub_1,cells=200,balanced=T,dims=37:45)

# UMAP
eng.subset.integ.prox_sub_1 = RunUMAP(eng.subset.integ.prox_sub_1, dims = c(1:18))
DimPlot(eng.subset.integ.prox_sub_1, label = TRUE)
# Cluster data (creating nearest neighbor graph, not clustering) 
eng.subset.integ.prox_sub_1 = FindNeighbors(eng.subset.integ.prox_sub_1, dims = c(1:18))
# Defining clusters
eng.subset.integ.prox_sub_1 = FindClusters(eng.subset.integ.prox_sub_1, res = 0.7)
DimPlot(eng.subset.integ.prox_sub_1, label = TRUE)
DimPlot(eng.subset.integ.prox_sub_1, label = TRUE, split.by = "Condition")
FeaturePlot(eng.subset.integ.prox_sub_1, features = c("Top2a","Krt5","Tp63","Krt13","Evpl","Sprr1a"), order = T)

p1 = FeaturePlot(eng.subset.integ.prox_sub_1, features = c("Top2a","Krt5","Krt13","Evpl","Sprr1a","Serpinb2"), order = T)
p1
p2 = DimPlot(eng.subset.integ.prox_sub_1, label = TRUE, split.by = "Orig_ID", ncol = 3)
p2
p1 | p2

# Lets cluster and annotate without annotations
eng.subset.prox = eng.subset.integ.prox_sub_1
FeaturePlot(eng.subset.prox, features = c("Top2a","Krt5","Krt13","Evpl","Sprr1a","Tp63"), order = T)
VlnPlot(eng.subset.prox, features = c("Top2a","Krt5","Krt13","Evpl","Sprr1a","Tp63"))
# Okay so breakdown as follows
# Cycling_Proximal_Epi = 4,7
# Basal = 0,10,11,5
# Hillock_Basal = 3,9,12,13
# Hillock_Luminal = 1,2,6,8

eng.subset.prox[["CellType_ProxSubset_no.int"]] <- rep(NA, ncol(eng.subset.prox))
c0 <- WhichCells(eng.subset.prox, idents = "0")
c1 <- WhichCells(eng.subset.prox, idents = "1")
c2 <- WhichCells(eng.subset.prox, idents = "2")
c3 <- WhichCells(eng.subset.prox, idents = "3")
c4 <- WhichCells(eng.subset.prox, idents = "4")
c5 <- WhichCells(eng.subset.prox, idents = "5")
c6 <- WhichCells(eng.subset.prox, idents = "6")
c7 <- WhichCells(eng.subset.prox, idents = "7")
c8 <- WhichCells(eng.subset.prox, idents = "8")
c9 <- WhichCells(eng.subset.prox, idents = "9")
c10 <- WhichCells(eng.subset.prox, idents = "10")
c11 <- WhichCells(eng.subset.prox, idents = "11")
c12 <- WhichCells(eng.subset.prox, idents = "12")
c13 <- WhichCells(eng.subset.prox, idents = "13")
eng.subset.prox$CellType_ProxSubset_no.int[c0] <- 'Basal_Like'
eng.subset.prox$CellType_ProxSubset_no.int[c1] <- 'Hillock_Luminal'
eng.subset.prox$CellType_ProxSubset_no.int[c2] <- 'Hillock_Luminal'
eng.subset.prox$CellType_ProxSubset_no.int[c3] <- 'Hillock_Basal'
eng.subset.prox$CellType_ProxSubset_no.int[c4] <- 'Cycling_Proximal_Epi'
eng.subset.prox$CellType_ProxSubset_no.int[c5] <- 'Basal_Like'
eng.subset.prox$CellType_ProxSubset_no.int[c6] <- 'Hillock_Luminal'
eng.subset.prox$CellType_ProxSubset_no.int[c7] <- 'Cycling_Proximal_Epi'
eng.subset.prox$CellType_ProxSubset_no.int[c8] <- 'Hillock_Luminal' 
eng.subset.prox$CellType_ProxSubset_no.int[c9] <- 'Hillock_Basal'
eng.subset.prox$CellType_ProxSubset_no.int[c10] <- 'Basal_Like'
eng.subset.prox$CellType_ProxSubset_no.int[c11] <- 'Basal_Like'
eng.subset.prox$CellType_ProxSubset_no.int[c12] <- 'Hillock_Basal'
eng.subset.prox$CellType_ProxSubset_no.int[c13] <- 'Hillock_Basal'

# Check the distribution to confirm all annotations were applied
table(eng.subset.prox$CellType_ProxSubset_no.int)
View(eng.subset.prox@meta.data)
sum(is.na(eng.subset.prox$CellType_ProxSubset_no.int))
# Looks good

# Marker list
Idents(eng.subset.prox) = eng.subset.prox$Condition
mark.sub_cond = FindAllMarkers(eng.subset.prox, min.pct = 0.1,logfc.threshold = 0.1)
mark.sub_cond$ratio = mark.sub_cond$pct.1/mark.sub_cond$pct.2
mark.sub_cond$power = mark.sub_cond$ratio*mark.sub_cond$avg_log2FC
View(mark.sub_cond)

Idents(eng.subset.prox) = eng.subset.prox$Condition
mark.sub.pol = FindMarkers(eng.subset.prox, ident.1 = "BAL_3D", ident.2 = c("Mixed_3D", "PD_3D"), min.pct = 0.1,logfc.threshold = 0.1)
mark.sub.pol$ratio = mark.sub.pol$pct.1/mark.sub.pol$pct.2
mark.sub.pol$power = mark.sub.pol$ratio*mark.sub.pol$avg_log2FC
View(mark.sub.pol)


# Plot with annotations
DimPlot(eng.subset.prox, reduction = "umap", label = TRUE, group.by = "CellType_ProxSubset_no.int") + ggtitle("Proximal Subset w/Annotations") +
  theme(plot.title = element_text(hjust = 0.5))

setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets")
setwd("~/Desktop/Single Cell")
save(eng.subset.prox, file = "eng.subset.prox.Robj")


# Integrate?
# Copy before integration
eng.subset.prox = eng.subset.integ.prox_sub_1

###### RE-INTEGRATE #######
# Split object first
eng.subset.prox[["RNA"]] <- split(eng.subset.prox[["RNA"]], f = eng.subset.prox$Orig_ID, layers = c("data", "counts", "scale.data")) # need to split all layers for IntegrateLayers to work

# Embed data, via individual 'layers'
eng.subset.prox <- NormalizeData(eng.subset.prox)
eng.subset.prox <- FindVariableFeatures(eng.subset.prox)
eng.subset.prox <- ScaleData(eng.subset.prox)
eng.subset.prox <- RunPCA(eng.subset.prox)

# Now we can run integration
eng.subset.prox.int <- IntegrateLayers(object = eng.subset.prox, method = RPCAIntegration,
                                          orig.reduction ="pca", new.reduction = "integrated.rpca", verbose = F, k.weight = 50)
eng.subset.prox.int[["RNA"]] <- JoinLayers(eng.subset.prox.int[["RNA"]])

# Choose number of PCs
ElbowPlot(eng.subset.prox.int,ndims=50)+ggtitle('Engineered Proximal Subset Integrated')
PCHeatmap(eng.subset.prox.int,cells=200,balanced=T,dims=1:9)
PCHeatmap(eng.subset.prox.int,cells=200,balanced=T,dims=10:18) # 15 PCs
PCHeatmap(eng.subset.prox.int,cells=200,balanced=T,dims=19:27)
PCHeatmap(eng.subset.prox.int,cells=200,balanced=T,dims=28:36)
PCHeatmap(eng.subset.prox.int,cells=200,balanced=T,dims=37:45) 

# Check final results
eng.subset.prox.int <- FindNeighbors(eng.subset.prox.int, reduction='integrated.rpca', dims=1:16)
eng.subset.prox.int <- FindClusters(eng.subset.prox.int, resolution = 1.2, cluster.name ='rpca_clusters')
eng.subset.prox.int <- RunUMAP(eng.subset.prox.int, reduction='integrated.rpca', dims=1:16, reduction.name = 'umap.rpca')

a = DimPlot(eng.subset.prox.int, reduction = "umap.rpca", label = TRUE)
a
b = FeaturePlot(eng.subset.prox.int, features = c("Top2a","Krt5","Krt13","Evpl","Sprr1a","Serpinb2"), reduction = "umap.rpca", ncol = 2, label = T, order = T)
a | b
FeaturePlot(eng.subset.prox.int, features = c("Top2a","Krt5","Krt13","Evpl","Sprr1a","Serpinb2"), reduction = "umap.rpca", ncol = 2, label = T)
VlnPlot(eng.subset.prox.int, features = c("Top2a","Krt5","Krt13","Evpl","Sprr1a","Tp63"))

VlnPlot(eng.subset.integ.prox_sub_1, features = c("Top2a","Krt5","Krt13","Evpl","Sprr1a","Tp63"))
DimPlot(eng.subset.integ.prox_sub_1, label = T)
# Cycling = 8,10
# Basal = 1,3,7,10,11,14
# Basal Hillock (Krt13+/Krt5+/Tp63+) = 5,6,9
# Luminal Hillock (Krt13+/Krt5-) = 0,2,4,12,13
DimPlot(eng.subset.prox.int, reduction = "umap.rpca", label = TRUE, split.by = "Condition")

# Let's look at our subset object actually
Idents(eng.subset.prox.int) = eng.subset.prox.int$Condition
mark.sub = FindAllMarkers(eng.subset.prox.int, min.pct = 0.1,logfc.threshold = 0.1)
mark.sub$ratio = mark.sub$pct.1/mark.sub$pct.2
mark.sub$power = mark.sub$ratio*mark.sub$avg_log2FC
View(mark.sub)


# Okay so breakdown as follows
# Clusters 8 + 10: Cycling_Proximal_Epi
# Clusters 9,6,3,4,12 = Basal
# Clusters 1: Hillock_Basal
# Clusters 11,2,0,5,7: Hillock_Luminal

eng.subset.prox.int[["CellType_Prox_Subset"]] <- rep(NA, ncol(eng.subset.prox.int))
c0 <- WhichCells(eng.subset.prox.int, idents = "0")
c1 <- WhichCells(eng.subset.prox.int, idents = "1")
c2 <- WhichCells(eng.subset.prox.int, idents = "2")
c3 <- WhichCells(eng.subset.prox.int, idents = "3")
c4 <- WhichCells(eng.subset.prox.int, idents = "4")
c5 <- WhichCells(eng.subset.prox.int, idents = "5")
c6 <- WhichCells(eng.subset.prox.int, idents = "6")
c7 <- WhichCells(eng.subset.prox.int, idents = "7")
c8 <- WhichCells(eng.subset.prox.int, idents = "8")
c9 <- WhichCells(eng.subset.prox.int, idents = "9")
c10 <- WhichCells(eng.subset.prox.int, idents = "10")
c11 <- WhichCells(eng.subset.prox.int, idents = "11")
c12 <- WhichCells(eng.subset.prox.int, idents = "12")
eng.subset.prox.int$CellType_Prox_Subset[c0] <- 'Hillock_Luminal'
eng.subset.prox.int$CellType_Prox_Subset[c1] <- 'Hillock_Basal'
eng.subset.prox.int$CellType_Prox_Subset[c2] <- 'Hillock_Luminal'
eng.subset.prox.int$CellType_Prox_Subset[c3] <- 'Basal_Like'
eng.subset.prox.int$CellType_Prox_Subset[c4] <- 'Basal_Like'
eng.subset.prox.int$CellType_Prox_Subset[c5] <- 'Hillock_Luminal'
eng.subset.prox.int$CellType_Prox_Subset[c6] <- 'Basal_Like'
eng.subset.prox.int$CellType_Prox_Subset[c7] <- 'Hillock_Luminal'
eng.subset.prox.int$CellType_Prox_Subset[c8] <- 'Cycling_Proximal_Epi' 
eng.subset.prox.int$CellType_Prox_Subset[c9] <- 'Basal_Like'
eng.subset.prox.int$CellType_Prox_Subset[c10] <- 'Cycling_Proximal_Epi'
eng.subset.prox.int$CellType_Prox_Subset[c11] <- 'Hillock_Luminal'
eng.subset.prox.int$CellType_Prox_Subset[c12] <- 'Basal_Like'

# Check the distribution to confirm all annotations were applied
table(eng.subset.prox.int$CellType_Prox_Subset)
View(eng.subset.prox.int@meta.data)
sum(is.na(eng.subset.prox.int$CellType_Prox_Subset))
# Looks good

# Plot with annotations
DimPlot(eng.subset.prox.int, reduction = "umap.rpca", label = TRUE, group.by = "CellType_Prox_Subset") + ggtitle("Proximal Subset w/Annotations") +
  theme(plot.title = element_text(hjust = 0.5))

setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets")
setwd("~/Desktop/Single Cell")
# Save
save(eng.subset.prox.int, file = "eng.subset.prox.int.Robj")

### Making ModuleScore for Hillock_Basal + Hillock_Luminal
hillock_basal_genes <- c("Krt5","Aqp3","Il33","Snai2","Ces2g",
                         "Lgals7","Tppp3","Jag2","Mt2A","Tp63","Aldh3a1")
hillock_luminal_genes <- c("Krt13", "Sprr1a","Ppl","Mal","S100a7a","Evpl","Cnfn","Wnt11","Alox5ap","S100a8",
                           "Krt7","Plat","Aldh1a3","Krt80","Slpi")
eng.subset.integrated_HK <- AddModuleScore(
  object = eng.subset.integrated_HK,
  features = list(hillock_basal_genes, hillock_luminal_genes),
  name = c("Hillock_Basal_Score", "Hillock_Luminal_Score"))
VlnPlot(eng.subset.integrated_HK, features = c("Hillock_Basal_Score1", "Hillock_Luminal_Score2"), group.by = "Condition", pt.size = 0) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"))

VlnPlot(eng.subset.prox.int, features = c("Krt13", "Krt5","Tp63","Sprr1a","Slpi", "Ppl","Evpl","Smad7"),
        split.by = "Condition", pt.size = 0) +
  coord_cartesian(ylim = c(0, 3))  # crop y-axis visually

# Run marker list by cell type
Idents(eng.subset.prox.int) = eng.subset.prox.int$CellType_Prox_Subset
mark2 = FindAllMarkers(eng.subset.prox.int, min.pct = 0.1,logfc.threshold = 0.1)
mark2$ratio = mark2$pct.1/mark2$pct.2
mark2$power = mark2$ratio*mark2$avg_log2FC
View(mark2)

a = FeaturePlot(PD_3D.integrated_HK, features = c("Krt13","Krt5","Ppl","Sprr1a"), reduction = "umap.rpca", ncol = 4, max.cutoff = 4)
b = FeaturePlot(Mixed_3D.integrated_HK, features = c("Krt13","Tp63","Krt5","Ppl","Sprr1a"), reduction = "umap.rpca", ncol = 4, max.cutoff = 4)
c = FeaturePlot(BAL_3D.integrated_HK, features = c("Krt13","Tp63","Krt5","Ppl","Sprr1a"), reduction = "umap.rpca", ncol = 4, max.cutoff = 4)
a / b / c

PD_3D.integrated_HK
Mixed_3D.integrated_HK
BAL_3D.integrated_HK

# Run marker list by condition for exploration
Idents(eng.subset.prox.int) = eng.subset.prox.int$Condition
dge_results_cond = FindAllMarkers(eng.subset.prox.int, min.pct = 0.1,logfc.threshold = 0.1)
dge_results_cond$ratio = dge_results_cond$pct.1/dge_results_cond$pct.2
dge_results_cond$power = dge_results_cond$ratio*dge_results_cond$avg_log2FC
View(dge_results_cond)

# Merged plots
mergedDim = DimPlot(eng.subset.integ.prox_sub, label = TRUE) + ggtitle("Engineered Subset Merged (HK)") +
  theme(plot.title = element_text(hjust = 0.5))
mergedDim_Split = DimPlot(eng.subset.integ.prox_sub, label = TRUE, split.by = "Orig_ID", ncol = 3)
merged_eng.HK = mergedDim + mergedDim_Split
merged_eng.HK

# Lineage
FeaturePlot(eng.subset_HK.cleaned, features = c("Epcam", "Cdh5", "Col1a1", "Ptprc"), label = TRUE)
