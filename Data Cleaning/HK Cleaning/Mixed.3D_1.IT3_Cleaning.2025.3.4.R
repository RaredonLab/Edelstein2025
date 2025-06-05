#####Mixed_3D_1 Cleaning


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/Mixed_3D_1")

library(dplyr)
library(Seurat)
library(patchwork)

# Create a starting Seurat object, loading and using the provided raw data.  
Mixed.3D_1.data <- Read10X(data.dir = "/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/Mixed_3D_1/Raw_Data")

Mixed.3D_1 <- CreateSeuratObject(counts = Mixed.3D_1.data)

Mixed.3D_1[["percent.mt"]] <- PercentageFeatureSet(Mixed.3D_1, pattern = "^Mt-")

# An object of class Seurat 
# 32883 features across 2008437 samples within 1 assay 
# Active assay: RNA (32883 features, 0 variable features)
# 1 layer present: counts

# We will use this to determine where to set our thresholds for nCount_RNA and percent.mt
setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/Mixed_3D_1/Plots")
png("Mixed.3D_1_Pre_Thresholding.VlnPlot.png", res=300, unit="in", height=8, width=11)
VlnPlot(Mixed.3D_1, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), raster = FALSE)
dev.off()

png("Mixed.3D_1_Pre_Thresholding.VlnPlot_Log.png", res=300, unit="in", height=8, width=11)
VlnPlot(Mixed.3D_1, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), log = TRUE, raster = FALSE)
dev.off()

VlnPlot(Mixed.3D_1, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), raster = FALSE, pt.size = 0)

Mixed.3D_1<-subset(Mixed.3D_1, subset = nCount_RNA > 500 & percent.mt < 25)

# An object of class Seurat 
# 32883 features across 8399 samples within 1 assay 
# Active assay: RNA (32883 features, 0 variable features)
# 1 layer present: counts

# Again, create a separate set of Vln Plot images for our newly thresholded object
png("Mixed.3D_1_VlnPlotPostSubset.png", res=300, unit="in", height=8, width=11)
VlnPlot(Mixed.3D_1, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
dev.off()

png("Mixed.3D_1_VlnPlot_LogPostSubset.png", res=300, unit="in", height=8, width=11)
VlnPlot(Mixed.3D_1, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), log = TRUE)
dev.off()


Mixed.3D_1<- NormalizeData(Mixed.3D_1, verbose=FALSE)

Mixed.3D_1<- FindVariableFeatures(Mixed.3D_1, verbose=FALSE, nfeatures=2000)
Mixed.3D_1<- ScaleData(Mixed.3D_1, features = VariableFeatures(Mixed.3D_1))
Mixed.3D_1<- RunPCA(Mixed.3D_1, features = VariableFeatures(Mixed.3D_1), npcs =100)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/Mixed_3D_1/PCA/Iteration_1")
png("Mixed.3D_1_PCA1.png", res=300, unit="in", height=8, width=11)
DimHeatmap(Mixed.3D_1, dims = 1:9, cells = 500, balanced = TRUE)
dev.off() 

png("Mixed.3D_1_PCA2.png", res=300, unit="in", height=8, width=11)
DimHeatmap(Mixed.3D_1, dims = 10:18, cells = 500, balanced = TRUE)
dev.off()

png("Mixed.3D_1_PCA3.png", res=300, unit="in", height=8, width=11)
DimHeatmap(Mixed.3D_1, dims = 19:27, cells = 500, balanced = TRUE)
dev.off()

png("Mixed.3D_1_PCA4.png", res=300, unit="in", height=8, width=11)
DimHeatmap(Mixed.3D_1, dims = 28:36, cells = 500, balanced = TRUE)
dev.off()

png("Mixed.3D_1_PCA5.png", res=300, unit="in", height=8, width=11)
DimHeatmap(Mixed.3D_1, dims = 37:45, cells = 500, balanced = TRUE)
dev.off()

png("Mixed.3D_1_PCA6.png", res=300, unit="in", height=8, width=11)
DimHeatmap(Mixed.3D_1, dims = 46:54, cells = 500, balanced = TRUE)
dev.off()

png(filename = 'Mixed.3D_1.ElbowPlot.png', res=300, unit="in", height=8, width=11)
ElbowPlot(Mixed.3D_1, ndims = 100)
dev.off()

## We will use PCs 1:36 for our embedding.

dims.use <- c(1:37, 38:43, 45)



Mixed.3D_1 <- FindNeighbors(Mixed.3D_1,
                         reduction="pca",
                         k.param = 20,
                         dims=dims.use
                         
)

Mixed.3D_1 <- FindClusters(Mixed.3D_1,
                        resolution= 2.5
                        
)

Mixed.3D_1<- RunUMAP(Mixed.3D_1, 
                  dims=dims.use,
                  min.dist = .2
)
# We create a DimPlot to show our resulting UMAP.
DimPlot(Mixed.3D_1, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/Mixed_3D_1/Plots/UMAP/Iteration_1")

png("Mixed.3D_1_UMAP_1.png", res=300, unit="in", height=8, width=11)
DimPlot(Mixed.3D_1, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE, repel = TRUE) + NoLegend()
dev.off()

png("Mixed.3D_1_nCount_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(Mixed.3D_1, features = "nCount_RNA", label = TRUE, repel = TRUE)
dev.off()

png("Mixed.3D_1_nFeature_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(Mixed.3D_1, features = "nFeature_RNA", label = TRUE, repel = TRUE)
dev.off()

png("Mixed.3D_1_percent.mt.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(Mixed.3D_1, features = "percent.mt", label = TRUE, repel = TRUE)
dev.off()

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/Mixed_3D_1/Plots/UMAP/Iteration_1")

png("Mixed.3D_1_cluster_VlnPlot.nCount_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(Mixed.3D_1, features = "nCount_RNA", split.by = "seurat_clusters")
dev.off()

png("Mixed.3D_1_cluster_VlnPlot.nFeature_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(Mixed.3D_1, features = "nFeature_RNA", split.by = "seurat_clusters")
dev.off()

png("Mixed.3D_1_cluster_VlnPlot.percent.mt.png", res=300, unit="in", height=8, width=11)
VlnPlot(Mixed.3D_1, features = "percent.mt", split.by = "seurat_clusters")
dev.off()


Mixed.3D_1.Markers <- FindAllMarkers(Mixed.3D_1, logfc.threshold = 0.25, only.pos = T)
Mixed.3D_1.Markers$ratio <- Mixed.3D_1.Markers$pct.1/Mixed.3D_1.Markers$pct.2
Mixed.3D_1.Markers$power <- Mixed.3D_1.Markers$ratio*Mixed.3D_1.Markers$avg_log2FC


Mixed.3D_1.IT2<- subset(Mixed.3D_1, idents = c(0, 1, 2, 4, 5, 10), invert = T)

# An object of class Seurat 
# 32883 features across 5099 samples within 1 assay 
# Active assay: RNA (32883 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap


Mixed.3D_1.IT2<- FindVariableFeatures(Mixed.3D_1.IT2, verbose=FALSE, nfeatures=2000)
Mixed.3D_1.IT2<- ScaleData(Mixed.3D_1.IT2, features = VariableFeatures(Mixed.3D_1.IT2))
Mixed.3D_1.IT2<- RunPCA(Mixed.3D_1.IT2, features = VariableFeatures(Mixed.3D_1.IT2), npcs =100)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/Mixed_3D_1/PCA/Iteration_2")
png("Mixed.3D_1_PCA1.png", res=300, unit="in", height=8, width=11)
DimHeatmap(Mixed.3D_1.IT2, dims = 1:9, cells = 500, balanced = TRUE)
dev.off() 

png("Mixed.3D_1_PCA2.png", res=300, unit="in", height=8, width=11)
DimHeatmap(Mixed.3D_1.IT2, dims = 10:18, cells = 500, balanced = TRUE)
dev.off()

png("Mixed.3D_1_PCA3.png", res=300, unit="in", height=8, width=11)
DimHeatmap(Mixed.3D_1.IT2, dims = 19:27, cells = 500, balanced = TRUE)
dev.off()

png("Mixed.3D_1_PCA4.png", res=300, unit="in", height=8, width=11)
DimHeatmap(Mixed.3D_1.IT2, dims = 28:36, cells = 500, balanced = TRUE)
dev.off()

png("Mixed.3D_1_PCA5.png", res=300, unit="in", height=8, width=11)
DimHeatmap(Mixed.3D_1.IT2, dims = 37:45, cells = 500, balanced = TRUE)
dev.off()

png("Mixed.3D_1_PCA6.png", res=300, unit="in", height=8, width=11)
DimHeatmap(Mixed.3D_1.IT2, dims = 46:54, cells = 500, balanced = TRUE)
dev.off()

png(filename = 'Mixed.3D_1.ElbowPlot.png', res=300, unit="in", height=8, width=11)
ElbowPlot(Mixed.3D_1.IT2, ndims = 100)
dev.off()

## We will use PCs 1:46 for our embedding.

dims.use <- c(1:46)



Mixed.3D_1.IT2 <- FindNeighbors(Mixed.3D_1.IT2,
                             reduction="pca",
                             k.param = 20,
                             dims=dims.use
                             
)

Mixed.3D_1.IT2 <- FindClusters(Mixed.3D_1.IT2,
                            resolution= 2.5
                            
)

Mixed.3D_1.IT2<- RunUMAP(Mixed.3D_1.IT2, 
                      dims=dims.use,
                      min.dist = .2
)
# We create a DimPlot to show our resulting UMAP.
DimPlot(Mixed.3D_1.IT2, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/Mixed_3D_1/Plots/UMAP/Iteration_2")

png("Mixed.3D_1_UMAP_1.png", res=300, unit="in", height=8, width=11)
DimPlot(Mixed.3D_1.IT2, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE, repel = TRUE) + NoLegend()
dev.off()

png("Mixed.3D_1_nCount_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(Mixed.3D_1.IT2, features = "nCount_RNA", label = TRUE, repel = TRUE)
dev.off()

png("Mixed.3D_1_nFeature_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(Mixed.3D_1.IT2, features = "nFeature_RNA", label = TRUE, repel = TRUE)
dev.off()

png("Mixed.3D_1_percent.mt.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(Mixed.3D_1.IT2, features = "percent.mt", label = TRUE, repel = TRUE)
dev.off()

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/Mixed_3D_1/Plots/UMAP/Iteration_2")

png("Mixed.3D_1_cluster_VlnPlot.nCount_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(Mixed.3D_1.IT2, features = "nCount_RNA", split.by = "seurat_clusters")
dev.off()

png("Mixed.3D_1_cluster_VlnPlot.nFeature_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(Mixed.3D_1.IT2, features = "nFeature_RNA", split.by = "seurat_clusters")
dev.off()

png("Mixed.3D_1_cluster_VlnPlot.percent.mt.png", res=300, unit="in", height=8, width=11)
VlnPlot(Mixed.3D_1.IT2, features = "percent.mt", split.by = "seurat_clusters")
dev.off()


Mixed.3D_1.IT2.Markers <- FindAllMarkers(Mixed.3D_1.IT2, logfc.threshold = 0.25, only.pos = T)
Mixed.3D_1.IT2.Markers$ratio <- Mixed.3D_1.IT2.Markers$pct.1/Mixed.3D_1.IT2.Markers$pct.2
Mixed.3D_1.IT2.Markers$power <- Mixed.3D_1.IT2.Markers$ratio*Mixed.3D_1.IT2.Markers$avg_log2FC

Mixed.3D_1.IT3<- subset(Mixed.3D_1.IT2, idents = c(1), invert = T)

# An object of class Seurat 
# 32883 features across 4794 samples within 1 assay 
# Active assay: RNA (32883 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap


Mixed.3D_1.IT3<- FindVariableFeatures(Mixed.3D_1.IT3, verbose=FALSE, nfeatures=2000)
Mixed.3D_1.IT3<- ScaleData(Mixed.3D_1.IT3, features = VariableFeatures(Mixed.3D_1.IT3))
Mixed.3D_1.IT3<- RunPCA(Mixed.3D_1.IT3, features = VariableFeatures(Mixed.3D_1.IT3), npcs =100)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/Mixed_3D_1/PCA/Iteration_3")
png("BAL3D_3_PCA1.png", res=300, unit="in", height=8, width=11)
DimHeatmap(Mixed.3D_1.IT3, dims = 1:9, cells = 500, balanced = TRUE)
dev.off() 

png("BAL3D_3_PCA2.png", res=300, unit="in", height=8, width=11)
DimHeatmap(Mixed.3D_1.IT3, dims = 10:18, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_3_PCA3.png", res=300, unit="in", height=8, width=11)
DimHeatmap(Mixed.3D_1.IT3, dims = 19:27, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_3_PCA4.png", res=300, unit="in", height=8, width=11)
DimHeatmap(Mixed.3D_1.IT3, dims = 28:36, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_3_PCA5.png", res=300, unit="in", height=8, width=11)
DimHeatmap(Mixed.3D_1.IT3, dims = 37:45, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_3_PCA6.png", res=300, unit="in", height=8, width=11)
DimHeatmap(Mixed.3D_1.IT3, dims = 46:54, cells = 500, balanced = TRUE)
dev.off()

png(filename = 'BAL3D_3.ElbowPlot.png', res=300, unit="in", height=8, width=11)
ElbowPlot(Mixed.3D_1.IT3, ndims = 100)
dev.off()

## We will use PCs 1:45 for our embedding.

dims.use <- c(1:45)



Mixed.3D_1.IT3 <- FindNeighbors(Mixed.3D_1.IT3,
                             reduction="pca",
                             k.param = 20,
                             dims=dims.use
                             
)

Mixed.3D_1.IT3 <- FindClusters(Mixed.3D_1.IT3,
                            resolution= 2.5
                            
)

Mixed.3D_1.IT3<- RunUMAP(Mixed.3D_1.IT3, 
                      dims=dims.use,
                      min.dist = .2
)
# We create a DimPlot to show our resulting UMAP.
DimPlot(Mixed.3D_1.IT3, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/Mixed_3D_1/Plots/UMAP/Iteration_3")

png("BAL3D_3_UMAP_1.png", res=300, unit="in", height=8, width=11)
DimPlot(Mixed.3D_1.IT3, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE, repel = TRUE) + NoLegend()
dev.off()

png("BAL3D_3_nCount_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(Mixed.3D_1.IT3, features = "nCount_RNA", label = TRUE, repel = TRUE)
dev.off()

png("BAL3D_3_nFeature_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(Mixed.3D_1.IT3, features = "nFeature_RNA", label = TRUE, repel = TRUE)
dev.off()

png("BAL3D_3_percent.mt.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(Mixed.3D_1.IT3, features = "percent.mt", label = TRUE, repel = TRUE)
dev.off()

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/Mixed_3D_1//Plots/UMAP/Iteration_3")

png("BAL3D_3_cluster_VlnPlot.nCount_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(Mixed.3D_1.IT3, features = "nCount_RNA", split.by = "seurat_clusters")
dev.off()

png("BAL3D_3_cluster_VlnPlot.nFeature_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(Mixed.3D_1.IT3, features = "nFeature_RNA", split.by = "seurat_clusters")
dev.off()

png("BAL3D_3_cluster_VlnPlot.percent.mt.png", res=300, unit="in", height=8, width=11)
VlnPlot(Mixed.3D_1.IT3, features = "percent.mt", split.by = "seurat_clusters")
dev.off()


Mixed.3D_1.IT3.Markers <- FindAllMarkers(Mixed.3D_1.IT3, logfc.threshold = 0.25, only.pos = T)
Mixed.3D_1.IT3.Markers$ratio <- Mixed.3D_1.IT3.Markers$pct.1/Mixed.3D_1.IT3.Markers$pct.2
Mixed.3D_1.IT3.Markers$power <- Mixed.3D_1.IT3.Markers$ratio*Mixed.3D_1.IT3.Markers$avg_log2FC






Mixed.3D_1.IT3 <- FindClusters(Mixed.3D_1.IT3,
                            resolution= 1
                            
)

Mixed.3D_1.IT3<- RunUMAP(Mixed.3D_1.IT3, 
                      dims=dims.use,
                      min.dist = .2
)

DimPlot(Mixed.3D_1.IT3, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE)
save(Mixed.3D_1.IT3, file = "Mixed.3D_1.IT3.2025.3.4.Robj")
