#####BAL_3D_2 Cleaning


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_2")

library(dplyr)
library(Seurat)
library(patchwork)

# Create a starting Seurat object, loading and using the provided raw data.  
BAL3D_2.data <- Read10X(data.dir = "/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_2/Raw_Data")

BAL3D_2 <- CreateSeuratObject(counts = BAL3D_2.data)

BAL3D_2[["percent.mt"]] <- PercentageFeatureSet(BAL3D_2, pattern = "^Mt-")

# An object of class Seurat 
# 32883 features across 1094517 samples within 1 assay 
# Active assay: RNA (32883 features, 0 variable features)
# 1 layer present: counts

# We will use this to determine where to set our thresholds for nCount_RNA and percent.mt
setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_2/Plots")
png("BAL3D_2_Pre_Thresholding.VlnPlot.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_2, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), raster = FALSE)
dev.off()

png("BAL3D_2_Pre_Thresholding.VlnPlot_Log.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_2, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), log = TRUE, raster = FALSE)
dev.off()


BAL3D_2<-subset(BAL3D_2, subset = nCount_RNA > 500 & percent.mt < 25)

# An object of class Seurat 
# 32883 features across 3968 samples within 1 assay 
# Active assay: RNA (32883 features, 0 variable features)
# 1 layer present: counts

# Again, create a separate set of Vln Plot images for our newly thresholded object
png("BAL3D_2_VlnPlotPostSubset.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_2, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
dev.off()

png("BAL3D_2_VlnPlot_LogPostSubset.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_2, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), log = TRUE)
dev.off()


BAL3D_2<- NormalizeData(BAL3D_2, verbose=FALSE)

BAL3D_2<- FindVariableFeatures(BAL3D_2, verbose=FALSE, nfeatures=2000)
BAL3D_2<- ScaleData(BAL3D_2, features = VariableFeatures(BAL3D_2))
BAL3D_2<- RunPCA(BAL3D_2, features = VariableFeatures(BAL3D_2), npcs =100)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_2/PCA/Iteration_1")
png("BAL3D_2_PCA1.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_2, dims = 1:9, cells = 500, balanced = TRUE)
dev.off() 

png("BAL3D_2_PCA2.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_2, dims = 10:18, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_2_PCA3.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_2, dims = 19:27, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_2_PCA4.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_2, dims = 28:36, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_2_PCA5.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_2, dims = 37:45, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_2_PCA6.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_2, dims = 46:54, cells = 500, balanced = TRUE)
dev.off()

png(filename = 'BAL3D_2.ElbowPlot.png', res=300, unit="in", height=8, width=11)
ElbowPlot(BAL3D_2, ndims = 100)
dev.off()

## We will use PCs 1:22 through 24:27 for our first embedding.

dims.use <- c(1:30)



BAL3D_2 <- FindNeighbors(BAL3D_2,
                         reduction="pca",
                         k.param = 20,
                         dims=dims.use
                         
)

BAL3D_2 <- FindClusters(BAL3D_2,
                        resolution= 2.5
                        
)

BAL3D_2<- RunUMAP(BAL3D_2, 
                  dims=dims.use,
                  min.dist = .2
)
# We create a DimPlot to show our resulting UMAP.
DimPlot(BAL3D_2, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_2/Plots/UMAP/Iteration_1")

png("BAL3D_2_UMAP_1.png", res=300, unit="in", height=8, width=11)
DimPlot(BAL3D_2, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE, repel = TRUE) + NoLegend()
dev.off()

png("BAL3D_2_nCount_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BAL3D_2, features = "nCount_RNA", label = TRUE, repel = TRUE)
dev.off()

png("BAL3D_2_nFeature_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BAL3D_2, features = "nFeature_RNA", label = TRUE, repel = TRUE)
dev.off()

png("BAL3D_2_percent.mt.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BAL3D_2, features = "percent.mt", label = TRUE, repel = TRUE)
dev.off()

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_2/Plots/UMAP/Iteration_1")

png("BAL3D_2_cluster_VlnPlot.nCount_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_2, features = "nCount_RNA", split.by = "seurat_clusters")
dev.off()

png("BAL3D_2_cluster_VlnPlot.nFeature_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_2, features = "nFeature_RNA", split.by = "seurat_clusters")
dev.off()

png("BAL3D_2_cluster_VlnPlot.percent.mt.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_2, features = "percent.mt", split.by = "seurat_clusters")
dev.off()


BAL3D_2.Markers <- FindAllMarkers(BAL3D_2, logfc.threshold = 0.25, only.pos = T)
BAL3D_2.Markers$ratio <- BAL3D_2.Markers$pct.1/BAL3D_2.Markers$pct.2
BAL3D_2.Markers$power <- BAL3D_2.Markers$ratio*BAL3D_2.Markers$avg_log2FC


BAL3D_2.IT2<- subset(BAL3D_2, idents = c(10), invert = T)

# An object of class Seurat 
# 32883 features across 3807 samples within 1 assay 
# Active assay: RNA (32883 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap


BAL3D_2.IT2<- FindVariableFeatures(BAL3D_2.IT2, verbose=FALSE, nfeatures=2000)
BAL3D_2.IT2<- ScaleData(BAL3D_2.IT2, features = VariableFeatures(BAL3D_2.IT2))
BAL3D_2.IT2<- RunPCA(BAL3D_2.IT2, features = VariableFeatures(BAL3D_2.IT2), npcs =100)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_2/PCA/Iteration_2")
png("BAL3D_2_PCA1.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_2.IT2, dims = 1:9, cells = 500, balanced = TRUE)
dev.off() 

png("BAL3D_2_PCA2.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_2.IT2, dims = 10:18, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_2_PCA3.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_2.IT2, dims = 19:27, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_2_PCA4.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_2.IT2, dims = 28:36, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_2_PCA5.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_2.IT2, dims = 37:45, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_2_PCA6.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_2.IT2, dims = 46:54, cells = 500, balanced = TRUE)
dev.off()

png(filename = 'BAL3D_2.ElbowPlot.png', res=300, unit="in", height=8, width=11)
ElbowPlot(BAL3D_2.IT2, ndims = 100)
dev.off()

## We will use PCs 1:36 for our embedding.

dims.use <- c(1:36)



BAL3D_2.IT2 <- FindNeighbors(BAL3D_2.IT2,
                             reduction="pca",
                             k.param = 20,
                             dims=dims.use
                             
)

BAL3D_2.IT2 <- FindClusters(BAL3D_2.IT2,
                            resolution= 2.5
                            
)

BAL3D_2.IT2<- RunUMAP(BAL3D_2.IT2, 
                      dims=dims.use,
                      min.dist = .2
)
# We create a DimPlot to show our resulting UMAP.
DimPlot(BAL3D_2.IT2, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_2/Plots/UMAP/Iteration_2")

png("BAL3D_2_UMAP_1.png", res=300, unit="in", height=8, width=11)
DimPlot(BAL3D_2.IT2, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE, repel = TRUE) + NoLegend()
dev.off()

png("BAL3D_2_nCount_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BAL3D_2.IT2, features = "nCount_RNA", label = TRUE, repel = TRUE)
dev.off()

png("BAL3D_2_nFeature_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BAL3D_2.IT2, features = "nFeature_RNA", label = TRUE, repel = TRUE)
dev.off()

png("BAL3D_2_percent.mt.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BAL3D_2.IT2, features = "percent.mt", label = TRUE, repel = TRUE)
dev.off()

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_2/Plots/UMAP/Iteration_2")

png("BAL3D_2_cluster_VlnPlot.nCount_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_2.IT2, features = "nCount_RNA", split.by = "seurat_clusters")
dev.off()

png("BAL3D_2_cluster_VlnPlot.nFeature_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_2.IT2, features = "nFeature_RNA", split.by = "seurat_clusters")
dev.off()

png("BAL3D_2_cluster_VlnPlot.percent.mt.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_2.IT2, features = "percent.mt", split.by = "seurat_clusters")
dev.off()


BAL3D_2.IT2.Markers <- FindAllMarkers(BAL3D_2.IT2, logfc.threshold = 0.25, only.pos = T)
BAL3D_2.IT2.Markers$ratio <- BAL3D_2.IT2.Markers$pct.1/BAL3D_2.IT2.Markers$pct.2
BAL3D_2.IT2.Markers$power <- BAL3D_2.IT2.Markers$ratio*BAL3D_2.IT2.Markers$avg_log2FC

BAL3D_2.IT2 <- FindClusters(BAL3D_2.IT2,
                            resolution= 1
                            
)

BAL3D_2.IT2<- RunUMAP(BAL3D_2.IT2, 
                      dims=dims.use,
                      min.dist = .2
)

DimPlot(BAL3D_2.IT2, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE)
save(BAL3D_2.IT2, file = "BAL3D_2.IT2.2025.2.26.Robj")
