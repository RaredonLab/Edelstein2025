#####BAL_3D_1 Cleaning


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_1")

library(dplyr)
library(Seurat)
library(patchwork)

# Create a starting Seurat object, loading and using the provided raw data.  
BAL3D_1.data <- Read10X(data.dir = "/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_1/Raw_Data")

BAL3D_1 <- CreateSeuratObject(counts = BAL3D_1.data)

BAL3D_1[["percent.mt"]] <- PercentageFeatureSet(BAL3D_1, pattern = "^Mt-")

# An object of class Seurat 
# 32883 features across 892208 samples within 1 assay 
# Active assay: RNA (32883 features, 0 variable features)
# 1 layer present: counts

# We will use this to determine where to set our thresholds for nCount_RNA and percent.mt
setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_1/Plots")
png("BAL3D_1_Pre_Thresholding.VlnPlot.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_1, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), raster = FALSE)
dev.off()

png("BAL3D_1_Pre_Thresholding.VlnPlot_Log.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_1, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), log = TRUE, raster = FALSE)
dev.off()


BAL3D_1<-subset(BAL3D_1, subset = nCount_RNA > 500 & percent.mt < 25)

# An object of class Seurat 
# 32883 features across 2886 samples within 1 assay 
# Active assay: RNA (32883 features, 0 variable features)
# 1 layer present: counts

# Again, create a separate set of Vln Plot images for our newly thresholded object
png("BAL3D_1_VlnPlotPostSubset.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_1, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
dev.off()

png("BAL3D_1_VlnPlot_LogPostSubset.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_1, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), log = TRUE)
dev.off()


BAL3D_1<- NormalizeData(BAL3D_1, verbose=FALSE)

BAL3D_1<- FindVariableFeatures(BAL3D_1, verbose=FALSE, nfeatures=2000)
BAL3D_1<- ScaleData(BAL3D_1, features = VariableFeatures(BAL3D_1))
BAL3D_1<- RunPCA(BAL3D_1, features = VariableFeatures(BAL3D_1), npcs =100)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_1/PCA/Iteration_1")
png("BAL3D_1_PCA1.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_1, dims = 1:9, cells = 500, balanced = TRUE)
dev.off() 

png("BAL3D_1_PCA2.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_1, dims = 10:18, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_1_PCA3.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_1, dims = 19:27, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_1_PCA4.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_1, dims = 28:36, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_1_PCA5.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_1, dims = 37:45, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_1_PCA6.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_1, dims = 46:54, cells = 500, balanced = TRUE)
dev.off()

png(filename = 'BAL3D_1.ElbowPlot.png', res=300, unit="in", height=8, width=11)
ElbowPlot(BAL3D_1, ndims = 100)
dev.off()

## We will use PCs 1:22 through 24:27 for our first embedding.

dims.use <- c(1:22, 24:27)



BAL3D_1 <- FindNeighbors(BAL3D_1,
                                reduction="pca",
                                k.param = 20,
                                dims=dims.use
                                
)

BAL3D_1 <- FindClusters(BAL3D_1,
                               resolution= 2.5
                               
)

BAL3D_1<- RunUMAP(BAL3D_1, 
                         dims=dims.use,
                         min.dist = .2
)
# We create a DimPlot to show our resulting UMAP.
DimPlot(BAL3D_1, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_1/Plots/UMAP/Iteration_1")

png("BAL3D_1_UMAP_1.png", res=300, unit="in", height=8, width=11)
DimPlot(BAL3D_1, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE, repel = TRUE) + NoLegend()
dev.off()

png("BAL3D_1_nCount_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BAL3D_1, features = "nCount_RNA", label = TRUE, repel = TRUE)
dev.off()

png("BAL3D_1_nFeature_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BAL3D_1, features = "nFeature_RNA", label = TRUE, repel = TRUE)
dev.off()

png("BAL3D_1_percent.mt.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BAL3D_1, features = "percent.mt", label = TRUE, repel = TRUE)
dev.off()

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_1/Plots/UMAP/Iteration_1")

png("BAL3D_1_cluster_VlnPlot.nCount_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_1, features = "nCount_RNA", split.by = "seurat_clusters")
dev.off()

png("BAL3D_1_cluster_VlnPlot.nFeature_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_1, features = "nFeature_RNA", split.by = "seurat_clusters")
dev.off()

png("BAL3D_1_cluster_VlnPlot.percent.mt.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_1, features = "percent.mt", split.by = "seurat_clusters")
dev.off()


BAL3D_1.Markers <- FindAllMarkers(BAL3D_1, logfc.threshold = 0.25, only.pos = T)
BAL3D_1.Markers$ratio <- BAL3D_1.Markers$pct.1/BAL3D_1.Markers$pct.2
BAL3D_1.Markers$power <- BAL3D_1.Markers$ratio*BAL3D_1.Markers$avg_log2FC


BAL3D_1.IT2<- subset(BAL3D_1, idents = c(17), invert = T)

# An object of class Seurat 
# 32883 features across 2818 samples within 1 assay 
# Active assay: RNA (32883 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap


BAL3D_1.IT2<- FindVariableFeatures(BAL3D_1.IT2, verbose=FALSE, nfeatures=2000)
BAL3D_1.IT2<- ScaleData(BAL3D_1.IT2, features = VariableFeatures(BAL3D_1.IT2))
BAL3D_1.IT2<- RunPCA(BAL3D_1.IT2, features = VariableFeatures(BAL3D_1.IT2), npcs =100)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_1/PCA/Iteration_2")
png("BAL3D_1_PCA1.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_1.IT2, dims = 1:9, cells = 500, balanced = TRUE)
dev.off() 

png("BAL3D_1_PCA2.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_1.IT2, dims = 10:18, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_1_PCA3.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_1.IT2, dims = 19:27, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_1_PCA4.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_1.IT2, dims = 28:36, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_1_PCA5.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_1.IT2, dims = 37:45, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_1_PCA6.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_1.IT2, dims = 46:54, cells = 500, balanced = TRUE)
dev.off()

png(filename = 'BAL3D_1.ElbowPlot.png', res=300, unit="in", height=8, width=11)
ElbowPlot(BAL3D_1.IT2, ndims = 100)
dev.off()

## We will use PCs 1:22 through 24:27 for our first embedding.

dims.use <- c(1:45)



BAL3D_1.IT2 <- FindNeighbors(BAL3D_1.IT2,
                         reduction="pca",
                         k.param = 20,
                         dims=dims.use
                         
)

BAL3D_1.IT2 <- FindClusters(BAL3D_1.IT2,
                        resolution= 2.5
                        
)

BAL3D_1.IT2<- RunUMAP(BAL3D_1.IT2, 
                  dims=dims.use,
                  min.dist = .2
)
# We create a DimPlot to show our resulting UMAP.
DimPlot(BAL3D_1.IT2, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_1/Plots/UMAP/Iteration_2")

png("BAL3D_1_UMAP_1.png", res=300, unit="in", height=8, width=11)
DimPlot(BAL3D_1.IT2, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE, repel = TRUE) + NoLegend()
dev.off()

png("BAL3D_1_nCount_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BAL3D_1.IT2, features = "nCount_RNA", label = TRUE, repel = TRUE)
dev.off()

png("BAL3D_1_nFeature_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BAL3D_1.IT2, features = "nFeature_RNA", label = TRUE, repel = TRUE)
dev.off()

png("BAL3D_1_percent.mt.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BAL3D_1.IT2, features = "percent.mt", label = TRUE, repel = TRUE)
dev.off()

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_1/Plots/UMAP/Iteration_2")

png("BAL3D_1_cluster_VlnPlot.nCount_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_1.IT2, features = "nCount_RNA", split.by = "seurat_clusters")
dev.off()

png("BAL3D_1_cluster_VlnPlot.nFeature_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_1.IT2, features = "nFeature_RNA", split.by = "seurat_clusters")
dev.off()

png("BAL3D_1_cluster_VlnPlot.percent.mt.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_1.IT2, features = "percent.mt", split.by = "seurat_clusters")
dev.off()


BAL3D_1.IT2.Markers <- FindAllMarkers(BAL3D_1.IT2, logfc.threshold = 0.25, only.pos = T)
BAL3D_1.IT2.Markers$ratio <- BAL3D_1.IT2.Markers$pct.1/BAL3D_1.IT2.Markers$pct.2
BAL3D_1.IT2.Markers$power <- BAL3D_1.IT2.Markers$ratio*BAL3D_1.IT2.Markers$avg_log2FC


BAL3D_1.IT3<- subset(BAL3D_1.IT2, idents = c(3), invert = T)

# An object of class Seurat 
# 32883 features across 2621 samples within 1 assay 
# Active assay: RNA (32883 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap

BAL3D_1.IT3<- FindVariableFeatures(BAL3D_1.IT3, verbose=FALSE, nfeatures=2000)
BAL3D_1.IT3<- ScaleData(BAL3D_1.IT3, features = VariableFeatures(BAL3D_1.IT3))
BAL3D_1.IT3<- RunPCA(BAL3D_1.IT3, features = VariableFeatures(BAL3D_1.IT3), npcs =100)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_1/PCA/Iteration_3")
png("BAL3D_1_PCA1.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_1.IT3, dims = 1:9, cells = 500, balanced = TRUE)
dev.off() 

png("BAL3D_1_PCA2.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_1.IT3, dims = 10:18, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_1_PCA3.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_1.IT3, dims = 19:27, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_1_PCA4.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_1.IT3, dims = 28:36, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_1_PCA5.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_1.IT3, dims = 37:45, cells = 500, balanced = TRUE)
dev.off()

png("BAL3D_1_PCA6.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BAL3D_1.IT3, dims = 46:54, cells = 500, balanced = TRUE)
dev.off()

png(filename = 'BAL3D_1.ElbowPlot.png', res=300, unit="in", height=8, width=11)
ElbowPlot(BAL3D_1.IT3, ndims = 100)
dev.off()

## We will use PCs 1:3 for our embedding.

dims.use <- c(1:36)



BAL3D_1.IT3 <- FindNeighbors(BAL3D_1.IT3,
                         reduction="pca",
                         k.param = 20,
                         dims=dims.use
                         
)

BAL3D_1.IT3 <- FindClusters(BAL3D_1.IT3,
                        resolution= 2.5
                        
)

BAL3D_1.IT3<- RunUMAP(BAL3D_1.IT3, 
                  dims=dims.use,
                  min.dist = .2
)
# We create a DimPlot to show our resulting UMAP.
DimPlot(BAL3D_1.IT3, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_1/Plots/UMAP/Iteration_3")

png("BAL3D_1_UMAP_1.png", res=300, unit="in", height=8, width=11)
DimPlot(BAL3D_1.IT3, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE, repel = TRUE) + NoLegend()
dev.off()

png("BAL3D_1_nCount_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BAL3D_1.IT3, features = "nCount_RNA", label = TRUE, repel = TRUE)
dev.off()

png("BAL3D_1_nFeature_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BAL3D_1.IT3, features = "nFeature_RNA", label = TRUE, repel = TRUE)
dev.off()

png("BAL3D_1_percent.mt.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BAL3D_1.IT3, features = "percent.mt", label = TRUE, repel = TRUE)
dev.off()

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/BAL_3D_1/Plots/UMAP/Iteration_3")

png("BAL3D_1_cluster_VlnPlot.nCount_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_1.IT3, features = "nCount_RNA", split.by = "seurat_clusters")
dev.off()

png("BAL3D_1_cluster_VlnPlot.nFeature_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_1.IT3, features = "nFeature_RNA", split.by = "seurat_clusters")
dev.off()

png("BAL3D_1_cluster_VlnPlot.percent.mt.png", res=300, unit="in", height=8, width=11)
VlnPlot(BAL3D_1.IT3, features = "percent.mt", split.by = "seurat_clusters")
dev.off()


BAL3D_1.IT3.Markers <- FindAllMarkers(BAL3D_1.IT3, logfc.threshold = 0.25, only.pos = T)
BAL3D_1.IT3.Markers$ratio <- BAL3D_1.IT3.Markers$pct.1/BAL3D_1.IT3.Markers$pct.2
BAL3D_1.IT3.Markers$power <- BAL3D_1.IT3.Markers$ratio*BAL3D_1.IT3.Markers$avg_log2FC
save(BAL3D_1.IT3, file = "BAL3D_1.IT3.2025.2.26.Robj")


BAL3D_1.IT3 <- FindClusters(BAL3D_1.IT3,
                            resolution= 1
                            
)

BAL3D_1.IT3<- RunUMAP(BAL3D_1.IT3, 
                      dims=dims.use,
                      min.dist = .2
)
DimPlot(BAL3D_1.IT3)
FeaturePlot(BAL3D_1.IT3, features = "Il33")

BAL3D_1.IT3.Markers <- FindAllMarkers(BAL3D_1.IT3, logfc.threshold = 0.25, only.pos = T)
BAL3D_1.IT3.Markers$ratio <- BAL3D_1.IT3.Markers$pct.1/BAL3D_1.IT3.Markers$pct.2
BAL3D_1.IT3.Markers$power <- BAL3D_1.IT3.Markers$ratio*BAL3D_1.IT3.Markers$avg_log2FC