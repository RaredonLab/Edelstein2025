#####PD_3D_3 Cleaning


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/PD_3D_3")

library(dplyr)
library(Seurat)
library(patchwork)

# Create a starting Seurat object, loading and using the provided raw data.  
PD.3D_3.data <- Read10X(data.dir = "/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/PD_3D_3/Raw_Data")

PD.3D_3 <- CreateSeuratObject(counts = PD.3D_3.data)

PD.3D_3[["percent.mt"]] <- PercentageFeatureSet(PD.3D_3, pattern = "^Mt-")

# An object of class Seurat 
# 32883 features across 1073533 samples within 1 assay 
# Active assay: RNA (32883 features, 0 variable features)
# 1 layer present: counts

# We will use this to determine where to set our thresholds for nCount_RNA and percent.mt
setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/PD_3D_3/Plots")
png("PD.3D_3_Pre_Thresholding.VlnPlot.png", res=300, unit="in", height=8, width=11)
VlnPlot(PD.3D_3, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), raster = FALSE)
dev.off()

png("PD.3D_3_Pre_Thresholding.VlnPlot_Log.png", res=300, unit="in", height=8, width=11)
VlnPlot(PD.3D_3, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), log = TRUE, raster = FALSE)
dev.off()


PD.3D_3<-subset(PD.3D_3, subset = nCount_RNA > 500 & percent.mt < 25)

# An object of class Seurat 
# 32883 features across 5429 samples within 1 assay 
# Active assay: RNA (32883 features, 0 variable features)
# 1 layer present: counts

# Again, create a separate set of Vln Plot images for our newly thresholded object
png("PD.3D_3_VlnPlotPostSubset.png", res=300, unit="in", height=8, width=11)
VlnPlot(PD.3D_3, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
dev.off()

png("PD.3D_3_VlnPlot_LogPostSubset.png", res=300, unit="in", height=8, width=11)
VlnPlot(PD.3D_3, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), log = TRUE)
dev.off()


PD.3D_3<- NormalizeData(PD.3D_3, verbose=FALSE)

PD.3D_3<- FindVariableFeatures(PD.3D_3, verbose=FALSE, nfeatures=2000)
PD.3D_3<- ScaleData(PD.3D_3, features = VariableFeatures(PD.3D_3))
PD.3D_3<- RunPCA(PD.3D_3, features = VariableFeatures(PD.3D_3), npcs =100)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/PD_3D_3/PCA/Iteration_1")
png("PD.3D_3_PCA1.png", res=300, unit="in", height=8, width=11)
DimHeatmap(PD.3D_3, dims = 1:9, cells = 500, balanced = TRUE)
dev.off() 

png("PD.3D_3_PCA2.png", res=300, unit="in", height=8, width=11)
DimHeatmap(PD.3D_3, dims = 10:18, cells = 500, balanced = TRUE)
dev.off()

png("PD.3D_3_PCA3.png", res=300, unit="in", height=8, width=11)
DimHeatmap(PD.3D_3, dims = 19:27, cells = 500, balanced = TRUE)
dev.off()

png("PD.3D_3_PCA4.png", res=300, unit="in", height=8, width=11)
DimHeatmap(PD.3D_3, dims = 28:36, cells = 500, balanced = TRUE)
dev.off()

png("PD.3D_3_PCA5.png", res=300, unit="in", height=8, width=11)
DimHeatmap(PD.3D_3, dims = 37:45, cells = 500, balanced = TRUE)
dev.off()

png("PD.3D_3_PCA6.png", res=300, unit="in", height=8, width=11)
DimHeatmap(PD.3D_3, dims = 46:54, cells = 500, balanced = TRUE)
dev.off()

png(filename = 'PD.3D_3.ElbowPlot.png', res=300, unit="in", height=8, width=11)
ElbowPlot(PD.3D_3, ndims = 100)
dev.off()

## We will use PCs 1:45 for our first embedding.

dims.use <- c(1:45)



PD.3D_3 <- FindNeighbors(PD.3D_3,
                         reduction="pca",
                         k.param = 20,
                         dims=dims.use
                         
)

PD.3D_3 <- FindClusters(PD.3D_3,
                        resolution= 2.5
                        
)

PD.3D_3<- RunUMAP(PD.3D_3, 
                  dims=dims.use,
                  min.dist = .2
)
# We create a DimPlot to show our resulting UMAP.
DimPlot(PD.3D_3, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/PD_3D_3/Plots/UMAP/Iteration_1")

png("PD.3D_3_UMAP_1.png", res=300, unit="in", height=8, width=11)
DimPlot(PD.3D_3, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE, repel = TRUE) + NoLegend()
dev.off()

png("PD.3D_3_nCount_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(PD.3D_3, features = "nCount_RNA", label = TRUE, repel = TRUE)
dev.off()

png("PD.3D_3_nFeature_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(PD.3D_3, features = "nFeature_RNA", label = TRUE, repel = TRUE)
dev.off()

png("PD.3D_3_percent.mt.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(PD.3D_3, features = "percent.mt", label = TRUE, repel = TRUE)
dev.off()

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/PD_3D_3/Plots/UMAP/Iteration_1")

png("PD.3D_3_cluster_VlnPlot.nCount_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(PD.3D_3, features = "nCount_RNA", split.by = "seurat_clusters")
dev.off()

png("PD.3D_3_cluster_VlnPlot.nFeature_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(PD.3D_3, features = "nFeature_RNA", split.by = "seurat_clusters")
dev.off()

png("PD.3D_3_cluster_VlnPlot.percent.mt.png", res=300, unit="in", height=8, width=11)
VlnPlot(PD.3D_3, features = "percent.mt", split.by = "seurat_clusters")
dev.off()


PD.3D_3.Markers <- FindAllMarkers(PD.3D_3, logfc.threshold = 0.25, only.pos = T)
PD.3D_3.Markers$ratio <- PD.3D_3.Markers$pct.1/PD.3D_3.Markers$pct.2
PD.3D_3.Markers$power <- PD.3D_3.Markers$ratio*PD.3D_3.Markers$avg_log2FC


PD.3D_3.IT2<- subset(PD.3D_3, idents = c(4, 7, 20, 24), invert = T)

# An object of class Seurat 
# 32883 features across 4724 samples within 1 assay 
# Active assay: RNA (32883 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap

PD.3D_3.IT2<- FindVariableFeatures(PD.3D_3.IT2, verbose=FALSE, nfeatures=2000)
PD.3D_3.IT2<- ScaleData(PD.3D_3.IT2, features = VariableFeatures(PD.3D_3.IT2))
PD.3D_3.IT2<- RunPCA(PD.3D_3.IT2, features = VariableFeatures(PD.3D_3.IT2), npcs =100)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/PD_3D_3/PCA/Iteration_2")
png("PD.3D_3_PCA1.png", res=300, unit="in", height=8, width=11)
DimHeatmap(PD.3D_3.IT2, dims = 1:9, cells = 500, balanced = TRUE)
dev.off() 

png("PD.3D_3_PCA2.png", res=300, unit="in", height=8, width=11)
DimHeatmap(PD.3D_3.IT2, dims = 10:18, cells = 500, balanced = TRUE)
dev.off()

png("PD.3D_3_PCA3.png", res=300, unit="in", height=8, width=11)
DimHeatmap(PD.3D_3.IT2, dims = 19:27, cells = 500, balanced = TRUE)
dev.off()

png("PD.3D_3_PCA4.png", res=300, unit="in", height=8, width=11)
DimHeatmap(PD.3D_3.IT2, dims = 28:36, cells = 500, balanced = TRUE)
dev.off()

png("PD.3D_3_PCA5.png", res=300, unit="in", height=8, width=11)
DimHeatmap(PD.3D_3.IT2, dims = 37:45, cells = 500, balanced = TRUE)
dev.off()

png("PD.3D_3_PCA6.png", res=300, unit="in", height=8, width=11)
DimHeatmap(PD.3D_3.IT2, dims = 46:54, cells = 500, balanced = TRUE)
dev.off()

png(filename = 'PD.3D_3.ElbowPlot.png', res=300, unit="in", height=8, width=11)
ElbowPlot(PD.3D_3.IT2, ndims = 100)
dev.off()

## We will use PCs 1:40 for our  embedding.

dims.use <- c(1:40)



PD.3D_3.IT2 <- FindNeighbors(PD.3D_3.IT2,
                             reduction="pca",
                             k.param = 20,
                             dims=dims.use
                             
)

PD.3D_3.IT2 <- FindClusters(PD.3D_3.IT2,
                            resolution= 2.5
                            
)

PD.3D_3.IT2<- RunUMAP(PD.3D_3.IT2, 
                      dims=dims.use,
                      min.dist = .2
)
# We create a DimPlot to show our resulting UMAP.
DimPlot(PD.3D_3.IT2, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/PD_3D_3/Plots/UMAP/Iteration_2")

png("PD.3D_3_UMAP_1.png", res=300, unit="in", height=8, width=11)
DimPlot(PD.3D_3.IT2, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE, repel = TRUE) + NoLegend()
dev.off()

png("PD.3D_3_nCount_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(PD.3D_3.IT2, features = "nCount_RNA", label = TRUE, repel = TRUE)
dev.off()

png("PD.3D_3_nFeature_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(PD.3D_3.IT2, features = "nFeature_RNA", label = TRUE, repel = TRUE)
dev.off()

png("PD.3D_3_percent.mt.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(PD.3D_3.IT2, features = "percent.mt", label = TRUE, repel = TRUE)
dev.off()

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/Sophie_Second_Opinion/PD_3D_3/Plots/UMAP/Iteration_2")

png("PD.3D_3_cluster_VlnPlot.nCount_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(PD.3D_3.IT2, features = "nCount_RNA", split.by = "seurat_clusters")
dev.off()

png("PD.3D_3_cluster_VlnPlot.nFeature_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(PD.3D_3.IT2, features = "nFeature_RNA", split.by = "seurat_clusters")
dev.off()

png("PD.3D_3_cluster_VlnPlot.percent.mt.png", res=300, unit="in", height=8, width=11)
VlnPlot(PD.3D_3.IT2, features = "percent.mt", split.by = "seurat_clusters")
dev.off()


PD.3D_3.IT2.Markers <- FindAllMarkers(PD.3D_3.IT2, logfc.threshold = 0.25, only.pos = T)
PD.3D_3.IT2.Markers$ratio <- PD.3D_3.IT2.Markers$pct.1/PD.3D_3.IT2.Markers$pct.2
PD.3D_3.IT2.Markers$power <- PD.3D_3.IT2.Markers$ratio*PD.3D_3.IT2.Markers$avg_log2FC





PD.3D_3.IT2 <- FindClusters(PD.3D_3.IT2,
                            resolution= 1
                            
)

PD.3D_3.IT2<- RunUMAP(PD.3D_3.IT2, 
                      dims=dims.use,
                      min.dist = .2
)
DimPlot(PD.3D_3.IT2)
FeaturePlot(PD.3D_3.IT2, features = "Il33")


save(PD.3D_3.IT2, file = "PD.3D_3.IT2.2025.3.3.Robj")


