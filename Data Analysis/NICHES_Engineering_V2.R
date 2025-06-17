options(future.globals.maxSize = 16000 * 1024^2)
# Set WD

library(Seurat)
library(SeuratWrappers)
library(NICHES)
library(ggplot2)

setwd("/Volumes/T7/Sophie/NICHES_Engineering")

cellclass.col.pal <- c("#01ba8b","#d9afed","#f7a4a3","#fff720")

condition.col.pal <- c("#f5d7d7","#E9c500","#a0c63d","#b2f1e9","#854AF2")

celltype.col.pal <- c("#ce378f","#b0e3f7","#cac7ff","#f7a4a3","#834AD9",
                      "#d3f5b0","#89f98f","#097539","#ffda9e","#f0a442","#ff141a",
                      "#fff410","#17c1dd","#1179fa","#ae8674") 

full.col.pal <- c("#f5d7d7","#fa9730","#00c3ff","#f9b3ef","#cac7ff","#d3f5b0","#2bff88","#1dfffb",
                  "#1167dd","#01ba8b","#f7a4a3","#d486b2","#f9f797","#d0a10a","#d0ff00","#cd378f",
                  "#f06361","#8734F0","#3cd500","#ffda9e","#c9522a","#b0e3f7","#f7d042","#a0c63d",
                  "#e331f7","#ff045c","#c0c4bc","#fff720","#605acf","#309500","#ad8674","#d9afed",
                  "#097538","#151dc2","#a0009b","#b22b29","#cac426")

load('eng.subset.integrated_HK.Robj')

table(eng.subset.integrated_HK$Orig_ID)
table(eng.subset.integrated_HK$Condition)
table(eng.subset.integrated_HK$CellType.combined.Integrated)

BAL_3D.integrated_HK <- subset(eng.subset.integrated_HK, subset = Condition == 'BAL_3D')
Mixed_3D.integrated_HK <- subset(eng.subset.integrated_HK, subset = Condition == 'Mixed_3D')
PD_3D.integrated_HK <- subset(eng.subset.integrated_HK, subset = Condition == 'PD_3D')

### NICHES: based on conditions, no imputation
# Condition: BAL_3D
BAL_3D.list <- SplitObject(BAL_3D.integrated_HK,split.by='Orig_ID')

NICHES.list <- list()

for(i in 1:length(BAL_3D.list)){
  print(i)
  CellType_Dist <- table(BAL_3D.list[[i]]$CellType.combined.Integrated)
  to.delete <- names(CellType_Dist[CellType_Dist<=1])
  BAL_3D.list[[i]] <- subset(BAL_3D.list[[i]],subset = CellType.combined.Integrated %in% to.delete,invert=T)
  
  NICHES.list[[i]] <- RunNICHES(BAL_3D.list[[i]],
                                LR.database = "fantom5",
                                species = "rat",
                                assay = "RNA",
                                cell_types = "CellType.combined.Integrated",
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

# Condition: PD_3D
PD_3D.list <- SplitObject(PD_3D.integrated_HK,split.by='Orig_ID')

NICHES.list <- list()

for(i in 1:length(PD_3D.list)){
  print(i)
  CellType_Dist <- table(PD_3D.list[[i]]$CellType.combined.Integrated)
  to.delete <- names(CellType_Dist[CellType_Dist<=1])
  PD_3D.list[[i]] <- subset(PD_3D.list[[i]],subset = CellType.combined.Integrated %in% to.delete,invert=T)
  
  NICHES.list[[i]] <- RunNICHES(PD_3D.list[[i]],
                                LR.database = "fantom5",
                                species = "rat",
                                assay = "RNA",
                                cell_types = "CellType.combined.Integrated",
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

# Condition: Mixed_3D
Mixed_3D.list <- SplitObject(Mixed_3D.integrated_HK,split.by='Orig_ID')

NICHES.list <- list()

for(i in 1:length(Mixed_3D.list)){
  print(i)
  CellType_Dist <- table(Mixed_3D.list[[i]]$CellType.combined.Integrated)
  to.delete <- names(CellType_Dist[CellType_Dist<=1])
  Mixed_3D.list[[i]] <- subset(Mixed_3D.list[[i]],subset = CellType.combined.Integrated %in% to.delete,invert=T)
  
  NICHES.list[[i]] <- RunNICHES(Mixed_3D.list[[i]],
                                LR.database = "fantom5",
                                species = "rat",
                                assay = "RNA",
                                cell_types = "CellType.combined.Integrated",
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

eng.CTC.bySample <- merge(BAL_3D.CTC.bySample, c(PD_3D.CTC.bySample,Mixed_3D.CTC.bySample))
eng.CTC.bySample <- JoinLayers(eng.CTC.bySample)

# Filter out low-information signals
VlnPlot(eng.CTC.bySample,c('nFeature_CellToCell'),group.by = 'Orig_ID.Sending',raster=F,pt.size = 0.1)+ggtitle('Before Filtration')+ylab('nFeature_CellToCell')

eng.CTC.sub.bySample <- subset(eng.CTC.bySample,nFeature_CellToCell>60)

VlnPlot(eng.CTC.sub.bySample,c('nFeature_CellToCell'),group.by = 'Orig_ID.Sending',raster=F,pt.size = 0.1)+ggtitle('After Filtration')+ylab('nFeature_CellToCell')

# Embedding and clustering
tmp <- ScaleData(eng.CTC.sub.bySample)
tmp <- FindVariableFeatures(tmp)
tmp <- RunPCA(tmp,npcs = 100)
pdf(file='Eng.CTC.bySample.PCs.pdf',width=10,height=8)
ElbowPlot(tmp,ndims = 100)
PCHeatmap(tmp,cells=200,balanced=T,dims=1:9)
PCHeatmap(tmp,cells=200,balanced=T,dims=10:18)
PCHeatmap(tmp,cells=200,balanced=T,dims=19:27)
PCHeatmap(tmp,cells=200,balanced=T,dims=28:36)
PCHeatmap(tmp,cells=200,balanced=T,dims=37:45)
PCHeatmap(tmp,cells=200,balanced=T,dims=46:54)
PCHeatmap(tmp,cells=200,balanced=T,dims=55:63)
PCHeatmap(tmp,cells=200,balanced=T,dims=64:72)
PCHeatmap(tmp,cells=200,balanced=T,dims=73:81)
PCHeatmap(tmp,cells=200,balanced=T,dims=82:90)
PCHeatmap(tmp,cells=200,balanced=T,dims=91:99)
dev.off()

tmp <- RunUMAP(tmp,dims = c(1:18))
tmp <- FindNeighbors(tmp,dims = c(1:18))
tmp <- FindClusters(tmp,resolution=0.4)

# Plot
DimPlot(tmp,label=T,raster = F,cols=full.col.pal)
DimPlot(tmp,label=T,raster = F, cols=full.col.pal,split.by = 'Condition.Sending',ncol=3)
DimPlot(tmp,group.by = 'CellType.Sending',raster=F,shuffle = T,cols=celltype.col.pal)
DimPlot(tmp,group.by = 'CellType.Receiving',raster=F,shuffle = T,cols=celltype.col.pal)

# Save regular clusters, in column "clusters"
tmp$clusters <- tmp$seurat_clusters

# Try integration by condition
# Split first
tmp[["CellToCell"]] <- split(tmp[["CellToCell"]], f = tmp$Condition.Sending)

# Embed data, via individual 'layers'
tmp <- NormalizeData(tmp)
tmp <- FindVariableFeatures(tmp)
tmp <- ScaleData(tmp)
tmp <- RunPCA(tmp)

tmp.integrated <- IntegrateLayers(object = tmp, method=RPCAIntegration,
                                  orig.reduction="pca",new.reduction="integrated.rpca",verbose=F)
tmp.integrated[["CellToCell"]] <- JoinLayers(tmp.integrated[["CellToCell"]])

# Choose number of PCs
ElbowPlot(tmp.integrated,ndims=50)+ggtitle('Integration by Condition ElbowPlot')
PCHeatmap(tmp.integrated,cells=200,balanced=T,dims=1:9)
PCHeatmap(tmp.integrated,cells=200,balanced=T,dims=10:18)
PCHeatmap(tmp.integrated,cells=200,balanced=T,dims=19:27)
PCHeatmap(tmp.integrated,cells=200,balanced=T,dims=28:36)
PCHeatmap(tmp.integrated,cells=200,balanced=T,dims=37:45) # number of PCs = 1:16, 18:19 is appropriate!

# Check final results
tmp.integrated <- RunUMAP(tmp.integrated,reduction='integrated.rpca',dims=c(1:16,18:19),reduction.name='umap.rpca')
tmp.integrated <- FindNeighbors(tmp.integrated,reduction='integrated.rpca',dims=c(1:16,18:19))
tmp.integrated <- FindClusters(tmp.integrated,resolution=0.3,cluster.name='rpca_clusters')

DimPlot(tmp.integrated,label=T,raster = F,cols=full.col.pal,reduction = 'umap.rpca')
DimPlot(tmp.integrated,label=T,raster = F, cols=full.col.pal,split.by = 'Condition.Sending',ncol=3,reduction='umap.rpca')
DimPlot(tmp.integrated,group.by = 'CellType.Sending',raster=F,shuffle = T,cols=celltype.col.pal,reduction='umap.rpca')
DimPlot(tmp.integrated,group.by = 'CellType.Receiving',raster=F,shuffle = T,cols=celltype.col.pal,reduction='umap.rpca')

# Save integrated rpca clusters, in column "integrated.clusters"
tmp.integrated$integrated.clusters <- tmp.integrated$rpca_clusters
# Rename
eng.CTC.bySample.byCondition <- tmp.integrated


## Add integration by Sample information
tmp <- eng.CTC.bySample.byCondition

## Try integration by condition
# Split first
tmp[["CellToCell"]] <- split(tmp[["CellToCell"]], f = tmp$Orig_ID.Sending)

# Embed data, via individual 'layers'
tmp <- NormalizeData(tmp)
tmp <- FindVariableFeatures(tmp)
tmp <- ScaleData(tmp)
tmp <- RunPCA(tmp)

tmp.integrated <- IntegrateLayers(object = tmp, method=RPCAIntegration,
                                  orig.reduction="pca",new.reduction="integrated.rpca",verbose=F)
tmp.integrated[["CellToCell"]] <- JoinLayers(tmp.integrated[["CellToCell"]])

# Choose number of PCs
ElbowPlot(tmp.integrated,ndims=50)+ggtitle('Integration by Sample ElbowPlot')
PCHeatmap(tmp.integrated,cells=200,balanced=T,dims=1:9)
PCHeatmap(tmp.integrated,cells=200,balanced=T,dims=10:18)
PCHeatmap(tmp.integrated,cells=200,balanced=T,dims=19:27)
PCHeatmap(tmp.integrated,cells=200,balanced=T,dims=28:36)
PCHeatmap(tmp.integrated,cells=200,balanced=T,dims=37:45) # number of PCs = 1:14, 16:18 is appropriate!

# Check final results
tmp.integrated <- RunUMAP(tmp.integrated,reduction='integrated.rpca',dims=c(1:14,16:18),reduction.name='umap.rpca.2')
tmp.integrated <- FindNeighbors(tmp.integrated,reduction='integrated.rpca',dims=c(1:14,16:18))
tmp.integrated <- FindClusters(tmp.integrated,resolution=0.6,cluster.name='rpca_clusters')
# Plots
a = DimPlot(eng.CTC.final,label=T,raster = F,cols=full.col.pal,reduction = 'umap.rpca.2',shuffle=T)
DimPlot(eng.CTC.final,label=T,raster = F, cols=full.col.pal,split.by = 'Condition.Sending',ncol=3,reduction='umap.rpca.2')
b = DimPlot(eng.CTC.final,group.by = 'CellType.Sending',raster=F,shuffle = T,cols=celltype.col.pal,reduction='umap.rpca.2')
c = DimPlot(eng.CTC.final,group.by = 'CellType.Receiving',raster=F,shuffle = T,cols=celltype.col.pal,reduction='umap.rpca.2')
d = DimPlot(eng.CTC.final,group.by = 'CellType.combined.Integrated.Receiving',raster=F,shuffle = T,cols=celltype.col.pal,reduction='umap.rpca.2')
e = DimPlot(eng.CTC.final,group.by = 'CellType.combined.Integrated.Receiving',raster=F,shuffle = T,cols=celltype.col.pal,reduction='umap.rpca.2')
f = DimPlot(eng.CTC.final,group.by = 'CellClass.combined.Integrated.Receiving',raster=F,shuffle = T,cols=celltype.col.pal,reduction='umap.rpca.2')
g = DimPlot(eng.CTC.final,group.by = 'CellClass.combined.Integrated.Receiving',raster=F,shuffle = T,cols=celltype.col.pal,reduction='umap.rpca.2')
a | b | c
a | f | g

tmp.integrated$integrated.clusters.2 <- tmp.integrated$rpca_clusters
eng.CTC.final <- tmp.integrated

save(eng.CTC.final, file = 'eng.CTC.final.Robj')

# Plot example
DimPlot(eng.CTC.final,group.by = 'clusters',shuffle=T, raster = F, cols = full.col.pal,label = T, reduction = 'umap') # CellToCell regular embedding and clustering
DimPlot(eng.CTC.final,group.by = 'integrated.clusters',shuffle=T, raster = F, cols = full.col.pal,label = T, reduction = 'umap.rpca') # Embedding integrated by condition
DimPlot(eng.CTC.final,group.by = 'integrated.clusters.2',shuffle=T, raster = F, cols = full.col.pal,label = T, reduction = 'umap.rpca.2') # Embedding integrated by sample

p1 = DimPlot(eng.CTC.final,group.by = 'clusters',shuffle=T, raster = F, cols = full.col.pal,label = F, reduction = 'umap') # CellToCell regular embedding and clustering
p2 = DimPlot(eng.CTC.final,group.by = 'CellType.Sending',shuffle=T, raster = F, cols = full.col.pal,label = F, reduction = 'umap') # CellToCell regular embedding and clustering
p3 = DimPlot(eng.CTC.final,group.by = 'CellType.Receiving',shuffle=T, raster = F, cols = full.col.pal,label = F, reduction = 'umap') # CellToCell regular embedding and clustering
p1 | p2 | p3

table(eng.CTC.final$CellType.Sending)
table(eng.CTC.final$CellType.Receiving)

p1 = DimPlot(eng.CTC.final,group.by = 'clusters',shuffle=T, raster = F, cols = full.col.pal,label = F, reduction = 'umap') # CellToCell regular embedding and clustering
p2 = DimPlot(eng.CTC.final,group.by = 'CellType.combined.Integrated.Sending',shuffle=T, raster = F, cols = full.col.pal,label = F, reduction = 'umap') # CellToCell regular embedding and clustering
p3 = DimPlot(eng.CTC.final,group.by = 'CellType.combined.Integrated.Receiving',shuffle=T, raster = F, cols = full.col.pal,label = F, reduction = 'umap') # CellToCell regular embedding and clustering
p1 | p2 | p3
