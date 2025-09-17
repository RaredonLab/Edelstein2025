

# Attempting to look at just the mesenchyme <-> immune signaling paradigms
# Pull out: Rspo3+_Mes, Pdgfrb+_Pericyte, Anti_Inflamm_Mac, Pro_Inflamm_Mac (BUT FROM THE ORIGINAL SINGLE CELL OBJECT NOT NICHES)

# Recall, we only want to compare PD_3D and Mixed_3D here
selected_types <- c("Rspo3+_Mes","Pdgfrb+_Pericyte", 
                    "Anti_Inflamm_Mac", "Pro_Inflamm_Mac")
mes.imm.obj <- subset(
  eng.subset.integrated_HK,
  subset = CellType %in% selected_types & 
    Condition %in% c("PD_3D", "Mixed_3D"))

# Now re-scale, etc
mes.imm.obj <- NormalizeData(mes.imm.obj) # Don't have to do but doing it anyway
# Identify highly variable features
mes.imm.obj <- FindVariableFeatures(mes.imm.obj)
# Scale the data
mes.imm.obj <- ScaleData(mes.imm.obj)

# Perform PCA
mes.imm.obj <- RunPCA(mes.imm.obj, features = VariableFeatures(object = mes.imm.obj))
# Visualize PCA with an Elbow Plot
ElbowPlot(mes.imm.obj, ndims = 200)
PCHeatmap(mes.imm.obj,cells=200,balanced=T,dims=1:9) # Skip 8 b/c mostly epithelial contamination; 
PCHeatmap(mes.imm.obj,cells=200,balanced=T,dims=10:18) 
PCHeatmap(mes.imm.obj,cells=200,balanced=T,dims=19:27)

# Run UMAP and clustering - SE picking specific PCs based on what we know about structure
mes.imm.obj <- RunUMAP(mes.imm.obj, dims = c(1:7,9:10))
mes.imm.obj <- FindNeighbors(mes.imm.obj, dims = c(1:7,9:10))
mes.imm.obj <- FindClusters(mes.imm.obj, resolution = 0.4)

Idents(mes.imm.obj) = mes.imm.obj$seurat_clusters
DimPlot(mes.imm.obj, label = T) + ggtitle("Clusters")
# Now look at lineage markers
FeaturePlot(mes.imm.obj, features = c("Epcam","Ptprc","Col1a1"), order = T, label = T)

# Okay, let's subset again to get rid of cluster 5 (epi contamination)
mes.imm.obj1 = subset(mes.imm.obj, idents = c("5"),invert=TRUE)
# Rescale the data
mes.imm.obj1 = ScaleData(mes.imm.obj1)
# Normalizing the data
mes.imm.obj1 = NormalizeData(mes.imm.obj1)
# Finding variable features
mes.imm.obj1 = FindVariableFeatures(mes.imm.obj1)

# Run PCA
mes.imm.obj1 = RunPCA(mes.imm.obj1, npcs = 100)
# Visualizing the PCAs
setwd("~/Desktop/iScience Transfer Submission/Additional Analyses/Mesenchyme_Immune_CircuitAnalysis")
pdf(file="mes.imm.obj1_PCs.pdf",width=10,height=8)
ElbowPlot(mes.imm.obj1,ndims = 50)
PCHeatmap(mes.imm.obj1,cells=200,balanced=T,dims=1:9) # Skip 3
PCHeatmap(mes.imm.obj1,cells=200,balanced=T,dims=10:18)
PCHeatmap(mes.imm.obj1,cells=200,balanced=T,dims=19:27)
PCHeatmap(mes.imm.obj1,cells=200,balanced=T,dims=28:36)
dev.off()
# GO with PCs 1,2,4,5,6,7,9,10

mes.imm.obj1 = RunUMAP(mes.imm.obj1, dims = c(1,2,4:7,9,10))
DimPlot(mes.imm.obj1, label = TRUE)
# Creating nearest neighbor graph, not clustering.
mes.imm.obj1 = FindNeighbors(mes.imm.obj1, dims =c(1,2,4:7,9,10))
# Defining clusters - very high res. There is purpose in doing this.
mes.imm.obj1 = FindClusters(mes.imm.obj1, res = 0.2)
DimPlot(mes.imm.obj1, label = TRUE)
DimPlot(mes.imm.obj1, label = TRUE, group.by = "CellType")
# Look at how old cell type annotations play out
DimPlot(mes.imm.obj1, group.by = "CellType", label = T) # These are old annotations
# Look at condition level distribution (PD_3D and Mixed_3D)
DimPlot(mes.imm.obj1, group.by = "Condition", label = F, cols = Condition.cols) # Condition; cool, cool, this checks out


Idents(mes.imm.obj1) = mes.imm.obj1$seurat_clusters
FeaturePlot(mes.imm.obj1, features = c("Ptprc","Col1a1"), order = T, label = T)
cell.immu = WhichCells(mes.imm.obj1, idents = c(0,1,2))
cell.mes = WhichCells(mes.imm.obj1, idents = c(3,4))
mes.imm.obj1 = SetIdent(mes.imm.obj1, cells = cell.mes, value = 'Mesenchyme')
mes.imm.obj1 = SetIdent(mes.imm.obj1, cells = cell.immu, value = 'Immune')
# Stash the cell class
mes.imm.obj1$CellClass.sub = Idents(mes.imm.obj1)
table(mes.imm.obj1$CellClass.sub)
# confirming that everything is labeled
sum(is.na(mes.imm.obj1$CellClass.sub))
DimPlot(mes.imm.obj1, label = T)

Idents(mes.imm.obj) = mes.imm.obj$seurat_clusters
# Cell Type annos - just adding another meta-data layer
mes.imm.obj1$CellType.sub <- NA
c0 <- WhichCells(mes.imm.obj1, idents = "0") # Polarized_Mac
c1 <- WhichCells(mes.imm.obj1, idents = "1") # Polarized_Mac
c2 <- WhichCells(mes.imm.obj1, idents = "2") # Polarized_Mac
c3 <- WhichCells(mes.imm.obj1, idents = "3") # Polarized_Mac
c4 <- WhichCells(mes.imm.obj1, idents = "4") # Polarized_Mac
mes.imm.obj1$CellType.sub[c0] <- 'Polarized_Mac'
mes.imm.obj1$CellType.sub[c1] <- 'Polarized_Mac'
mes.imm.obj1$CellType.sub[c2] <- 'Polarized_Mac'
mes.imm.obj1$CellType.sub[c3] <- 'Rspo3+_Mes'
mes.imm.obj1$CellType.sub[c4] <- 'Pdgfrb+_Pericyte'

# Check that all applied properly
table(mes.imm.obj1$CellType.sub)
DimPlot(mes.imm.obj1, group.by = 'CellType.sub', label = T)

# Save this subset transcriptomic object
setwd("~/Desktop/iScience Transfer Submission/Additional Analyses/Mesenchyme_Immune_CircuitAnalysis")
save(mes.imm.obj1, file = "mes.imm.obj.Robj")

####### Now, let's run NICHES ##########
# Check data structure
table(mes.imm.obj1$Orig_ID)
table(mes.imm.obj1$Condition)
table(mes.imm.obj1$CellType)
table(mes.imm.obj1$CellClass.sub)
table(mes.imm.obj1$CellType.sub)

# Split by conditions first
Mixed_3D.mes.imm <- subset(mes.imm.obj1, subset = Condition == 'Mixed_3D')
PD_3D.mes.imm <- subset(mes.imm.obj1, subset = Condition == 'PD_3D')

### NICHES: based on conditions, no imputation

##### Condition: PD_3D
PD_3D.list <- SplitObject(Mixed_3D.mes.imm,split.by='Orig_ID')
NICHES.list <- list()

for(i in 1:length(PD_3D.list)){
  print(i)
  CellType_Dist <- table(PD_3D.list[[i]]$CellType.sub)
  to.delete <- names(CellType_Dist[CellType_Dist<=1])
  PD_3D.list[[i]] <- subset(PD_3D.list[[i]],subset = CellType.sub %in% to.delete,invert=T)
  
  NICHES.list[[i]] <- RunNICHES(PD_3D.list[[i]],
                                LR.database = "fantom5",
                                species = "rat",
                                assay = "RNA",
                                cell_types = "CellType.sub",
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
Mixed_3D.list <- SplitObject(PD_3D.mes.imm,split.by='Orig_ID')

NICHES.list <- list()

for(i in 1:length(Mixed_3D.list)){
  print(i)
  CellType_Dist <- table(Mixed_3D.list[[i]]$CellType.sub)
  to.delete <- names(CellType_Dist[CellType_Dist<=1])
  Mixed_3D.list[[i]] <- subset(Mixed_3D.list[[i]],subset = CellType.sub %in% to.delete,invert=T)
  
  NICHES.list[[i]] <- RunNICHES(Mixed_3D.list[[i]],
                                LR.database = "fantom5",
                                species = "rat",
                                assay = "RNA",
                                cell_types = "CellType.sub",
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
eng.CTC.bySample_mes.imm <- merge(PD_3D.CTC.bySample, c(Mixed_3D.CTC.bySample))
eng.CTC.bySample_mes.imm <- JoinLayers(eng.CTC.bySample_mes.imm)

### Filter out low-information signals
VlnPlot(eng.CTC.bySample_mes.imm,c('nFeature_CellToCell'),group.by = 'Orig_ID.Sending',raster=F,pt.size = 0.1) + ggtitle('Before Filtration')+ylab('nFeature_CellToCell')
eng.CTC.bySample_mes.imm <- subset(eng.CTC.bySample_mes.imm,nFeature_CellToCell>20)
# Check again after filtering
VlnPlot(eng.CTC.bySample_mes.imm,c('nFeature_CellToCell'),group.by = 'Orig_ID.Sending',raster=F,pt.size = 0.1)+ggtitle('After Filtration')+ylab('nFeature_CellToCell')

# Make copy just in case
mes.imm.eng.CTC.bySample <- eng.CTC.bySample_mes.imm

# Embedding and clustering
mes.imm.eng.CTC.bySample <- ScaleData(mes.imm.eng.CTC.bySample)
mes.imm.eng.CTC.bySample <- FindVariableFeatures(mes.imm.eng.CTC.bySample)

mes.imm.eng.CTC.bySample <- RunPCA(mes.imm.eng.CTC.bySample,npcs = 100)
setwd("~/Desktop/iScience Transfer Submission/Additional Analyses/Mesenchyme_Immune_CircuitAnalysis")
pdf(file='Mes.Imm_Subset.CTC.bySample.PCs.pdf',width=10,height=8)
ElbowPlot(mes.imm.eng.CTC.bySample,ndims = 100)
PCHeatmap(mes.imm.eng.CTC.bySample,cells=200,balanced=T,dims=1:9) # 1-9
PCHeatmap(mes.imm.eng.CTC.bySample,cells=200,balanced=T,dims=10:18)
PCHeatmap(mes.imm.eng.CTC.bySample,cells=200,balanced=T,dims=19:27)
PCHeatmap(mes.imm.eng.CTC.bySample,cells=200,balanced=T,dims=28:36)
PCHeatmap(mes.imm.eng.CTC.bySample,cells=200,balanced=T,dims=37:45)
PCHeatmap(mes.imm.eng.CTC.bySample,cells=200,balanced=T,dims=46:54)
PCHeatmap(mes.imm.eng.CTC.bySample,cells=200,balanced=T,dims=55:63)
dev.off()

# Embedding, finding nearest neighbors, and finding clusters (handpicking PCs)
mes.imm.eng.CTC.bySample <- RunUMAP(mes.imm.eng.CTC.bySample,dims = c(1:9))
mes.imm.eng.CTC.bySample <- FindNeighbors(mes.imm.eng.CTC.bySample,dims = c(1:9))
mes.imm.eng.CTC.bySample <- FindClusters(mes.imm.eng.CTC.bySample,resolution=1.3) # Deliberately choosing this high resolution so that the circular cluster in the middle can be 
DimPlot(mes.imm.eng.CTC.bySample,label=T,raster = F,cols=full.col.pal)

# Plot specifics for condition sending, cell types (sending and receiving)
DimPlot(mes.imm.eng.CTC.bySample,label=T,raster = F, cols=full.col.pal,split.by = 'Condition.Sending',ncol=3)
DimPlot(mes.imm.eng.CTC.bySample,group.by = 'CellType.sub.Sending',raster=F,shuffle = T,cols=full.col.pal)
DimPlot(mes.imm.eng.CTC.bySample,group.by = 'CellType.sub.Receiving',raster=F,shuffle = T,cols=full.col.pal)

# Want to see three plots side-by-side
cols.cellclass = c(
  "Epithelium" = "#D982C6",
  "Immune" = "#87B37A",
  "Mesenchyme" = "#F4A261")
bb <- DimPlot(
  mes.imm.eng.CTC.bySample,
  group.by = 'CellClass.sub.Sending',
  raster = FALSE,
  shuffle = TRUE,
  cols = cols.cellclass) +
  ggtitle("Sending Lineage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()
bb
cc <- DimPlot(
  mes.imm.eng.CTC.bySample,
  group.by = 'CellClass.sub.Receiving',
  raster = FALSE,
  shuffle = TRUE,
  cols = cols.cellclass) +
  ggtitle("Receiving Lineage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()
cc
send.rec.subset = bb | cc
print(send.rec.subset)

DimPlot(
  mes.imm.eng.CTC.bySample,
  group.by = 'CellClass.sub.Sending',
  raster = FALSE,
  shuffle = TRUE,
  cols = cols.cellclass, split.by = "Condition.Sending") +
  ggtitle("Sending Lineage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()

# Object
mes.imm.eng.CTC.bySample
setwd("~/Desktop/iScience Transfer Submission/Additional Analyses/Mesenchyme_Immune_CircuitAnalysis")
save(mes.imm.eng.CTC.bySample, file = "mes.imm.eng.CTC.bySample.Robj")

Idents(mes.imm.eng.CTC.bySample) = mes.imm.eng.CTC.bySample$CellType.sub.Sending # Below adjusting min.pct and logfc.threshold, but 0.1 is going to give us the most specific markers/genes
mes_imm.mark <- FindAllMarkers(mes.imm.eng.CTC.bySample, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1, assay = "CellToCell") # Can change the min.pct as needed
mes_imm.mark$ratio = mes_imm.mark$pct.1/mes_imm.mark$pct.2 # Creating our ratio slot by deviding pct.1 by pct.2
mes_imm.mark$power = mes_imm.mark$ratio*mes_imm.mark$avg_log2FC # Creating our power slow by dividing ratio by avg_log2FC
View(mes_imm.mark)

Idents(mes.imm.eng.CTC.bySample) = mes.imm.eng.CTC.bySample$CellType.sub.Joint # Below adjusting min.pct and logfc.threshold, but 0.1 is going to give us the most specific markers/genes
mes_imm.mark <- FindAllMarkers(mes.imm.eng.CTC.bySample, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1, assay = "CellToCell") # Can change the min.pct as needed
mes_imm.mark$ratio = mes_imm.mark$pct.1/mes_imm.mark$pct.2 # Creating our ratio slot by deviding pct.1 by pct.2
mes_imm.mark$power = mes_imm.mark$ratio*mes_imm.mark$avg_log2FC # Creating our power slow by dividing ratio by avg_log2FC
View(mes_imm.mark)

FeaturePlot(mes.imm.eng.CTC.bySample, features = c("Tnc—Itgav"))
FeaturePlot(mes.imm.obj, features = c("Bdnf","Ddr1"))
FeaturePlot(mes.imm.obj, features = c("Ncam1","Ptpra"))
FeaturePlot(mes.imm.obj, features = c("Mmp2","Sdc2"))
FeaturePlot(mes.imm.obj, features = c("Cxcl12","Itgb1"))
FeaturePlot(mes.imm.obj, features = c("Cxcl12","Sdc4"))
FeaturePlot(mes.imm.eng.CTC.bySample, features = c("Cxcl12—Sdc4"))
FeaturePlot(mes.imm.obj, features = c("Ngf","Kidins220"), order = T)
FeaturePlot(mes.imm.eng.CTC.bySample, features = c("Ngf—Kidins220"), split.by = "Condition.Sending")

Ngf—Kidins220


Idents(mes.imm.obj1) = mes.imm.obj1$CellType.sub # Below adjusting min.pct and logfc.threshold, but 0.1 is going to give us the most specific markers/genes
mes_imm.mark.trans <- FindAllMarkers(mes.imm.obj1, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1, assay = "RNA") # Can change the min.pct as needed
mes_imm.mark.trans$ratio = mes_imm.mark.trans$pct.1/mes_imm.mark.trans$pct.2 # Creating our ratio slot by deviding pct.1 by pct.2
mes_imm.mark.trans$power = mes_imm.mark.trans$ratio*mes_imm.mark.trans$avg_log2FC # Creating our power slow by dividing ratio by avg_log2FC
View(mes_imm.mark.trans)
