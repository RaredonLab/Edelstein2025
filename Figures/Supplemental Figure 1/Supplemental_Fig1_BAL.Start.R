
# Supplemental Data: BAL Starting Population - Subclustering Epithelium
# close all, clear all
graphics.off()  
rm(list = ls())

# Load packages
library(Seurat)
library(ggplot2)
library(SeuratObject)
library(dplyr)
library(tidyr)

# Set working directory to be true and load data
setwd("~/Desktop/Datasets/Starting Populations")
load("BAL.integrated_03-29-2025.Robj")

# Check out the structure of the data
str(BAL.integrated@meta.data)
table(BAL.integrated$CellClass_Integrated)
# Take a look at the integrated embedding
DimPlot(BAL.integrated, reduction = 'umap.rpca', group.by = "CellType")
# Check the current cell type meta-data designations
table(BAL.integrated$CellType)

# Rename "LEPs" to "Pan_Epithelium" in CellType
BAL.integrated$CellType <- ifelse(
  BAL.integrated$CellType == "LEPs",
  "Pan_Epithelium",
  BAL.integrated$CellType)
# Check updated table
table(BAL.integrated$CellType)
# Make sure you save this update
setwd("~/Desktop/Datasets/Starting Populations")
save(BAL.integrated, file = "BAL.integrated_03-29-2025.Robj")

# Subset out just the epithelium
Idents(BAL.integrated) = BAL.integrated$CellClass_Integrated
BAL.epi.sub = subset(BAL.integrated, idents = "Epithelium")

# Scale
BAL.epi.sub = ScaleData(BAL.epi.sub)
# Normalizing the data
BAL.epi.sub = NormalizeData(BAL.epi.sub)
# Finding variable features
BAL.epi.sub = FindVariableFeatures(BAL.epi.sub)

# Principle component analysis
BAL.epi.sub = RunPCA(BAL.epi.sub, npcs = 100)
# Visualizing the PCAs
ElbowPlot(BAL.epi.sub,ndims = 100)
PCHeatmap(BAL.epi.sub,cells=200,balanced=T,dims=1:9)
PCHeatmap(BAL.epi.sub,cells=200,balanced=T,dims=10:18)
PCHeatmap(BAL.epi.sub,cells=200,balanced=T,dims=19:27)
PCHeatmap(BAL.epi.sub,cells=200,balanced=T,dims=28:36)
PCHeatmap(BAL.epi.sub,cells=200,balanced=T,dims=37:45)
dev.off()

# UMAP
BAL.epi.sub = RunUMAP(BAL.epi.sub, dims = 1:5)
DimPlot(BAL.epi.sub, label = TRUE)
# Cluster data (creating nearest neighbor graph, not clustering)
BAL.epi.sub = FindNeighbors(BAL.epi.sub, dims = 1:5)
# Defining clusters
BAL.epi.sub = FindClusters(BAL.epi.sub, res = 0.3)
# See what this looks like
DimPlot(BAL.epi.sub)

# Get a sense of a few common features
FeaturePlot(BAL.epi.sub, features = c("Epcam","Ccdc153","Krt8","Sox2","Sox9"))

# Plot original object, grouped by CellType
DimPlot(BAL.integrated, reduction = 'umap.rpca', group.by = "CellType")
# Check cell type labels in the original
table(BAL.integrated$CellType)

# Subset again for no immune contamination
BAL.epi.sub1 = subset(BAL.epi.sub, idents = c("0"),invert = TRUE)

# Re-scale
BAL.epi.sub1 = ScaleData(BAL.epi.sub1)
# Normalizing the data
BAL.epi.sub1 = NormalizeData(BAL.epi.sub1)
# Finding variable features
BAL.epi.sub1 = FindVariableFeatures(BAL.epi.sub1)

# Principle component analysis
BAL.epi.sub1 = RunPCA(BAL.epi.sub1, npcs = 100)
# Visualizing the PCAs
ElbowPlot(BAL.epi.sub1,ndims = 100)
PCHeatmap(BAL.epi.sub1,cells=200,balanced=T,dims=1:9)
PCHeatmap(BAL.epi.sub1,cells=200,balanced=T,dims=10:18)
PCHeatmap(BAL.epi.sub1,cells=200,balanced=T,dims=19:27)
PCHeatmap(BAL.epi.sub1,cells=200,balanced=T,dims=28:36)
PCHeatmap(BAL.epi.sub1,cells=200,balanced=T,dims=37:45)
dev.off()

# UMAP; cherry-picking PCs 
BAL.epi.sub1 = RunUMAP(BAL.epi.sub1, dims =c(1:3, 5:6, 9))
# Look at embedding before finding neighbors
DimPlot(BAL.epi.sub1, label = TRUE)
# cluster data (creating nearest neighbor graph, not clustering)
BAL.epi.sub1 = FindNeighbors(BAL.epi.sub1, dims =c(1:3, 5:6, 9))
# Defining clusters
BAL.epi.sub1 = FindClusters(BAL.epi.sub1, res = 1.2)
# See what this looks like
DimPlot(BAL.epi.sub1)

# Plot some common markers
FeaturePlot(BAL.epi.sub1, features = c("Epcam","Ptprc","Ccdc153","Krt8","Sox2","Sox9"))

# Subset again (2nd time) - 1 is contamination w/ immune
BAL.epi.sub2 = subset(BAL.epi.sub1, idents = c("1"),invert = TRUE)

# Re-scale
BAL.epi.sub2 = ScaleData(BAL.epi.sub2)
# Normalizing the data
BAL.epi.sub2 = NormalizeData(BAL.epi.sub2)
# Finding variable features
BAL.epi.sub2 = FindVariableFeatures(BAL.epi.sub2)

# Principle component analysis
BAL.epi.sub2 = RunPCA(BAL.epi.sub2, npcs = 100)
# Visualizing the PCAs
ElbowPlot(BAL.epi.sub2,ndims = 100)
PCHeatmap(BAL.epi.sub2,cells=200,balanced=T,dims=1:9)
PCHeatmap(BAL.epi.sub2,cells=200,balanced=T,dims=10:18)
PCHeatmap(BAL.epi.sub2,cells=200,balanced=T,dims=19:27)
PCHeatmap(BAL.epi.sub2,cells=200,balanced=T,dims=28:36)
PCHeatmap(BAL.epi.sub2,cells=200,balanced=T,dims=37:45)
dev.off()

# UMAP
BAL.epi.sub2 = RunUMAP(BAL.epi.sub2, dims = 1:7)
DimPlot(BAL.epi.sub2, label = TRUE)
# cluster data (creating nearest neighbor graph, not clustering)
BAL.epi.sub2 = FindNeighbors(BAL.epi.sub2, dims = 1:7)
# Defining clusters
BAL.epi.sub2 = FindClusters(BAL.epi.sub2, res = 4.0)
# See what this looks like
DimPlot(BAL.epi.sub2, label = T)

# Get a sense of features
FeaturePlot(BAL.epi.sub2, features = c("Ccdc153","Pifo","Epcam",
                                      "Krt8","Krt18","Pou2f3","Trpm5",
                                      "Aqp5","Napsa","Lamp3","Rtkn2","Pdpn"), order = T)
FeaturePlot(BAL.epi.sub2, features = c("Sox9","Scgb3a2","Aqp5","Sftpc","Ager","Sema3e","Sox2"), label = T)
FeaturePlot(BAL.epi.sub2, features = c("Sox9","Dclk1","Krt18","Trpm5","Scgb1a1","Msln"), label = T)
FeaturePlot(BAL.epi.sub2, features = c("Epcam","Ptprc"), label = T)

# 3 is also contamination - remove
# Subset again
BAL.epi.sub3 = subset(BAL.epi.sub2, idents = c("3"),invert = TRUE)

# Re-scale
BAL.epi.sub3 = ScaleData(BAL.epi.sub3)
# Normalizing the data
BAL.epi.sub3 = NormalizeData(BAL.epi.sub3)
# Finding variable features
BAL.epi.sub3 = FindVariableFeatures(BAL.epi.sub3)

# Principle component analysis
BAL.epi.sub3 = RunPCA(BAL.epi.sub3, npcs = 100)
# Visualizing the PCAs
ElbowPlot(BAL.epi.sub3,ndims = 100)
PCHeatmap(BAL.epi.sub3,cells=200,balanced=T,dims=1:9)
PCHeatmap(BAL.epi.sub3,cells=200,balanced=T,dims=10:18)
PCHeatmap(BAL.epi.sub3,cells=200,balanced=T,dims=19:27)
PCHeatmap(BAL.epi.sub3,cells=200,balanced=T,dims=28:36)
PCHeatmap(BAL.epi.sub3,cells=200,balanced=T,dims=37:45)
dev.off()

# UMAP
BAL.epi.sub3 = RunUMAP(BAL.epi.sub3, dims =c(1:5))
DimPlot(BAL.epi.sub3, label = TRUE)
# cluster data (creating nearest neighbor graph, not clustering)
BAL.epi.sub3 = FindNeighbors(BAL.epi.sub3, dims =c(1:5))
# Defining clusters
BAL.epi.sub3 = FindClusters(BAL.epi.sub3, res = 2.8)
# See what this looks like
DimPlot(BAL.epi.sub3, label = T)

# Cluster 0 = Krt15+_Epi
# Cluster 1 = Ciliated
# Cluster 2 = Ciliated
# Cluster 3 = Secretory
# Cluster 4 = Secretory
# Cluster 5 = Ciliated
# Cluster 6 = Secretory
# Cluster 7 = Sox9+_Epi

# FeaturePlots for common markers in marker list
FeaturePlot(BAL.epi.sub3, features = c("Aqp5","Scgb3a2","Ccdc153","Sema3e","Sox9","Krt15","Sox2"), order = T)
FeaturePlot(BAL.epi.sub3)
# Developing marker list
mark1 = FindAllMarkers(BAL.epi.sub3, min.pct = 0.1,logfc.threshold = 0.1)
mark1$ratio = mark1$pct.1/mark1$pct.2
mark1$power = mark1$ratio*mark1$avg_log2FC
# View markers for clusters
View(mark1)

## Annotating
# Create empty meta-data column
BAL.epi.sub3$CellType.epi <- NA
# Annotate cells for other clusters in the global object
c0 <- WhichCells(BAL.epi.sub3, idents = "0")
c1 <- WhichCells(BAL.epi.sub3, idents = "1")
c2 <- WhichCells(BAL.epi.sub3, idents = "2")
c3 <- WhichCells(BAL.epi.sub3, idents = "3")
c4 <- WhichCells(BAL.epi.sub3, idents = "4")
c5 <- WhichCells(BAL.epi.sub3, idents = "5")
c6 <- WhichCells(BAL.epi.sub3, idents = "6")
c7 <- WhichCells(BAL.epi.sub3, idents = "7")

BAL.epi.sub3$CellType.epi[c0] <- 'Krt15+_Epi'
BAL.epi.sub3$CellType.epi[c1] <- 'Ciliated'
BAL.epi.sub3$CellType.epi[c2] <- 'Ciliated'
BAL.epi.sub3$CellType.epi[c3] <- 'Secretory'
BAL.epi.sub3$CellType.epi[c4] <- 'Secretory'
BAL.epi.sub3$CellType.epi[c5] <- 'Ciliated'
BAL.epi.sub3$CellType.epi[c6] <- 'Secretory'
BAL.epi.sub3$CellType.epi[c7] <- 'Tuft'

# Check the distribution to confirm all annotations were applied
table(BAL.epi.sub3$CellType.epi)
View(BAL.epi.sub3@meta.data)
DimPlot(BAL.epi.sub3, group.by = 'CellType.epi', label = T, repel = T)

# Rename and save
BAL.epi.subset = BAL.epi.sub3
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 1")
save(BAL.epi.subset, file = "BAL.epi.subset.Robj")
# load("BAL.epi.subset.Robj")

# Define CellType color palette
CellType.cols <- c(
  "ATI_Like" = "#b22222","ATI" = "#b22222",
  "ATII" = "#4682B4","ATII_Like" = "#4682B4",
  "Secretory_Like" = "#2E8B57", "Secretory" = "#2E8B57",
  "Ciliated" = "#9B489B","Ciliated_Like" = "#9B489B",
  "Stressed_Progenitor" = "#c71585",
  "Hillock_Luminal" = "#ff0000","Hillock_Like" = "#ff0000",
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
  "Pro_Inflamm_Mac" = "#9497fd",
  "Anti_Inflamm_Mac" = "#ff9e80","Polarized_Mac" = "#ff9e80",
  "Rspo3+_Mes" = "#89a5a5","Actc1_Mural" = "#1179fa",
  "Cycling_Distal_Epi" = "#93d0ff",
  "Cycling_Proximal_Epi" = "#c4c2ff",
  "Pdgfrb+_Pericyte" = "#F4AFB4",
  "Cell_Cycle" = "#e0bf1b","Cycling_Immune" = "#e0bf1b",
  "Fzd7+_Stressed" = "#2F004F",
  "Pan_Epithelium" = "#E2AEDD") # Add LEPs explicitly if needed

# So now let's pull in our whole object first and get the colors right
DimPlot(BAL.integrated, reduction = "umap.rpca", group.by = "CellType", cols = CellType.cols)

# Make cleaner plot
DimPlot(BAL.integrated, reduction = "umap.rpca", group.by = "CellType", cols = CellType.cols, label = T, repel = T) +
  theme_void() +
  theme(legend.position = "right", legend.background = element_blank(),
    legend.key = element_blank(), legend.text = element_text(size = 12), legend.key.size = unit(0.9, "cm"),     
    legend.title = element_text(size = 14, face = "bold"),  # bold and larger title
    plot.title = element_blank(), plot.background = element_blank(), panel.background = element_blank(),
    panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) +
  labs(color = "Cell Type") +  # this sets the legend title
  guides(color = guide_legend(override.aes = list(size = 5)))  # controls circle size in legend

## No legend
BAL_UMAP = DimPlot(BAL.integrated, reduction = "umap.rpca", group.by = "CellType", cols = CellType.cols, label = TRUE, repel = TRUE) +
  theme_void() +
  theme(legend.position = "none", plot.title = element_blank(), plot.background = element_blank(),
    panel.background = element_blank(), panel.grid = element_blank(),axis.text = element_blank(),
    axis.ticks = element_blank(),axis.title = element_blank())
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 1")
ggsave("BAL_UMAP.glbl.png", plot = BAL_UMAP, width = 7, height = 6, dpi = 600)

## Now make a umap of just the epi subset
DimPlot(BAL.epi.subset, group.by = "CellType.epi")

# Need to make color palette
# Conserving colors from broader palette, where possible
BAL.epi.cols = c("Ciliated" = "#9B489B",
                 "Krt15+_Epi" = "#FF531A",
                 "Secretory" = "#2E8B57",
                 "Tuft" = "#15DDB5")

BAL.epi.sub_UMAP = DimPlot(BAL.epi.subset, group.by = "CellType.epi", cols = BAL.epi.cols, label = TRUE, repel = TRUE) +
  theme_void() +
  theme(legend.position = "none", plot.title = element_blank(), plot.background = element_blank(),
    panel.background = element_blank(), panel.grid = element_blank(),
    axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
BAL.epi.sub_UMAP
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 1")
ggsave("BAL.epi.sub_UMAP.png", plot = BAL.epi.sub_UMAP, width = 4, height = 4, dpi = 600)

Idents(BAL.integrated) = BAL.integrated$CellType
BAL.mark <- FindAllMarkers(
  BAL.integrated,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  assay = "RNA")
View(BAL.mark)

# Now dot plot to show expression by cell type
FeaturePlot(BAL.integrated, reduction = "umap.rpca",features = c("Abcd2"))

BAL.integrated$CellType <- factor(
  BAL.integrated$CellType,
  levels = c("B", "Mac_Alv", "Monocytes", "Neutrophils", "T", "Pan_Epithelium", "Cell_Cycle"))

genes.to.plot <- c("Cd79b", "Bank1","Pax5","Cr2", # B cells
                   "Mrc1","Abcd2","Siglec10","Spic", # Mac_Alv
                   "Clec4d","Itgam","Cd99", # Mono
                   "Ifit3","Mx1","Cd84","Csf3r", # Neutro
                   "Cd3g", "Cd3e", "Cd3d","Skap1", # T
                   "Scgb3a2","Krt18","Krt8","Ccdc153","Sox9","Dclk1", # Epi
                   "Mki67","Top2a","Ube2c","Ccna2") #Cell cycle
# Make comprehensive dot plot for whole object
DotPlot(
  BAL.integrated,
  features = genes.to.plot,
  group.by = "CellType") +
  scale_color_gradient(low = "#f0f0f0", high = "#08306b", name = "Avg. Expression") +
  scale_size(range = c(1, 5), name = "% Expressing") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, face = "plain", color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, face = "plain", color = "black"),
    axis.title.x = element_text(size = 14, face = "bold", color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, face = "bold", color = "black", margin = margin(r = 10)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold", color = "black"),
    legend.text = element_text(size = 12, color = "black")) +
  labs(x = "Gene", y = "Cell Type")

# With axis lines
BAL.dotplot = DotPlot(
  BAL.integrated,
  features = genes.to.plot,
  group.by = "CellType") +
  scale_color_gradientn(
    colors = c("white", "#d80032"),
    name = "Avg. Expression") +
  scale_size(range = c(1, 5), name = "% Expressing") +
  theme_minimal(base_size = 12) +
  theme(
    # Axis lines
    axis.line.x = element_line(color = "black", size = 0.3),
    axis.line.y = element_line(color = "black", size = 0.3),
    # Axis tick labels (plain, not bold)
    axis.text.x = element_text(size = 12, face = "plain", color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, face = "plain", color = "black"),
    # Axis titles (bold, black)
    axis.title.x = element_text(size = 14, face = "bold", color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, face = "bold", color = "black", margin = margin(r = 10)),
    # Remove grid and background
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    # Legend
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold", color = "black"),
    legend.text = element_text(size = 12, color = "black")) +
  labs(x = "Gene", y = "Cell Type")
BAL.dotplot
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 1")
ggsave("BAL.dotplot.celltype.png", plot = BAL.dotplot, width = 10, height = 5, dpi = 600)

## Make BAL epi dot plot (this is just on the subset epithelium)
BAL.epi.subset$CellType.epi <- factor(
  BAL.epi.subset$CellType.epi,
  levels = c("Ciliated","Secretory","Krt15+_Epi","Sox9+_Epi"))

# Cherry-picking the candidate markers for each cell type
genes.to.plot1 <- c("Ccdc153","Pifo","Dnah6","Scgb3a2","Chia","Bpifa5","Krt15","Aqp5","Tgfb2","Nkain3",
                    "Sox9","Dclk1","Pou2f3","Trpm5") 

# Make dot plot
BAL.epi.dotplot = DotPlot(BAL.epi.subset, features = genes.to.plot1, group.by = "CellType.epi") +
  scale_color_gradientn(colors = c("white", "#2B6F1B"), name = "Avg. Expression") +
  scale_size(range = c(1, 5), name = "% Expressing") +
  theme_minimal(base_size = 14) +
  theme(axis.line.x = element_line(color = "black", size = 0.3),
    axis.line.y = element_line(color = "black", size = 0.3),
    # Axis tick labels (plain, not bold)
    axis.text.x = element_text(size = 12, face = "plain", color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, face = "plain", color = "black"),
    # Axis titles (bold, black)
    axis.title.x = element_text(size = 14, face = "bold", color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, face = "bold", color = "black", margin = margin(r = 10)),
    # Remove grid and background
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    # Legend
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold", color = "black"),
    legend.text = element_text(size = 12, color = "black")) +
  labs(x = "Gene", y = "Cell Type")
BAL.epi.dotplot
# Save the figure (make sure working directory is set correctly)
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 1")
ggsave("BAL.epi.sub_dotplot.celltype.png", plot = BAL.epi.dotplot, width = 6.5, height = 5, dpi = 600)

### Making select feature plots based on text references (section 1 of results)
FeaturePlot(BAL.integrated, reduction = "umap.rpca", features = c("Pparg","Prodh2","Cd79b","Cd3e","Itgam","Krt8"))
# Define genes to plot
features_to_plot <- c("Pparg", "Prodh2", "Cd79b", "Cd3e", "Itgam", "Krt8")
# Function to create individual FeaturePlot
feature_umap <- function(gene) {
  FeaturePlot(
    BAL.integrated,
    features = gene,
    reduction = "umap.rpca",
    pt.size = 0.6) +
    ggtitle(gene) +
    theme_void() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 10),           # increased font size
      legend.key.height = unit(0.65, "cm"),             # taller colorbar
      legend.key.width = unit(0.4, "cm"),              # wider colorbar
      legend.margin = margin(t = 2, b = 2, unit = "mm"),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
    )
}

# Generate and arrange the plots
plots_list <- lapply(features_to_plot, feature_umap)
gene_plot_grid <- wrap_plots(plots_list, ncol = 3)
gene_plot_grid
# Save figure
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 1")
ggsave("FeaturePlot_Panel.png", gene_plot_grid, width = 9, height = 7.75, dpi = 600)



