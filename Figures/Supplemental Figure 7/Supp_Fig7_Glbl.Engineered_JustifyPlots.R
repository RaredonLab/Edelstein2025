
# close all, clear all
graphics.off()  
rm(list = ls())

library(Seurat)
library(ggplot2)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(cowplot)
library(grid)
library(forcats)
library(patchwork)
library(viridis)  # If not installed: install.packages("viridis")

## MEMO: Supplemental Figure 7: Engineered Object (Global) Cluster Justifications

# Set working directory to be true
setwd("~/Desktop/Datasets/Engineered Global Objects")
# Load object
load("eng.subset.integrated_NodeAligned.Robj")
# Inspect the cell type meta-data
table(eng.subset.integrated_HK$CellType.NodeAligned)

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
  "Pdgfrb+_Pericytes" = "#F4AFB4",
  "Cell_Cycle" = "#e0bf1b","Cycling_Immune" = "#e0bf1b",
  "Fzd7+_Stressed" = "#2F004F",
  "Pan_Epithelium" = "#E2AEDD") 

# Check cell type annotations
table(eng.subset.integrated_HK$CellType.NodeAligned)

# Ensure CellType.NodeAligned is a factor with correct levels
eng.subset.integrated_HK$CellType.NodeAligned <- factor(
  eng.subset.integrated_HK$CellType.NodeAligned,
  levels = names(CellType.cols))

global.dimplot = DimPlot(eng.subset.integrated_HK,
  group.by = "CellType.NodeAligned", reduction = "umap.rpca",
  cols = CellType.cols) +
  theme_void() +  # Removes axes, grid, background
  theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm"), plot.title = element_blank())
global.dimplot
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 7")
ggsave("global.dimplot.png", plot = global.dimplot, width = 7, height = 7, dpi = 600)

## Make legend that is ordered (with a dummy plot)
# Set order
celltype_levels <- c(
  "ATI_Like", "ATII_Like", "Basal_Like", "Ciliated", "Hillock_Like", "RAS_Like", "Secretory",
  "Stressed_Progenitor", "Cycling_Distal_Epi", "Cycling_Proximal_Epi", 
  "Polarized_Mac", "Cycling_Immune",
  "Pdgfrb+_Pericytes", "Rspo3+_Mes")
# Create dummy data-frame
legend_df <- data.frame(
  CellType = factor(celltype_levels, levels = celltype_levels),
  dummy_x = 1,
  dummy_y = 1)

# Dummy plot for legend
legend_plot <- ggplot(legend_df, aes(x = dummy_x, y = dummy_y, color = CellType)) +
  geom_point(size = 5) +
  scale_color_manual(values = CellType.cols, breaks = celltype_levels) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
  theme_void() +
  theme(legend.position = "right", legend.text = element_text(size = 10),
    legend.key.size = unit(0.7, "cm"), legend.margin = margin(5, 5, 5, 5),
    legend.spacing.y = unit(0.4, "cm"))
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 7")
ggsave("GlblObject_Legend.Standalone.png", plot = legend_plot, width = 4, height = 10, dpi = 600)

## Now plots for CellClass and Phase
table(eng.subset.integrated_HK$CellClass.NodeAligned)
# Set order for lineage (for legend)
eng.subset.integrated_HK$CellClass.NodeAligned <- factor(
  eng.subset.integrated_HK$CellClass.NodeAligned,
  levels = c("Epithelium", "Immune", "Mesenchyme"))

# Colors from project color palette - lineage
CellClass.cols <- c(
  "Epithelium" = "#D982C6",
  "Immune" = "#87B37A",
  "Mesenchyme" = "#F4A261",
  "Endothelium" = "#2A9D8F")

# Make plot without axes, but legend and labels
umap_cellclass_plot <- DimPlot(
  eng.subset.integrated_HK,
  group.by = "CellClass.NodeAligned",
  reduction = "umap.rpca",
  cols = CellClass.cols,
  label = FALSE,
  repel = TRUE) +
  theme_void() +
  theme(legend.position = "right", legend.text = element_text(size = 12),
    legend.key.size = unit(0.8, "cm"), plot.margin = unit(c(0, 0, 0, 0), "cm"),
    plot.title = element_blank()) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 6), order = 1))
umap_cellclass_plot
# Save
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 7")
ggsave("GlblObject_CellClass_UMAP.png", plot = umap_cellclass_plot, width = 7, height = 6, dpi = 600)

## Now for phase
# Define Phase colors
phase_colors <- c(
  "S" = "#16CEA6",
  "G1" = "#D2691E",
  "G2M" = "#AEC5EB")
# Set order; not necessary but doing it
eng.subset.integrated_HK$Phase <- factor(
  eng.subset.integrated_HK$Phase,
  levels = c("G1", "S", "G2M"))
# Make UMAP
phase_umap_plot <- DimPlot(eng.subset.integrated_HK, group.by = "Phase", reduction = "umap.rpca", cols = phase_colors,
  label = FALSE) +
  theme_void() +
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.key.size = unit(0.8, "cm"),
    plot.margin = unit(c(0, 0, 0, 0), "cm"), plot.title = element_blank()) +
  guides(color = guide_legend(
    title = NULL,
    override.aes = list(size = 6)))
phase_umap_plot
# Save
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 7")
ggsave("GlblObject_Phase_UMAP.png", plot = phase_umap_plot, width = 7, height = 6, dpi = 600)

### Dot Plot for expression
# Run quick marker list just to double check (I know for the most part what each consists of but...)
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType.NodeAligned
Glbl.mark <- FindAllMarkers(
  eng.subset.integrated_HK,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.1,
  assay = "RNA")
View(Glbl.mark)

## Set order
eng.subset.integrated_HK$CellType.NodeAligned <- forcats::fct_relevel(
  eng.subset.integrated_HK$CellType.NodeAligned,
  "ATI_Like", "ATII_Like", "Basal_Like", "Ciliated", "Hillock_Like", "RAS_Like", 
  "Secretory", "Stressed_Progenitor", "Cycling_Distal_Epi", "Cycling_Proximal_Epi", 
  "Polarized_Mac", "Cycling_Immune", "Pdgfrb+_Pericytes", "Rspo3+_Mes")

# Genes by cell type we want to plot
genes.to.plot <- c("Pdpn","Ager","Wnt7a","Akap5", #ATI
                   "Defb3","S100g","Hhip","Npw", # ATII
                   "Krt5","Krt14","Wnt10a","Col17a1", # Basal
                   "Ccdc153","Pifo","Dnah6","Ak9", # Ciliated
                   "Krt13","Sprr1a","Krt16","Cnfn", # Hillock_Like
                   "Scgb3a2","Sftpc","Scgb1a1","Sftpb", # RAS_Like
                   "Scgb3a1","Bpifa1","Agr2","Muc5b", # Secretory
                   "Sox9","Cox4i2","Stc1","Ptges", # Stressed_Progenitor
                   "Lamp3", "Sema3e","Top2a","Mki67", # Cycling_Distal_Epi
                   "Tp63","Krt17","Ube2c","Ccna2", #  Cycling_Prox_Epi
                   "Spic","Siglec1","Pparg","Prodh2", # Polarized_Mac
                   "C1qb","Il1b","Itgam","Bub1", # Cycling_Immune
                   "Pdgfrb","Gucy1b1","Postn","Fgf10", # Pericytes
                   "Rspo3","Des","Pdgfra","Twist1") # Rspo3+_Mes
# Dot plot with axis lines
eng.glbl.dotplot = DotPlot(
  eng.subset.integrated_HK,
  features = genes.to.plot,
  group.by = "CellType.NodeAligned") +
  scale_color_gradientn(colors = c("white", "#E04B00"), name = "Avg. Expression") +
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
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank(),
    # Legend
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold", color = "black"),
    legend.text = element_text(size = 12, color = "black")) +
  labs(x = "Gene", y = "Cell Type")
eng.glbl.dotplot
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 7")
ggsave("eng.glbl.dotplot_celltype.png", plot = eng.glbl.dotplot, width = 20, height = 8, dpi = 600)

## Now some feature plots
# Checking out the features we want
FeaturePlot(eng.subset.integrated_HK, reduction = "umap.rpca", 
            features = c("Rtkn2","Col4a4","Napsa","Sftpc","Krt5",
             "Ccdc153","Krt13","Scgb3a2","Sox9","Cox4i2"), ncol = 5)
## EPI FIRST
features_to_plot <- c("Rtkn2","Akap5","Napsa","Sftpc","Krt5",
                      "Ccdc153","Krt13","Scgb3a2","Sox9","Cox4i2")
# Create plotting function
feature_umap <- function(gene) {
  FeaturePlot(
    eng.subset.integrated_HK,
    features = gene,
    reduction = "umap.rpca",
    pt.size = 0.6) +
    scale_color_viridis_c(option = "D") +  # Viridis color scale
    ggtitle(gene) +
    theme_void() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.title = element_blank(), legend.text = element_text(size = 10),
      legend.key.height = unit(0.65, "cm"), legend.key.width = unit(0.4, "cm"),
      legend.margin = margin(t = 2, b = 2, unit = "mm"), plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
}
# Generate and arrange the plots
plots_list <- lapply(features_to_plot, feature_umap)
gene_plot_grid_epi <- wrap_plots(plots_list, ncol = 5) # Two rows of 5 each
gene_plot_grid_epi
# Save plot
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 7")
ggsave("Glbl_Features_Epi.png", gene_plot_grid_epi, width = 18, height = 6.5, dpi = 600)

## Now Immune
FeaturePlot(eng.subset.integrated_HK, reduction = "umap.rpca", 
            features = c("Prodh2","Spic","Siglec10","Il1b"), ncol = 2)
## Immune genes
features_to_plot <- c("Prodh2","Spic","Siglec10","Il1b")
# Create plotting function
feature_umap <- function(gene) {
  FeaturePlot(
    eng.subset.integrated_HK,
    features = gene,
    reduction = "umap.rpca",
    pt.size = 0.6, order = T) +
    scale_color_viridis_c(option = "plasma") +  # Viridis color scale
    ggtitle(gene) +
    theme_void() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.title = element_blank(), legend.text = element_text(size = 10),
          legend.key.height = unit(0.65, "cm"), legend.key.width = unit(0.4, "cm"),
          legend.margin = margin(t = 2, b = 2, unit = "mm"), plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
}
# Generate and arrange the plots
plots_list <- lapply(features_to_plot, feature_umap)
gene_plot_grid_imm <- wrap_plots(plots_list, ncol = 2) # Two rows of 5 each
gene_plot_grid_imm
# Save plot
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 7")
ggsave("Glbl_Features_Imm.png", gene_plot_grid_imm, width = 8, height = 6.5, dpi = 600)

## Now Mesenchyme
FeaturePlot(eng.subset.integrated_HK, reduction = "umap.rpca", 
            features = c("Pdgfrb","Gucy1b1","Postn","Rspo3","Des","Pdgfra"), ncol = 3, order = T)
## Mesenchyme genes
features_to_plot <- c("Pdgfrb","Gucy1b1","Fgf10","Rspo3","Des","Pdgfra")
# Create plotting function
feature_umap <- function(gene) {
  FeaturePlot(
    eng.subset.integrated_HK,
    features = gene,
    reduction = "umap.rpca",
    pt.size = 0.6, order = T) +
    scale_color_viridis_c(option = "cividis") +  # Change color palette for mesenchyme
    ggtitle(gene) +
    theme_void() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.title = element_blank(), legend.text = element_text(size = 10),
          legend.key.height = unit(0.65, "cm"), legend.key.width = unit(0.4, "cm"),
          legend.margin = margin(t = 2, b = 2, unit = "mm"), plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
}
# Generate and arrange the plots
plots_list <- lapply(features_to_plot, feature_umap)
gene_plot_grid_mes <- wrap_plots(plots_list, ncol = 3) # Two rows of 5 each
gene_plot_grid_mes
# Save plot
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 7")
ggsave("Glbl_Features_Mes.png", gene_plot_grid_mes, width = 11, height = 6.5, dpi = 600)
