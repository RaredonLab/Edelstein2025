
### MEMO: Supplemental Data, (Figure 2) PD Starting Population 

# close all, clear all
graphics.off()  
rm(list = ls())

# Load libraries
library(Seurat)
library(ggplot2)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)

# Set working directory to be true and load the object
setwd("~/Desktop/Datasets/Starting Populations")
load("PD.integrated_03-29-2025.Robj")

# Take a look at structure of data/meta-data
str(PD.integrated@meta.data)

# Define colors for CellClass
CellClass.cols <- c(
  "Epithelium" = "#D982C6",
  "Immune" = "#87B37A",
  "Mesenchyme" = "#F4A261",
  "Endothelium" = "#2A9D8F")

# Define CellType color palette
CellType.cols <- c(
  "ATI" = "#b22222",
  "ATII" = "#4682B4",
  "Secretory" = "#2E8B57",
  "Ciliated" = "#9B489B",
  "Basal" = "#c6e308",
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
  "Cell_Cycle" = "#e0bf1b","Actc1_Mural" = "#1179fa") 

# Look at the cell type meta-data
table(PD.integrated$CellType)
# Get a sense of the integrated embedding 
DimPlot(PD.integrated, reduction = "umap.rpca", cols = CellType.cols, group.by = "CellType")

## UMAP by Lineage
PD.umap_lineage <- DimPlot(PD.integrated, reduction = "umap.rpca", group.by = "CellClass_Merged", 
                           cols = CellClass.cols) +
  theme_void() +
  theme(legend.position = "none", plot.title = element_blank(), plot.background = element_blank(),
    panel.background = element_blank(), panel.grid = element_blank(),
    axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
# Print and take a look
PD.umap_lineage
# Save to file (set WD first)
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 2")
ggsave("PD.umap_lineage.png", plot = PD.umap_lineage, width = 5, height = 4, dpi = 600)

## UMAP by CellType
PD.umap_celltype = DimPlot(
  PD.integrated, reduction = "umap.rpca",
  group.by = "CellType", cols = CellType.cols, label = T) +
  theme_void() +
  theme(legend.position = "none", plot.title = element_blank(),
    plot.background = element_blank(), panel.background = element_blank(),
    panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank())
# Take a look at it
PD.umap_celltype
# Before 
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 2")
ggsave("PD.umap_celltype.png", plot = PD.umap_celltype, width = 8.5, height = 7, dpi = 600)

## UMAP by Cell Cycle Score (Phase)
# Define phase colors (from our project color palette)
phase_colors <- c(
  "S" = "#16CEA6",
  "G1" = "#D2691E",
  "G2M" = "#AEC5EB")

# Create UMAP, grouped by PHase
PD.umap_phase <- DimPlot(PD.integrated, reduction = "umap.rpca", group.by = "Phase", cols = phase_colors) +
  theme_void() +
  theme(legend.position = "none", plot.title = element_blank(), plot.background = element_blank(),
    panel.background = element_blank(), panel.grid = element_blank(), axis.text = element_blank(),
    axis.ticks = element_blank(), axis.title = element_blank())
PD.umap_phase
# Save to file
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 2")
ggsave("PD.umap_phase.png", plot = PD.umap_phase, width = 5, height = 4, dpi = 600)

## Now we want to make a custom legend (Separate from plots)
# Already defined cell class cols above

# Create dummy data for color swatches
cellclass_df <- data.frame(Label = factor(names(CellClass.cols), levels = names(CellClass.cols)))
# Generate dummy plot for legend
cellclass_legend_plot <- ggplot(cellclass_df, aes(x = 1, y = Label, color = Label)) +
  geom_point(size = 8) +  # larger dummy points
  scale_color_manual(values = CellClass.cols) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 7)))  # larger dots in legend
# Extract and display legend
cellclass_legend <- cowplot::get_legend(cellclass_legend_plot)
cowplot::plot_grid(cellclass_legend)
# Save it
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 2")
ggsave("legend_CellClass.png", plot = cowplot::plot_grid(cellclass_legend), width = 2.5, height = 2, dpi = 600)

# Now for cell type
# Define CellType colors
CellType.cols <- c("ATI" = "#b22222", "ATII" = "#4682B4","Basal" = "#c6e308",
  "Ciliated" = "#9B489B", "Secretory" = "#2E8B57","Tuft" = "#15DDB5",
  "B" = "#41571b","Mac_Alv" = "#00FBFF", "Mac_Inter" = "#FF481B","Monocytes" = "#3cd500",
  "Neutrophils" = "#fcfd1d","NK" = "#F51BEA","pDCs" = "#E9165D","T" = "#95B8D1"
  ,"Arterial" = "#7E5109","gCaps" = "#A0522D", "Venous" = "#660708",
  "Actc1_Mural" = "#1179fa","Mesothelium" = "#4A314D",
  "Cell_Cycle" = "#e0bf1b")

# Dummy data
celltype_df <- data.frame(Label = factor(names(CellType.cols), levels = names(CellType.cols)))
# Legend plot
celltype_legend_plot <- ggplot(celltype_df, aes(x = 1, y = Label, color = Label)) +
  geom_point(size = 8) +
  scale_color_manual(values = CellType.cols) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 7)))

# Extract and save
celltype_legend <- cowplot::get_legend(celltype_legend_plot)
cowplot::plot_grid(celltype_legend)
ggsave("legend_CellType.png", plot = cowplot::plot_grid(celltype_legend), width = 3.5, height = 7, dpi = 600)

# Phase
# Define Phase colors
phase_colors <- c(
  "S" = "#16CEA6",
  "G1" = "#D2691E",
  "G2M" = "#AEC5EB")

# Dummy data
phase_df <- data.frame(Label = factor(names(phase_colors), levels = names(phase_colors)))
# Legend plot
phase_legend_plot <- ggplot(phase_df, aes(x = 1, y = Label, color = Label)) +
  geom_point(size = 8) +
  scale_color_manual(values = phase_colors) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 7)))
# Extract and save
phase_legend <- cowplot::get_legend(phase_legend_plot)
cowplot::plot_grid(phase_legend)
ggsave("legend_Phase.png", plot = cowplot::plot_grid(phase_legend), width = 2.5, height = 2, dpi = 600)

### Dot Plot for expression

# Run quick marker list just to double check (I know for the most part what each consists of but...)
Idents(PD.integrated) = PD.integrated$CellType
PD.mark <- FindAllMarkers(
  PD.integrated,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.1,
  assay = "RNA")
View(PD.mark)

PD.integrated$CellType <- factor(
  PD.integrated$CellType,
  levels = c( "ATI","ATII","Basal","Ciliated","Secretory","Tuft",
               "B", "Mac_Alv","Mac_Inter","Monocytes", "Neutrophils","NK","pDCs", "T", 
               "Arterial","gCaps","Venous",
               "Actc1_Mural","Mesothelium",
               "Cell_Cycle"))

genes.to.plot <- c("Pdpn","Ager","Rtkn2","Akap5", #ATI
                   "Napsa","Lamp3","S100g","Defb4", # ATII
                   "Krt5","Tp63","Krt14","Wnt10a", # Basal
                   "Ccdc153","Pifo","Dnah6","Ak9", # Ciliated
                   "Chia","Agr2","Lypd2","Muc5b", # Secretory
                   "Pou2f3","Dclk1","Trpm5","Avil", # Tuft
                   "Cd79b", "Bank1","Pax5","Cr2", # B cells
                   "Prodh2","Spic","Krt79","Abcd2", # Mac_Alv
                   "Cd163","C1qc","C1qa","Folr2", # Mac_inter
                   "Mal","Eno3","Spn","Tgfbi", # Mono
                   "S100a9","S100a8","Ifit3","Stfa2", # Neutro
                   "Gzma","Gzmk","Il2rb","Klrg1", # NK
                   "Siglech","Kmo","Jaml","Blnk", # pDC
                   "Cd3g", "Cd3e", "Cd3d","Skap1", # T
                   "Gja5","Plat","Ecm1","Sox17", # Arterial
                   "Acer2","Wif1","Apln","Nostrin", # gCaps
                   "Amigo2","Jam2","Slc6a2","Csrp2", # Venous
                   "Actc1","Des","Pdgfrb","Col1a1", # Mural
                   "Msln","Rspo1","Aldh1a2","Upk3b", # Meso
                   "Mki67","Top2a","Ube2c","Ccna2") # Cell cycle

## Dot plot of top markers by cell type
# With axis lines
PD.dotplot = DotPlot(
  PD.integrated, features = genes.to.plot, group.by = "CellType") +
  scale_color_gradientn(colors = c("white", "#2d00f7"), name = "Avg. Expression") +
  scale_size(range = c(1, 5), name = "% Expressing") +
  theme_minimal(base_size = 14) +
  theme(axis.line.x = element_line(color = "black", size = 0.3), # Axis lines
    axis.line.y = element_line(color = "black", size = 0.3),
    axis.text.x = element_text(size = 12, face = "plain", color = "black", angle = 45, hjust = 1), # Axis tick labels (plain, not bold)
    axis.text.y = element_text(size = 12, face = "plain", color = "black"),
    axis.title.x = element_text(size = 14, face = "bold", color = "black", margin = margin(t = 10)), # Axis titles (bold, black)
    axis.title.y = element_text(size = 14, face = "bold", color = "black", margin = margin(r = 10)),
    # Remove grid and background
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank(),
    # Legend at right 
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold", color = "black"),
    legend.text = element_text(size = 12, color = "black")) +
  labs(x = "Gene", y = "Cell Type")
# Print to make sure we like it
PD.dotplot
# Set wd to correct location and save
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 2")
ggsave("PD.dotplot.celltype.png", plot = PD.dotplot, width = 20, height = 8, dpi = 600)

### Feature Plots of interest
# Making select feature plots based on text references
FeaturePlot(PD.integrated, reduction = "umap.rpca", features = c("Ager","Lamp3","Krt5","Ccdc153","Agr2",
                                                                 "Pou2f3","Trpm5"), ncol = 7)
## Epithelium FIRST
# Define genes to plot
features_to_plot <- c("Ager","Lamp3","Krt5","Ccdc153","Agr2",
                      "Pou2f3","Trpm5")
# Function to create individual FeaturePlot
feature_umap <- function(gene) {
  FeaturePlot(
    PD.integrated,
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
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
}

# Generate and arrange the plots
plots_list <- lapply(features_to_plot, feature_umap)
gene_plot_grid_epi <- wrap_plots(plots_list, ncol = 7)
gene_plot_grid_epi
# Save figure (with correct WD)
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 2")
ggsave("gene_plot_grid_epi.png", gene_plot_grid_epi, width = 18, height = 3.5, dpi = 600)

## IMMUNE
# Define genes to plot
features_to_plot <- c("Cd79b","Prodh2","C1qc","Mal","S100a8",
                      "Jaml","Cd3e")
# Function to create individual FeaturePlot
feature_umap <- function(gene) {
  FeaturePlot(
    PD.integrated,
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
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
}

# Generate and arrange the plots
plots_list <- lapply(features_to_plot, feature_umap)
gene_plot_grid_imm <- wrap_plots(plots_list, ncol = 7)
gene_plot_grid_imm
# Save figure (set WD to be true first)
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 2")
ggsave("gene_plot_grid_imm.png", gene_plot_grid_imm, width = 18, height = 3.5, dpi = 600)

## Endothelium
# Define genes to plot
features_to_plot <- c("Gja5","Wif1","Slc6a2")
# Function to create individual FeaturePlot
feature_umap <- function(gene) {
  FeaturePlot(
    PD.integrated,
    features = gene,
    reduction = "umap.rpca",
    pt.size = 0.6) +
    ggtitle(gene) +
    theme_void() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 10),         
      legend.key.height = unit(0.65, "cm"),           
      legend.key.width = unit(0.4, "cm"),          
      legend.margin = margin(t = 2, b = 2, unit = "mm"),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
}

# Generate and arrange the plots
plots_list <- lapply(features_to_plot, feature_umap)
gene_plot_grid_endo <- wrap_plots(plots_list, ncol = 3)
gene_plot_grid_endo
# Save figure
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 2")
ggsave("gene_plot_grid_endo.png", gene_plot_grid_endo, width = 9, height = 3, dpi = 600)

## Mesenchyme
# Define genes to plot
features_to_plot <- c("Col1a1","Rspo1","Pdgfrb")
# Function to create individual FeaturePlot
feature_umap <- function(gene) {
  FeaturePlot(
    PD.integrated,
    features = gene,
    reduction = "umap.rpca",
    pt.size = 0.6) +
    ggtitle(gene) +
    theme_void() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 10),          
      legend.key.height = unit(0.65, "cm"),           
      legend.key.width = unit(0.4, "cm"),            
      legend.margin = margin(t = 2, b = 2, unit = "mm"),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
}

# Generate and arrange the plots
plots_list <- lapply(features_to_plot, feature_umap)
gene_plot_grid_mes <- wrap_plots(plots_list, ncol = 3)
# Print plot; all checks out
gene_plot_grid_mes
# Save figure (with correct WD)
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 2")
ggsave("gene_plot_grid_mes.png", gene_plot_grid_mes, width = 9, height = 3, dpi = 600)



