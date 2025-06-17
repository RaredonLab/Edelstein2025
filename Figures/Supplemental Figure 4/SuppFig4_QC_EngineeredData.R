
# close all, clear all
graphics.off()  
rm(list = ls())

# Load packages
library(Seurat)
library(ggplot2)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(tidyverse)
library(patchwork)

## MEMO: QC Figure for Supplement (Supp Fig 4.)

# Set working directory to be true and load the data
setwd("~/Desktop/Datasets/Engineered Global Objects")
load("eng.subset.integrated_NodeAligned.Robj")

# Color palette for replicates
rep.cols = c(
# Blues for (PD_3D)
"PD_3D_1" = "#205975",
"PD_3D_2" = "#62bed9",
"PD_3D_3" = "#c8f4ff",
# Greens for  (Mixed_3D)
"Mixed_3D_1" = "#536e0a",
"Mixed_3D_2" = "#7da50f",
"Mixed_3D_3" = "#bdda0f",
# Reds for BAL3D
"BAL_3D_1" = "#b6042a",
"BAL_3D_2" = "#f50538",
"BAL_3D_3" = "#ff8ca5")

# Vln plots for general info by Orig_ID
VlnPlot(eng.subset.integrated_HK, feature = "nFeature_RNA", group.by = "Orig_ID", cols = rep.cols)
VlnPlot(eng.subset.integrated_HK, feature = "nCount_RNA", group.by = "Orig_ID", cols = rep.cols)
VlnPlot(eng.subset.integrated_HK, feature = "percent.mt", group.by = "Orig_ID", cols = rep.cols)

## Making a nicer looking plot
# Extract metadata and reshaping
qc_df <- eng.subset.integrated_HK@meta.data %>%
  select(Orig_ID, nFeature_RNA, nCount_RNA, percent.mt) %>%
  pivot_longer(cols = c(nFeature_RNA, nCount_RNA, percent.mt),
               names_to = "Metric", values_to = "Value")
qc_df$Orig_ID <- factor(qc_df$Orig_ID, levels = names(rep.cols))
qc_df$Metric <- factor(qc_df$Metric, levels = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

# PLotting the data we just pulled
QC.engineered_vlns = ggplot(qc_df, aes(x = Orig_ID, y = Value, fill = Orig_ID)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black", alpha = 0.6) +
  geom_jitter(aes(color = Orig_ID), size = 0.7, width = 0.15, alpha = 0.6, stroke = 0) +
  scale_fill_manual(values = rep.cols) +
  scale_color_manual(values = rep.cols) +
  facet_wrap(~ Metric, ncol = 1, scales = "free_y", strip.position = "left") +
  theme_minimal(base_size = 14) +
  theme(
    strip.placement = "outside",
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    legend.position = "none",
    panel.grid = element_blank())
QC.engineered_vlns
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 4")
ggsave("QC.engineered_vlns.png", plot = QC.engineered_vlns, width = 8, height = 7, dpi = 600)
# Make plot of whole umap embedding, but grouped by sample
DimPlot(eng.subset.integrated_HK, reduction = "umap.rpca", group.by = "Orig_ID", cols = rep.cols) + ggtitle("")

# Without axes, cleaned up
glbl.embedding_UMAP_sample = DimPlot(
  eng.subset.integrated_HK, reduction = "umap.rpca", group.by = "Orig_ID", cols = rep.cols) +
  ggtitle("") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
    axis.title = element_blank(), panel.grid = element_blank(), panel.border = element_blank(), panel.background = element_blank())
print(glbl.embedding_UMAP_sample)
# Save (WD set to be correct)
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 4")
ggsave("glbl.embedding_UMAP_bysample.png", plot = glbl.embedding_UMAP_sample, width = 7, height = 6, dpi = 600)

# Make QC Metrics split by Sample on embedding
nFeature.eng = FeaturePlot(eng.subset.integrated_HK, features = "nFeature_RNA", reduction = "umap.rpca", split.by = "Orig_ID", cols = c("lightgrey", "#004FFF"), order = T)
nCount.eng = FeaturePlot(eng.subset.integrated_HK, features = "nCount_RNA", reduction = "umap.rpca", split.by = "Orig_ID", cols = c("lightgrey", "#F75C03"), order = T)
percent.mt.eng = FeaturePlot(eng.subset.integrated_HK, features = "percent.mt", reduction = "umap.rpca", split.by = "Orig_ID", cols = c("lightgrey", "#FF007F"))
# Stack plots
eng.QC.embeddings = nFeature.eng / nCount.eng / percent.mt.eng
# Take a look
eng.QC.embeddings
# Save
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 4")
ggsave("Embedding_QC.Engineered.png", plot = eng.QC.embeddings, width = 20, height = 8, dpi = 600)

## Making these simply to get the legend, but not using for the plots (using plots above)
# nFeature_RNA
p1 <- FeaturePlot(
  eng.subset.integrated_HK,
  features = "nFeature_RNA",
  reduction = "umap.rpca",
  split.by = "Orig_ID",
  pt.size = 0.1,
  cols = c("lightgrey", "#004FFF"), order = T) +
  theme_void() +
  theme(
    plot.title = element_blank(),
    legend.position = "right")
# nCount_RNA
p2 <- FeaturePlot(
  eng.subset.integrated_HK,
  features = "nCount_RNA",
  reduction = "umap.rpca",
  split.by = "Orig_ID",
  pt.size = 0.1,
  cols = c("lightgrey", "#F75C03"), order = T) +
  theme_void() +
  theme(
    plot.title = element_blank(),
    legend.position = "right")
# percent.mt
p3 <- FeaturePlot(
  eng.subset.integrated_HK,
  features = "percent.mt",
  reduction = "umap.rpca",
  split.by = "Orig_ID",
  pt.size = 0.1,
  cols = c("lightgrey", "#FF007F"), order = T) +
  theme_void() +
  theme(
    plot.title = element_blank(),
    legend.position = "right")
# Combining 
legend_QC.eng_plot <- p1 / p2 / p3 +
  plot_layout(ncol = 1, heights = c(1, 1, 1))
legend_QC.eng_plot
# Save 
ggsave("Legend_forEng.QC.png",
       legend_QC.eng_plot,
       width = 20, height = 8, dpi = 600)

## Viewing data quality with diagonal trends (From Seurat Vingette)
# This shows us a general quality distribution by replicate; where as total nCounts increase, the nFeatures also increase
# OR at least this is what we want for good quality single cell data
# Red dotted lines mark cells with low gene count (<300) or low total UMI count (<1000).

# Pull out metadata
qc_df <- eng.subset.integrated_HK@meta.data %>%
  mutate(high.mito = percent.mt > 10) %>%
  mutate(Orig_ID = factor(Orig_ID, levels = names(rep.cols)))

# Plot as scatter
qc_plot_scatter <- ggplot(qc_df, aes(x = nCount_RNA, y = nFeature_RNA, color = high.mito)) +
  geom_point(alpha = 0.6, size = 1) +
  facet_wrap(~Orig_ID, nrow = 1) +
  scale_x_log10() +
  scale_color_manual(
    values = c("FALSE" = "#e41a1c", "TRUE" = "#377eb8"),
    name = "percent.mt > 10",  # make title
    labels = c("FALSE" = "No", "TRUE" = "Yes"),
    guide = guide_legend(
      override.aes = list(size = 3))) +
  geom_hline(yintercept = 300, linetype = "dotted", color = "darkred") +
  geom_vline(xintercept = 1000, linetype = "dotted", color = "darkred") +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11, face = "bold"),  # bold legend title
    legend.text = element_text(size = 10)) +
  labs(x = "nCount_RNA", y = "nFeature_RNA")
# Print
qc_plot_scatter
# Save
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 4")
ggsave("qc_plot_scatter.png", plot = qc_plot_scatter, width = 20, height = 6, dpi = 600)


