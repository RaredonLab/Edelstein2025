
# close all, clear all
graphics.off()
rm(list = ls())

library(Seurat)
library(ggplot2)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(tidyverse)
library(patchwork)

## MEMO: QC Supplement (Supp Fig 5) for BAL Starting Population

# Set working directory to be true
setwd("~/Desktop/Datasets/Starting Populations")
# Load the BAL starting population
load("BAL.integrated_03-29-2025.Robj")

# Check structure of data
str(BAL.integrated@meta.data)

# Color palette for replicates
rep.cols = c(
  # Oranges for BAL
  "BAL_1" = "#FF5733",     
  "BAL_2" = "#FFA07A",     
  "BAL_3" = "#FFD1BA")

# Getting a sense of main QC metrics, split by replicate (Orig_ID)
VlnPlot(BAL.integrated, feature = "nFeature_RNA", group.by = "Orig_ID", cols = rep.cols)
VlnPlot(BAL.integrated, feature = "nCount_RNA", group.by = "Orig_ID", cols = rep.cols)
VlnPlot(BAL.integrated, feature = "percent.mt", group.by = "Orig_ID", cols = rep.cols)

# Extract metadata and reshape
qc_df <- BAL.integrated@meta.data %>%
  select(Orig_ID, nFeature_RNA, nCount_RNA, percent.mt) %>%
  pivot_longer(cols = c(nFeature_RNA, nCount_RNA, percent.mt),
               names_to = "Metric", values_to = "Value")
qc_df$Orig_ID <- factor(qc_df$Orig_ID, levels = names(rep.cols))
qc_df$Metric <- factor(qc_df$Metric, levels = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

# Plot (as violins)
QC.BAL_vlns = ggplot(qc_df, aes(x = Orig_ID, y = Value, fill = Orig_ID)) +
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
# View
QC.engineered_vlns
# Set wd and save 
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 5")
ggsave("QC.BAL_vlns.png", plot = QC.BAL_vlns, width = 8, height = 8, dpi = 600)

## Make plot of whole umap embedding, but grouped by sample
DimPlot(BAL.integrated, reduction = "umap.rpca", group.by = "Orig_ID", cols = rep.cols) + ggtitle("")

## Without axes, cleaned up
glbl.embedding_UMAP_sample_BAL = DimPlot(
  BAL.integrated,
  reduction = "umap.rpca",
  group.by = "Orig_ID",
  cols = rep.cols) +
  ggtitle("") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(), panel.border = element_blank(), panel.background = element_blank())
print(glbl.embedding_UMAP_sample_BAL)
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 5")
ggsave("glbl.embedding_UMAP_BAL.St_bysample.png", plot = glbl.embedding_UMAP_sample_BAL, width = 7, height = 6, dpi = 600)

# Make QC Metrics split by Sample on embedding
nFeature.BAL = FeaturePlot(BAL.integrated, features = "nFeature_RNA", reduction = "umap.rpca", split.by = "Orig_ID", cols = c("lightgrey", "#004FFF"), order = T, pt.size = 0.1)
nCount.BAL = FeaturePlot(BAL.integrated, features = "nCount_RNA", reduction = "umap.rpca", split.by = "Orig_ID", cols = c("lightgrey", "#F75C03"), order = T, pt.size = 0.1)
percent.mt.BAL = FeaturePlot(BAL.integrated, features = "percent.mt", reduction = "umap.rpca", split.by = "Orig_ID", cols = c("lightgrey", "#FF007F"), pt.size = 0.1)
# Stack plots
BAL.QC.embeddings = nFeature.BAL / nCount.BAL / percent.mt.BAL
# Take a look
BAL.QC.embeddings
# Save
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 5")
ggsave("BAL.QC.embeddings.png", plot = BAL.QC.embeddings, width = 9, height = 8, dpi = 600)

## Making these simply to get the legend, but not using for the plots (using plots above)
# nFeature_RNA
p1 <- FeaturePlot(
  BAL.integrated,
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
  BAL.integrated,
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
  BAL.integrated,
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
legend_QC.BAL_plot <- p1 / p2 / p3 +
  plot_layout(ncol = 1, heights = c(1, 1, 1))
legend_QC.BAL_plot
# Save 
ggsave("Legend_forBAL.QC.png",
       legend_QC.BAL_plot,
       width = 9, height = 8, dpi = 600)

## Viewing data quality with diagonal trends (From Seurat Vingette)
# This shows us a general quality distribution by replicate; where as total nCounts increase, the nFeatures also increase
# OR at least this is what we want for good quality single cell data
# Red dotted lines mark cells with low gene count (<300) or low total UMI count (<1000).

# Pull out metadata
qc_df <- BAL.integrated@meta.data %>%
  mutate(high.mito = percent.mt > 10) %>%
  mutate(Orig_ID = factor(Orig_ID, levels = names(rep.cols)))

# Plot as scatter
qc_plot_scatter_BAL <- ggplot(qc_df, aes(x = nCount_RNA, y = nFeature_RNA, color = high.mito)) +
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
# View the plot
qc_plot_scatter_BAL
# Save the plot
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 5")
ggsave("qc_plot_scatter_BAL.png", plot = qc_plot_scatter_BAL, width = 14, height = 4, dpi = 600)


