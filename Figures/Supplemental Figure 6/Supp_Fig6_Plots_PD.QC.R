
# close all, clear all
graphics.off()
rm(list = ls())

# Load libraries
library(Seurat)
library(ggplot2)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(tidyverse)
library(patchwork)

## MEMO: QC Supplement (Supp Fig 6) for PD Starting Population

# Set working directory to be true
setwd("~/Desktop/Datasets/Starting Populations")
# Load population (object): starting pulmonary dissociation (merged, integrated)
load("PD.integrated_03-29-2025.Robj")

# Check structure of data
str(PD.integrated@meta.data)
# Fix Orig_ID
PD.integrated$Orig_ID <- recode(
  PD.integrated$Orig_ID,
  "iBASC_1" = "PD_1",
  "iBASC_2" = "PD_2",
  "iBASC_3" = "PD_3")
# Re-save
setwd("~/Desktop/Datasets/Starting Populations")
save(PD.integrated, file = "PD.integrated_03-29-2025.Robj")

# Color palette for replicates
# Purples for PD
rep.cols = c(
  "PD_1" = "#731963",   
  "PD_2" = "#9A4FB9",   
  "PD_3" = "#B79CCF")

# Get a sense of general QC metrics, split by Orig_ID (replicate)
VlnPlot(PD.integrated, feature = "nFeature_RNA", group.by = "Orig_ID", cols = rep.cols)
VlnPlot(PD.integrated, feature = "nCount_RNA", group.by = "Orig_ID", cols = rep.cols)
VlnPlot(PD.integrated, feature = "percent.mt", group.by = "Orig_ID", cols = rep.cols)

# Extract metadata and reshape
qc_df <- PD.integrated@meta.data %>%
  select(Orig_ID, nFeature_RNA, nCount_RNA, percent.mt) %>%
  pivot_longer(cols = c(nFeature_RNA, nCount_RNA, percent.mt),
               names_to = "Metric", values_to = "Value")
qc_df$Orig_ID <- factor(qc_df$Orig_ID, levels = names(rep.cols))
qc_df$Metric <- factor(qc_df$Metric, levels = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

# Plot as violin
QC.PD_vlns = ggplot(qc_df, aes(x = Orig_ID, y = Value, fill = Orig_ID)) +
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
QC.PD_vlns
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 6")
ggsave("QC.PD_vlns.png", plot = QC.PD_vlns, width = 8, height = 8, dpi = 600)

## Make plot of whole umap embedding, but grouped by sample
DimPlot(PD.integrated, reduction = "umap.rpca", group.by = "Orig_ID", cols = rep.cols, shuffle = T) + ggtitle("")

## Without axes, cleaned up
glbl.embedding_UMAP_sample_PD = DimPlot(
  PD.integrated,
  reduction = "umap.rpca",
  group.by = "Orig_ID",
  cols = rep.cols, shuffle = T) +
  ggtitle("") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(), panel.border = element_blank(), panel.background = element_blank())
print(glbl.embedding_UMAP_sample_PD)
# Set working directory to be true and save
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 6")
ggsave("glbl.embedding_UMAP_PD.St_bysample.png", plot = glbl.embedding_UMAP_sample_PD, width = 7, height = 6, dpi = 600)

# Make QC Metrics split by Sample on embedding
nFeature.PD = FeaturePlot(PD.integrated, features = "nFeature_RNA", reduction = "umap.rpca", split.by = "Orig_ID", cols = c("lightgrey", "#004FFF"), order = T, pt.size = 0.1)
nCount.PD = FeaturePlot(PD.integrated, features = "nCount_RNA", reduction = "umap.rpca", split.by = "Orig_ID", cols = c("lightgrey", "#F75C03"), order = T, pt.size = 0.1)
percent.mt.PD = FeaturePlot(PD.integrated, features = "percent.mt", reduction = "umap.rpca", split.by = "Orig_ID", cols = c("lightgrey", "#FF007F"), pt.size = 0.1)
# Stack plots
PD.QC.embeddings = nFeature.PD / nCount.PD / percent.mt.PD
# Take a look
PD.QC.embeddings
# Save
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 6")
ggsave("PD.QC.embeddings.png", plot = PD.QC.embeddings, width = 9, height = 8, dpi = 600)

## Making these simply to get the legend, but not using for the plots (using plots above)
# nFeature_RNA
p1 <- FeaturePlot(
  PD.integrated,
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
  PD.integrated,
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
  PD.integrated,
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
legend_QC.PD_plot <- p1 / p2 / p3 +
  plot_layout(ncol = 1, heights = c(1, 1, 1))
legend_QC.PD_plot
# Save 
ggsave("Legend_forPD.QC.png",
       legend_QC.PD_plot,
       width = 9, height = 8, dpi = 600)

## Viewing data quality with diagonal trends (From Seurat Vingette)
# This shows us a general quality distribution by replicate; where as total nCounts increase, the nFeatures also increase
# OR at least this is what we want for good quality single cell data
# Red dotted lines mark cells with low gene count (<300) or low total UMI count (<1000).

# Pull out metadata
qc_df <- PD.integrated@meta.data %>%
  mutate(high.mito = percent.mt > 10) %>%
  mutate(Orig_ID = factor(Orig_ID, levels = names(rep.cols)))

# Plot as scatter
qc_plot_scatter_PD <- ggplot(qc_df, aes(x = nCount_RNA, y = nFeature_RNA, color = high.mito)) +
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
# View
qc_plot_scatter_PD
# Save
setwd("~/Desktop/Manuscript Figures/Supplemental Figure 6")
ggsave("qc_plot_scatter_PD.png", plot = qc_plot_scatter_PD, width = 12, height = 4, dpi = 600)
