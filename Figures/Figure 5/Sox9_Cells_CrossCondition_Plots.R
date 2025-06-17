
# close all, clear all
graphics.off()  
rm(list = ls())

# Load packages
library(dplyr)
library(viridis)
library(patchwork)
library(Seurat)
library(cowplot)

# Set working directory to be true 
setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets/SE Legacy Objects")
# Load starting population Seurat objects
load("PD.integrated_03-29-2025.Robj")
load("BAL.integrated_03-29-2025.Robj")

# Load global object (have to change WD) (organoids)
setwd("~/Desktop/Datasets/Engineered Global Objects")
load("eng.subset.integrated_NodeAligned.Robj")

# Inspect meta-data structure
str(PD.integrated@meta.data)
str(BAL.integrated@meta.data)
str(eng.subset.integrated_HK@meta.data)

# Fix meta-data for PD starting "Condition" label; standarddize condition
PD.integrated$Condition <- "PD"

# Set cell identities and visualize Sox9 expression across datasets
Idents(PD.integrated) = PD.integrated$Condition
VlnPlot(PD.integrated, features = "Sox9")
# Do the same for BAL
Idents(BAL.integrated) = BAL.integrated$Condition
VlnPlot(BAL.integrated, features = "Sox9")
# Do the same for our engineered samples
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$Condition
VlnPlot(eng.subset.integrated_HK, features = "Sox9")

# Create binary column indicating Sox9+ cells (threshold = 0.20 expression)
eng.subset.integrated_HK$SOX9_Positive <- FetchData(eng.subset.integrated_HK, vars = "Sox9") > 0.20
PD.integrated$SOX9_Positive <- FetchData(PD.integrated, vars = "Sox9") > 0.20
BAL.integrated$SOX9_Positive <- FetchData(BAL.integrated, vars = "Sox9") > 0.20

## Calculate % Sox9+ cells per sample
# Group and summarize by sample (Orig_ID)
PD_summary <- PD.integrated@meta.data %>%
  group_by(Orig_ID) %>%
  summarize(
    percent_SOX9 = 100 * sum(SOX9_Positive, na.rm = TRUE) / n(),
    Condition = unique(Condition))
BAL_summary <- BAL.integrated@meta.data %>%
  group_by(Orig_ID) %>%
  summarize(
    percent_SOX9 = 100 * sum(SOX9_Positive, na.rm = TRUE) / n(),
    Condition = unique(Condition))
ENG_summary <- eng.subset.integrated_HK@meta.data %>%
  group_by(Orig_ID, Condition) %>%
  summarize(
    percent_SOX9 = 100 * sum(SOX9_Positive, na.rm = TRUE) / n())
# Combine all summaries
combined_sox9_summary <- bind_rows(PD_summary, BAL_summary, ENG_summary)

# Reorder condition factor levels for plotting
combined_sox9_summary$Condition <- factor(
  combined_sox9_summary$Condition,
  levels = c("PD", "PD_3D", "Mixed_3D", "BAL", "BAL_3D"))

### Plot Boxplot of % Sox9+ cells
# Custom color palette
condition_colors <- c(
  "PD" = "#29339B",
  "BAL" = "#e9c716",
  "PD_3D" = "#4a2377",
  "Mixed_3D" = "#ff585e",
  "BAL_3D" = "#00b0be")

# Generate boxplot with jittered points
sox9_pos = ggplot(combined_sox9_summary, aes(x = Condition, y = percent_SOX9, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 1.0) +
  geom_jitter(width = 0.2, shape = 21, color = "black", size = 2) +
  scale_fill_manual(values = condition_colors) +
  labs(
    title = "",
    x = NULL,
    y = "% Sox9+ Cells") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, NA))
print(sox9_pos)
# Save plot
setwd("~/Desktop/Manuscript Figures/Figure 5")
ggsave("sox9_pos_box.png", plot = sox9_pos, width = 7, height = 4, dpi = 600)

## Proportion test for Sox9+ cell frequencies
# Prepare raw counts by sample
# Select only the columns needed to avoid type mismatches
pd_meta <- PD.integrated@meta.data %>%
  select(Orig_ID, Condition, SOX9_Positive) %>%
  mutate(dataset = "PD")
bal_meta <- BAL.integrated@meta.data %>%
  select(Orig_ID, Condition, SOX9_Positive) %>%
  mutate(dataset = "BAL")
eng_meta <- eng.subset.integrated_HK@meta.data %>%
  select(Orig_ID, Condition, SOX9_Positive) %>%
  mutate(dataset = "3D")

# Combine metadata and summarize counts
raw_counts <- bind_rows(pd_meta, bal_meta, eng_meta) %>%
  group_by(Orig_ID, Condition) %>%
  summarize(
    total_cells = n(),
    sox9_pos = sum(SOX9_Positive, na.rm = TRUE),
    .groups = "drop")

# Summarize total Sox9+ and total cells by condition
grouped_counts <- raw_counts %>%
  group_by(Condition) %>%
  summarize(
    total_pos = sum(sox9_pos),
    total_n = sum(total_cells),
    .groups = "drop")

# Pairwise comparisons using prop.test
condition_pairs <- combn(unique(grouped_counts$Condition), 2, simplify = FALSE)
pairwise_tests <- lapply(condition_pairs, function(pair) {
  group1 <- grouped_counts %>% filter(Condition == pair[1])
  group2 <- grouped_counts %>% filter(Condition == pair[2])
  test_result <- prop.test(
    x = c(group1$total_pos, group2$total_pos),
    n = c(group1$total_n, group2$total_n))
  tibble(
    group1 = pair[1],
    group2 = pair[2],
    p.value = test_result$p.value
  )
})

# Combine and adjust
pairwise_results <- bind_rows(pairwise_tests) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))

# View statistical results
pairwise_results

########## Making Plots ########

# Make dot plot for stressed_progenitor and RAS signature 
ras_genes <- c("Scgb1a1", "Scgb3a2", "Sftpc", "Muc5b", "Muc5ac","Lamp3","Napsa")
stress_genes <- c("Sox9", "Cox4i2", "Ptges", "Rasd2", "Stc1")
all_markers <- c(ras_genes, stress_genes)

# Set identities and plot for each 3D condition
Idents(PD_3D.integrated_HK) <- "CellType"
Idents(Mixed_3D.integrated_HK) <- "CellType"
Idents(BAL_3D.integrated_HK) <- "CellType"

# Generate plots by condition
pd3d = DotPlot(
  PD_3D.integrated_HK,
  features = all_markers,
  idents = c("RAS_Like", "Stressed_Progenitor", "Secretory", "ATII_Like", "ATI_Like")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distinct RAS and Stressed Progenitor Profiles in PD_3D",
    x = "Gene", y = "Cell Type")
# Print
pd3d

mixed3d = DotPlot(
  Mixed_3D.integrated_HK,
  features = all_markers,
  idents = c("RAS_Like", "Stressed_Progenitor", "Secretory", "ATII_Like", "ATI_Like")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distinct RAS and Stressed Progenitor Profiles in PD_3D",
    x = "Gene", y = "Cell Type")
# Print
mixed3d

bal3d = DotPlot(
  BAL_3D.integrated_HK,
  features = all_markers,
  idents = c("RAS_Like", "Stressed_Progenitor", "Secretory", "ATII_Like", "ATI_Like")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distinct RAS and Stressed Progenitor Profiles in PD_3D",
    x = "Gene", y = "Cell Type")
# Print
bal3d

## FeaturePlots of RAS Signature genes
# Set identity and plot RAS markers
Idents(BAL_3D.integrated_HK) = BAL_3D.integrated_HK$Condition
# RAS Signature Feature Plots
Idents(BAL_3D.integrated_HK) = BAL_3D.integrated_HK$CellType
p1 <- FeaturePlot(
  BAL_3D.integrated_HK,
  features = "Scgb3a2",
  reduction = "umap.rpca",
  order = TRUE, label = T, repel = T) + scale_color_viridis_c(option = "C")
p2 <- FeaturePlot(
  BAL_3D.integrated_HK,
  features = "Scgb1a1",
  reduction = "umap.rpca",
  order = TRUE, label = T, repel = T) + scale_color_viridis_c(option = "C")
p3 <- FeaturePlot(
  BAL_3D.integrated_HK,
  features = "Sftpc",
  reduction = "umap.rpca",
  order = TRUE, label = T, repel = T) + scale_color_viridis_c(option = "C")
RAS_signature = p1 | p2 | p3
setwd("~/Desktop/Manuscript Figures/Figure 5")
ggsave("RAS_signature.png", plot = RAS_signature, width = 8, height = 4, dpi = 600)

## Dimensional plots to highlight signatures
# Highlight specific epithelial cell types
highlighted <- c("Stressed_Progenitor", "RAS_Like", "Secretory")
# Set color scheme
highlight_colors <- c(
  "Stressed_Progenitor" = "#c71585",
  "RAS_Like" = "#595db0",
  "Secretory" = "#2E8B57")
# Create a new column for plotting
eng.subset.integrated_HK$plot_highlight <- ifelse(
  eng.subset.integrated_HK$CellType.combined.Integrated %in% highlighted,
  eng.subset.integrated_HK$CellType.combined.Integrated,
  "Other")

# # Define colors and plot
celltype_colors <- c(highlight_colors, "Other" = "lightgray")
# Final plot
fig5_signatures = DimPlot(
  eng.subset.integrated_HK,
  group.by = "plot_highlight",
  cols = celltype_colors,
  reduction = "umap.rpca") +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_blank(),
    legend.text = element_text(size = 10), legend.key = element_blank(),
    plot.title = element_blank()) +
  scale_color_manual(values = celltype_colors, breaks = names(highlight_colors))
print(fig5_signatures)
# Save plot
setwd("~/Desktop/Manuscript Figures/Figure 5")
ggsave("fig5_signatures.png", plot = fig5_signatures, width = 5, height = 5, dpi = 600)
