

## MEMO: Trying different AddModuleScore() approaches for the RAS cells

# close all, clear all
graphics.off()  
rm(list = ls())

library(Seurat)
library(ggplot2)
library(SeuratObject)
library(dplyr)
library(tidyr)

# Set wd
setwd("~/Desktop/Datasets/Engineered Global Objects")
load("eng.subset.integrated_NodeAligned.Robj")

# Check structure
str(eng.subset.integrated_HK@meta.data)
table(eng.subset.integrated_HK$CellType.NodeAligned)

#### Thresholding approach
goi <- c('Sftpc','Scgb1a1','Scgb3a2','Sftpb')
temp <- eng.subset.integrated_HK[goi,]
rownames(temp@assays$RNA@layers$scale.data) <- rownames(temp)
colnames(temp@assays$RNA@layers$scale.data) <- colnames(temp)

# Trying working from the scaled data
scale.data <- temp@assays$RNA@layers$scale.data

# we want to set a threshold for negative vs positive for each individual gene
# Scgb1a1
rownames(scale.data)[1]
test <- data.frame(Scgb1a1 = scale.data[1,])
View(test)
test$cell <- rownames(test)
ggplot(test,
       aes(x = Scgb1a1))+geom_density()+theme_classic()
FeaturePlot(eng.subset.integrated_HK,'Scgb1a1',min.cutoff = 0,split.by = 'Condition',reduction = 'umap.rpca',slot = 'scale.data')
# Scgb1a1 1.5 (threshold)
eng.subset.integrated_HK$Scgb1a1.pos <- NA
thresh.scgb1a1 <- 2
eng.subset.integrated_HK$Scgb1a1.pos <- eng.subset.integrated_HK['Scgb1a1',]@assays$RNA@layers$scale.data > thresh.scgb1a1
View(eng.subset.integrated_HK@meta.data)

# Sftpb
rownames(scale.data)[2]
test <- data.frame(Sftpb = scale.data[2,])
#View(test)
test$cell <- rownames(test)
ggplot(test,
       aes(x = Sftpb))+geom_density()+theme_classic()
FeaturePlot(eng.subset.integrated_HK,'Sftpb',split.by = 'Condition',reduction = 'umap.rpca',slot = 'scale.data',min.cutoff = 1.0, order = T)
# Sftpb 0 (threshold)
eng.subset.integrated_HK$Sftpb.pos <- NA
thresh.sftpb <- 1.0
eng.subset.integrated_HK$Sftpb.pos <- eng.subset.integrated_HK['Sftpb',]@assays$RNA@layers$scale.data > thresh.sftpb
View(eng.subset.integrated_HK@meta.data)

# Sftpc
rownames(scale.data)[3]
test <- data.frame(Sftpc = scale.data[3,])
#View(test)
test$cell <- rownames(test)
ggplot(test,
       aes(x = Sftpc))+geom_density()+theme_classic()
FeaturePlot(eng.subset.integrated_HK,'Sftpc',split.by = 'Condition',reduction = 'umap.rpca',slot = 'scale.data', min.cutoff = -1.5, order = T)
# Sftpc 1 (threshold)
eng.subset.integrated_HK$Sftpc.pos <- NA
thresh.sftpc <- -1.5
eng.subset.integrated_HK$Sftpc.pos <- eng.subset.integrated_HK['Sftpc',]@assays$RNA@layers$scale.data > thresh.sftpc
View(eng.subset.integrated_HK@meta.data)

# Scgb3a2
rownames(scale.data)[4]
test <- data.frame(Scgb3a2 = scale.data[4,])
#View(test)
test$cell <- rownames(test)
ggplot(test,
       aes(x = Scgb3a2))+geom_density()+theme_classic()
FeaturePlot(eng.subset.integrated_HK,'Scgb3a2',split.by = 'Condition',reduction = 'umap.rpca',slot = 'scale.data',min.cutoff = - 2.5, order = T)
# Sftpc 1 (threshold) - 0 okay because even in the proximal part of PD_3D, those shouldn't be double or triple pos
eng.subset.integrated_HK$Scgb3a2.pos <- NA
thresh.scgb3a2 <- -2.5
eng.subset.integrated_HK$Scgb3a2.pos <- eng.subset.integrated_HK['Scgb3a2',]@assays$RNA@layers$scale.data > thresh.scgb3a2
View(eng.subset.integrated_HK@meta.data)

# Conditional logic
eng.subset.integrated_HK$Double.pos <- NA
eng.subset.integrated_HK$Double.pos <- eng.subset.integrated_HK$Scgb1a1.pos & eng.subset.integrated_HK$Sftpb.pos

# Conditional logic for RAS signature
eng.subset.integrated_HK$Triple.pos_RAS <- NA
eng.subset.integrated_HK$Triple.pos_RAS <- eng.subset.integrated_HK$Scgb1a1.pos & eng.subset.integrated_HK$Sftpc.pos & eng.subset.integrated_HK$Scgb3a2.pos

# Conditional logic for AT0 signature
eng.subset.integrated_HK$Triple.pos_AT0 <- NA
eng.subset.integrated_HK$Triple.pos_AT0 <- eng.subset.integrated_HK$Scgb1a1.pos & eng.subset.integrated_HK$Sftpb.pos & eng.subset.integrated_HK$Scgb3a2.pos

# Quad-positive (AT0 + RAS)
eng.subset.integrated_HK$Quad.pos_AT0_RAS <- NA
eng.subset.integrated_HK$Quad.pos_AT0_RAS <- eng.subset.integrated_HK$Scgb1a1.pos & eng.subset.integrated_HK$Sftpb.pos & eng.subset.integrated_HK$Scgb3a2.pos & eng.subset.integrated_HK$Scgb1a1.pos

## Create data-frame and 
# Filter epi
epi_meta <- eng.subset.integrated_HK@meta.data %>%
  filter(CellClass.combined.Integrated == "Epithelium")
# Summarize counts per Orig_ID and Condition
epi_summary <- epi_meta %>%
  group_by(Condition, Orig_ID) %>%
  summarise(
    total_epi = n(),
    n_Triple.pos_RAS = sum(Triple.pos_RAS, na.rm = TRUE),
    n_Triple.pos_AT0 = sum(Triple.pos_AT0, na.rm = TRUE),
    n_Quad.pos_AT0_RAS = sum(Quad.pos_AT0_RAS, na.rm = TRUE)) %>%
  mutate(
    pct_Triple.pos_RAS = 100 * n_Triple.pos_RAS / total_epi,
    pct_Triple.pos_AT0 = 100 * n_Triple.pos_AT0 / total_epi,
    pct_Quad.pos_AT0_RAS = 100 * n_Quad.pos_AT0_RAS / total_epi)
# Now convert to long format
epi_summary_long <- epi_summary %>%
  select(Condition, Orig_ID,
         pct_Triple.pos_RAS, pct_Triple.pos_AT0, pct_Quad.pos_AT0_RAS) %>%
  pivot_longer(
    cols = starts_with("pct_"),
    names_to = "Signature",
    values_to = "Percent_Positive") %>%
  mutate(Signature = recode(Signature,
                            pct_Triple.pos_RAS = "RAS Signature",
                            pct_Triple.pos_AT0 = "AT0 Signature",
                            pct_Quad.pos_AT0_RAS = "AT0 + RAS Signature"))

ggplot(epi_summary_long, aes(x = Condition, y = Percent_Positive, fill = Signature)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(color = Signature), size = 2, shape = 21,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
  labs(
    x = "Condition",
    y = "% Positive of Epithelium",
    title = "Percent Signature-Positive Epithelial Cells per Replicate",
    fill = "Signature",
    color = "Signature") +
  theme_minimal(base_size = 14)

## What about percent positive of all RAS_Like + Secretory?
target_meta <- eng.subset.integrated_HK@meta.data %>%
  filter(CellType.NodeAligned %in% c("RAS_Like", "Secretory"))
type_summary <- target_meta %>%
  group_by(CellType.NodeAligned, Condition, Orig_ID) %>%
  summarise(
    total_cells = n(),
    n_Triple.pos_RAS = sum(Triple.pos_RAS, na.rm = TRUE),
    n_Triple.pos_AT0 = sum(Triple.pos_AT0, na.rm = TRUE),
    n_Quad.pos_AT0_RAS = sum(Quad.pos_AT0_RAS, na.rm = TRUE)) %>%
  mutate(
    pct_Triple.pos_RAS = 100 * n_Triple.pos_RAS / total_cells,
    pct_Triple.pos_AT0 = 100 * n_Triple.pos_AT0 / total_cells,
    pct_Quad.pos_AT0_RAS = 100 * n_Quad.pos_AT0_RAS / total_cells)
# Change format for plotting
type_summary_long <- type_summary %>%
  select(CellType.NodeAligned, Condition, Orig_ID,
         pct_Triple.pos_RAS, pct_Triple.pos_AT0, pct_Quad.pos_AT0_RAS) %>%
  pivot_longer(
    cols = starts_with("pct_"),
    names_to = "Signature",
    values_to = "Percent_Positive") %>%
  mutate(Signature = recode(Signature,
                            pct_Triple.pos_RAS = "RAS Signature",
                            pct_Triple.pos_AT0 = "AT0 Signature",
                            pct_Quad.pos_AT0_RAS = "AT0 + RAS Signature"))
# Changing axes
ggplot(type_summary_long, aes(x = Condition, y = Percent_Positive, fill = Signature)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(color = Signature), size = 2, shape = 21,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
  facet_wrap(~ CellType.NodeAligned, nrow = 1) +
  labs(
    x = "Condition",
    y = "% Positive of Cell Type",
    title = "Percent Signature-Positive Cells per Replicate (RAS_Like & Secretory)",
    fill = "Signature",
    color = "Signature") +
  theme_minimal(base_size = 14)


FeaturePlot(eng.subset.integrated_HK, reduction = "umap.rpca", split.by = "Orig_ID", feature = "Triple.pos_RAS")

######################################## CONTINUATION ######################################################### 
# Checking structure of data for refresher; ensuring scoring is stored
str(eng.subset.integrated_HK@meta.data)
# Pull meta-data
meta <- eng.subset.integrated_HK@meta.data

### Plot 1: % Triple.pos_RAS of all cells in each replicate ###
plot1_df <- meta %>%
  group_by(Orig_ID, Condition) %>%
  summarise(
    Triple_Pos = sum(Triple.pos_RAS, na.rm = TRUE),
    Total_Cells = n(),
    Percent_Tissue = 100 * Triple_Pos / Total_Cells,
    .groups = "drop")

plot1 <- ggplot(plot1_df, aes(x = Condition, y = Percent_Tissue)) +
  geom_jitter(width = 0.2, size = 3, color = "#E76F51") +
  geom_boxplot(aes(group = Condition), alpha = 0.2, fill = "#E76F51", outlier.shape = NA) +
  labs(
    title = "% of All Cells that are Triple.pos_RAS",
    x = "Condition",
    y = "Percent of Tissue") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot1

### Plot 2: % Triple.pos_RAS of Epithelium only (CellClass.NodeAligned) ###
# This is the plot we want****
condition_colors <- c(
  "PD" = "#29339B",
  "BAL" = "#e9c716",
  "PD_3D" = "#4a2377",
  "Mixed_3D" = "#ff585e",
  "BAL_3D" = "#00b0be")
# Load required libraries
library(dplyr)
library(ggplot2)
library(ggbeeswarm)

# Define condition color palette
condition_colors <- c(
  "PD" = "#29339B",
  "BAL" = "#e9c716",
  "PD_3D" = "#4a2377",
  "Mixed_3D" = "#ff585e",
  "BAL_3D" = "#00b0be")

# Subset metadata
meta <- eng.subset.integrated_HK@meta.data
# Calculate % Triple.pos_RAS (Sftpc+, Scgb3a2+, Scgb1a1+) per replicate, normalized by epithelial cells
plot_df_triple.pos.ras <- meta %>%
  filter(CellClass.NodeAligned == "Epithelium") %>%
  group_by(Orig_ID, Condition) %>%
  summarise(
    Triple_Pos = sum(Triple.pos_RAS, na.rm = TRUE),
    Total_Epi = n(),
    Percent_Epi_TriplePos = 100 * Triple_Pos / Total_Epi,
    .groups = "drop")

# Make the plot
plot_ras_epi <- ggplot(plot_df_triple.pos.ras, aes(x = Condition, y = Percent_Epi_TriplePos, fill = Condition)) +
  geom_boxplot(width = 0.5, size = 0.6, outlier.shape = NA, color = "black") +
  geom_quasirandom(
    method = "smiley", width = 0.2, 
    shape = 21, size = 2.5, stroke = 0.3,
    fill = "black", color = "black") +
  scale_fill_manual(values = condition_colors) +
  labs(title = "", x = NULL,
    y = "% Positive for RAS Signature (of Epithelium)") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 13),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none")
# print the plot
print(plot_ras_epi)
setwd("~/Desktop/Manuscript Figures/Figure 5")
ggsave("RAS_Signature_by_Condition_PrctPos.png", plot = plot_ras_epi, width = 5, height = 6, dpi = 600)

# Make the horizontal plot
plot_ras_epi <- ggplot(plot_df_triple.pos.ras, aes(x = Condition, y = Percent_Epi_TriplePos, fill = Condition)) +
  geom_boxplot(width = 0.5, size = 0.6, outlier.shape = NA, color = "black") +
  geom_quasirandom(
    method = "smiley", width = 0.2, 
    shape = 21, size = 2.5, stroke = 0.3,
    fill = "black", color = "black") +
  scale_fill_manual(values = condition_colors) +
  labs(title = "",x = NULL,
    y = "% Positive for RAS Signature (of Epithelium)") +
  coord_flip() +  # Makes the plot horizontal
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none")

# Print the plot
print(plot_ras_epi)
# Save
setwd("~/Desktop/Manuscript Figures/Figure 5")
ggsave("RAS_Signature_by_Condition_PrctPos_horizontal.png", plot = plot_ras_epi, width = 7, height = 4, dpi = 600)

## Run stats
pairwise_result <- pairwise.wilcox.test(
  x = plot_df_triple.pos.ras$Percent_Epi_TriplePos,
  g = plot_df_triple.pos.ras$Condition,
  p.adjust.method = "BH")
# Look at results
print(pairwise_result)

# Kruskal-Wallis; non-parametric (less than 5 samples)
kruskal_result <- kruskal.test(Percent_Epi_TriplePos ~ Condition, data = plot_df_triple.pos.ras)
print(kruskal_result)

### Doing individual comparison statistical tests (t-tests); just checking
# Subset to two conditions
df_ttest <- plot_df_triple.pos.ras %>%
  filter(Condition %in% c("PD_3D", "BAL_3D"))
df_ttest <- plot_df_triple.pos.ras %>%
  filter(Condition %in% c("Mixed_3D", "BAL_3D"))
# Run Welch’s t-test
t.test(Percent_Epi_TriplePos ~ Condition, data = df_ttest, var.equal = FALSE)
# Auto prints results

### Plotting expression for signature genes
# Subset to RAS_Like cells
ras_like_obj <- subset(eng.subset.integrated_HK, subset = CellType.NodeAligned == "RAS_Like")
# Define genes of interest
ras_genes <- c("Sftpc", "Scgb3a2", "Scgb1a1")
# Extract expression + condition
expr <- FetchData(ras_like_obj, vars = c(ras_genes, "Condition"))
expr$Cell <- rownames(expr)

# Convert to long format for proper plotting
expr_long <- expr %>%
  pivot_longer(cols = all_of(ras_genes), names_to = "Gene", values_to = "Expression")
# Plot
RAS_sig_Vlns = ggplot(expr_long, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = TRUE, color = "black") +
  geom_jitter(width = 0.2, size = 0.8, color = "black") +  # Larger, solid black points
  scale_fill_manual(values = condition_colors) +
  facet_grid(rows = vars(Gene), switch = "y", scales = "free_y") +
  theme_classic(base_size = 14) +
  theme(
    strip.placement = "outside",
    strip.position = "left",  # Label on left
    strip.text.y.left = element_text(angle = 90, size = 13, face = "bold"),  # Bold, vertical
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    strip.background = element_blank(),
    legend.position = "none") +
  labs(y = "Normalized Expression Level")
setwd("~/Desktop/Manuscript Figures/Figure 5")
ggsave("RAS_sig_Vlns.png", plot = RAS_sig_Vlns, width = 8, height = 7, dpi = 600)

