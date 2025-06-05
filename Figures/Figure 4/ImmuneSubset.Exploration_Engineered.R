
# MEMO: Exploration for "M1" v "M2" Mac Polarization

# Packages
library(grid)  # For unit()
library(ggplot2)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(ggalluvial)

# Subset out immune
# Check the metadata for available cell class annotations
table(eng.subset.integrated_HK$CellClass)
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$seurat_clusters
DimPlot(eng.subset.integrated_HK, reduction = "umap", label = T)
DimPlot(eng.subset.integrated_HK, reduction = "umap", label = T, split.by = "Orig_ID")
FeaturePlot(eng.subset.integrated_HK, features = c("Epcam","Ptprc","Col1a1","Cdh5"), reduction = "umap")

# Look at object
eng.subset.integrated_HK

# Define the immune cell classes (modify based on your dataset)
immune_classes <- c("Immune")  
# Subset immune cells
immune_subset <- subset(eng.subset.integrated_HK, subset = CellClass.combined.Integrated %in% immune_classes)
# Print confirmation
print(immune_subset)
# Rescale the data
immune_subset = ScaleData(immune_subset)
# Normalizing the data
immune_subset = NormalizeData(immune_subset)
# Finding variable features
immune_subset = FindVariableFeatures(immune_subset)

# PCA
immune_subset = RunPCA(immune_subset, npcs = 100)
ElbowPlot(immune_subset,ndims = 100)
PCHeatmap(immune_subset,cells=200,balanced=T,dims=1:9)
PCHeatmap(immune_subset,cells=200,balanced=T,dims=10:18)
PCHeatmap(immune_subset,cells=200,balanced=T,dims=19:27) # handpick

# Cluster + embed
immune_subset = RunUMAP(immune_subset, dims = c(1:3, 5:7, 9))
DimPlot(immune_subset, label = TRUE)
# Creating nearest neighbor graph, not clustering.
immune_subset = FindNeighbors(immune_subset, dims = c(1:3, 5:7, 9))
# Defining clusters - very high res. There is purpose in doing this.
immune_subset = FindClusters(immune_subset, res = 0.1)
DimPlot(immune_subset, label = TRUE, split.by = "Orig_ID", cols = mac.cols)

## Get rid of 3 
immune_subset.clean <- subset(immune_subset, subset = seurat_clusters != 3)
# Rescale the data
immune_subset.clean = ScaleData(immune_subset.clean)
# Normalizing the data
immune_subset.clean = NormalizeData(immune_subset.clean)
# Finding variable features
immune_subset.clean = FindVariableFeatures(immune_subset.clean)

# PCA
immune_subset.clean = RunPCA(immune_subset.clean, npcs = 100)
ElbowPlot(immune_subset.clean,ndims = 100)
PCHeatmap(immune_subset.clean,cells=200,balanced=T,dims=1:9)
PCHeatmap(immune_subset.clean,cells=200,balanced=T,dims=10:18)
PCHeatmap(immune_subset.clean,cells=200,balanced=T,dims=19:27) # Liking 19 PCs

immune_subset.clean = RunUMAP(immune_subset.clean, dims = c(1:8))
DimPlot(immune_subset.clean, label = TRUE)
# Creating nearest neighbor graph, not clustering.
immune_subset.clean = FindNeighbors(immune_subset.clean, dims = c(1:8))
# Defining clusters - very high res. There is purpose in doing this.
immune_subset.clean = FindClusters(immune_subset.clean, res = 0.05)
DimPlot(immune_subset.clean, label = TRUE)
DimPlot(immune_subset.clean, label = TRUE, split.by = "Condition")

# BAL_3D is M2 - Look at FindAllMarkers test output first
Idents(immune_subset.clean) = immune_subset.clean$seurat_clusters
mac.pol_DGE = FindAllMarkers(immune_subset.clean, 
                                        min.pct = 0.1,logfc.threshold = 0.1)
mac.pol_DGE$ratio = mac.pol_DGE$pct.1/mac.pol_DGE$pct.2
mac.pol_DGE$power = mac.pol_DGE$ratio*mac.pol_DGE$avg_log2FC
View(mac.pol_DGE)

# Run a marker test with cluster 0 vs. 1 and 2
Idents(immune_subset.clean) <- immune_subset.clean$seurat_clusters
mac.pol_DGE <- FindMarkers(immune_subset.clean, ident.1 = "0", ident.2 = c("1","2"),
                           min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE, test.use = "wilcox")
mac.pol_DGE$ratio <- mac.pol_DGE$pct.1 / mac.pol_DGE$pct.2
mac.pol_DGE$power <- mac.pol_DGE$ratio * mac.pol_DGE$avg_log2FC
View(mac.pol_DGE)

# Icam1/Itgam/Timp1/Socs1/Il1b/Il4r/Lcn2/Ccl12/Junb/Vegfa/Myc/Mmp9/Il1a/Ccnd1/Birc5/Socs3/Osm/Nos2/Fn1
FeaturePlot(immune_subset.clean, features = c("Icam1","Itgam","Il1b","Il1a","Nos2","Fn1","Osm"), order = T)
FeaturePlot(immune_subset.clean, features = c("Stat1","Jak2","Jak1","Ifng"), order = T)
FeaturePlot(immune_subset.clean, features = c("Il1a","Il1b","Tlr4","Pparg","Irf5","Tlr2","Stat6","Stat3","Klf4"), order = T)

# Things to promote my M1/M2 like ideas; they're all Jak1/Jak2 positive
# Clusters 1 + 2 = M1-like (Nos2)
# M1 shown to be stimulated by IFN-gamma pathway (look for genes in this pathway upregulated then)

FeaturePlot(immune_subset.clean, features = c("Cadm1","Pparg"), order = T, group.by = "CellType.sub")
FeaturePlot(immune_subset.clean, features = c("Cd14","Cxcl1","Vcan","Cxcl3","Trem3","Arg1"))
FeaturePlot(immune_subset.clean, features = c("Cadm1","Rcn3","Atp1b1","Prrt1","Gpr155","Prrt1"))
FeaturePlot(immune_subset.clean, features = c("Nos2","Il1b","Arg1","Il4","Tgfb1","Tgfbi"))

# M1 Profile (mixed): Cd14, Cxcl1, Vcan+ (CALLING CXL1^high/Vcan^high)
# M2 profile (bal): Cadm1, Rcn3 (CALLING Atp1b1^High/Cadm1^high)
# Create empty column first
immune_subset.clean$CellType.sub <- NA
# Subset cells by identity
temp0 <- WhichCells(immune_subset.clean, idents = "0")
temp1 <- WhichCells(immune_subset.clean, idents = "1")
temp2 <- WhichCells(immune_subset.clean, idents = "2")
# Assign cell type labels
immune_subset.clean$CellType.sub[temp0] <- "Anti_Inflamm_Mac"
immune_subset.clean$CellType.sub[temp1] <- "Pro_Inflamm_Mac"
immune_subset.clean$CellType.sub[temp2] <- "Pro_Inflamm_Mac"

str(immune_subset.clean@meta.data)
# Clean up a bit before saving
cols_to_remove <- c("RNA_snn_res.0.32", "RNA_snn_res.0.1", "RNA_snn_res.0.05", "RNA_snn_res.0.025")
# Drop them from the metadata
immune_subset.clean@meta.data <- immune_subset.clean@meta.data[, !(colnames(immune_subset.clean@meta.data) %in% cols_to_remove)]
# Save to server
setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets")
save(immune_subset.clean, file = "immune_subset.clean.Robj")

# Now plot
mac.cols = c("Pro_Inflamm_Mac" = "#9497fd",
"Anti_Inflamm_Mac" = "#ff9e80")
DimPlot(immune_subset.clean, group.by = "CellType.sub", cols = mac.cols)
DimPlot(immune_subset.clean, group.by = "Phase")
Idents(immune_subset.clean) = immune_subset.clean$Phase

# Subset plot for figure 4
mac.cols <- c(
  "Pro_Inflamm_Mac" = "#9497fd",
  "Anti_Inflamm_Mac" = "#ff9e80")
mac.subset_UMAP = DimPlot(immune_subset.clean, group.by = "CellType.sub", cols = mac.cols) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14),  # Adjust size as needed
    plot.title = element_blank())
setwd("~/Desktop/Manuscript Figures/Figure 4")
ggsave("mac.subset_UMAP.png", plot = mac.subset_UMAP, width = 5, height = 5, dpi = 600)
mac.subset_UMAP

## Now we want to plot the distribution of each type across conditions
# Compute proportion of cell types per condition
cell_counts_condition = immune_subset.clean@meta.data %>%
  filter(CellType.sub %in% c("Pro_Inflamm_Mac", "Anti_Inflamm_Mac")) %>%
  group_by(Condition, CellType.sub) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(proportion = n / sum(n))

# Plot as horizontal stacked bar
stackedbar_mac_pol_condition <- ggplot(cell_counts_condition, aes(x = Condition, y = proportion, fill = CellType.sub)) +
  geom_col(width = 0.7, color = "black") +
  coord_flip(clip = "off") +  # allowing labels to extend past plot panel
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c(
    "Pro_Inflamm_Mac" = "#9497fd",
    "Anti_Inflamm_Mac" = "#ff9e80")) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(color = "black"),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.position = "bottom",
    plot.margin = unit(c(1, 1.5, 1, 1), "lines")  # top, right, bottom, left prevents cutoff
  ) +
  labs(y = "Contribution to Macrophage Population",
    fill = "Polarized State")
print(stackedbar_mac_pol_condition)
setwd("~/Desktop/Manuscript Figures/Figure 4")
ggsave("stackedbar_mac_pol_condition.png", plot = stackedbar_mac_pol_condition, width = 7, height = 5, dpi = 600)

# Split by replicate (Orig_ID)
cell_counts_rep <- immune_subset.clean@meta.data %>%
  filter(CellType.sub %in% c("Pro_Inflamm_Mac", "Anti_Inflamm_Mac")) %>%
  group_by(Condition, CellType.sub, Orig_ID) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Orig_ID) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()


cell_counts_rep$Orig_ID <- factor(cell_counts_rep$Orig_ID, levels = c(
  "PD_3D_1", "PD_3D_3",
  "Mixed_3D_1", "Mixed_3D_2", "Mixed_3D_3",
  "BAL_3D_1", "BAL_3D_2", "BAL_3D_3"))

stackedbar_mac_pol <- ggplot(cell_counts, aes(x = Orig_ID, y = proportion, fill = CellType.sub)) +
  geom_col(width = 0.7, color = "black") +
  coord_flip(clip = "off") +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(
    "Pro_Inflamm_Mac" = "#9497fd",
    "Anti_Inflamm_Mac" = "#ff9e80")) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(color = "black"),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.position = "bottom",
    plot.margin = unit(c(1, 1.5, 1, 1), "lines")) +
  labs(
    y = "Proportion of Macrophage Polarization",
    fill = "Polarized State")
print(stackedbar_mac_pol)



str(immune_subset.clean@meta.data)

## Can we show how condition --> polarization --> cell cycle score may be related?
# More color stuff - for conditions
Condition = c(
  "PD" = "#29339B",
  "BAL" = "#e9c716",
  "PD_3D" = "#4a2377",
  "Mixed_3D" = "#ff585e",
  "BAL_3D" = "#00b0be")

# Alluvial of Polar --> Cell Cycle ## DON'T REALLY LIKE BUT SAVING BC
# Phase color 
phase_colors <- c(
  "S" = "#16CEA6",
  "G1" = "#D2691E",
  "G2M" = "#AEC5EB")

mac_colors <- c(
  "Pro_Inflamm_Mac" = "#9497fd",
  "Anti_Inflamm_Mac" = "#ff9e80")
alluvial_data_simple <- immune_subset.clean@meta.data %>%
  filter(CellType.sub %in% c("Pro_Inflamm_Mac", "Anti_Inflamm_Mac")) %>%
  count(CellType.sub, Phase)
# Alluvial plot - just a test
alluv_polar_phase = ggplot(alluvial_data_simple,
       aes(axis1 = CellType.sub, axis2 = Phase, y = n)) +
  geom_alluvium(aes(fill = CellType.sub), width = 0.25, alpha = 0.9, color = "gray30") +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.25, color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5, color = "black") +
  scale_x_discrete(limits = c("CellType.sub", "Phase"), expand = c(.05, .05)) +
  scale_fill_manual(values = c(mac_colors, phase_colors)) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 20, 10, 20)) +
  labs(title = NULL, y = NULL, x = NULL)
print(alluv_polar_phase)
setwd("~/Desktop/Manuscript Figures/Figure 4")
ggsave("alluv_polar_phase.png", plot = alluv_polar_phase, width = 7, height = 5, dpi = 600)

## UMAP by Phase
mac.sub_umap_Phase = DimPlot(immune_subset.clean, group.by = "Phase", cols = phase_colors) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 14),
    legend.key.size = unit(1.2, "lines"),  # Increase key (color box) size
    plot.title = element_blank())
print(mac.sub_umap_Phase)
setwd("~/Desktop/Manuscript Figures/Figure 4")
ggsave("mac.sub_umap_Phase.png", plot = mac.sub_umap_Phase, width = 5, height = 5, dpi = 600)

# Volcano or heatmap??? To be determined...
