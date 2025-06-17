

## Memo: Slingshot on macrophage polarization?

# Set working directory to be true and load data
setwd("~/Desktop/Manuscript Figures/Figure 4/Objects")
load("immune_subset.clean.Robj")

library(slingshot)
# library(tradeSeq) # not compatible with this version of R
library(Seurat)
library(ggplot2)
library(dplyr)
library(SingleCellExperiment)
library(grid)
library(viridis)
library(patchwork)

# Plot general object
Idents(immune_subset.clean) = immune_subset.clean$CellType.sub
DimPlot(immune_subset.clean, reduction = "umap")
Idents(immune_subset.clean) = immune_subset.clean$seurat_clusters
DimPlot(immune_subset.clean, reduction = "umap")

# Get a sense of where cycling cells are
FeaturePlot(immune_subset.clean, features = c("Top2a"))
# Check embeddings
names(immune_subset.clean@reductions)
# pca, umap, integrated.rpca, umap.rpca

# Convert Seurat object to SingleCellExperiment
sce.mac <- as.SingleCellExperiment(immune_subset.clean)

# Assign Seurat's integrated UMAP to the reducedDims slot
reducedDims(sce.mac)$UMAP <- immune_subset.clean@reductions$umap@cell.embeddings

# Make sure our cell type annotation data is part of the single cell object (sce.hillock)
colData(sce.mac)$CellType.sub <- immune_subset.clean$CellType.sub
colData(sce.mac)$clusters <- immune_subset.clean$seurat_clusters

# Run Slingshot on the integrated embedding
sce.mac <- slingshot(
  sce.mac,
  clusterLabels = sce.mac$CellType.sub,
  reducedDim = "UMAP", # which has been assifned as umap.rpca
  #start.clus = "Pro_Inflamm_Mac",
  omega = T)
summary(sce.mac$slingPseudotime_1)

# Basic plot
summary(sce.mac$slingPseudotime_1)
plot(reducedDim(sce.mac, 'UMAP'), col = colorRampPalette(c("blue", "yellow", "red"))(100)[cut(sce.mac$slingPseudotime_1, breaks=100)], pch = 16, asp = 1)
lines(SlingshotDataSet(sce.mac), lwd = 2, col = 'black')

##### More refined #####
# Extract UMAP coordinates and pseudotime values
umap_coords <- as.data.frame(reducedDim(sce.mac, 'UMAP'))
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$pseudotime <- sce.mac$slingPseudotime_1
# Extract lineage curve coordinates from Slingshot
sling_curves <- SlingshotDataSet(sce.mac)
# Need to convert curves to a data frame
curve_df <- do.call(rbind, lapply(sling_curves@curves, function(curve) {
  data.frame(UMAP1 = curve$s[curve$ord, 1],
             UMAP2 = curve$s[curve$ord, 2])
}))

p0 <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = pseudotime)) +
  geom_point(size = 1.0, alpha = 0.6) +
  scale_color_viridis(option = "viridis", name = "Pseudotime") +  
  geom_path(data = curve_df, aes(x = UMAP1, y = UMAP2),
            color = "black", size = 1) +  # Add trajectory lines
  theme_minimal() +
  labs(title = "Slingshot Pseudotime Trajectory", x = "umap_1", y = "umap_2") +
  theme(
    axis.title = element_text(size = 12),  # Bold axis labels
    axis.text = element_text(size = 12),  # Restore axis text
    axis.ticks = element_line(),  # Restore axis ticks
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove full panel border
    axis.line = element_line(color = "black", size = 1),
    legend.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
# Display the plot
print(p0)

# Plot by cell type labels
umap_coords <- as.data.frame(reducedDim(sce.mac, 'UMAP'))
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$CellType.sub <- sce.mac$CellType.sub
# Define custom color palette for cell types
celltype_colors <- c("Pro_Inflamm_Mac" = "#9497fd","Anti_Inflamm_Mac" = "#ff9e80")
# Extract lineage curve coordinates from Slingshot
sling_curves <- SlingshotDataSet(sce.mac)
# Convert curves to a data frame
curve_df <- do.call(rbind, lapply(sling_curves@curves, function(curve) {
  data.frame(UMAP1 = curve$s[curve$ord, 1],
             UMAP2 = curve$s[curve$ord, 2])
}))

# Plot by CellType.sub with custom colors
p1 <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = CellType.sub)) +
  geom_point(size = 0.75, alpha = 0.8) +
  scale_color_manual(values = celltype_colors, name = "Cell Type") +  
  geom_path(data = curve_df, aes(x = UMAP1, y = UMAP2),
            color = "black", size = 1) +  # Add trajectory lines
  theme_minimal() +
  labs(title = "Slingshot Trajectory by Cell Type",
       x = "umap_1",  
       y = "umap_2") +
  theme(
    axis.title = element_text(size = 12),  # Bold axis labels
    axis.text = element_text(size = 12),  # Restore axis text
    axis.ticks = element_line(),  # Restore axis ticks
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove full panel border
    axis.line = element_line(color = "black", size = 1),  # Add X and Y axis lines
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center title
    legend.title = element_text(size = 12))
# Display the plot
print(p1)

# Extract UMAP coordinates and pseudotime values
umap_coords <- as.data.frame(reducedDim(sce.mac, 'UMAP'))
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$pseudotime <- sce.mac$slingPseudotime_1
# Extract lineage curve coordinates from Slingshot
sling_curves <- SlingshotDataSet(sce.mac)
# Need to convert curves to a data frame
curve_df <- do.call(rbind, lapply(sling_curves@curves, function(curve) {
  data.frame(UMAP1 = curve$s[curve$ord, 1],
             UMAP2 = curve$s[curve$ord, 2])
}))
p0 <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = pseudotime)) +
  geom_point(size = 1.0, alpha = 0.6) +
  scale_color_viridis(option = "inferno", name = "Pseudotime") +  
  geom_path(data = curve_df, aes(x = UMAP1, y = UMAP2),
            color = "black", size = 1) +
  theme_void() +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = "bottom",  # keep it clean and centered
    plot.title = element_blank())
# Show the plot
print(p0)
# Save it
setwd("~/Desktop/Manuscript Figures/Figure 4")
ggsave(filename = "slingshot_mac.pol_pseudotime.png", plot = p0, width = 7, height = 6, dpi = 600)

# Plot by cell type labels
umap_coords <- as.data.frame(reducedDim(sce.mac, 'UMAP'))
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$CellType.sub <- sce.mac$CellType.sub
# Define custom color palette for cell types
celltype_colors <- c("Pro_Inflamm_Mac" = "#9497fd","Anti_Inflamm_Mac" = "#ff9e80")
# Extract lineage curve coordinates from Slingshot
sling_curves <- SlingshotDataSet(sce.mac)
# Convert curves to a data frame
curve_df <- do.call(rbind, lapply(sling_curves@curves, function(curve) {
  data.frame(UMAP1 = curve$s[curve$ord, 1],
             UMAP2 = curve$s[curve$ord, 2])
}))
p1 <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = CellType.sub)) +
  geom_point(size = 0.75, alpha = 0.8) +
  scale_color_manual(values = celltype_colors, name = "Cell Type") +
  geom_path(data = curve_df, aes(x = UMAP1, y = UMAP2),
            color = "black", size = 1) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_blank()
  )
print(p1)
setwd("~/Desktop/Manuscript Figures/Figure 4")
ggsave(filename = "slingshot_mac.pol_curve.celltype.png", plot = p1, width = 7, height = 4, dpi = 600)


