
# Date: 03-26-2025
# Memo: Re-running slingshot on my subset proximal object

# close all, clear all always :) 
graphics.off()
rm(list = ls())

library(slingshot)
# library(tradeSeq) # not compatible with this version of R
library(Seurat)
library(ggplot2)
library(dplyr)
library(SingleCellExperiment)
library(grid)
library(viridis)
library(patchwork)
library(tidyr)
library(purrr)


# Load object
setwd("~/Desktop/Single Cell/Engineered Subset (Proximal) for Pseudotime")
load("eng.subset.prox.int.Robj")

# Check structure
str(eng.subset.prox.int@meta.data)

# Look at embedding
DimPlot(eng.subset.prox.int, reduction = "umap.rpca", group.by = "Condition")
DimPlot(eng.subset.prox.int, reduction = "umap.rpca", group.by = "Orig_ID")
DimPlot(eng.subset.prox.int, reduction = "umap.rpca", group.by = "Phase")
DimPlot(eng.subset.prox.int, reduction = "umap.rpca", group.by = "CellType") # These are annotations from our old condition level ojects
DimPlot(eng.subset.prox.int, reduction = "umap.rpca", group.by = "CellType_Prox_Subset")

DimPlot(eng.subset.prox.int, reduction = "umap.rpca", group.by = "CellType_Prox_Subset", cols = CellType.cols)
# Set desired order for the legend
eng.subset.prox.int$CellType_Prox_Subset <- factor(
  eng.subset.prox.int$CellType_Prox_Subset,
  levels = c("Cycling_Proximal_Epi", "Basal_Like", "Hillock_Basal", "Hillock_Luminal"))

# Then make the plot
subset.umap <- DimPlot(
  eng.subset.prox.int,
  reduction = "umap.rpca",
  group.by = "CellType_Prox_Subset",
  cols = CellType.cols,
  label = TRUE,
  repel = TRUE) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),  # 👈 bigger legend text
    plot.margin = margin(10, 10, 10, 10)) +
  ggtitle(NULL)
print(subset.umap)
setwd("~/Desktop/Manuscript Figures/Figure 3")
ggsave("subset.umap.png", plot = subset.umap, width = 7.5, height = 6, dpi = 600)

# DimPlot but split by condition
DimPlot(eng.subset.prox.int, reduction = "umap.rpca", split.by = "Condition")
# Plot general object
DimPlot(eng.subset.prox.int, reduction = "umap.rpca")

# Check embeddings
names(eng.subset.prox.int@reductions)
# pca, umap, integrated.rpca, umap.rpca

# Convert Seurat object to SingleCellExperiment
sce.hillock <- as.SingleCellExperiment(eng.subset.prox.int)

# Assign Seurat's integrated UMAP to the reducedDims slot
reducedDims(sce.hillock)$UMAP <- eng.subset.prox.int@reductions$umap.rpca@cell.embeddings

# Make sure our cell type annotation data is part of the single cell object (sce.hillock)
colData(sce.hillock)$CellType_Prox_Subset <- eng.subset.prox.int$CellType_Prox_Subset

# Run Slingshot on the integrated embedding
sce.hillock <- slingshot(
  sce.hillock,
  clusterLabels = sce.hillock$CellType.combined.Integrated,
  reducedDim = "UMAP", # which has been assifned as umap.rpca
  start.clus = "Cycling_Proximal_Epi",
  omega = T)
summary(sce.hillock$slingPseudotime_1)

# Basic plot
summary(sce.hillock$slingPseudotime_1)
plot(reducedDim(sce.hillock, 'UMAP'), col = colorRampPalette(c("blue", "yellow", "red"))(100)[cut(sce.hillock$slingPseudotime_1, breaks=100)], pch = 16, asp = 1)
lines(SlingshotDataSet(sce.hillock), lwd = 2, col = 'black')

##### More refined #####
# Extract UMAP coordinates and pseudotime values
umap_coords <- as.data.frame(reducedDim(sce.hillock, 'UMAP'))
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$pseudotime <- sce.hillock$slingPseudotime_1
# Extract lineage curve coordinates from Slingshot
sling_curves <- SlingshotDataSet(sce.hillock)
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
umap_coords <- as.data.frame(reducedDim(sce.hillock, 'UMAP'))
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$CellType_Prox_Subset <- sce.hillock$CellType_Prox_Subset
# Define custom color palette for cell types
celltype_colors <- c("Hillock_Basal" = "#2800c7","Hillock_Luminal" = "#ff0000", "Basal_Like" = "#c6e308",
                     "Cycling_Proximal_Epi" = "#c4c2ff")
# Extract lineage curve coordinates from Slingshot
sling_curves <- SlingshotDataSet(sce.hillock)
# Convert curves to a data frame
curve_df <- do.call(rbind, lapply(sling_curves@curves, function(curve) {
  data.frame(UMAP1 = curve$s[curve$ord, 1],
             UMAP2 = curve$s[curve$ord, 2])
}))

# Plot by CellType.sub with custom colors
p1 <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = CellType_Prox_Subset)) +
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
ggsave(filename = "Slingshot_hillock.sub_CellType.png", plot = p1, width = 7, height = 5, dpi = 600)

# Add pseudotime to meta-data
eng.subset.prox.int$Pseudotime <- sce.hillock$slingPseudotime_1
# Add to Seurat metadata
eng.subset.prox.int <- AddMetaData(eng.subset.prox.int, metadata = sce.hillock$slingPseudotime_1, col.name = "Pseudotime_Slingshot")
# Check to make sure it was added
str(eng.subset.prox.int@meta.data)
eng.subset.prox.int@meta.data <- eng.subset.prox.int@meta.data[, !colnames(eng.subset.prox.int@meta.data) %in% c("Pseudotime")]
save(eng.subset.prox.int, file = "eng.subset.prox.int_sling_03-17-2025.Robj")


########### Now that it is saved, make plots for SE Manuscript Figure 3 ###########

# Recreate umap_coords with pseudotime
umap_coords <- as.data.frame(reducedDim(sce.hillock, 'UMAP'))
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$pseudotime <- sce.hillock$slingPseudotime_1

# Pseudotime plot without curve or axes
p0_clean <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = pseudotime)) +
  geom_point(size = 1.0, alpha = 0.6) +
  scale_color_viridis(
    option = "inferno",
    name = "Pseudotime",
    guide = guide_colorbar(
      title.position = "top",    # Places title above color bar
      title.hjust = 0.5,         # Centers the title
      barwidth = 5,             # Adjust width of color bar
      barheight = 1.0            # Adjust height of color bar
    )) +
  theme_void() +
  labs(title = "") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12, face = "bold"),  # Make legend title bold
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
print(p0_clean)
ggsave(filename = "Slingshot_hillock.sub_pseudotime.png", plot = p0_clean, width = 7, height = 5, dpi = 600)

# Recreate umap_coords with cell type info
umap_coords <- as.data.frame(reducedDim(sce.hillock, 'UMAP'))
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$CellType_Prox_Subset <- sce.hillock$CellType_Prox_Subset
p1_clean <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = CellType_Prox_Subset)) +
  geom_point(size = 0.75, alpha = 0.8) +
  scale_color_manual(values = celltype_colors) +  
  geom_path(data = curve_df, aes(x = UMAP1, y = UMAP2),
            color = "black", size = 0.8) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_blank())
print(p1_clean)
ggsave(filename = "Slingshot_hillock.sub_celltype.curve.png", plot = p1_clean, width = 7, height = 5, dpi = 600)

####################################################### BREAK ##############################################################

# Add UMAP embedding from Seurat object to SingleCellExperiment object
reducedDims(sce.hillock)$UMAP <- eng.subset.prox.int@reductions$umap.rpca@cell.embeddings
# Verify that UMAP is now available
print(reducedDimNames(sce.hillock))

plotGeneExpressionUMAP <- function(gene, sce, dim_name = "UMAP") {
  # Check if the gene exists
  if (!(gene %in% rownames(assays(sce)$counts))) {
    stop(paste("Gene", gene, "not found in the dataset"))
  }
  
  # Extract UMAP coordinates using the correct dimension name
  umap_coords <- reducedDim(sce, dim_name)
  
  df <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    Expression = log1p(assays(sce)$counts[gene, ]))
  
  # Order by Expression (low to high) to plot higher expression on top
  df <- df[order(df$Expression, decreasing = FALSE), ]  
  
  # Plot UMAP with gene expression overlay
  ggplot(df, aes(x = UMAP1, y = UMAP2, color = Expression)) +
    geom_point(size = 1.5, alpha = 0.8) +
    scale_color_viridis(option = "plasma", name = "Expression") +
    labs(
      title = paste(gene),
      x = "umap_1",
      y = "umap_2") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(size = 12),  # Bold axis labels
      axis.text = element_text(size = 12),  # Restore axis text
      axis.ticks = element_line(color = "black"),  # Restore axis ticks
      panel.grid = element_blank(),  # Remove grid lines
      axis.line = element_line(color = "black", size = 0.5))
}

Idents(eng.subset.prox.int) = eng.subset.prox.int$CellType_Prox_Subset
mark = FindAllMarkers(eng.subset.prox.int, min.pct = 0.1,logfc.threshold = 0.1)
mark$ratio = mark$pct.1/mark$pct.2
mark$power = mark$ratio*mark$avg_log2FC
View(mark)

# Example usage after adding UMAP
p2 = plotGeneExpressionUMAP("Krt5", sce.hillock, dim_name = "UMAP")
p3 = plotGeneExpressionUMAP("Il33", sce.hillock, dim_name = "UMAP")
p4 = plotGeneExpressionUMAP("Aqp3", sce.hillock, dim_name = "UMAP")
p5 = plotGeneExpressionUMAP("Tp63", sce.hillock, dim_name = "UMAP")
p6 = plotGeneExpressionUMAP("Krt13", sce.hillock, dim_name = "UMAP")
p7 = plotGeneExpressionUMAP("Serpinb2", sce.hillock, dim_name = "UMAP")
# p8 = plotGeneExpressionUMAP("Dsc2", sce.hillock, dim_name = "UMAP")
p9 = plotGeneExpressionUMAP("Evpl", sce.hillock, dim_name = "UMAP")
p10 = plotGeneExpressionUMAP("Mal", sce.hillock, dim_name = "UMAP")
p11 = plotGeneExpressionUMAP("Ppl", sce.hillock, dim_name = "UMAP")
p12 = plotGeneExpressionUMAP("Tgfb2", sce.hillock, dim_name = "UMAP")

# Random experimenting
plotGeneExpressionUMAP("Cd68", sce.hillock, dim_name = "UMAP")
plotGeneExpressionUMAP("Cldn17", sce.hillock, dim_name = "UMAP")
plotGeneExpressionUMAP("Il17b", sce.hillock, dim_name = "UMAP")
plotGeneExpressionUMAP("Cited2", sce.hillock, dim_name = "UMAP")
plotGeneExpressionUMAP("Tgfbi", sce.hillock, dim_name = "UMAP")
plotGeneExpressionUMAP("Snai2", sce.hillock, dim_name = "UMAP")
plotGeneExpressionUMAP("Cavin3", sce.hillock, dim_name = "UMAP")
plotGeneExpressionUMAP("Itgb6", sce.hillock, dim_name = "UMAP")

# Plotting immune things for MSBR
VlnPlot(eng.subset.prox.int, features = c("Krt13","Krt5","Sprr1a","Il1a","S100a8","Mal","Slpi","Cd68"))
FeaturePlot(eng.subset.prox.int, features = c("Krt13","Krt5","Sprr1a","Il1a","S100a8","Mal","Slpi","Cd68"))

####################################################### BREAK ##############################################################

## Figure 3 Smooth-Spline Curves of Gene Expression Over Pseudotime
# Packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

# Select genes we want to plot
goi <- c("Krt13", "Sprr1a","Evpl",
         "Krt5", "Tp63","Il33","Notch1","Aldh3a1",
         "Prss22","Prdm1","Krt7","Maff","Krt80","Krt78","Adm","Cnfn")

# Pull gene expression, pseudotime, and grouping variable
df <- FetchData(eng.subset.prox.int, 
                vars = c("Pseudotime_Slingshot", "CellType_Prox_Subset", goi, "Condition"))

# Tidy to long format
df_long <- df %>%
  pivot_longer(cols = all_of(goi), names_to = "Gene", values_to = "Expression")

# By Condition
ggplot(df_long, aes(x = Pseudotime_Slingshot, y = Expression, color = Condition)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, size = 1.2) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("BAL_3D" = "#00b0be", 
                                "Mixed_3D" = "#ff585e", 
                                "PD_3D" = "#4a2377")) +
  theme_minimal(base_size = 14) +
  labs(x = "Pseudotime", y = "Expression", color = "Condition") +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    legend.position = "right")


## Filtering
library(dplyr)
library(tidyr)
library(ggplot2)

# Define gene groups
all_conditions <- c("Krt13", "Sprr1a", "Evpl")
bal_only <- c("Krt5", "Tp63", "Il33", "Notch1", "Aldh3a1","Col17a1")
mixed_pd_only <- c("Prss22", "Prdm1", "Krt7", "Maff", "Krt80", "Adm")

# Pull data
goi <- c(all_conditions, bal_only, mixed_pd_only)

df <- FetchData(eng.subset.prox.int, 
                vars = c("Pseudotime_Slingshot", "Condition", goi))

# Reshape to long format
df_long <- df %>%
  pivot_longer(cols = all_of(goi), names_to = "Gene", values_to = "Expression")

# Filter based on gene + condition rules
df_filtered <- df_long %>%
  filter(
    (Gene %in% all_conditions) |
      (Gene %in% bal_only & Condition == "BAL_3D") |
      (Gene %in% mixed_pd_only & Condition %in% c("Mixed_3D", "PD_3D"))
  )

df_filtered$Gene <- factor(df_filtered$Gene, levels = c(
  "Krt13", "Sprr1a", "Evpl",
  "Krt5", "Tp63", "Il33", "Notch1", "Aldh3a1","Col17a1",
  "Prss22", "Prdm1", "Krt7", "Maff", "Krt80", "Adm"))

pseudotime_DGE.curves = ggplot(df_filtered, aes(x = Pseudotime_Slingshot, y = Expression, color = Condition)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, size = 1.2) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +  # adjust ncol to fit your figure
  scale_color_manual(values = c("BAL_3D" = "#00b0be", 
                                "Mixed_3D" = "#ff585e", 
                                "PD_3D" = "#4a2377")) +
  theme_minimal(base_size = 14) +
  labs(x = "Pseudotime", y = "Expression", color = "Condition") +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    legend.position = "right")
print(pseudotime_DGE.curves)

#### Now if we want to shade under the curves
# Define pseudotime grid
pseudotime_grid <- seq(0, 20, length.out = 200)
# Create smooth curves using regularized grid
smoothed_curves <- df_filtered %>%
  group_by(Gene, Condition) %>%
  nest() %>%
  mutate(
    fitted = map(data, ~ {
      # Fit loess model
      model <- loess(Expression ~ Pseudotime_Slingshot, data = .x, span = 0.5)
      
      # Predict on a consistent grid
      data.frame(
        Pseudotime_Slingshot = pseudotime_grid,
        Smoothed = predict(model, newdata = data.frame(Pseudotime_Slingshot = pseudotime_grid))
      )
    })
  ) %>%
  select(Gene, Condition, fitted) %>%
  unnest(fitted)

# Now plot shading
pseudotime_DGE_shaded = ggplot(smoothed_curves, aes(x = Pseudotime_Slingshot, y = Smoothed, color = Condition, fill = Condition)) +
  geom_ribbon(aes(ymin = 0, ymax = Smoothed), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("BAL_3D" = "#00b0be", "Mixed_3D" = "#ff585e", "PD_3D" = "#4a2377")) +
  scale_fill_manual(values = c("BAL_3D" = "#00b0be", "Mixed_3D" = "#ff585e", "PD_3D" = "#4a2377")) +
  labs(x = "Pseudotime", y = "Expression", color = "Condition", fill = "Condition") +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.3, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"))
print(pseudotime_DGE_shaded)
ggsave("pseudotime_DGE_shaded.png", pseudotime_DGE_shaded, width = 12, height = 10)

### Shading 2.0 - narrower (more selective for figure 3 - limited space)
# Define gene groups
all_conditions <- c("Krt13", "Sprr1a")
bal_only <- c("Krt5", "Tp63", "Notch1","Col17a1")
mixed_pd_only <- c("Prss22", "Krt7", "Krt80", "Adm")

# Pull data
goi <- c(all_conditions, bal_only, mixed_pd_only)

df <- FetchData(eng.subset.prox.int, 
                vars = c("Pseudotime_Slingshot", "Condition", goi))

# Reshape to long format
df_long <- df %>%
  pivot_longer(cols = all_of(goi), names_to = "Gene", values_to = "Expression")

# Filter based on gene + condition rules
df_filtered <- df_long %>%
  filter(
    (Gene %in% all_conditions) |
      (Gene %in% bal_only & Condition == "BAL_3D") |
      (Gene %in% mixed_pd_only & Condition %in% c("Mixed_3D", "PD_3D"))
  )

df_filtered$Gene <- factor(df_filtered$Gene, levels = c(
  "Krt13", "Sprr1a",
  "Krt5", "Tp63", "Notch1","Col17a1",
  "Prss22", "Krt7", "Krt80", "Adm"))

pseudotime_DGE.curves = ggplot(df_filtered, aes(x = Pseudotime_Slingshot, y = Expression, color = Condition)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, size = 1.2) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 2) +  # adjust ncol to fit your figure
  scale_color_manual(values = c("BAL_3D" = "#00b0be", 
                                "Mixed_3D" = "#ff585e", 
                                "PD_3D" = "#4a2377")) +
  theme_minimal(base_size = 14) +
  labs(x = "Pseudotime", y = "Expression", color = "Condition") +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    legend.position = "right")
print(pseudotime_DGE.curves)

#### Now if we want to shade under the curves
# Define pseudotime grid
pseudotime_grid <- seq(0, 21, length.out = 200)
# Create smooth curves using regularized grid
smoothed_curves <- df_filtered %>%
  group_by(Gene, Condition) %>%
  nest() %>%
  mutate(
    fitted = map(data, ~ {
      # Fit loess model
      model <- loess(Expression ~ Pseudotime_Slingshot, data = .x, span = 0.5)
      
      # Predict on a consistent grid
      data.frame(
        Pseudotime_Slingshot = pseudotime_grid,
        Smoothed = predict(model, newdata = data.frame(Pseudotime_Slingshot = pseudotime_grid))
      )
    })
  ) %>%
  select(Gene, Condition, fitted) %>%
  unnest(fitted)

# Now plot shading
pseudotime_DGE_shaded = ggplot(smoothed_curves, aes(x = Pseudotime_Slingshot, y = Smoothed, color = Condition, fill = Condition)) +
  geom_ribbon(aes(ymin = 0, ymax = Smoothed), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("BAL_3D" = "#00b0be", "Mixed_3D" = "#ff585e", "PD_3D" = "#4a2377")) +
  scale_fill_manual(values = c("BAL_3D" = "#00b0be", "Mixed_3D" = "#ff585e", "PD_3D" = "#4a2377")) +
  labs(x = "Pseudotime", y = "Expression", color = "Condition", fill = "Condition") +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.3, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9))
print(pseudotime_DGE_shaded)
ggsave("pseudotime_DGE_shaded.png", pseudotime_DGE_shaded, width = 12, height = 10)

# If we want jittered, too
pseudotime_DGE_shaded.jittered <- ggplot(smoothed_curves, aes(x = Pseudotime_Slingshot, y = Smoothed, color = Condition, fill = Condition)) +
  # Add jittered points from raw data
  geom_point(data = df_filtered,
             aes(x = Pseudotime_Slingshot, y = Expression, color = Condition),
             size = 0.4, alpha = 0.4, position = position_jitter(width = 0.1)) +
  # Shaded curve underlay
  geom_ribbon(aes(ymin = 0, ymax = Smoothed), alpha = 0.2, color = NA) +
  # Smooth line
  geom_line(size = 1.2) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("BAL_3D" = "#00b0be", "Mixed_3D" = "#ff585e", "PD_3D" = "#4a2377")) +
  scale_fill_manual(values = c("BAL_3D" = "#00b0be", "Mixed_3D" = "#ff585e", "PD_3D" = "#4a2377")) +
  labs(x = "Pseudotime", y = "Expression", color = "Condition", fill = "Condition") +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.3, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9))
print(pseudotime_DGE_shaded.jittered)
ggsave("pseudotime_DGE_shaded.jittered.png", pseudotime_DGE_shaded.jittered, width = 5.25, height = 8, dpi = 600)


## Make pseudotime bar to add manually in illustrator
# Create pseudotime gradient from 0 to 20
# Define pseudotime range
pseudotime_gradient <- data.frame(
  x = seq(0, 20, length.out = 500),  # smoothness
  y = 1)

# Control the y height (bar thickness) here
pseudotime_bar <- ggplot(pseudotime_gradient, aes(x = x, y = y, fill = x)) +
  geom_tile(height = 0.05) +  # Change this to control bar thickness (e.g., 0.2 for thinner)
  scale_fill_viridis_c(option = "magma", direction = -1) +
  labs(x = "Pseudotime", y = NULL) +
  coord_cartesian(ylim = c(0.5, 1.5)) +  # vertically center it and control bar range
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    legend.position = "none",
    plot.margin = margin(5, 20, 5, 20))
print(pseudotime_bar)
ggsave("pseudotime_gradient_magma.png", pseudotime_bar, width = 3, height = 3, dpi = 600)








