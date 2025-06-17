
# close all, clear all
graphics.off()  
rm(list = ls())

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Load global object (have to change WD) (organoids)
setwd("~/Desktop/Datasets/Engineered Global Objects")
load("eng.subset.integrated_NodeAligned.Robj")

# Define cell types to highlight (from global object)
highlighted_types <- c("Stressed_Progenitor", "RAS_Like", "Secretory")

# Set Seurat cell identities to use node-aligned annotations (cell types)
Idents(eng.subset.integrated_HK) <- eng.subset.integrated_HK$CellType.NodeAligned

# Create a new metadata column categorizing cells into:
#   - highlighted epithelial subtypes (of interest)
#   - all others ("Other")
eng.subset.integrated_HK$Highlight_Group <- ifelse(
  eng.subset.integrated_HK$CellType.NodeAligned %in% highlighted_types,
  eng.subset.integrated_HK$CellType.NodeAligned,
  "Other")

# Define colors for selected subtypes and a muted gray for all others
highlight_colors <- c(
  "Stressed_Progenitor" = "#D6368A",
  "RAS_Like" = "#484FA1",
  "Secretory" = "#208B3A",
  "Other" = "lightgray")

# Generate a UMAP plot using the custom Highlight_Group and color palette
fig5_selec = DimPlot(
  eng.subset.integrated_HK,
  group.by = "Highlight_Group",
  pt.size = 0.5,
  cols = highlight_colors, reduction = "umap.rpca") +
  theme_void(base_size = 14) +  # No axes, labels, or grid
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_blank())
# Take a look at the plot
print(fig5_selec)
# Save plot as high-resolution PNG for publication
setwd("~/Desktop/Manuscript Figures/Figure 5")
ggsave("EpithelialFig5_Select.png", plot = fig5_selec, width = 8, height = 8, dpi = 600)


