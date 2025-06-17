
# Load packages
library(dplyr)
library(tibble)
library(EnhancedVolcano)
library(ggplot2)
library(Seurat)

# load in mac subset
setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets")
# Load data - this is the object 
load("immune_subset.clean.Robj")

# Look at structure
str(immune_subset.clean@meta.data)

# Get expression matrix (binary: expressed or not)
expression_matrix <- GetAssayData(immune_subset.clean, assay = "RNA", slot = "counts") > 0
Idents(immune_subset.clean) = immune_subset.clean$CellType.sub
Idents(immune_subset.clean) = immune_subset.clean$Condition

# Perform DGE analysis to be used for volcano plot
Idents(immune_subset.clean) <- immune_subset.clean$seurat_clusters
mac.pol_DGE <- FindMarkers(immune_subset.clean, ident.1 = "0", ident.2 = c("1","2"),
                           min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE, test.use = "wilcox")
mac.pol_DGE$ratio <- mac.pol_DGE$pct.1 / mac.pol_DGE$pct.2
mac.pol_DGE$power <- mac.pol_DGE$ratio * mac.pol_DGE$avg_log2FC
View(mac.pol_DGE)

##### We should be able to just run EnhancedVolcano on these dge_results ######
# Picking out genes from dge
highlighted_genes <- c('Tlr7','Tlr2','Tlr4','Il1b', # Broad
         'Fjx1','Cadm1','Prrt1','Gcgr','Mrgprx2','Atp1b1','Crip2','Ms4a4a','Ereg','Ccl4','Ccl2',
         'Rnf150','Gpr155','Rcn3','Pparg','Cd84','Cela1','Bmp1','Il18bp','Scarb1','Tnfsf9',
         'Vcan','Lsp1','Cxcl3','Cxcl1','Cxcl2','Trem3','Arg1','Cd40','S100a4','S100a6','Tnfsf9','Cd14',
         'Osm','Rps15al4','Slpi','Lsp1','Nos2','Il1a','Cxcl6','Fxyd2','Pclaf','Ereg','Prg4')
# Store EnhancedVolcano plot as a variable
volcano_plot_1 <- EnhancedVolcano(mac.pol_DGE,
                                  lab = rownames(mac.pol_DGE),
                                  x = 'avg_log2FC',
                                  y = 'p_val_adj',
                                  selectLab = highlighted_genes,
                                  title = "",
                                  subtitle = NULL,
                                  xlab = bquote(~Log[2]~ 'fold change'),
                                  ylab = bquote(~"-"~Log[10]~ 'Adjusted P-Value'),
                                  pCutoff = 0.05,
                                  col = c('black', '#432371', '#0D41E1', '#A2D638'),
                                  FCcutoff = 0.5,
                                  pointSize = 2.0,
                                  labSize = 4.0,
                                  boxedLabels = FALSE,
                                  drawConnectors = TRUE,
                                  widthConnectors = 0.5,
                                  colAlpha = 0.6,
                                  legendPosition = "none")
# Modify the theme to center the title
volcano_plot_1 + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(volcano_plot_1)
ggsave("Hillock.Polarization_Volcano.png", volcano_plot_1,
       width = 12, height = 8, dpi = 600, bg = "white")
# Cleaning up
volcano_plot_1_clean <- volcano_plot_1 +
  theme_minimal(base_size = 12) +  # Start from clean theme
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none")
print(volcano_plot_1_clean)
setwd("~/Desktop/Manuscript Figures/Figure 4")
ggsave("Macrophage.Polarization_Volcano_clean.png", volcano_plot_1_clean,
       width = 8.5, height = 6, dpi = 600, bg = "white")
