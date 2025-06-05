# close all, clear all
graphics.off()  
rm(list = ls())
# Load packages
require(Seurat)
require(dplyr)

# Load ComplexHeatmap_3 (file is ComplexHeatmap_3.continuous) function
# Also load ComplexHeatmap_2 (for two vars, one being pseudotime continuous)

# Define object as whatever object we are working with
# Load local functions
# Or you can source it locally

# Set working directory to be true to get object
setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/Datasets")

# Load data - this is the object (our proximal subset with pseudotime_slingshot as meta-data)
load("immune_subset.clean.Robj")

# Define object as whatever object we are working with
object = immune_subset.clean
# Inspect data to get a handle on the metadata slot names
names(immune_subset.clean@meta.data)

# Look at the ordering of different metadata slots
table(immune_subset.clean$CellClass.combined.Integrated)
table(immune_subset.clean$CellType.sub)
table(immune_subset.clean$Condition)

# Re-order the metadata for plotting
immune_subset.clean$CellClass.combined.Integrated <- factor(immune_subset.clean$CellClass.combined.Integrated,
                                                            levels = c('Immune')) # Modify the order here if needed
immune_subset.clean$CellType_Prox_Subset <- factor(immune_subset.clean$CellType_Prox_Subset,
                                                   levels = c('Anti_Inflamm_Mac','Pro_Inflamm_Mac')) # Modify the levels here
immune_subset.clean$Condition <- factor(immune_subset.clean$Condition,
                                        levels = c('BAL_3D','Mixed_3D','PD_3D')) # Modify the levels here
# Look at the effect of re-ordering
table(immune_subset.clean$CellType.sub, immune_subset.clean$CellClass.combined.Integrated, immune_subset.clean$Condition)

# Color palette - pulling these from ColorPalette_OrganoidProj.R
col.pal <- list()
col.pal$CellClass.combined.Integrated <- c('#87B37A') #SE Custom for immune
names(col.pal$CellClass.combined.Integrated) <- c('Immune')
col.pal$CellType.sub <- c('#ff9e80','#9497fd')
names(col.pal$CellType.sub) <- c('Anti_Inflamm_Mac','Pro_Inflamm_Mac')
col.pal$Condition <- c('#00b0be','#ff585e','#4a2377') #SE Custom
names(col.pal$Condition) <- c('BAL_3D','Mixed_3D','PD_3D')

# Define a continuous color palette for Pseudotime

# Not necessary for this sample because it is already so small, lol
# But keeping in code for future merged objects with greater number of cells perhaps
# Create a downsampled object
# downsampled <- subset(immune_subset.clean, downsample = 20000)
# table(downsampled$orig.ident)
# Scale the object
# downsampled <- ScaleData(downsampled, features = rownames(downsampled))

# Create a marker list (as requested by Sam)
Idents(immune_subset.clean) <- immune_subset.clean$seurat_clusters
mac.pol_DGE <- FindMarkers(immune_subset.clean, ident.1 = "0", ident.2 = c("1","2"),
                           min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE, test.use = "wilcox")
mac.pol_DGE$ratio <- mac.pol_DGE$pct.1 / mac.pol_DGE$pct.2
mac.pol_DGE$power <- mac.pol_DGE$ratio * mac.pol_DGE$avg_log2FC
View(mac.pol_DGE)

# Define the features to plot based on DEGs
# top.marker.list <- mark %>% group_by(cluster) %>% top_n(10, power)
# GOI <- top.marker.list$gene

# Make a PNG output plot at 300 dpi for publication, width and height in real-world inches (feel free to change width and height as desired)
# With 3 vars, one being cont. pseudotime
# Always remember to scale data or code will prompt a bug
immune_subset.clean <- ScaleData(immune_subset.clean,features = rownames(immune_subset.clean))

# Maybe redundant from above, but these are the features I want in this specific plot
GOI <- c('Tgfb1','Irf5','Tlr2','Stat6','Stat3','Tlr4','Il1b','Cdkn1a', # Broad
         'Fjx1','Cadm1','Prrt1','Gcgr','Mrgprx2','Atp1b1',
         'Rnf150','Gpr155','Rcn3','Pparg','Cd84','Cela1',
         'Vcan','Lsp1','Cxcl3','Cxcl1','Cxcl2','Trem3','Arg1','Cd40',
         'Osm','Rps15al4','Slpi','Lsp1','Nos2','Il1a','Cxcl6') 
# Define row labels (what I want labeled)
row.labels <- c('Tgfb1','Irf5','Tlr2','Stat6','Stat3','Tlr4','Il1b','Cdkn1a', # Broad
                'Fjx1','Cadm1','Prrt1','Gcgr','Mrgprx2','Atp1b1',
                'Rnf150','Gpr155','Rcn3','Pparg','Cd84','Cela1',
                'Vcan','Lsp1','Cxcl3','Cxcl1','Cxcl2','Trem3','Arg1','Cd40',
                'Osm','Rps15al4','Slpi','Lsp1','Nos2','Il1a','Cxcl6')

# With only two vars (pseudotime and cell type)
png(file = 'plot2.png', width = 6.5, height = 6, units = 'in', res = 300)
ComplexHeatMap_2.inferno(object = immune_subset.clean,
                 data.type = 'RNA',
                 primary = 'CellType.sub',
                 secondary = 'Condition',
                 primary.cols = col.pal$CellType.sub,  # Continuous color scale
                 secondary.cols = col.pal$Condition,      # Categorical color scale
                 features = GOI,
                 labels = c('Condition','Cell Type'),
                 selected.row.anotations = row.labels,
                 selected.label.size = 8,
                 use.scale.data = TRUE,
                 range.frac = 0.5,
                 row.dendrogram = FALSE)
dev.off()

png(file = 'plot4.png', width = 7.5, height = 6, units = 'in', res = 300)
ComplexHeatMap_3.inferno(object = immune_subset.clean,
                 data.type = 'RNA',
                 primary = 'CellClass.combined.Integrated',
                 secondary = 'CellType.sub',
                 tertiary = 'Condition',  # ← Add the Condition tier
                 primary.cols = col.pal$CellClass.combined.Integrated,
                 secondary.cols = col.pal$CellType.sub,
                 tertiary.cols = col.pal$Condition,  # ← Pass the palette
                 features = GOI,
                 labels = c('Cell Class', 'Cell Type', 'Condition'),
                 selected.row.anotations = row.labels,
                 selected.label.size = 8,
                 use.scale.data = TRUE,
                 range.frac = 0.5,
                 row.dendrogram = FALSE)
dev.off()
