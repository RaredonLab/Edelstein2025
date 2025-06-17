

## MEMO: Making condition level umaps by CellType for figure showing BF/HE/UMAP embeddings
# Updating to use NodeAligned annotations (from our global object) that have been mapped back to condition level objects

# Load in objects from folder on local machine (SEE)
setwd("~/Desktop/Datasets/Condition Level Objects")
load("BAL_3D.integrated_HK.NodeAligned.Robj") # BAL_3D
load("PD_3D.integrated_HK.NodeAligned.Robj") # PD_3D
load("Mixed_3D.integrated_HK.NodeAligned.Robj") # Mixed_3D

table(BAL_3D.integrated_HK$CellType.NodeAligned)
table(PD_3D.integrated_HK$CellType.NodeAligned)
table(Mixed_3D.integrated_HK$CellType.NodeAligned)

# Make specific node-aligned meta-data slots for the condition level that mostly stick with the global node aligned

## PD_3D ##
table(PD_3D.integrated_HK$CellType.NodeAligned)
# Copy the original metadata column
PD_3D.integrated_HK$CellType.NodeAligned.Condition <- PD_3D.integrated_HK$CellType.NodeAligned
# Rename specific entries
PD_3D.integrated_HK$CellType.NodeAligned.Condition[
  PD_3D.integrated_HK$CellType.NodeAligned.Condition == "Polarized_Mac"] <- "Pro_Inflamm_Mac"
PD_3D.integrated_HK$CellType.NodeAligned.Condition[
  PD_3D.integrated_HK$CellType.NodeAligned.Condition == "Hillock_Like"] <- "Hillock_Luminal"
# Check the result
table(PD_3D.integrated_HK$CellType.NodeAligned.Condition)
str(PD_3D.integrated_HK@meta.data)
save(PD_3D.integrated_HK, file = "PD_3D.integrated_HK.NodeAligned.Robj")

## BAL_3D ##
# Copy the original metadata column
BAL_3D.integrated_HK$CellType.NodeAligned.Condition <- BAL_3D.integrated_HK$CellType.NodeAligned
# Rename specific entries
BAL_3D.integrated_HK$CellType.NodeAligned.Condition[
  BAL_3D.integrated_HK$CellType.NodeAligned.Condition == "Polarized_Mac"] <- "Anti_Inflamm_Mac"
BAL_3D.integrated_HK$CellType.NodeAligned.Condition[
  BAL_3D.integrated_HK$CellType.NodeAligned.Condition == "Hillock_Like"] <- "Hillock_Basal"
# Reassign Pdgfrb+_Pericytes to Basal_Like
BAL_3D.integrated_HK$CellType.NodeAligned.Condition[
  BAL_3D.integrated_HK$CellType.NodeAligned.Condition == "Pdgfrb+_Pericytes"] <- "Basal_Like"
# Check the updated counts
table(BAL_3D.integrated_HK$CellType.NodeAligned.Condition)
str(BAL_3D.integrated_HK@meta.data)
save(BAL_3D.integrated_HK, file = "BAL_3D.integrated_HK.NodeAligned.Robj")

## Mixed_3D ##
# Copy the original metadata column
Mixed_3D.integrated_HK$CellType.NodeAligned.Condition <- Mixed_3D.integrated_HK$CellType.NodeAligned
# Rename polarized mac to their anti inflamm character
Mixed_3D.integrated_HK$CellType.NodeAligned.Condition[
  Mixed_3D.integrated_HK$CellType.NodeAligned.Condition == "Polarized_Mac"] <- "Pro_Inflamm_Mac"
# Rename hillock_like to the specific hillock flavo
Mixed_3D.integrated_HK$CellType.NodeAligned.Condition[
  Mixed_3D.integrated_HK$CellType.NodeAligned.Condition == "Hillock_Like"] <- "Hillock_Luminal"
# Check the updated result
table(Mixed_3D.integrated_HK$CellType.NodeAligned.Condition)
str(Mixed_3D.integrated_HK@meta.data)
save(Mixed_3D.integrated_HK, file = "Mixed_3D.integrated_HK.NodeAligned.Robj")

# Define colors for CellClass
CellClass.cols <- c(
  "Epithelium" = "#D982C6",
  "Immune" = "#87B37A",
  "Mesenchyme" = "#F4A261",
  "Endothelium" = "#2A9D8F",
  "General" = "gray70")

# Define CellType color palette
CellType.cols <- c(
  "ATI_Like" = "#b22222","ATI" = "#b22222",
  "ATII" = "#4682B4","ATII_Like" = "#4682B4",
  "Secretory_Like" = "#2E8B57", "Secretory" = "#2E8B57",
  "Ciliated" = "#9B489B","Ciliated_Like" = "#9B489B",
  "Stressed_Progenitor" = "#c71585",
  "Hillock_Luminal" = "#ff0000","Hillock_Like" = "#ff0000",
  "Hillock_Basal" = "#2800c7",
  "Basal_Like" = "#c6e308", "Basal" = "#c6e308",
  "B" = "#41571b",
  "Mac_Alv" = "#00FBFF",
  "Mac_Inter" = "#FF481B",
  "gCaps" = "#A0522D",
  "Tuft" = "#15DDB5",
  "Mesothelium" = "#4A314D",
  "Neutrophils" = "#fcfd1d",
  "pDCs" = "#E9165D",
  "Monocytes" = "#3cd500",
  "NK" = "#F51BEA",
  "T" = "#95B8D1",
  "Arterial" = "#7E5109",
  "Venous" = "#660708",
  "RAS_Like" = "#595db0",
  "Cycling_Epithelium" = "#DB7093",
  "Pro_Inflamm_Mac" = "#9497fd",
  "Anti_Inflamm_Mac" = "#ff9e80","Polarized_Mac" = "#ff9e80",
  "Rspo3+_Mes" = "#89a5a5",
  "Cycling_Distal_Epi" = "#93d0ff",
  "Cycling_Proximal_Epi" = "#c4c2ff",
  "Pdgfrb+_Pericytes" = "#F4AFB4",
  "Cell_Cycle" = "#e0bf1b","Cycling_Immune" = "#e0bf1b",
  "Fzd7+_Stressed" = "#2F004F",
  "Epithelium" = "#D982C6") # Adding for BAL start epithelium (using same color as cellclass palette to avoid confusion)

# Set wd to where our AI file is (for linking)
setwd("~/Desktop/Manuscript Figures/Figure 2")
# Define the desired legend order
celltype_order <- c("ATI_Like", "ATII_Like","Basal", "Basal_Like", "Ciliated","Cycling_Distal_Epi",
                    "Cycling_Proximal_Epi","Hillock_Basal", "Hillock_Luminal", "RAS_Like", "Secretory", 
                    "Stressed_Progenitor", "Cycling_Immune", 
                    "Pro_Inflamm_Mac", "Anti_Inflamm_Mac", "Rspo3+_Mes", "Pdgfrb+_Pericytes")

table(BAL_3D.integrated_HK$CellType.NodeAligned.Condition)
# BAL_3D UMAP - Filter legend order to match only present cell types
BAL_present <- intersect(celltype_order, unique(BAL_3D.integrated_HK$CellType.NodeAligned.Condition))
BAL_3D.umap <- DimPlot(BAL_3D.integrated_HK, reduction = "umap.rpca", label = TRUE, 
                       cols = CellType.cols, group.by = "CellType.NodeAligned.Condition", repel = TRUE) + 
  scale_color_manual(values = CellType.cols, limits = BAL_present) +  
  theme(axis.title = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), 
        plot.title = element_blank(),
        legend.position = "none")  # <-- hide legend
BAL_3D.umap
ggsave("BAL_3D.umap.png", plot = BAL_3D.umap, width = 6, height = 6, dpi = 600)

# Mixed_3D UMAP - Filter legend order to match only present cell types
Mixed_present <- intersect(celltype_order, unique(Mixed_3D.integrated_HK$CellType.NodeAligned.Condition))
Mixed_3D.umap <- DimPlot(Mixed_3D.integrated_HK, reduction = "umap.rpca", label = TRUE, 
                         cols = CellType.cols, group.by = "CellType.NodeAligned.Condition", repel = TRUE) + 
  scale_color_manual(values = CellType.cols, limits = Mixed_present) +  
  theme(axis.title = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), 
        plot.title = element_blank(),
        legend.position = "none")  # <-- hide legend
Mixed_3D.umap
ggsave("Mixed_3D.umap.png", plot = Mixed_3D.umap, width = 6, height = 6, dpi = 600)

# PD_3D UMAP - Filter legend order to match only present cell types
PD_present <- intersect(celltype_order, unique(PD_3D.integrated_HK$CellType.NodeAligned.Condition))
PD_3D.umap <- DimPlot(PD_3D.integrated_HK, reduction = "umap.rpca", label = TRUE, 
                      cols = CellType.cols, group.by = "CellType.NodeAligned.Condition", repel = TRUE) + 
  scale_color_manual(values = CellType.cols, limits = PD_present) +  
  theme(axis.title = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), 
        plot.title = element_blank(),
        legend.position = "none")  # <-- hide legend
print(PD_3D.umap)
ggsave("PD_3D.umap.png", plot = PD_3D.umap, width = 6, height = 6, dpi = 600)

BAL_3D.umap
Mixed_3D.umap
PD_3D.umap
table(PD_3D.integrated_HK$CellType.NodeAligned.Condition)
table(Mixed_3D.integrated_HK$CellType.NodeAligned.Condition)
table(BAL_3D.integrated_HK$CellType.NodeAligned.Condition)

### Making custom global legend
observed_celltypes <- unique(c(
  levels(factor(BAL_3D.integrated_HK$CellType.NodeAligned.Condition)),
  levels(factor(PD_3D.integrated_HK$CellType.NodeAligned.Condition)),
  levels(factor(Mixed_3D.integrated_HK$CellType.NodeAligned.Condition))))
celltype_to_class <- c(
  # Epithelium
  "ATI_Like" = "Epithelium", "ATII_Like" = "Epithelium", "Secretory" = "Epithelium",
  "Ciliated" = "Epithelium", "Cycling_Proximal_Epi" = "Epithelium", "Cycling_Distal_Epi" = "Epithelium",
  "Stressed_Progenitor" = "Epithelium", "Hillock_Luminal" = "Epithelium",
  "Hillock_Basal" = "Epithelium", "Basal_Like" = "Epithelium", "RAS_Like" = "Epithelium",
  # Immune
  "Anti_Inflamm_Mac" = "Immune", "Pro_Inflamm_Mac" = "Immune", "Cycling_Immune" = "Immune",
  # Mesenchyme
  "Rspo3+_Mes" = "Mesenchyme", "Pdgfrb+_Pericytes" = "Mesenchyme")

# Build data drame
legend_df <- data.frame(
  CellType = observed_celltypes,
  Class = celltype_to_class[observed_celltypes],
  Color = CellType.cols[observed_celltypes]) %>%
  mutate(
    CellType = factor(CellType, levels = observed_celltypes),
    Class = factor(Class, levels = c("Epithelium", "Immune", "Mesenchyme"))) %>%
  arrange(Class, CellType)

# Order cells manually
epi_order <- c("Cycling_Distal_Epi","Cycling_Proximal_Epi","RAS_Like","Hillock_Luminal","Hillock_Basal","Stressed_Progenitor","Secretory","Ciliated","Basal_Like","ATII_Like","ATI_Like")
immune_order <- c("Cycling_Immune","Anti_Inflamm_Mac", "Pro_Inflamm_Mac")
mes_order <- c("Rspo3+_Mes", "Pdgfrb+_Pericytes")


epi_order <- c("Stressed_Progenitor", "Secretory", "RAS_Like", "Hillock_Luminal",
               "Hillock_Basal", "Basal_Like", "Ciliated", "Cycling_Proximal_Epi",
               "Cycling_Distal_Epi", "ATII_Like", "ATI_Like")
immune_order <- c("Cycling_Immune","Anti_Inflamm_Mac", "Pro_Inflamm_Mac")
mes_order <- c("Rspo3+_Mes", "Pdgfrb+_Pericytes")

celltype_order <- c(epi_order, immune_order, mes_order)
legend_df$CellType <- factor(legend_df$CellType, levels = celltype_order)

# plot
glbl.legend.eng <- ggplot(legend_df, aes(x = 1, y = CellType)) +
  geom_point(aes(color = Color), size = 5) +
  scale_color_identity() +
  facet_grid(Class ~ ., scales = "free_y", space = "free_y") +
  geom_text(aes(label = CellType), hjust = 0, nudge_x = 0.25, size = 4) +
  theme_void(base_size = 14) +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", size = 12),
    strip.background = element_blank(),
    plot.margin = margin(10, 40, 10, 10)) +
  xlim(1, 4)
glbl.legend.eng
ggsave("glbl.legend.eng.png", plot = glbl.legend.eng, width = 4, height = 5, dpi = 600)

