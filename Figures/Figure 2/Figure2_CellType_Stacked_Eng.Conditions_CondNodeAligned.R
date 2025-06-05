
# Stacked bar plots for engineered conditions for Figure 2

library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

# Load in objects from folder on local machine (SEE)
setwd("~/Desktop/Datasets/Condition Level Objects")
load("BAL_3D.integrated_HK.NodeAligned.Robj") # BAL_3D
load("PD_3D.integrated_HK.NodeAligned.Robj") # PD_3D
load("Mixed_3D.integrated_HK.NodeAligned.Robj") # Mixed_3D

table(BAL_3D.integrated_HK$CellType.NodeAligned)
table(PD_3D.integrated_HK$CellType.NodeAligned)
table(Mixed_3D.integrated_HK$CellType.NodeAligned)

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
  "Epithelium" = "#D982C6") 

table(PD_3D.integrated_HK$CellType.NodeAligned.Condition)
table(Mixed_3D.integrated_HK$CellType.NodeAligned.Condition)
table(BAL_3D.integrated_HK$CellType.NodeAligned.Condition)

# Define the mapping of cell types to lineages
# Define the mapping of cell types to lineages
lineage_map <- list(
  Immune = c("Pro_Inflamm_Mac", "Anti_Inflamm_Mac", "Cycling_Immune"),
  Mesenchyme = c("Pdgfrb+_Pericytes", "Rspo3+_Mes"),
  Epithelium = c("ATI_Like", "ATII_Like", "Stressed_Progenitor", "Hillock_Basal", 
                 "Hillock_Luminal", "RAS_Like", "Cycling_Distal_Epi", "Cycling_Proximal_Epi",
                 "Basal_Like", "Ciliated", "Secretory"))

# Combine into one named vector for mapping
celltype_to_lineage <- unlist(lapply(names(lineage_map), function(lin) {
  setNames(rep(lin, length(lineage_map[[lin]])), lineage_map[[lin]])
}))

# Function to extract and annotate counts using CellType.NodeAligned.Condition
extract_counts <- function(seurat_obj, object_name) {
  celltypes <- seurat_obj$CellType.NodeAligned.Condition
  df <- as.data.frame(table(celltypes))
  colnames(df) <- c("CellType", "Count")
  df$Object <- object_name
  df$Lineage <- celltype_to_lineage[as.character(df$CellType)]
  df <- df %>% filter(!is.na(Lineage))
  return(df)
}

# Apply function to all three Seurat objects
df_BAL <- extract_counts(BAL_3D.integrated_HK, "BAL_3D")
df_PD <- extract_counts(PD_3D.integrated_HK, "PD_3D")
df_Mixed <- extract_counts(Mixed_3D.integrated_HK, "Mixed_3D")

# Combine into one dataframe
combined_df <- bind_rows(df_BAL, df_PD, df_Mixed)

# Normalize to get percentage within each lineage per object
final_df <- combined_df %>%
  group_by(Object, Lineage) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

# Set factor levels for consistent plotting
final_df$Object <- factor(final_df$Object, levels = c("BAL_3D", "PD_3D", "Mixed_3D"))
final_df$Lineage <- factor(final_df$Lineage, levels = c("Epithelium", "Immune", "Mesenchyme"))

# Order CellType by overall abundance
final_df$CellType <- factor(final_df$CellType, levels = final_df %>%
                              group_by(CellType) %>%
                              summarise(n = sum(Count)) %>%
                              arrange(-n) %>%
                              pull(CellType))

# Plotting function by lineage
plot_lineage <- function(lineage_name) {
  df_subset <- final_df %>% filter(Lineage == lineage_name)
  
  ggplot(df_subset, aes(x = Object, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    scale_fill_manual(values = CellType.cols) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste(lineage_name),
      x = NULL,
      y = "Cell Type Contribution to Class",
      fill = "Cell Type") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 12),
      plot.title = element_text(face = "bold", size = 16),
      legend.position = "bottom",
      legend.box.margin = margin(t = 10)
    )
}

# Generate and print plots
plot_epi <- plot_lineage("Epithelium")
plot_imm <- plot_lineage("Immune")
plot_mes <- plot_lineage("Mesenchyme")
print(plot_epi)
print(plot_imm)
print(plot_mes)

# Lineage contribution
lineage_summary <- final_df %>%
  group_by(Object, Lineage) %>%
  summarise(LineageCount = sum(Count), .groups = "drop") %>%
  group_by(Object) %>%
  mutate(LineageProportion = LineageCount / sum(LineageCount)) %>%
  ungroup()

lineage.contr.system <- ggplot(lineage_summary, aes(x = Object, y = LineageProportion, fill = Lineage)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +  # no percent_format!
  scale_fill_manual(values = c(
    "Epithelium" = "#D982C6",
    "Immune" = "#87B37A",
    "Mesenchyme" = "#F4A261")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Cell Class",
    x = NULL,
    y = "Cell Class Contribution to Organoid System",
    fill = "Lineage") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1))
print(lineage.contr.system)

eng.stacked = lineage.contr.system | plot_epi | plot_imm | plot_mes
eng.stacked
setwd("~/Desktop/Manuscript Figures/Figure 2")
ggsave("EngineeredObjCellTypeBreakdown_byLineage.png", plot = eng.stacked, width = 20, height = 7, dpi = 600)
