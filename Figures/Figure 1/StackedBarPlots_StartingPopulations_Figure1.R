
## MEMO: Stacked bar plots for engineered conditions for Figure 1 (starting populations)

# close all, clear all
graphics.off()  
rm(list = ls())

# Load packages
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

# Set wd to be true and load datasets
setwd("~/Desktop/Datasets/Starting Populations")
load("BAL.integrated_03-29-2025.Robj") # BAL Start
load("Mixed_3D.start.combined.Robj") # Mixed Pseudo Start
load("PD.integrated_03-29-2025.Robj") # PD Start

# Get a sense of the objects
DimPlot(Mixed_3D.start.combined)
table(Mixed_3D.start.combined$CellType)
table(BAL.integrated$CellType)
table(PD.integrated$CellType)

## Fix LEP annotation in BAL_3D
# Change "LEPs" to "Epithelium" in CellType metadata
BAL.integrated$CellType[BAL.integrated$CellType == "LEPs"] <- "Pan_Epithelium"
# Confirm the change
table(BAL.integrated$CellType)

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
  "Rspo3+_Mes" = "#89a5a5","Actc1_Mural" = "#1179fa",
  "Cycling_Distal_Epi" = "#93d0ff",
  "Cycling_Proximal_Epi" = "#c4c2ff",
  "Pdgfrb+_Pericyte" = "#F4AFB4",
  "Cell_Cycle" = "#e0bf1b","Cycling_Immune" = "#e0bf1b",
  "Fzd7+_Stressed" = "#2F004F",
  "Pan_Epithelium" = "#E2AEDD")

table(BAL.integrated$CellType)
table(Mixed_3D.start.combined$CellType)
table(PD.integrated$CellType)

# Define the mapping of cell types to lineages
lineage_map <- list(
  Immune = c("B","T","Mac_Alv","Mac_Inter","Monocytes","Neutrophils","NK","pDCs"),
  Mesenchyme = c("Actc1_Mural", "Mesothelium"),
  Epithelium = c("Secretory","ATII","Tuft","Basal","Ciliated","ATI","Pan_Epithelium"),
  Endothelium = c("gCaps","Arterial","Venous"),
  General = c("Cell_Cycle"))

# Combine into one named vector for mapping
celltype_to_lineage <- unlist(lapply(names(lineage_map), function(lin) {
  setNames(rep(lin, length(lineage_map[[lin]])), lineage_map[[lin]])
}))

# Function to extract and annotate counts from a Seurat object
extract_counts <- function(seurat_obj, object_name) {
  celltypes <- seurat_obj$CellType
  df <- as.data.frame(table(celltypes))
  colnames(df) <- c("CellType", "Count")
  df$Object <- object_name
  df$Lineage <- celltype_to_lineage[as.character(df$CellType)]
  df <- df %>% filter(!is.na(Lineage))
  df
}

# Apply function to all three Seurat objects
df_BAL <- extract_counts(BAL.integrated, "BAL")
df_PD <- extract_counts(PD.integrated, "PD")
df_Mixed <- extract_counts(Mixed_3D.start.combined, "BAL + PD")

# Combine into one dataframe
combined_df <- bind_rows(df_BAL, df_PD, df_Mixed)

# Normalize to get percentage within each lineage per object
final_df <- combined_df %>%
  group_by(Object, Lineage) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()
# T0 help with normalization
final_df$CellType <- factor(final_df$CellType, levels = final_df %>%
                              group_by(CellType) %>%
                              summarise(n = sum(Count)) %>%
                              arrange(-n) %>%
                              pull(CellType))
# Preview
print(head(final_df))
View(final_df)

# First, let's make three plots by lineage where each bar is the condition (stacked breakdown)
# Ensure factor levels for consistent order (optional)
final_df$Object <- factor(final_df$Object, levels = c("BAL", "BAL + PD","PD"))
final_df$Lineage <- factor(final_df$Lineage, levels = c("Epithelium", "Immune", "Endothelium","Mesenchyme","General"))

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
      legend.box.margin = margin(t = 10))}

# Create and print each plot
plot_epi <- plot_lineage("Epithelium")
plot_imm <- plot_lineage("Immune")
plot_mes <- plot_lineage("Mesenchyme")
plot_endo <- plot_lineage("Endothelium")

print(plot_epi)
print(plot_imm)
print(plot_mes)
print(plot_endo)

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
    "Mesenchyme" = "#F4A261",
    "Endothelium" = "#2A9D8F",
    "General" = "gray70")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Cell Class",
    x = NULL,
    y = "Cell Class Contribution to Starting Population",
    fill = "Lineage") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1))
print(lineage.contr.system)

# Put plots together
start.stacked = lineage.contr.system | plot_epi | plot_imm | plot_mes
# View aggregated plots
start.stacked
# Save
setwd("~/Desktop/Manuscript Figures/Figure 1")
ggsave("StartingPopulationsContributions_byLineage.png", plot = start.stacked, width = 20, height = 7, dpi = 600)
ggsave("plot_epi.start.png", plot = plot_epi, width = 6, height = 6, dpi = 600)

