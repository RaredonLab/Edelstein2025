

# close all, clear all
graphics.off()  
rm(list = ls())

options(future.globals.maxSize = 16000 * 1024^2)

# Set seed
set.seed(2)

# Packages
library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(scales)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(cowplot)
library(patchwork)
library(SeuratWrappers)
library(reticulate)
library(stringr)
library(NICHES)

# Look at structure of data
str(eng.subset.integrated_HK@meta.data)

# Set idents to the desired meta-data label (we want this iteration of cell class labels)
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellClass.combined.Integrated
# Run a marker test using Seurat FindAllMarkers
eng.subset.int_lineage = FindAllMarkers(eng.subset.integrated_HK, 
                                  min.pct = 0.1,logfc.threshold = 0.1)
eng.subset.int_lineage$ratio = eng.subset.int_lineage$pct.1/eng.subset.int_lineage$pct.2
eng.subset.int_lineage$power = eng.subset.int_lineage$ratio*eng.subset.int_lineage$avg_log2FC
View(eng.subset.int_lineage)

# Take a look at the Ligands and Receptors that are built into NICHES
head(NICHES::ncomms8866_rat$Ligand.ApprovedSymbol)
head(NICHES::ncomms8866_rat$Receptor.ApprovedSymbol)

# Grab your ground truth ligand and receptor lists
ligands <- unique(NICHES::ncomms8866_rat$Ligand.ApprovedSymbol)
receptors <- unique(NICHES::ncomms8866_rat$Receptor.ApprovedSymbol)

## Filter marker list
# Filter ligands
lineage_ligands <- eng.subset.int_lineage[eng.subset.int_lineage$gene %in% ligands, ]
# Filter receptors
lineage_receptors <- eng.subset.int_lineage[eng.subset.int_lineage$gene %in% receptors, ]

# Split the ligand and receptor markers by cluster (lineage)
ligands_by_lineage <- split(lineage_ligands, lineage_ligands$cluster)
receptors_by_lineage <- split(lineage_receptors, lineage_receptors$cluster)
View(ligands_by_lineage)
View(receptors_by_lineage)

## Combine ligands and receptors into tidy data frames (so we can view more easily); also keeps all of the stats data
# Ligands
ligands_df <- bind_rows(ligands_by_lineage, .id = "Lineage")
# Receptors
receptors_df <- bind_rows(receptors_by_lineage, .id = "Lineage")
View(ligands_df)
View(receptors_df)

# Look at differences between basal and hillock; one ident set to Basal_like; the other to hillock
Idents(eng.subset.integrated_HK) = eng.subset.integrated_HK$CellType.combined.Integrated
eng.subset.int_bas.hill = FindMarkers(eng.subset.integrated_HK, ident.1 = "Basal_Like", ident.2 = "Hillock_Like",
                                        min.pct = 0.1,logfc.threshold = 0.1)
eng.subset.int_bas.hill$ratio = eng.subset.int_bas.hill$pct.1/eng.subset.int_bas.hill$pct.2
eng.subset.int_bas.hill$power = eng.subset.int_bas.hill$ratio*eng.subset.int_bas.hill$avg_log2FC
eng.subset.int_bas.hill$gene <- rownames(eng.subset.int_bas.hill)
View(eng.subset.int_bas.hill)
# High power = high in basal; low = high in hillock

# Filter by ligand/receptor
ligands <- unique(NICHES::ncomms8866_rat$Ligand.ApprovedSymbol)
receptors <- unique(NICHES::ncomms8866_rat$Receptor.ApprovedSymbol)

# Filter for ligand matches
ligand_markers_bas.hill <- eng.subset.int_bas.hill %>%
  filter(gene %in% ligands)
# Filter for receptor matches
receptor_markers_bas.hill <- eng.subset.int_bas.hill %>%
  filter(gene %in% receptors)
View(ligand_markers_bas.hill)
# Annotating for clarity
# Ligands
ligand_markers_bas.hill <- ligand_markers_bas.hill %>%
  mutate(Profile = if_else(avg_log2FC > 0, "Basal", "Hillock"))
# Receptors
receptor_markers_bas.hill <- receptor_markers_bas.hill %>%
  mutate(Profile = if_else(avg_log2FC > 0, "Basal", "Hillock"))
# Take a look at all
View(ligand_markers_bas.hill)
View(receptor_markers_bas.hill)
head(ligand_markers_bas.hill)
head(receptor_markers_bas.hill)

# Split one of our NICHES outputs - ligand receptor separate
eng.CTC.final.mark <- eng.CTC.final.mark %>%
  separate(gene, into = c("Ligand", "Receptor"), sep = "—")
# Ligands
basal_ligands <- ligand_markers_bas.hill %>%
  filter(Profile == "Basal") %>%
  pull(gene)
hillock_ligands <- ligand_markers_bas.hill %>%
  filter(Profile == "Hillock") %>%
  pull(gene)
# Receptors
basal_receptors <- receptor_markers_bas.hill %>%
  filter(Profile == "Basal") %>%
  pull(gene)
hillock_receptors <- receptor_markers_bas.hill %>%
  filter(Profile == "Hillock") %>%
  pull(gene)

### Filter NICHES
# Basal-enriched interactions: Ligand or Receptor is Basal-specific
basal_lig.rec.interactions.filtered <- hill.mes.imm.sub.mark.int_joint %>%
  filter(Ligand %in% basal_ligands | Receptor %in% basal_receptors) %>%
  mutate(Profile = "Basal")
# Hillock-enriched interactions
hillock_lig.rec.interactions.filtered <- hill.mes.imm.sub.mark.int_joint %>%
  filter(Ligand %in% hillock_ligands | Receptor %in% hillock_receptors) %>%
  mutate(Profile = "Hillock")
View(basal_lig.rec.interactions.filtered)
View(hillock_lig.rec.interactions.filtered)

# Bring in broader NICHES object (whole thing)
str(eng.CTC.final@meta.data)
# Plot example
DimPlot(eng.CTC.final,group.by = 'integrated.clusters.2',shuffle=T, raster = F, cols = full.col.pal,label = T, reduction = 'umap.rpca.2') # Embedding integrated by sample
DimPlot(eng.CTC.final,group.by = 'CellClass.combined.Integrated.Sending',shuffle=T, raster = F, cols = full.col.pal,label = T, reduction = 'umap.rpca.2') # Embedding integrated by sample
DimPlot(eng.CTC.final,group.by = 'CellClass.combined.Integrated.Receiving',shuffle=T, raster = F, cols = full.col.pal,label = T, reduction = 'umap.rpca.2') # Embedding integrated by sample

DimPlot(eng.CTC.final,group.by = 'CellType.combined.Integrated.Sending',shuffle=T, raster = F, cols = full.col.pal,label = T, reduction = 'umap.rpca.2') # Embedding integrated by sample
DimPlot(eng.CTC.final,group.by = 'CellType.combined.Integrated.Receiving',shuffle=T, raster = F, cols = full.col.pal,label = T, reduction = 'umap.rpca.2') # Embedding integrated by sample

# Run marker test on the whole NICHES object
Idents(eng.CTC.final) = eng.CTC.final$integrated.clusters.2
eng.CTC.final.mark = FindAllMarkers(eng.CTC.final, 
                                  min.pct = 0.1,logfc.threshold = 0.1)
eng.CTC.final.mark$ratio = eng.CTC.final.mark$pct.1/eng.CTC.final.mark$pct.2
eng.CTC.final.mark$power = eng.CTC.final.mark$ratio*eng.CTC.final.mark$avg_log2FC
View(eng.CTC.final.mark)

# Sorting based on our basal/hillock receptor ligand lists generated from the transcriptomic object
# Now filtering our NICHES marker list using it 
basal_lig.rec.interactions.filtered_glbl <- eng.CTC.final.mark %>%
  filter(Ligand %in% basal_ligands | Receptor %in% basal_receptors) %>%
  mutate(Profile = "Basal")
# Hillock-enriched interactions
hillock_lig.rec.interactions.filtered_glbl <- eng.CTC.final.mark %>%
  filter(Ligand %in% hillock_ligands | Receptor %in% hillock_receptors) %>%
  mutate(Profile = "Hillock")
View(basal_lig.rec.interactions.filtered_glbl)
View(hillock_lig.rec.interactions.filtered_glbl)

### Checking expression at transcriptomic levels and on NICHES level (eng.subset.integrated_HK = transcriptomic; eng.CTC.final = NICHES)
# Looking at cluster 8 (integrated.clusters.2) which seems to be specific to hillock sending, immune receiving (not whole cluster though)
FeaturePlot(eng.subset.integrated_HK, features = c("Gnai2","C5ar1"), reduction = "umap.rpca")
FeaturePlot(eng.CTC.final, features = c(""), reduction = "umap.rpca.2")

Idents(eng.CTC.final) = eng.CTC.final$integrated.clusters.2
VlnPlot(eng.CTC.final, features = c("Sfrp1—Fzd6"))
VlnPlot(eng.CTC.final, features = c("Plau—St14"))

# Basal vs. ALL
Idents(eng.subset.integrated_HK) <- eng.subset.integrated_HK$CellType.combined.Integrated
# Basal vs all other cell types
basal_vs_all <- FindMarkers(
  eng.subset.integrated_HK,
  ident.1 = "Basal_Like",
  ident.2 = NULL,  # NULL means vs. all other cells
  min.pct = 0.1,
  logfc.threshold = 0.1)

# Hillock vs all
hillock_vs_all <- FindMarkers(
  eng.subset.integrated_HK,
  ident.1 = "Hillock_Like",
  ident.2 = NULL,
  min.pct = 0.1,
  logfc.threshold = 0.1)

# Add gene names as column for filtering purposes
basal_vs_all$gene <- rownames(basal_vs_all)
hillock_vs_all$gene <- rownames(hillock_vs_all)
ligands <- unique(NICHES::ncomms8866_rat$Ligand.ApprovedSymbol)
receptors <- unique(NICHES::ncomms8866_rat$Receptor.ApprovedSymbol)
# Ligands specific to Basal
basal_ligands_specific <- basal_vs_all %>%
  filter(gene %in% ligands, avg_log2FC > 0)
# Receptors specific to Basal
basal_receptors_specific <- basal_vs_all %>%
  filter(gene %in% receptors, avg_log2FC > 0)
# Ligands specific to Hillock
hillock_ligands_specific <- hillock_vs_all %>%
  filter(gene %in% ligands, avg_log2FC > 0)
# Receptors specific to Hillock
hillock_receptors_specific <- hillock_vs_all %>%
  filter(gene %in% receptors, avg_log2FC > 0)
View(basal_ligands_specific)
View(basal_receptors_specific)
View(hillock_ligands_specific)
View(hillock_receptors_specific)

# Save outputs
setwd("~/Desktop/Manuscript Figures/Basal vs. Hillock LR Analysis for MSBR")
write_xlsx(basal_ligands_specific, "basal_ligands_specific.xlsx")
write_xlsx(basal_receptors_specific, "basal_receptors_specific.xlsx")
write_xlsx(hillock_ligands_specific, "hillock_ligands_specific.xlsx")
write_xlsx(hillock_receptors_specific, "hillock_receptors_specific.xlsx")

Idents(eng.CTC.final) = eng.CTC.final$CellType.combined.Integrated.Joint
eng.subset.int_joint_CTC_hill.mac = FindMarkers(eng.CTC.final, ident.1 = "Hillock_Like - Polarized_Mac",
                                       ident.2 = NULL,
                                        min.pct = 0.1,logfc.threshold = 0.1)
eng.subset.int_joint_CTC_hill.mac$ratio = eng.subset.int_joint_CTC_hill.mac$pct.1/eng.subset.int_joint_CTC_hill.mac$pct.2
eng.subset.int_joint_CTC_hill.mac$power = eng.subset.int_joint_CTC_hill.mac$ratio*eng.subset.int_joint_CTC_hill.mac$avg_log2FC
View(eng.subset.int_joint_CTC_hill.mac)


