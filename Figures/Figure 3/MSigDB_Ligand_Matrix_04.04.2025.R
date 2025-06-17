
# Gene set enrichment analysis for the Proximal Subset (SEE's data, by SEE)
# Date: 04.04.2025

library(Seurat)
library(ggplot2)
library(purrr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(circlize)
library(viridis)
library(devtools)
library(msigdbr) # version is 7.5.1 which is MSigDB Broad 2022
library(clusterProfiler)
library(forcats)
library(writexl)
library(NICHES)
library(janitor)

# Load data
# Can load both integrated and non-integrated objects
setwd("~/Desktop/Single Cell/Engineered Subset (Proximal) for Pseudotime")
load("eng.subset.prox.int.Robj")
load("eng.subset.prox.Robj")

# Load all MSigDB gene sets for rat
msigdb_rat <- msigdbr(species = "Rattus norvegicus")
# Filter for desired categories
msig_h <- msigdb_rat %>% filter(gs_collection == "H") # Hallmark paths
msig_c2 <- msigdb_rat %>% filter(gs_collection == "C2") # Curated pathways
msig_reactome <- msigdb_rat %>%
  filter(gs_collection == "C2", gs_subcollection == "CP:REACTOME") # reactome
msig_kegg <- msigdb_rat %>%
  filter(gs_collection == "C2", gs_subcollection == "CP:KEGG") # KEGG pathways

# Maybe subset out JUST the hillock lum/basal from unintegrated object
hillock_cells <- subset(eng.subset.prox, 
                        subset = CellType_ProxSubset_no.int %in% c("Hillock_Luminal", "Hillock_Basal"))

# 1. Normalize the data
hillock_cells <- NormalizeData(hillock_cells)
# 2. Find variable features
hillock_cells <- FindVariableFeatures(hillock_cells)
# 3. Scale the data (use all genes or just variable ones)
hillock_cells <- ScaleData(hillock_cells)
# 4. Run PCA
hillock_cells <- RunPCA(hillock_cells)
DimHeatmap(hillock_cells, dims = 1:9, cells = 500, balanced = TRUE)
DimHeatmap(hillock_cells, dims = 10:18, cells = 500, balanced = TRUE)

# 5. Run UMAP (optional, but useful)
hillock_cells <- RunUMAP(hillock_cells, dims = 1:8)
# 6. Run clustering (if you want)
hillock_cells <- FindNeighbors(hillock_cells, dims = 1:8)
hillock_cells <- FindClusters(hillock_cells, resolution = 0.5)
DimPlot(hillock_cells)
FeaturePlot(hillock_cells, features = c("Krt13"), order = T) 
FeaturePlot(hillock_cells, features = c("Krt5"), order = T)

# Getting rid of 7 which is pure basal
hillock.sub = subset(hillock_cells, idents = c("7"),invert = TRUE)

# Stash clusters
hillock.sub$stash = Idents(hillock.sub)
# 1. Normalize the data
hillock.sub <- NormalizeData(hillock.sub)
# 2. Find variable features
hillock.sub <- FindVariableFeatures(hillock.sub)
# 3. Scale the data (use all genes or just variable ones)
hillock.sub <- ScaleData(hillock.sub)
# 4. Run PCA
hillock.sub <- RunPCA(hillock.sub)
DimHeatmap(hillock.sub, dims = 1:9, cells = 500, balanced = TRUE)
DimHeatmap(hillock.sub, dims = 10:18, cells = 500, balanced = TRUE)
hillock.sub <- RunUMAP(hillock.sub, dims = 1:12)
# 6. Run clustering (if you want)
hillock.sub <- FindNeighbors(hillock.sub, dims = 1:12)
hillock.sub <- FindClusters(hillock.sub, resolution = 0.5)
DimPlot(hillock.sub)
FeaturePlot(hillock.sub, features = c("Krt13","Sprr1a"), order = T) 
FeaturePlot(hillock.sub, features = c("Krt5","Tp63","Maff","Mal"))

# Stash clusters
hillock.sub$stash_hillock = Idents(hillock.sub)
# Create another lineage just for the hell of it
cell.epi = WhichCells(hillock.sub, idents = c(0,1,2,3,4,5,6,7))
hillock.sub = SetIdent(hillock.sub, cells = cell.epi, value = 'Epithelium')
# Stash the cell class
hillock.sub$CellClass.hillock = Idents(hillock.sub)
table(hillock.sub$CellClass.hillock)
# confirming that everything is labeled
sum(is.na(hillock.sub$CellClass.hillock))

# Restores the cell identities that were previously stashed or stored elsewhere.
Idents(hillock.sub) = hillock.sub$stash_hillock
table(Idents(hillock.sub))
FeaturePlot(hillock.sub, features = c("Krt5","Tp63"), split.by = "Condition")
DimPlot(hillock.sub, split.by = "Condition")
VlnPlot(hillock.sub, features = c("Krt5","Tp63","Maff","Mal","Krt13","Aqp3","Col17a1","Krt14"))
FeaturePlot(hillock.sub, features = c("Krt5","Tp63","Maff","Mal"), label = T) 
Idents(hillock.sub) = hillock.sub$Condition
VlnPlot(hillock.sub, features = c("Krt5","Tp63","Krt80","Mal","Slpi","Cnfn","S100a7a","S100a8"))

Idents(hillock.sub) = hillock.sub$stash_hillock
DimPlot(hillock.sub, label = T)
# Create empty meta-data slot
hillock.sub[["CellType.Hillock"]] <- rep(NA, ncol(hillock.sub))
# Cell typing
c0 <- WhichCells(hillock.sub, idents = "0")
c1 <- WhichCells(hillock.sub, idents = "1")
c2 <- WhichCells(hillock.sub, idents = "2")
c3 <- WhichCells(hillock.sub, idents = "3")
c4 <- WhichCells(hillock.sub, idents = "4")
c5 <- WhichCells(hillock.sub, idents = "5")
c6 <- WhichCells(hillock.sub, idents = "6")
c7 <- WhichCells(hillock.sub, idents = "7")
hillock.sub$CellType.Hillock[c0] <- 'Hillock_Basal'
hillock.sub$CellType.Hillock[c1] <- 'Hillock_Luminal'
hillock.sub$CellType.Hillock[c2] <- 'Hillock_Luminal'
hillock.sub$CellType.Hillock[c3] <- 'Hillock_Luminal'
hillock.sub$CellType.Hillock[c4] <- 'Hillock_Luminal'
hillock.sub$CellType.Hillock[c5] <- 'Hillock_Luminal'
hillock.sub$CellType.Hillock[c6] <- 'Hillock_Luminal'
hillock.sub$CellType.Hillock[c7] <- 'Hillock_Luminal'

# Check the distribution to confirm all annotations were applied
table(hillock.sub$CellType.Hillock)
View(hillock.sub@meta.data)
sum(is.na(hillock.sub$CellType.Hillock))
# Looks good
setwd("~/Desktop/Manuscript Figures/Figure 3")
save(hillock.sub, file = "hillock.sub.Robj")

# Create DGE mark
Idents(eng.subset.prox.int) = eng.subset.prox.int$Condition
dge_hillock.pol = FindMarkers(eng.subset.prox.int, ident.1 = "BAL_3D", 
                              ident.2 = c("Mixed_3D", "PD_3D"), 
                              min.pct = 0.1,logfc.threshold = 0.1)
dge_hillock.pol$ratio = dge_hillock.pol$pct.1/dge_hillock.pol$pct.2
dge_hillock.pol$power = dge_hillock.pol$ratio*dge_hillock.pol$avg_log2FC
View(dge_hillock.pol)

# Create DGE mark
Idents(hillock.sub) = hillock.sub$CellType.Hillock
hillock.pol_DGE = FindAllMarkers(hillock.sub, 
                              min.pct = 0.1,logfc.threshold = 0.1)
hillock.pol_DGE$ratio = hillock.pol_DGE$pct.1/hillock.pol_DGE$pct.2
hillock.pol_DGE$power = hillock.pol_DGE$ratio*hillock.pol_DGE$avg_log2FC
View(hillock.pol_DGE)

# Use log2FC as ranking metric - we are using our DGE list that we generated by doing Conditions (BAL_3D (Basal) vs. PD_3D and Mixed_3D (Luminal))
gene_list <- hillock.pol_DGE$avg_log2FC
# Set names to gene symbols
names(gene_list) <- rownames(hillock.pol_DGE)
# Remove NA values if present
gene_list <- gene_list[!is.na(gene_list)]
# Sort in decreasing order for GSEA
gene_list <- sort(gene_list, decreasing = TRUE)

# Use the msig_h gene sets (prepared earlier with gs_collection == "H")
msigdb_list <- msig_h %>%
  select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)
msigdb_df <- msig_h %>%
  dplyr::select(gs_name, gene_symbol) %>%
  as.data.frame()

# Run GSEA
gsea_msig_hillock <- GSEA(
  geneList = gene_list,
  TERM2GENE = msigdb_df,
  pvalueCutoff = 0.05,
  verbose = FALSE)
head(gsea_msig_hillock@result, 10)
View(gsea_msig_hillock@result)

# B/c dge_results came from Hillock Basal vs. Hillock Luminal
# A positive NES (normalized enrichment score) means the pathway is enriched in Hillock Basal
# A negative NES means the pathway is enriched in Hillock Luminal

# Filter data
gsea_df_hillock.hallmark_cell.type <- gsea_msig_hillock@result %>% # pull the GSEA results table (@result) from gsea_msig object and pipe into dplyr for filtering
  mutate(
    EnrichedIn = ifelse(NES > 0, "PD_Mixed.Luminal_Hillock", "BAL_3D.Basal_Hillock"),
    log10p = -log10(p.adjust) # Adds a new column called EnrichedIn based on the Normalized Enrichment Score (NES)
  ) %>% # NES > 0, the pathway is enriched in Hillock Basal, NES < 0, the pathway is enriched in Hillock Luminal
  arrange(desc(log10p)) %>% # Add a column log10p, which is the -log10 of the adjusted p-value
  filter(p.adjust < 0.05)  # sort so most significant (lowest p-adjust) are first
View(gsea_df_hillock.hallmark_cell.type)
# Save as excel sheet
setwd("~/Desktop/Single Cell/Hillock_Pathway_Analysis")
write_xlsx(gsea_df_hillock.hallmark_cell.type, "gsea_BAL3D_vs_PD_Mixed_filtered.celltype_HALLMARK_MSigDB.xlsx")

# Trying to sort this data based on what Sam wants me to do
# So here is our working data-frame
gsea_df_hillock.hallmark_cell.type

# 1. Load Ground Truth Gene Sets
## 1a. Ligands from NICHES
ligand_genes <- unique(ncomms8866_rat$Ligand.ApprovedSymbol)

## 1b. Matrix proteins from MatrisomeDB (mouse, used for rat)
# Replace with your actual path to the Excel file
setwd("~/Desktop/Single Cell/Hillock_Pathway_Analysis/Ground Truths")
matrisome_df <- read_excel("Mm_Matrisome_Masterlist_Naba et al_2012.xlsx")
# Just pull gene names
# Matrisome includes more than just matrix stuff - includes ECM regulators, secreted factors, etc. but we don't want those right now
ecm_categories <- c("ECM Glycoproteins", "Collagens", "Proteoglycans", "ECM Regulators")
# Filter matrisome_df
filtered_matrisome_df <- matrisome_df %>%
  filter(`Matrisome Category` %in% ecm_categories)
# Extract gene symbols
matrix_genes <- unique(filtered_matrisome_df$`Gene Symbol`)

## 1c. Transcription Factors from AnimalTFDB (rat-specific)
# Replace with your actual path
tf_df <- read.delim("Rattus_norvegicus_TF.txt")  # from AnimalTFDB
# Pull genes
tf_genes <- unique(tf_df$Symbol)

# === 2. Annotate GSEA Results ===
gsea_df_annotated <- gsea_df_hillock.hallmark_cell.type
# Initialize columns
gsea_df_annotated$Ligands_NICHES <- NA
gsea_df_annotated$Ligand_NICHES_Count <- NA
gsea_df_annotated$Matrix_proteins <- NA
gsea_df_annotated$Matrix_protein_Count <- NA
gsea_df_annotated$TFs_present <- NA
gsea_df_annotated$TF_Count <- NA

# Loop through each row
for (i in seq_len(nrow(gsea_df_annotated))) {
  core_genes <- unlist(strsplit(gsea_df_annotated$core_enrichment[i], "/"))
  
  # Ligands (excluding matrix)
  ligands_present <- setdiff(intersect(core_genes, ligand_genes), matrix_genes)
  gsea_df_annotated$Ligands_NICHES[i] <- paste(ligands_present, collapse = ", ")
  gsea_df_annotated$Ligand_NICHES_Count[i] <- length(ligands_present)
  
  # Matrix proteins
  matrix_present <- intersect(core_genes, matrix_genes)
  gsea_df_annotated$Matrix_proteins[i] <- paste(matrix_present, collapse = ", ")
  gsea_df_annotated$Matrix_protein_Count[i] <- length(matrix_present)
  
  # Transcription factors
  tf_present <- intersect(core_genes, tf_genes)
  gsea_df_annotated$TFs_present[i] <- paste(tf_present, collapse = ", ")
  gsea_df_annotated$TF_Count[i] <- length(tf_present)
}
View(gsea_df_annotated)
# Now save
write_xlsx(gsea_df_annotated, "Hillock_MSigDB_Hallmark_Annotated_with_Ligands_Matrix_TFs.xlsx")

# Filter for separate plots
top_n <- 20
gsea_basal <- gsea_df_annotated %>% 
  filter(EnrichedIn == "BAL_3D.Basal_Hillock") %>%
  slice_max(log10p, n = top_n)

gsea_luminal <- gsea_df_annotated %>% 
  filter(EnrichedIn == "PD_Mixed.Luminal_Hillock") %>%
  slice_max(log10p, n = top_n)

# Plotting 
ggplot(gsea_luminal_filtered, aes(x = log10p, y = reorder(Description, log10p))) +
  geom_col(fill = "#ff0000") +
  geom_text(aes(x = 0.05, label = Description), 
            hjust = 0, size = 3.5, color = "black", fontface = "bold") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    title = "Luminal_Hillock",
    x = expression(-log[10]~adjusted~P~value),
    y = "MSigDB_Hallmark_2022") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(size = 1, color = "black"),
    axis.line.y = element_line(size = 1, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(size = 10))

############### Plotting as mirrors ############### 

selected_pathways <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_UV_RESPONSE_DN",
  "HALLMARK_HYPOXIA",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_ADIPOGENESIS",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_APICAL_JUNCTION",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_XENOBIOTIC_METABOLISM")

# Fix the EnrichedIn labels if needed:
gsea_df_annotated$EnrichedIn <- case_when(
  gsea_df_annotated$NES < 0 ~ "Hillock Basal",
  gsea_df_annotated$NES > 0 ~ "Hillock Luminal",
  TRUE ~ "NA")

gsea_mirror_df <- gsea_df_annotated %>%
  filter(Description %in% selected_pathways) %>%
  mutate(
    log10p_signed = ifelse(EnrichedIn == "Hillock Luminal", -log10p, log10p),
    Group = EnrichedIn,
    Label_Pos = ifelse(log10p_signed > 0, 0.05, -0.05),
    Description = fct_reorder(Description, log10p_signed))

# Plot: Pathways on Y-axis, enrichment on X
hallmark_hillock.mirror = ggplot(gsea_mirror_df, aes(x = log10p_signed, y = Description, fill = EnrichedIn)) +
  geom_col() +
  geom_text(aes(x = Label_Pos, label = Description),
            hjust = ifelse(gsea_mirror_df$log10p_signed > 0, 0, 1),
            size = 2.5, fontface = "bold") +
  scale_fill_manual(values = c("Hillock Luminal" = "#ff0000", "Hillock Basal" = "#2800c7")) +
  labs(
    x = expression(-log[10]~adjusted~P~value),
    y = "MSigDB_Hallmark_2022",
    fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),  # Labels are inside the plot now
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "bottom")
print(hallmark_hillock.mirror)
setwd("~/Desktop/Manuscript Figures/Figure 3")
ggsave("Hillock.MSigDBHallmarks_Mirror_Plot.png", plot = hallmark_hillock.mirror,
       width = 6, height = 5, dpi = 600)

# Let's try and generate lists for upregulated ligands, receptors, matrix proteins, etc in our DGE list
hillock_TFs <- hillock.pol_DGE %>%
  filter(gene %in% tf_genes)
### Matrix
# Define ECM categories of interest
ecm_categories <- c("ECM Glycoproteins", "Collagens", "Proteoglycans", "ECM Regulators")
# Filter to only those
filtered_matrisome_df <- matrisome_df %>%
  filter(`Matrisome Category` %in% ecm_categories)
# Rename for easier joining
filtered_matrisome_df <- filtered_matrisome_df %>%
  rename(gene = `Gene Symbol`)
hillock_matrix <- hillock.pol_DGE %>%
  filter(gene %in% matrix_genes)
hillock_matrix_annotated <- hillock_matrix %>%
  left_join(filtered_matrisome_df, by = "gene")
write_xlsx(hillock_matrix_annotated, "hillock_matrix_annotated.xlsx")

# Ligands (that are not in matrix)
# Unique ligands
ncomms8866_rat <- ncomms8866_rat %>%
  mutate(Ligand.ApprovedSymbol = toupper(Ligand.ApprovedSymbol))
hillock.pol_DGE <- hillock.pol_DGE %>%
  mutate(gene = toupper(gene))
# Extract unique ligands
ligand_genes <- unique(ncomms8866_rat$Ligand.ApprovedSymbol)
# Annotate ligands in marker list
hillock_annotated <- hillock.pol_DGE %>%
  mutate(Is_Ligand = gene %in% ligand_genes)
ligands_only <- setdiff(ligand_genes, matrix_genes)
hillock_ligands <- hillock.pol_DGE %>%
  filter(gene %in% ligands_only)

View(hillock_TFs)
View(hillock_matrix)
View(hillock_ligands)


