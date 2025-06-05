
## MEMO: Exploring biological process enrichment with MSigDB in the immune (macrophage) subset

# close all, clear all
graphics.off()  
rm(list = ls())

# Load packages
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

# Load all MSigDB gene sets for rat
msigdb_rat <- msigdbr(species = "Rattus norvegicus")
# Filter for desired categories
msig_h <- msigdb_rat %>% filter(gs_collection == "H") # Hallmark paths
msig_c2 <- msigdb_rat %>% filter(gs_collection == "C2") # Curated pathways
msig_reactome <- msigdb_rat %>%
  filter(gs_collection == "C2", gs_subcollection == "CP:REACTOME") # reactome
msig_kegg <- msigdb_rat %>%
  filter(gs_collection == "C2", gs_subcollection == "CP:KEGG") # KEGG pathways

# Run marker list of pro vs. anti-inflammatory
Idents(immune_subset.clean) <- immune_subset.clean$seurat_clusters
mac.pol_DGE <- FindMarkers(immune_subset.clean, ident.1 = "0", ident.2 = c("1","2"),
                           min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE, test.use = "wilcox")
mac.pol_DGE$ratio <- mac.pol_DGE$pct.1 / mac.pol_DGE$pct.2
mac.pol_DGE$power <- mac.pol_DGE$ratio * mac.pol_DGE$avg_log2FC
View(mac.pol_DGE)

# Use log2FC as ranking metric
gene_list <- mac.pol_DGE$avg_log2FC
# Set names to gene symbols
names(gene_list) <- rownames(mac.pol_DGE)
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
gsea_msig_mac <- GSEA(
  geneList = gene_list,
  TERM2GENE = msigdb_df,
  pvalueCutoff = 0.5,
  verbose = FALSE)
head(gsea_msig_mac@result, 10)
View(gsea_msig_mac@result)

# Filter data
gsea_msig_mac_filtered <- gsea_msig_mac@result %>% # pull the GSEA results table (@result) from gsea_msig object and pipe into dplyr for filtering
  mutate(log10p = -log10(p.adjust) # Adds a new column called EnrichedIn based on the Normalized Enrichment Score (NES)
  ) %>% # NES > 0, the pathway is enriched in Hillock Basal, NES < 0, the pathway is enriched in Hillock Luminal
  arrange(desc(log10p)) %>% # Add a column log10p, which is the -log10 of the adjusted p-value
  filter(p.adjust < 0.5)  # sort so most significant (lowest p-adjust) are first
View(gsea_msig_mac_filtered)

# Load Ground Truth Gene Sets
## Ligands from NICHES
ligand_genes <- unique(ncomms8866_rat$Ligand.ApprovedSymbol)

##  Matrix proteins from MatrisomeDB (mouse, used for rat)
# Replace with your actual path to the Excel file
setwd("~/Desktop/Single Cell/Hillock_Pathway_Analysis/Ground Truths")
matrisome_df <- writexl::read_excel("Mm_Matrisome_Masterlist_Naba et al_2012.xlsx")
# Just pull gene names
# Matrisome includes more than just matrix stuff - includes ECM regulators, secreted factors, etc. but we don't want those right now
ecm_categories <- c("ECM Glycoproteins", "Collagens", "Proteoglycans", "ECM Regulators")
# Filter matrisome_df
filtered_matrisome_df <- matrisome_df %>%
  filter(`Matrisome Category` %in% ecm_categories)
# Extract gene symbols
matrix_genes <- unique(filtered_matrisome_df$`Gene Symbol`)

## Transcription Factors from AnimalTFDB (rat-specific)
# Replace with your actual path
tf_df <- read.delim("Rattus_norvegicus_TF.txt")  # from AnimalTFDB
# Pull genes
tf_genes <- unique(tf_df$Symbol)

## Annotate GSEA Results
gsea_df_annotated <- gsea_msig_mac_filtered
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
  ligands_present <- intersect(core_genes, ligand_genes)
  gsea_df_annotated$Ligands_NICHES[i] <- paste(ligands_present, collapse = ", ")
  gsea_df_annotated$Ligand_NICHES_Count[i] <- length(ligands_present)
  
  # Transcription factors
  tf_present <- intersect(core_genes, tf_genes)
  gsea_df_annotated$TFs_present[i] <- paste(tf_present, collapse = ", ")
  gsea_df_annotated$TF_Count[i] <- length(tf_present)
}
View(gsea_df_annotated)

# Filter for separate plots
top_n <- 40 # top 40
gsea_mac <- gsea_df_annotated %>% 
  slice_max(log10p, n = top_n)

# Plotting the results
ggplot(gsea_mac, aes(x = log10p, y = reorder(Description, log10p))) +
  geom_col(fill = "#9497fd") +
  geom_text(aes(x = 0.05, label = Description), 
            hjust = 0, size = 3.5, color = "black", fontface = "bold") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    title = "Pathways Upregulated in Pro_Inflamm_Mac State",
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

# Cherry-picking top (significant) pathways
selected_pathways <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_HYPOXIA",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE")

# Filtering
gsea.mac.df <- gsea_df_annotated %>%
  filter(Description %in% selected_pathways)
  
# Plotting again with filtration
ggplot(gsea.mac.df, aes(x = log10p, y = reorder(Description, log10p))) +
  geom_col(fill = "#9497fd") +
  geom_text(aes(x = 0.05, label = Description), 
            hjust = 0, size = 3.5, color = "black", fontface = "bold") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    title = "Pathways Upregulated in Pro_Inflamm_Mac State",
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

# Fix the EnrichedIn labels if needed:
gsea.mac.df$EnrichedIn <- case_when(
  gsea.mac.df$NES < 0 ~ "Anti_Inflamm_Mac",
  gsea.mac.df$NES > 0 ~ "Pro_Inflamm_Mac",
  TRUE ~ "NA")

# Make a df for a mirrored plot
gsea_mirror_df <- gsea_mac %>%
  filter(Description %in% selected_pathways) %>%
  mutate(
    log10p_signed = ifelse(EnrichedIn == "Pro_Inflamm_Mac", -log10p, log10p),
    Group = EnrichedIn,
    Label_Pos = ifelse(log10p_signed > 0, 0.05, -0.05),
    Description = fct_reorder(Description, log10p_signed))

# Plot: Pathways on Y-axis, enrichment on X
hallmark_hillock.mirror = ggplot(gsea_mirror_df, aes(x = log10p_signed, y = Description, fill = EnrichedIn)) +
  geom_col() +
  geom_text(aes(x = Label_Pos, label = Description),
            hjust = ifelse(gsea_mirror_df$log10p_signed > 0, 0, 1),
            size = 2.5, fontface = "bold") +
  scale_fill_manual(values = c("Pro_Inflamm_Mac" = "#9497fd", "Anti_Inflamm_Mac" = "#ff9e80")) +
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
