


# close all, clear all
graphics.off()  
rm(list = ls())

# Load the required libraries
library(CellChat)
library(Seurat)
library(dplyr)
library(readxl)
library(patchwork)
library(openxlsx)
options(stringsAsFactors = FALSE)
# reticulate::use_python("/Users/suoqinjin/anaconda3/bin/python", required=T) 

# Load data
setwd("~/Desktop/Datasets/Engineered Global Objects")
load("eng.subset.integrated_NodeAligned.Robj")

# Assign standard name to object 
seurat_object <- eng.subset.integrated_HK

# Extract necessary information from seurat object
data.input <- seurat_object[["RNA"]]$data # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
labels <- Idents(seurat_object) # idents being the CellType.NodeAligned
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels

# create the cell chat object
cellchat <- createCellChat(object = seurat_object, group.by = "ident", assay = "RNA")

# using mouse instead of rat because more convenient - rat upload causing too much trouble
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show structure
dplyr::glimpse(CellChatDB$interaction)
cellchat@DB <- CellChatDB

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
# Issue identified!! Please check the official Gene Symbol of the following genes:  H2-BI H2-Ea-ps 
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is ___

ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, nboot = 20) # will speed things up
# > triMean is used for calculating the average gene expression per cell group. 
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# All possible pathways
all_pathways_db <- sort(unique(cellchat@DB$interaction$pathway_name))
length(all_pathways_db); head(all_pathways_db, 20)
# By category (useful for scanning)
split_by_cat <- split(cellchat@DB$interaction$pathway_name,
                      cellchat@DB$interaction$annotation)  # e.g. Secreted/ECM/Cell-Cell
lapply(split_by_cat, function(v) sort(unique(v)))

# In our object
available_pathways <- sort(unique(cellchat@netP$pathways))
length(available_pathways); head(available_pathways, 20)
available_pathways

# Looking at just Cxcl pathway
pathways.show <- c("EPHB") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

### Attempt at making cord plots

## Adding color palette
# your palette
CellType.cols <- c(
  "ATI_Like" = "#b22222","ATI" = "#b22222",
  "ATII" = "#4682B4","ATII_Like" = "#4682B4",
  "Secretory_Like" = "#2E8B57", "Secretory" = "#2E8B57",
  "Ciliated" = "#9B489B","Ciliated_Like" = "#9B489B",
  "Stressed_Progenitor" = "#c71585",
  "Hillock_Like" = "#ff0000",
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

# align to our object’s identities
lvls <- levels(factor(cellchat@idents))
fallback <- "#B0B0B0"  # neutral gray for any label not in your palette
color.use <- setNames(rep(fallback, length(lvls)), lvls)
matchable <- intersect(names(CellType.cols), lvls)
color.use[matchable] <- CellType.cols[matchable]
# (optional) quick audit: who’s using fallback?
setdiff(lvls, names(CellType.cols))

## Plot for Hillock -> Immune
sources.use <- "Hillock_Like"
targets.use <- c("Cycling_Immune","Polarized_Mac")
sel <- subset(edges_pair, source %in% sources.use & target %in% targets.use)
# look at sel to pick
View(sel)
# genes we want to show - specific for Hillock -> Distal_Epi
want_ligs <- c("Cd44","C3","Cxcl3","Ccl6","Il33","Adm","Cxcl1")
want_recs <- c("Lamc1","ITGAM_ITGB2","Cxcr1","Ccr1","IL1RL1_IL1RAP","Calcrl")
pairLR_df <- unique(subset(sel, ligand %in% want_ligs & receptor %in% want_recs)[, "interaction_name", drop=FALSE])
stopifnot(nrow(pairLR_df) > 0)
library(circlize); circos.clear()
hill_imm = netVisual_chord_gene(
  object       = cellchat,
  sources.use  = sources.use,
  targets.use  = targets.use,
  pairLR.use   = pairLR_df,      # <- from A/B/C above
  color.use    = color.use,      # your CellType.cols-aligned vector
  small.gap    = 0.15,
  big.gap      = 0.4,
  lab.cex      = 0.45,
  transparency = 0.4,
  scale        = TRUE,
  reduce       = 0.001,
  title.name   = sprintf("Hillock → Immune", nrow(pairLR_df)),
  thresh       = 0.05) # only significant
print(hill_imm)
setwd("~/Desktop/iScience Transfer Submission/Additional Analyses/CellChat")
tiff("hill_imm_CellChat_chord.tiff", units = "in", width = 10, height = 10, res = 600, compression = "lzw")
replayPlot(hill_imm)
dev.off()

## Plot for Hillock -> Mes
sources.use <- "Hillock_Like"
targets.use <- c("Rspo3+_Mes")
sel <- subset(edges_pair, source %in% sources.use & target %in% targets.use)
# look at sel to pick
View(sel)
# genes we want to show - specific for Hillock -> Distal_Epi
want_ligs <- c("Jag1","Nampt","Efnb2","Fgf9","Wnt4","Cdh1")
want_recs <- c("Notch2","ITGA5_ITGB1","Ephb4","Ephb1","Fgfr1","Fgfr2","FZD1_LRP6","ITGA2_ITGB1")
pairLR_df <- unique(subset(sel, ligand %in% want_ligs & receptor %in% want_recs)[, "interaction_name", drop=FALSE])
stopifnot(nrow(pairLR_df) > 0)
library(circlize); circos.clear()
hill_mes = netVisual_chord_gene(
  object       = cellchat,
  sources.use  = sources.use,
  targets.use  = targets.use,
  pairLR.use   = pairLR_df,      # <- from A/B/C above
  color.use    = color.use,      # your CellType.cols-aligned vector
  small.gap    = 1.0,
  big.gap      = 1.2,
  lab.cex      = 0.45,
  transparency = 0.4,
  scale        = TRUE,
  reduce       = 0.01,
  title.name   = sprintf("Hillock → Mesenchyme + Immune", nrow(pairLR_df)),
  thresh       = 0.05) # only significant
print(hill_mes)
setwd("~/Desktop/iScience Transfer Submission/Additional Analyses/CellChat")
tiff("hill_mes_CellChat_chord.tiff", units = "in", width = 10, height = 10, res = 600, compression = "lzw")
replayPlot(hill_mes)
dev.off()

## Plot for Hillock -> Mes + Immune
sources.use <- "Hillock_Like"
targets.use <- c("Rspo3+_Mes","Polarized_Mac","Cycling_Immune")
sel <- subset(edges_pair, source %in% sources.use & target %in% targets.use)
# look at sel to pick
View(sel)
# genes we want to show - specific for Hillock ->  Mes + Immune
want_ligs <- c("Jag1","Nampt","Efnb2","Fgf9","Cd44","C3","Cxcl3","Ccl6","Il33","Adm","Wnt11","Wnt4","Cdh1")
want_recs <- c("Notch2","ITGA5_ITGB1","Ephb4","Ephb1","Fgfr1","Fgfr2","Lamc1","ITGAM_ITGB2","Cxcr1","Ccr1","IL1RL1_IL1RAP","Calcrl","Fzd1","ITGA2_ITGB1","FZD1_LRP6")
pairLR_df <- unique(subset(sel, ligand %in% want_ligs & receptor %in% want_recs)[, "interaction_name", drop=FALSE])
stopifnot(nrow(pairLR_df) > 0)
library(circlize); circos.clear()
hill_imm.mes = netVisual_chord_gene(
  object       = cellchat,
  sources.use  = sources.use,
  targets.use  = targets.use,
  pairLR.use   = pairLR_df,      # <- from A/B/C above
  color.use    = color.use,      # your CellType.cols-aligned vector
  small.gap    = 1.0,
  big.gap      = 1.2,
  lab.cex      = 0.45,
  transparency = 0.4,
  scale        = TRUE,
  reduce       = 0.01,
  title.name   = sprintf("Hillock → Mesenchyme + Immune", nrow(pairLR_df)),
  thresh       = 0.05) # only significant
print(hill_imm.mes)
setwd("~/Desktop/iScience Transfer Submission/Additional Analyses/CellChat")
tiff("hill_imm.mes_CellChat_chord.tiff", units = "in", width = 10, height = 10, res = 600, compression = "lzw")
replayPlot(hill_imm.mes)
dev.off()

# Plot for Hillock -> Distal Epi
sources.use <- "Hillock_Like"
targets.use <- c("ATII_Like","RAS_Like","Stressed_Progenitor","ATI_Like","Cycling_Distal_Epi")
sel <- subset(edges_pair, source %in% sources.use & target %in% targets.use)
View(sel)
# genes we want to show - specific for Hillock -> Distal_Epi
want_ligs <- c("Wnt4","Slurp1","Il1a","Tgfa","Hbegf","Wnt7b","Wnt5a")
want_recs <- c("FZD6_LRP6","Chrnb1","IL1R1_IL1RAP","EGFR_ERBB2","Fzd6")
pairLR_df <- unique(subset(sel, ligand %in% want_ligs & receptor %in% want_recs)[, "interaction_name", drop=FALSE])
stopifnot(nrow(pairLR_df) > 0)
library(circlize); circos.clear()
hill_dist.epi = netVisual_chord_gene(
  object       = cellchat,
  sources.use  = sources.use,
  targets.use  = targets.use,
  pairLR.use   = pairLR_df,      # <- from A/B/C above
  color.use    = color.use,      # your CellType.cols-aligned vector
  small.gap    = 1.0,
  big.gap      = 1.2,
  lab.cex      = 0.45,
  transparency = 0.4,
  scale        = TRUE,
  reduce       = 0.001,
  title.name   = sprintf("Hillock → Distal Epithelium", nrow(pairLR_df)),
  thresh       = 0.05) # only significant
print(hill_dist.epi)
setwd("~/Desktop/iScience Transfer Submission/Additional Analyses/CellChat")
tiff("hill_dist.epi_CellChat_chord.tiff", units = "in", width = 10, height = 10, res = 600, compression = "lzw")
replayPlot(hill_dist.epi)
dev.off()

# Plot for Hillock -> Proximal_Epi
sources.use <- "Hillock_Like"
targets.use <- c("Basal_Like","Secretory","Cycling_Proximal_Epi")
sel <- subset(edges_pair, source %in% sources.use & target %in% targets.use)
View(sel)
# genes we want to show - specific for Hillock -> Distal_Epi
want_ligs <- c("Bmp6","Slurp1","Cdh1","Gas6","Lamc2","Jag1")
want_recs <- c("BMPR1B_ACVR2A","Chrnb1","ITGA6_ITGB1","Axl","Cd44","Notch1")
pairLR_df <- unique(subset(sel, ligand %in% want_ligs & receptor %in% want_recs)[, "interaction_name", drop=FALSE])
stopifnot(nrow(pairLR_df) > 0)
library(circlize); circos.clear()
hill_prox.epi = netVisual_chord_gene(
  object       = cellchat,
  sources.use  = sources.use,
  targets.use  = targets.use,
  pairLR.use   = pairLR_df,      # <- from A/B/C above
  color.use    = color.use,      # your CellType.cols-aligned vector
  small.gap    = 1.0,
  big.gap      = 1.2,
  lab.cex      = 0.45,
  transparency = 0.4,
  scale        = TRUE,
  reduce       = 0.001,
  title.name   = sprintf("Hillock → Proximal Epithelium", nrow(pairLR_df)),
  thresh       = 0.05) # only significant
print(hill_prox.epi)
setwd("~/Desktop/iScience Transfer Submission/Additional Analyses/CellChat")
tiff("hill_prox.epi_CellChat_chord.tiff", units = "in", width = 10, height = 10, res = 600, compression = "lzw")
replayPlot(hill_prox.epi)
dev.off()

## Making a singular sending plot and then receiving for figure 3

# First sending
# Plot for Hillock -> Others
sources.use <- "Hillock_Like"
targets.use <- c("Basal_Like","Secretory","Stressed_Progenitor","RAS_Like","ATI_Like","ATII_Like",
                 "Rspo3+_Mes","Polarized_Mac")
sel <- subset(edges_pair, source %in% sources.use & target %in% targets.use)
View(sel)
# genes we want to show - specific for Hillock -> Distal_Epi
want_ligs <- c("Tgfa","Efna1","Il1a","Wnt11","Wnt4","Jag1","Lamc2","Wnt5a","Wnt7b","Cxcl3","Cxcl1","Ccl6","Fgf9")
want_recs <- c("Egfr","Epha2","IL1R1_IL1RAP","Fzd6","FZD1_LRP6","Cd44","Fzd1","Cxcr1","Ccr1","Fgfr2","Fgfr1")
pairLR_df <- unique(subset(sel, ligand %in% want_ligs & receptor %in% want_recs)[, "interaction_name", drop=FALSE])
stopifnot(nrow(pairLR_df) > 0)
library(circlize); circos.clear()
hill_send = netVisual_chord_gene(
  object       = cellchat,
  sources.use  = sources.use,
  targets.use  = targets.use,
  pairLR.use   = pairLR_df,      # <- from A/B/C above
  color.use    = color.use,      # your CellType.cols-aligned vector
  small.gap    = 1.0,
  big.gap      = 1.2,
  lab.cex      = 0.45,
  transparency = 0.4,
  scale        = TRUE,
  reduce       = 0.001,
  title.name   = sprintf("Hillock → Others", nrow(pairLR_df)),
  thresh       = 0.05) # only significant
print(hill_send)
setwd("~/Desktop/iScience Transfer Submission/Additional Analyses/CellChat")
tiff("hill_send_CellChat_chord.tiff", units = "in", width = 10, height = 10, res = 600, compression = "lzw")
replayPlot(hill_send)
dev.off()

# Now for receiving
# Plot for Others -> Hillock
targets.use <- "Hillock_Like"
sources.use <- c("Basal_Like","Secretory","Stressed_Progenitor","RAS_Like","ATI_Like","ATII_Like",
                 "Rspo3+_Mes","Polarized_Mac")
sel <- subset(edges_pair, source %in% sources.use & target %in% targets.use)
View(sel)
# genes we want to show - specific for Hillock -> Distal_Epi
want_ligs <- c("Fn1","Angptl4","Spp1","Areg","Sema4a","Dsg2","Nectin1")
want_recs <- c("Sdc1","Cd44","Egfr","Cd47","Plxnb2","Dsc3","Nectin4")
pairLR_df <- unique(subset(sel, ligand %in% want_ligs & receptor %in% want_recs)[, "interaction_name", drop=FALSE])
stopifnot(nrow(pairLR_df) > 0)
library(circlize); circos.clear()
hill_rec = netVisual_chord_gene(
  object       = cellchat,
  sources.use  = sources.use,
  targets.use  = targets.use,
  pairLR.use   = pairLR_df,      # <- from A/B/C above
  color.use    = color.use,      # your CellType.cols-aligned vector
  small.gap    = 1.0,
  big.gap      = 2.2,
  lab.cex      = 1.2, # increases font size
  transparency = 0.4,
  scale        = TRUE,
  reduce       = 0.001,
  title.name   = sprintf("Hillock → Others", nrow(pairLR_df)),
  thresh       = 0.05) # only significant
print(hill_rec)
setwd("~/Desktop/iScience Transfer Submission/Additional Analyses/CellChat")
tiff("hill_rec_CellChat_chord.tiff", units = "in", width = 10, height = 10, res = 600, compression = "lzw")
replayPlot(hill_rec)
dev.off()


########################## SAVING #####################################

# Make copy to give unique name
cellchat_org.m = cellchat 
# Before moving any further, let's make sure we are saving the object
# Save
saveRDS(cellchat_org.m, file = sprintf("cellchat_%s.rds", Sys.Date()), compress = "xz")
save(cellchat_org.m, file = "cellchat_org.m.Robj")
# Load (in a new session)
# cellchat <- readRDS(sprintf("cellchat_%s.rds", Sys.Date()))

# 1) Save key tables so you can replot without the heavy object
edges_pairs   <- subsetCommunication(cellchat_org.m, slot.name = "net")   # LR edges (all pairs)
edges_pathway <- subsetCommunication(cellchat_org.m, slot.name = "netP")  # pathway-level edges
write.csv(edges_pairs,   "cellchat_edges_pairs.csv",   row.names = FALSE)
write.csv(edges_pathway, "cellchat_edges_pathway.csv", row.names = FALSE)

# 2) Record parameters & versions for reproducibility
params <- list(
  compute_type     = "truncatedMean",
  trim             = 0.1,
  nboot            = 20,
  db_species       = "mouse",
  identities       = levels(factor(cellchat@idents)),
  date             = as.character(Sys.time()),
  CellChat_version = as.character(packageVersion("CellChat")))
saveRDS(params, "cellchat_run_params.rds")
utils::capture.output(sessionInfo(), file = "sessionInfo.txt")


