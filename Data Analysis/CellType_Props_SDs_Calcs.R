

# Packages
library(grid)  # For unit()
library(ggplot2)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(ggalluvial)
library(tidyr)
library(tibble)

# Set wd and load data
setwd("~/Desktop/Datasets/Starting Populations")
load("BAL.integrated_03-29-2025.Robj")
load("PD.integrated_03-29-2025.Robj")

# Check datasets (how many cells, etc)
BAL.integrated
PD.integrated

# Structure of PD_integrated
str(PD.integrated@meta.data)
table(PD.integrated$CellType)
table(PD.integrated$Orig_ID)
# We want to calculate the percentages by cell type and replicate to report in the text, as we do for BAL

# Define object and replicate variable
obj <- PD.integrated
rep_var <- if ("Orig_ID" %in% colnames(obj@meta.data)) "Orig_ID" else "Replicate.No."

# Clean tibble
md <- tibble(
  replicate = as.character(obj@meta.data[[rep_var]]),
  CellType  = as.character(obj@meta.data$CellType))

# Per-replicate totals
totals_by_rep <- md %>%
  dplyr::group_by(replicate) %>%
  dplyr::summarise(total_cells = dplyr::n(), .groups = "drop")

# per-replicate counts by CellType
byrep_counts <- md %>%
  dplyr::group_by(replicate, CellType) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::left_join(totals_by_rep, by = "replicate") %>%
  dplyr::mutate(pct = 100 * n / total_cells)

# zero-fill missing cell types per replicate
all_types <- sort(unique(byrep_counts$CellType))
all_reps  <- sort(unique(byrep_counts$replicate))

byrep <- tidyr::complete(
  byrep_counts,
  replicate = all_reps,
  CellType  = all_types,
  fill      = list(n = 0, pct = 0)) %>%
  dplyr::group_by(replicate) %>%
  tidyr::fill(total_cells, .direction = "downup") %>%
  dplyr::ungroup()

# Summarize across replicates (mean ± SD % of total cells)
summary_all <- byrep %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(
    n_reps   = dplyr::n_distinct(replicate),
    mean_n   = mean(n,   na.rm = TRUE),
    sd_n     = sd(n,     na.rm = TRUE),
    mean_pct = mean(pct, na.rm = TRUE),
    sd_pct   = sd(pct,   na.rm = TRUE),
    .groups  = "drop") %>%
  dplyr::arrange(dplyr::desc(mean_pct)) %>%
  dplyr::mutate(
    `%(mean)` = round(mean_pct, 1),
    `%(SD)`   = round(sd_pct, 1),
    report    = sprintf("%.1f ± %.1f", `%(mean)`, `%(SD)`)) %>%
  dplyr::select(CellType, n_reps, `%(mean)`, `%(SD)`, report, `n(mean)` = mean_n, `n(SD)` = sd_n)
# Print
summary_all

# Also want to check numbers on engineered data-sets
setwd("~/Desktop/Datasets/Condition Level Objects")
load("Mixed_3D.integrated_HK.NodeAligned.Robj")
load("BAL_3D.integrated_HK.NodeAligned.Robj")
load("PD_3D.integrated_HK.NodeAligned.Robj")

# numbers in each
BAL_3D.integrated_HK # 9676
PD_3D.integrated_HK # 8724
Mixed_3D.integrated_HK # 6344

str(eng.subset.integrated_HK@meta.data)
# For global object
library(dplyr)
library(tidyr)

# Define object and replicate/condition variables
obj <- eng.subset.integrated_HK
rep_var <- if ("Orig_ID" %in% colnames(obj@meta.data)) "Orig_ID" else "Replicate.No."
cond_var <- "Condition"

# Choose which cell-type field you want
ctype_col <- "CellType.NodeAligned"  # or "CellType.NodeAligned"

# Clean tibble
md <- tibble(
  condition = as.character(obj@meta.data[[cond_var]]),
  replicate = as.character(obj@meta.data[[rep_var]]),
  CellType.NodeAligned  = as.character(obj@meta.data[[ctype_col]]))

# Per-replicate totals (within condition)
totals_by_rep <- md %>%
  dplyr::group_by(condition, replicate) %>%
  dplyr::summarise(total_cells = dplyr::n(), .groups = "drop")

# Per-replicate counts by CellType (within condition)
byrep_counts <- md %>%
  dplyr::group_by(condition, replicate, CellType.NodeAligned) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::left_join(totals_by_rep, by = c("condition","replicate")) %>%
  dplyr::mutate(pct = 100 * n / total_cells)

# Zero-fill missing cell types per replicate (within each condition × replicate)
all_types <- sort(unique(byrep_counts$CellType))

byrep <- byrep_counts %>%
  dplyr::group_by(condition, replicate) %>%
  tidyr::complete(CellType.NodeAligned = all_types, fill = list(n = 0, pct = 0)) %>%
  tidyr::fill(total_cells, .direction = "downup") %>%
  dplyr::ungroup()

# Summarize across replicates: mean ± SD (% and counts) for each Condition × CellType
summary_all <- byrep %>%
  dplyr::group_by(condition, CellType.NodeAligned) %>%
  dplyr::summarise(
    n_reps   = dplyr::n_distinct(replicate),
    mean_n   = mean(n,   na.rm = TRUE),
    sd_n     = sd(n,     na.rm = TRUE),
    mean_pct = mean(pct, na.rm = TRUE),
    sd_pct   = sd(pct,   na.rm = TRUE),
    .groups  = "drop") %>%
  dplyr::arrange(condition, dplyr::desc(mean_pct)) %>%
  dplyr::mutate(
    `%(mean)` = round(mean_pct, 1),
    `%(SD)`   = round(sd_pct, 1),
    report    = sprintf("%.1f ± %.1f", `%(mean)`, `%(SD)`)) %>%
  dplyr::select(
    condition, CellType.NodeAligned, n_reps,
    `%(mean)`, `%(SD)`, report,
    `n(mean)` = mean_n, `n(SD)` = sd_n)

# Print
summary_all
View(summary_all)
