
# Close all, clear all
graphics.off()  
rm(list = ls())

# Load packages
library(Seurat)
library(ggplot2)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(patchwork)
library(stringr)

# load data
# Set wd to be true
setwd("~/Desktop/Datasets/Starting Populations")
# Load datasets needed
load("BAL.integrated_03-29-2025.Robj")
load("PD.integrated_03-29-2025.Robj")
load("Mixed_3D.start.combined.Robj")

# now global object (engineered)
setwd("~/Desktop/Datasets/Engineered Global Objects")
load("eng.subset.integrated_NodeAligned.Robj")

# check structures
str(Mixed_3D.start.combined@meta.data)
str(PD.integrated@meta.data)
str(BAL.integrated@meta.data)

# MEMO: Supplemental Analysis to test hypothesis in reference to reviewer 1s comment that epithelial trajectories and organoid morphology across conditions are primarily the result of initial epithelial cell number and physical crowding/fusion dynamics (space effects)
# What we are testing (working hypothesis): Differences in epithelial trajectories and organoid morphology, particularly in BAL_3D are driven by starting cellular composition
# Test 1: Epithelial Efficiency (Day 0 -> Day 10): If epithelial scarcity alone dictated outcmoes, BAL_3D should show lower epithelial fold expansion. But we predict it shows comparable or greater efficiency, suggesting a composition-driven effect
# Test 2: If space alone mattered, subtype distributions at D10 would scale up/down with epithelial number. Does BAL_3D truly enrich for specific epithelial fates?

# BAL doesn't just have fewer organoids, but rather more efficiently generates specific epithelial fates.
# Epithelial Fold Expansion (EFE), variable by condition (c): proportion of epithelial cells among all cells in condition (c), at timepoint t (day 0 or day 10)
# EFEc = (Pepi(c,D10)) / (Pepi(c, Day0)); quantifying how much the epithelial compartment expands relative to its starting fraction within each condition

## EFE (Epithelial Fold-Expansion)
# Uses Orig_ID as the replicate unit throughout

# Define epithelial labels
# Day 10 (eng.subset.integrated_HK$CellType): epithelial subtypes we want to tabulate
epi_day10 <- c(
  "ATI_Like","ATII_Like","Basal_Like","Secretory","Ciliated",
  "Hillock_Basal","Hillock_Luminal","RAS_Like",
  "Cycling_Epithelium","Stressed_Progenitor")
  # Adjust if you want stricter/looser definition

## Day 0 (start objects' $CellType):
epi_start_generic <- c(
  "ATI","ATII","Basal","Secretory","Ciliated","Tuft","Pan_Epithelium")

# This function calculates the epithelial fraction per replicate (Replicate.No.)
# and then averages across all replicates for that condition
p_epi_by_replicate <- function(meta, label_col = "CellType",
                               epi_levels, rep_col = "Orig_ID") {
  # make sure the metadata actually has the columns we need/want
  # now build a small data.frame with replicate IDs and whether each cell is epithelial
  stopifnot(all(c(label_col, rep_col) %in% names(meta)))
  d <- data.frame(
    Rep   = meta[[rep_col]],
    IsEpi = meta[[label_col]] %in% epi_levels)
  # compute the fraction of epithelial cells within each replicate
  per_rep <- aggregate(IsEpi ~ Rep, data = d, FUN = mean)
  # average the replicate fractions so each replicate counts equally
  mean(per_rep$IsEpi)
}

## Extract the metadata from the integrated Day 10 object
# D10 epithelial fraction by Condition (from eng.subset.integrated_HK - global obj)
meta10 <- eng.subset.integrated_HK@meta.data
# make sure the key columns are present
stopifnot(all(c("Condition","Orig_ID","CellType") %in% names(meta10)))

# split by Condition, compute per-Orig_ID means, then average replicates
cond_splits <- split(meta10[, c("Condition","Orig_ID","CellType")], meta10$Condition)
# For each condition:
# mark whether each cell is epithelial,
# compute per-Replicate.No. epithelial fraction,
# average across replicates
p10 <- sapply(cond_splits, function(d) {
  d$IsEpi <- d$CellType %in% epi_day10
  per_rep <- aggregate(IsEpi ~ Orig_ID, data = d, FUN = mean)
  mean(per_rep$IsEpi)
})
# we should get results where p10 is a named numeric vector with conditions as: "PD_3D","BAL_3D","Mixed_3D"

# D0 epithelial fraction by Condition (from separate start objects, b/c all of our starting objects are separate)
# empty vector to hold results
p0 <- c()

# PD start object -> map to PD_3D
# For those not familiar with this kind of conditional logic:
# We are saying, if PD.integrated exists,  pull out its metadata slot (@meta.data) and store it in a variable called m
# m is then a data frame containing all the cell-level metadata for the PD start condition
if (exists("PD.integrated")) {
  m <- PD.integrated@meta.data
  stopifnot(all(c("Orig_ID","CellType") %in% names(m)))
  p0["PD_3D"] <- p_epi_by_replicate(m, label_col = "CellType",
                                    epi_levels = epi_start_generic,
                                    rep_col = "Orig_ID")
}

# BAL start object -> map to BAL_3D
if (exists("BAL.integrated")) {
  m <- BAL.integrated@meta.data
  stopifnot(all(c("Orig_ID","CellType") %in% names(m)))
  p0["BAL_3D"] <- p_epi_by_replicate(m, label_col = "CellType",
                                     epi_levels = epi_start_generic,
                                     rep_col = "Orig_ID")
}

# Mixed start object -> map to Mixed_3D
if (exists("Mixed_3D.start.combined")) {
  m <- Mixed_3D.start.combined@meta.data
  stopifnot(all(c("Orig_ID","CellType") %in% names(m)))
  ## If you prefer well-level replicates instead of Orig_ID, switch to rep_col = "FinalReplicate"
  p0["Mixed_3D"] <- p_epi_by_replicate(m, label_col = "CellType",
                                       epi_levels = epi_start_generic,
                                       rep_col = "Orig_ID")
}

# Combine D0 and D10 and compute EFE
# only keep conditions present in both D0 and D10
conds <- intersect(names(p10), names(p0))
# build a tidy results data.frame
# sort conditions alphabetically and print with 3 significant digits
out <- data.frame(
  Condition   = conds,
  p_epi_Day0  = as.numeric(p0[conds]),
  p_epi_Day10 = as.numeric(p10[conds]),
  EFE         = as.numeric(p10[conds]) / as.numeric(p0[conds]),
  row.names   = NULL)

## Sort nicely and print
out <- out[order(out$Condition), ]
print(out, digits = 3)

# Checking without replicate averaging (raw)
# This shows the Day 10 epithelial fractions by Condition *without replicate averaging*
# i.e., treating every cell equally instead of every Orig_ID equally
prop.table(table(meta10$Condition, meta10$CellType %in% epi_day10), 1)

# Now if we want to get the 95% CIs
## Only way I got to work ##
## Make per-replicate data-frames
# Day 10 per-replicate epithelial fraction
meta10 <- eng.subset.integrated_HK@meta.data
meta10$IsEpi10 <- meta10$CellType %in% epi_day10
perrep_day10 <- aggregate(IsEpi10 ~ Condition + Orig_ID, data=meta10, mean)

## Day 0 per-replicate epithelial fraction
mk_day0_perrep <- function(seu_meta, cond_label, rep_col){
  stopifnot(all(c("CellType", rep_col) %in% names(seu_meta)))
  d <- data.frame(
    Condition = cond_label,
    Rep       = seu_meta[[rep_col]],
    IsEpi0    = seu_meta$CellType %in% epi_start_generic)
  aggregate(IsEpi0 ~ Condition + Rep, data=d, mean)
}

perrep_day0_list <- list()
if (exists("PD.integrated"))   perrep_day0_list[["PD_3D"]]    <- mk_day0_perrep(PD.integrated@meta.data, "PD_3D", "Orig_ID")
if (exists("BAL.integrated"))  perrep_day0_list[["BAL_3D"]]   <- mk_day0_perrep(BAL.integrated@meta.data, "BAL_3D", "Orig_ID")
if (exists("Mixed_3D.start.combined")) perrep_day0_list[["Mixed_3D"]] <- mk_day0_perrep(Mixed_3D.start.combined@meta.data, "Mixed_3D", "FinalReplicate")
perrep_day0 <- do.call(rbind, perrep_day0_list)

## Boostrapping helper function
boot_ci <- function(x, R=2000){
n <- length(x)
if (n < 2) return(c(mean=mean(x), lwr=NA, upr=NA))
b <- replicate(R, mean(sample(x, n, replace=TRUE)))
c(mean=mean(b), lwr=quantile(b,0.025), upr=quantile(b,0.975))
}

## Now run for each condition
conds <- intersect(unique(perrep_day0$Condition),
                   unique(perrep_day10$Condition))

set.seed(123)
efe_ci <- do.call(rbind, lapply(conds, function(cc){
  x0 <- perrep_day0$IsEpi0[perrep_day0$Condition==cc]
  x10 <- perrep_day10$IsEpi10[perrep_day10$Condition==cc]
  
  ci0  <- boot_ci(x0)
  ci10 <- boot_ci(x10)
  
  # bootstrap ratios by pairing means
  br <- replicate(2000, {
    mean(sample(x10, length(x10), replace=TRUE)) /
      mean(sample(x0, length(x0), replace=TRUE))
  })
  
  data.frame(
    Condition = cc,
    n_reps_day0 = length(x0),
    n_reps_day10 = length(x10),
    p0_mean=ci0["mean"], p0_lwr=ci0["lwr"], p0_upr=ci0["upr"],
    p10_mean=ci10["mean"], p10_lwr=ci10["lwr"], p10_upr=ci10["upr"],
    EFE_mean=mean(br), EFE_lwr=quantile(br,0.025), EFE_upr=quantile(br,0.975)
  )
}))
print(efe_ci, digits=3)

# Clean to remove NA columns
# We are getting NAs b/c within each condition we only have 3 replicates so the bootstrap helper function we create is faulting
# BECAUSE we set parameter such that return CIs when n < 2. But we still get the EFE info we want
efe_clean <- efe_ci %>%
  dplyr::select(Condition, n_reps_day0, n_reps_day10,
                p0_mean, p10_mean,
                EFE_mean, EFE_lwr, EFE_upr)
print(efe_clean, digits = 3)

### Making a better formatted table
# if efe_clean doesn't exist, derive it from efe_ci
if (!exists("efe_clean") && exists("efe_ci")) {
  efe_clean <- dplyr::select(
    efe_ci,
    Condition, n_reps_day0, n_reps_day10,
    p0_mean, p10_mean,
    EFE_mean, EFE_lwr, EFE_upr)
}

# small format helpers
fmt_pct <- function(x, acc = 0.1) scales::percent(as.numeric(x), accuracy = acc)
fmt_num <- function(x, k = 2) formatC(as.numeric(x), format = "f", digits = k)

# build the pretty table (ordered by effect size)
efe_table <- efe_clean |>
  dplyr::arrange(dplyr::desc(EFE_mean)) |>
  dplyr::mutate(
    `Replicates (Day0/Day10)` = sprintf("%d / %d", n_reps_day0, n_reps_day10),
    `Day 0 epithelial`        = fmt_pct(p0_mean, 0.1),
    `Day 10 epithelial`       = fmt_pct(p10_mean, 0.1),
    `EFE (fold)`              = sprintf(
      "%s (%s–%s)",
      fmt_num(EFE_mean, 2),
      fmt_num(EFE_lwr, 2),
      fmt_num(EFE_upr, 2))
  ) |>
  dplyr::select(
    Condition,
    `Replicates (Day0/Day10)`,
    `Day 0 epithelial`,
    `Day 10 epithelial`,
    `EFE (fold)`)
# print nicely in console
efe_table

### Making forest plot of EFE with 95% CIs
#  pick clean log10 breaks spanning the CIs
nice_breaks <- c(1, 2, 3, 5, 10, 20, 50, 100, 200, 500)
rng <- range(c(plot_df$EFE_lwr, plot_df$EFE_upr), na.rm = TRUE)
brks <- nice_breaks[nice_breaks >= max(1, floor(rng[1])) & nice_breaks <= ceiling(rng[2])]
if (length(brks) < 3) brks <- nice_breaks  # fallback for when range is tight

# final plot; still log-scaled, mention in caption
gg_efe_overlay_col <- ggplot() +
  # bootstrap cloud (semi-transparent), jittered horizontally
  geom_point(
    data = transform(boot_df, Condition = factor(Condition, levels = levels(plot_df$Condition))),
    aes(x = Condition, y = EFE_boot, color = Condition),
    alpha = 0.22, size = 1.2,
    position = position_jitter(width = 0.08, height = 0)) +
  # 95% CI whiskers
  geom_errorbar(data = plot_df,
    aes(x = Condition, ymin = EFE_lwr, ymax = EFE_upr, color = Condition),
    width = 0.25, linewidth = 0.9, lineend = "butt") +
  # mean point
  geom_point(data = plot_df,
    aes(x = Condition, y = EFE_mean, fill = Condition),
    shape = 21, size = 3.8, stroke = 0.9, color = "black") +
  # reference line at no change
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.4, color = "grey40") +
  # IMPORTANT: keep log scale but show plain numeric ticks
  scale_y_log10(breaks = brks, labels = brks, expand = expansion(mult = c(0.02, 0.15))) +
  scale_color_manual(values = cond.cols) +
  scale_fill_manual(values = cond.cols) +
  coord_flip() +
  labs(x = "Condition",
    # label has no "log scale" mention, put in caption
    y = "Epithelial Fold-Expansion (EFE)",
    caption = "Cloud: bootstrap draws (R=4000); point: bootstrap mean; whiskers: 95% CI; dashed line: 1× (no expansion)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
    plot.caption = element_text(hjust = 0, size = 10),
    axis.title.x = element_text(margin = margin(t = 6)),
    panel.grid.major.y = element_blank())
# print
gg_efe_overlay_col
# ggsave("FigSx_EFE_forest_bootstrap_colored.pdf", gg_efe_overlay_col, width=6.6, height=3.4, units="in")
setwd("~/Desktop/iScience Transfer Submission/Additional Analyses/Epithelial Expansion")
ggsave("Supp_EFE_forest_bootstrap_colored.png", gg_efe_overlay_col, width=6.0, height=4.4, units="in", dpi=600)
