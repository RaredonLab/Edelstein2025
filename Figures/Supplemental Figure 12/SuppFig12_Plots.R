
## MEMO: Plotting the organoid to macrophage distance data (Supplemental Fig 12)

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
library(readxl) # for loading in distance data

# Set wd
setwd("~/Library/CloudStorage/OneDrive-YaleUniversity/Raredon Lab/iScience Transfer Submission/Submission2_Revisions_2025-12-__/Figures/Supplemental Figure 12")

# read both sheets from excel spreadsheet where ImageJ/FIJI data is compiled - one for each sheet
# Treat "N/A" as NA at import 
tenx <- read_excel(
  "Organoid_Macrophage_Distances_ImageJ.xlsx",
  sheet = "10x Image",
  na = "N/A")
# 20x image data
twentyx <- read_excel(
  "Organoid_Macrophage_Distances_ImageJ.xlsx",
  sheet = "20x Image",
  na = "N/A")

# keep only rows with real distances
tenx_clean <- tenx %>%
  filter(!is.na(Distance)) %>%
  mutate(SourceImage = "10x")
twentyx_clean <- twentyx %>%
  filter(!is.na(Distance)) %>%
  mutate(SourceImage = "20x")

# combine into one data frame
df_all <- bind_rows(tenx_clean, twentyx_clean)

# make sure 'Distance' is numeric (it should be, but doing to be safe)
df_all <- df_all %>%
  mutate(Distance = as.numeric(Distance))

# plot all distances pooled
p <- ggplot(df_all, aes(x = "", y = Distance)) +
  geom_boxplot(width = 0.25, outlier.shape = NA,
    alpha = 0.15, color = "grey50", lwd = 0.4) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.85,
    color = "#1B58A6") + # color for dots/points
  theme_classic(base_size = 16) +
  theme(axis.line = element_line(size = 0.6),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "", y = "Macrophage-Organoid Distance (µm)")
# print
p
# Save
setwd("~/Library/CloudStorage/OneDrive-YaleUniversity/Raredon Lab/iScience Transfer Submission/Submission2_Revisions_2025-12-__/Figures/Supplemental Figure 12")
ggsave("mac_org.distance_compiled.png", p, width=3.0, height= 6, units="in", dpi=600)

# if we want to plot the 10x vs 20x as separate columns
dist.plot_splitby_source.img <- ggplot(df_all, aes(x = SourceImage, y = Distance, color = SourceImage)) +
  geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.15,
    color = "grey50", lwd = 0.5) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.85) +
  scale_color_manual(values = c(
    "10x" = "#4A90E2", "20x" = "#003F7D")) +
  theme_classic(base_size = 16) +
  theme(axis.line = element_line(size = 0.6), panel.border = element_blank(),
    legend.position = "none", axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 12)),
    plot.margin = margin(5, 15, 5, 5)) +
  labs(x = "Source Image", y = "Macrophage-Organoid Distance (µm)")
# print
dist.plot_splitby_source.img
# Save
setwd("~/Library/CloudStorage/OneDrive-YaleUniversity/Raredon Lab/iScience Transfer Submission/Submission2_Revisions_2025-12-__/Figures/Supplemental Figure 12")
ggsave("mac_org.distance_split.by_sourceimg.png", dist.plot_splitby_source.img, width=5.0, height= 6, units="in", dpi=600)



