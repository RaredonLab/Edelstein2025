# Everything is in this folder on the server:
'/Volumes/Home-1/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/NICHES CircuitPlot Troubleshooting'

## Sophie adjusting Sam's troubleshooting code to make new plots for figure 6

#### Replicable Rotation Bug from Sophie 2025-04-25 ####
# close all, clear all
graphics.off()
rm(list = ls())
options(future.globals.maxSize = 16000 * 1024^2)
# Set seed
set.seed(2)
# pack
require(ggplot2)

# Load data
setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/NICHES CircuitPlot Troubleshooting/Regen Circuit - for MSBR")
load("regen.circuit.object.Robj")
load("regen.circuit_CTC.Robj")
regen.circuit.object # transcriptomic object
regen.subset_CTC_byCondition # NICHES object
table(regen.circuit.object$CellType.regen.spec) # Good, good.
table(regen.subset_CTC_byCondition$VectorType) # good.

# Load niches package from github
library(devtools)
setwd("/Users/msbr/GitHub/NICHESMethods")
getwd() # confirm that in the right place
devtools::load_all() # load local package
#check() # check that it works

# Set wd back
setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/NICHES CircuitPlot Troubleshooting")

# Specific cols from color palette for this object
regen.cols = c("Hillock_Like" = "#ff0000", "Pdgfrb+_Pericyte" = "#F4AFB4",
               "Polarized_Mac" = "#ff9e80", "Rspo3+_Mes" = "#89a5a5")

# Check out structure of the data
str(regen.circuit.object@meta.data) # Cell types to be used here are CellType.regen.spc
str(regen.subset_CTC_byCondition@meta.data)
table(regen.circuit.object$CellType.regen.spec) # Good. checks out
table(regen.circuit.object$Condition)
table(regen.subset_CTC_byCondition$VectorType)
table(regen.subset_CTC_byCondition$Condition.Joint)

# Specifying our "global" node list
global.node.list <- c("Hillock_Like", "Polarized_Mac",
                      "Rspo3+_Mes","Pdgfrb+_Pericyte")

# Angle 0 through 360
for(i in 0:361){
  angle <- i
  plot1 = CircuitPlot(transcr.obj = regen.circuit.object,
                      connect.obj = regen.subset_CTC_byCondition,
                      feature = 'Il1a—Il1r2',
                      group.by = 'CellType.regen.spec',
                      graph.angle = angle,
                      h = 0.02,
                      offset = 0.05,
                      autocrine.offset = 0.02,
                      edge.scale.factor = 1, # If you reduce this, you can highlight more pronounced signals?
                      arrow.head.angle = 15,
                      arrow.head.length = 0.03,
                      autocrine.arrow.curvature = 10,
                      cols.use = regen.cols,
                      global.node.list = global.node.list,
                      edge.fixed = T,
                      min.edge.value = NULL,
                      max.edge.value = NULL)
  #print(plot1)
  ggsave(paste("Il1a—Il1r2_CircuitPlot.subset_",angle,".png",sep=''), plot = plot1, width = 6, height = 5, dpi = 600)
}

# Just for checking that the plot is representing things correctly
node.object <- DefineNodeObject(transcr.obj = regen.circuit.object, # transcriptomic object
                                feature = 'Il1a')

edge.object <- DefineEdgeObject(connect.obj = regen.subset_CTC_byCondition, # NICHES object
                                feature = 'Il1a—Il1r2',
                                assay = 'CellToCell')

node.aggregate <- AggregateNodeData(node.object = node.object,
                                    group.by = 'CellType.regen.spec',
                                    global.node.list = unique(node.object$CellType.regen.spec))

edge.aggregate <- AggregateEdgeData(edge.object = edge.object,
                                    group.by = 'CellType.regen.spec.Joint') # Warnings but should be fine
View(node.aggregate)
View(edge.aggregate)

# Serpina-Lrp1
plot1 = CircuitPlot(transcr.obj = regen.circuit.object,
                    connect.obj = regen.subset_CTC_byCondition,
                    feature = 'Serpina1—Lrp1',
                    group.by = 'CellType.regen.spec',
                    graph.angle = 45.25,
                    h = 0.02,
                    offset = 0.05,
                    autocrine.offset = 0.02,
                    edge.scale.factor = 1, # If you reduce this, you can highlight more pronounced signals?
                    arrow.head.angle = 15,
                    arrow.head.length = 0.03,
                    autocrine.arrow.curvature = 10,
                    cols.use = regen.cols,
                    global.node.list = global.node.list,
                    edge.fixed = T,
                    min.edge.value = NULL,
                    max.edge.value = NULL,
                    unity.normalize = T) +
  guides(color = guide_legend(override.aes = list(size = 5)))  # Changes the dot size in legend
print(plot1)
setwd("~/Desktop/Manuscript Figures/Figure 6/Correct Plots")
ggsave("Serpina1_Lrp1_CircuitPlot.subset.png", plot = plot1, width = 8, height = 6.5, dpi = 600)

# Cd24-Selp
plot2 = CircuitPlot(transcr.obj = regen.circuit.object,
                    connect.obj = regen.subset_CTC_byCondition,
                    feature = 'Cd24—Selp',
                    group.by = 'CellType.regen.spec',
                    graph.angle = 45.25, # To address bug
                    h = 0.02,
                    offset = 0.05,
                    autocrine.offset = 0.02,
                    edge.scale.factor = 1, # If you reduce this, you can highlight more pronounced signals?
                    arrow.head.angle = 15,
                    arrow.head.length = 0.03,
                    autocrine.arrow.curvature = 10,
                    cols.use = regen.cols,
                    global.node.list = global.node.list,
                    edge.fixed = T,
                    min.edge.value = NULL,
                    max.edge.value = NULL,
                    unity.normalize = T) +
  guides(color = guide_legend(override.aes = list(size = 5))) 
print(plot2)
setwd("~/Desktop/Manuscript Figures/Figure 6/Correct Plots")
ggsave("Cd24_Selp_CircuitPlot.subset.png", plot = plot2, width = 8, height = 6.5, dpi = 600)

# Combine these two 
library(patchwork)
# Combine the two CircuitPlots and collect a single shared legend
combined_plot <- (plot1 | plot2) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")
# Save to file
ggsave("Serpina1_Cd24_CircuitPlots.png", plot = combined_plot, width = 12, height = 6, dpi = 600)

### Now make the Interleukins
# Il1a-Il1r2
plot3 = CircuitPlot(transcr.obj = regen.circuit.object,
                    connect.obj = regen.subset_CTC_byCondition,
                    feature = 'Il1a—Il1r2',
                    group.by = 'CellType.regen.spec',
                    graph.angle = 45.25,
                    h = 0.02,
                    offset = 0.05,
                    autocrine.offset = 0.02,
                    edge.scale.factor = 1, # If you reduce this, you can highlight more pronounced signals?
                    arrow.head.angle = 15,
                    arrow.head.length = 0.03,
                    autocrine.arrow.curvature = 10,
                    cols.use = regen.cols,
                    global.node.list = global.node.list,
                    edge.fixed = T,
                    min.edge.value = NULL,
                    max.edge.value = NULL,
                    unity.normalize = T)
print(plot3)
setwd("~/Desktop/Manuscript Figures/Figure 6/Correct Plots")
ggsave("Il1a_Il1r2_CircuitPlot.subset.png", plot = plot3, width = 8, height = 6.5, dpi = 600)

# Il1rn-Il1r1
plot4 = CircuitPlot(transcr.obj = regen.circuit.object,
                    connect.obj = regen.subset_CTC_byCondition,
                    feature = 'Il1rn—Il1r1',
                    group.by = 'CellType.regen.spec',
                    graph.angle = 45.25,
                    h = 0.02,
                    offset = 0.05,
                    autocrine.offset = 0.02,
                    edge.scale.factor = 1, # If you reduce this, you can highlight more pronounced signals?
                    arrow.head.angle = 15,
                    arrow.head.length = 0.03,
                    autocrine.arrow.curvature = 10,
                    cols.use = regen.cols,
                    global.node.list = global.node.list,
                    edge.fixed = T,
                    min.edge.value = NULL,
                    max.edge.value = NULL,
                    unity.normalize = T)
print(plot4)
setwd("~/Desktop/Manuscript Figures/Figure 6/Correct Plots")
ggsave("Il1rn_Il1r1_CircuitPlot.subset.png", plot = plot4, width = 8, height = 6.5, dpi = 600)

# Il1rn-Il1rl2
plot5 = CircuitPlot(transcr.obj = regen.circuit.object,
                    connect.obj = regen.subset_CTC_byCondition,
                    feature = 'Il1rn—Il1rl2',
                    group.by = 'CellType.regen.spec',
                    graph.angle = 45.25,
                    h = 0.02,
                    offset = 0.05,
                    autocrine.offset = 0.02,
                    edge.scale.factor = 1, # If you reduce this, you can highlight more pronounced signals?
                    arrow.head.angle = 15,
                    arrow.head.length = 0.03,
                    autocrine.arrow.curvature = 10,
                    cols.use = regen.cols,
                    global.node.list = global.node.list,
                    edge.fixed = T,
                    min.edge.value = NULL,
                    max.edge.value = NULL,
                    unity.normalize = T)
print(plot5)
setwd("~/Desktop/Manuscript Figures/Figure 6/Correct Plots")
ggsave("Il1rn_Il1rl2_CircuitPlot.subset_larger.png", plot = plot5, width = 9, height = 7.5, dpi = 600)

# Combine the IL-1 family plots
combined_plot2 <- (plot3 | plot4 | plot5) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")
setwd("~/Desktop/Manuscript Figures/Figure 6/Correct Plots")
ggsave("regen.subset_Interleukins_larger.png", plot = combined_plot2, width = 12, height = 4, dpi = 600)


#### Now for the Mdk, EPhrins, Eregs, etc.
### Now make the Interleukins
# Mdk—Tspan1
plot6 = CircuitPlot(transcr.obj = regen.circuit.object,
                    connect.obj = regen.subset_CTC_byCondition,
                    feature = 'Mdk—Tspan1',
                    group.by = 'CellType.regen.spec',
                    graph.angle = 45.25,
                    h = 0.02,
                    offset = 0.05,
                    autocrine.offset = 0.02,
                    edge.scale.factor = 1, # If you reduce this, you can highlight more pronounced signals?
                    arrow.head.angle = 15,
                    arrow.head.length = 0.03,
                    autocrine.arrow.curvature = 10,
                    cols.use = regen.cols,
                    global.node.list = global.node.list,
                    edge.fixed = T,
                    min.edge.value = NULL,
                    max.edge.value = NULL,
                    unity.normalize = T)
print(plot6)
setwd("~/Desktop/Manuscript Figures/Figure 6/Correct Plots")
ggsave("Mdk—Tspan1_CircuitPlot.subset.png", plot = plot6, width = 8, height = 6.5, dpi = 600)

# Efnb2—Ephb3
plot7 = CircuitPlot(transcr.obj = regen.circuit.object,
                    connect.obj = regen.subset_CTC_byCondition,
                    feature = 'Efnb2—Ephb3',
                    group.by = 'CellType.regen.spec',
                    graph.angle = 45.25,
                    h = 0.02,
                    offset = 0.05,
                    autocrine.offset = 0.02,
                    edge.scale.factor = 1, # If you reduce this, you can highlight more pronounced signals?
                    arrow.head.angle = 15,
                    arrow.head.length = 0.03,
                    autocrine.arrow.curvature = 10,
                    cols.use = regen.cols,
                    global.node.list = global.node.list,
                    edge.fixed = T,
                    min.edge.value = NULL,
                    max.edge.value = NULL,
                    unity.normalize = T)
print(plot7)
setwd("~/Desktop/Manuscript Figures/Figure 6/Correct Plots")
ggsave("Efnb2—Ephb3_CircuitPlot.subset.png", plot = plot7, width = 8, height = 6.5, dpi = 600)

# Ereg—Erbb3
plot8 = CircuitPlot(transcr.obj = regen.circuit.object,
                    connect.obj = regen.subset_CTC_byCondition,
                    feature = 'Ereg—Erbb3',
                    group.by = 'CellType.regen.spec',
                    graph.angle = 45.25,
                    h = 0.02,
                    offset = 0.05,
                    autocrine.offset = 0.02,
                    edge.scale.factor = 1, # If you reduce this, you can highlight more pronounced signals?
                    arrow.head.angle = 15,
                    arrow.head.length = 0.03,
                    autocrine.arrow.curvature = 10,
                    cols.use = regen.cols,
                    global.node.list = global.node.list,
                    edge.fixed = T,
                    min.edge.value = NULL,
                    max.edge.value = NULL,
                    unity.normalize = T)
print(plot8)
setwd("~/Desktop/Manuscript Figures/Figure 6/Correct Plots")
ggsave("Ereg—Erbb3_CircuitPlot.subset.png", plot = plot8, width = 8, height = 6.5, dpi = 600)

# Combine the 3 previous plots
combined_plot3 <- (plot6 | plot7 | plot8) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")
setwd("~/Desktop/Manuscript Figures/Figure 6/Correct Plots")
ggsave("ereg.mdk.eph.combined.png", plot = combined_plot3, width = 12, height = 4, dpi = 600)

# Fn1-Itgav
plot9 = CircuitPlot(transcr.obj = regen.circuit.object,
                    connect.obj = regen.subset_CTC_byCondition,
                    feature = 'Fn1—Itgav',
                    group.by = 'CellType.regen.spec',
                    graph.angle = 45.25,
                    h = 0.02,
                    offset = 0.05,
                    autocrine.offset = 0.02,
                    edge.scale.factor = 1, # If you reduce this, you can highlight more pronounced signals?
                    arrow.head.angle = 15,
                    arrow.head.length = 0.03,
                    autocrine.arrow.curvature = 10,
                    cols.use = regen.cols,
                    global.node.list = global.node.list,
                    edge.fixed = T,
                    min.edge.value = NULL,
                    max.edge.value = NULL,
                    unity.normalize = T)
print(plot9)
setwd("~/Desktop/Manuscript Figures/Figure 6/Correct Plots")
ggsave("Fn1_itgav_subset.png", plot = plot9, width = 8, height = 6.5, dpi = 600)

# Fn1-Itgav
plot10 = CircuitPlot(transcr.obj = regen.circuit.object,
                    connect.obj = regen.subset_CTC_byCondition,
                    feature = 'Tgfa—Egfr',
                    group.by = 'CellType.regen.spec',
                    graph.angle = 45.25,
                    h = 0.02,
                    offset = 0.05,
                    autocrine.offset = 0.02,
                    edge.scale.factor = 1, # If you reduce this, you can highlight more pronounced signals?
                    arrow.head.angle = 15,
                    arrow.head.length = 0.03,
                    autocrine.arrow.curvature = 10,
                    cols.use = regen.cols,
                    global.node.list = global.node.list,
                    edge.fixed = T,
                    min.edge.value = NULL,
                    max.edge.value = NULL,
                    unity.normalize = T)
print(plot10)
setwd("~/Desktop/Manuscript Figures/Figure 6/Correct Plots")
ggsave("Tgfa.Egfr_subset.png", plot = plot10, width = 8, height = 6.5, dpi = 600)

# Combine
combined_plot4 <- (plot9 | plot10) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")
setwd("~/Desktop/Manuscript Figures/Figure 6/Correct Plots")
ggsave("tri_circuit_combined.png", plot = combined_plot4, width = 12, height = 6, dpi = 600)
