

# Ciruit Plot Demonstration for Mes <--> Imm Signaling (Reviewer 1 Requests)

setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/NICHES CircuitPlot Troubleshooting")

options(future.globals.maxSize = 16000 * 1024^2)
# Set seed
set.seed(2)
# pack
require(ggplot2)

# Load data
setwd("~/Desktop/iScience Transfer Submission/Additional Analyses/Mesenchyme_Immune_CircuitAnalysis")
load("mes.imm.obj.Robj")
load("mes.imm.eng.CTC.bySample.Robj")

# Load niches package from github
library(devtools)
# setwd("/Users/msbr/GitHub/NICHESMethods")
setwd("~/Documents/GitHub/NICHESMethods")
getwd() # confirm that in the right place
devtools::load_all() # load local package
#check() # check that it works

# Set wd back
setwd("/Volumes/Home/RaredonLab-CC1126-MEDANE/Raredon_Lab_Internal_Collaboration/Organoid Project/NICHES CircuitPlot Troubleshooting")

# Specific cols from color palette for this object
mes.imm.cols = c("Pdgfrb+_Pericyte" = "#F4AFB4",
               "Polarized_Mac" = "#ff9e80", "Rspo3+_Mes" = "#89a5a5")

# Check out structure of the data
str(mes.imm.obj1@meta.data) # Cell types to be used here are CellType.regen.spc
str(mes.imm.eng.CTC.bySample@meta.data)
table(mes.imm.obj1$CellType.sub) # Good. checks out
table(mes.imm.obj1$Condition)
table(mes.imm.eng.CTC.bySample$VectorType)
table(mes.imm.eng.CTC.bySample$Condition.Joint)

# Specifying our "global" node list
global.node.list <- c("Polarized_Mac",
                      "Rspo3+_Mes","Pdgfrb+_Pericyte")

# Angle 0 through 360
for(i in 0:361){
  angle <- i
  plot1 = CircuitPlot(transcr.obj = mes.imm.obj1,
                      connect.obj = mes.imm.eng.CTC.bySample,
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
node.object <- DefineNodeObject(transcr.obj = mes.imm.obj1, # transcriptomic object
                                feature = 'Ncam1')

edge.object <- DefineEdgeObject(connect.obj = mes.imm.eng.CTC.bySample, # NICHES object
                                feature = 'Ncam1—Ptpra',
                                assay = 'CellToCell')

node.aggregate <- AggregateNodeData(node.object = node.object,
                                    group.by = 'Condition',
                                    global.node.list = unique(node.object$CellType.sub))

edge.aggregate <- AggregateEdgeData(edge.object = edge.object,
                                    group.by = 'CellType.sub.Joint') # Warnings but should be fine
# View(node.aggregate)
# View(edge.aggregate)

plot1 = CircuitPlot(transcr.obj = mes.imm.obj1,
                    connect.obj = mes.imm.eng.CTC.bySample,
                    feature = 'Ncam1—Ptpra',
                    group.by = 'CellType.sub',
                    graph.angle = 30,
                    h = 0.02,
                    offset = 0.05,
                    autocrine.offset = 0.02,
                    edge.scale.factor = 1, # If you reduce this, you can highlight more pronounced signals?
                    arrow.head.angle = 15,
                    arrow.head.length = 0.03,
                    autocrine.arrow.curvature = 10,
                    cols.use = mes.imm.cols,
                    global.node.list = global.node.list,
                    edge.fixed = T,
                    min.edge.value = NULL,
                    max.edge.value = NULL,
                    unity.normalize = T) +
  guides(color = guide_legend(override.aes = list(size = 5)))  # Changes the dot size in legend
print(plot1)
# Save
setwd("~/Desktop/iScience Transfer Submission/Additional Analyses/Mesenchyme_Immune_CircuitAnalysis/Misc Circuit Plots")
ggsave("Ncam1—Ptpra_CircuitPlot.subset.png", plot = plot1, width = 8, height = 6.5, dpi = 600)

Idents(mes.imm.obj1) = mes.imm.obj1$CellType.sub
FeaturePlot(mes.imm.obj1, features = c("Vcam1","Itgb2"), label = T, split.by = "Condition")
FeaturePlot(mes.imm.eng.CTC.bySample, features = c("Rspo3—Sdc4"))
