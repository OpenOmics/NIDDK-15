setwd("/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/cellchat_BD/trimester1")

# ----Library Install ----
libloc = "/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/R4.3.0_cellchat/"

# devtools::install_github('immunogenomics/presto', lib = libloc)
# devtools::install_github("jinworks/CellChat", lib = libloc)

library(CellChat, lib.loc = libloc)
library(Seurat, lib.loc = libloc)
library(patchwork)
library(dplyr)
options(stringsAsFactors = FALSE)

# ----Seurat Object import and subset ----
# seur <- readRDS("/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/merge_processing/harmony/timepoint_patient/patient_timepoint-cr.multi_azimuth_sctransform_harmony_fix_crmetadata.rds")
seur <- readRDS("/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/merge_processing/harmony/timepoint_patient/analysis1/analysis1_subset.rds")

# rename B-cells
b_cells <- read.csv("/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/merge_processing/seurat_merge/B_cells/Bcells_1Tri_FCRL1_metadata.csv", row.names = 1)

seur <- AddMetaData(seur, b_cells)
seur@meta.data <- seur@meta.data %>%
  mutate(predicted.celltype.l1.split = coalesce(FCRL1, predicted.celltype.l1))

# head(seur@meta.data)
#
# table(seur@meta.data$flare) # BD edit --> changed to cr.flare
# after  before control      NA without
# 29487    5909   20558     850   14256
seur.a <- subset(seur, subset = cr.flare == "after")
seur.b <- subset(seur, subset = cr.flare == "before")
seur.c <- subset(seur, subset = cr.flare == "control")
seur.w <- subset(seur, subset = cr.flare == "without")
rm(seur)

# head(seur.a@meta.data)

# ----Define CellChat functions ----

create.cellchat.obj <- function(seurat_obj, CellChatDB){
  cellchat <- createCellChat(object = seurat_obj, group.by = "predicted.celltype.l1.split", assay = "RNA")
  groupSize <- as.numeric(table(cellchat@idents))

  cellchat@DB <- CellChatDB.use

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  cellchat <- computeCommunProb(cellchat, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)

  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  
  return(cellchat)
}

plotCellChat <- function(cellchat){
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
}

# ----Database ----
CellChatDB <- CellChatDB.human
# showDatabaseCategory(CellChatDB)
# dplyr::glimpse(CellChatDB$interaction)
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
CellChatDB.use <- CellChatDB

# ---- Execute and plot -----
cellchat.a <- create.cellchat.obj(seur.a, CellChatDB.use)
saveRDS(cellchat.a, file = "Analysis1_after_cellchat.rds")
png("Cellchat_after.png", width=3200, height=1800, res = 300)
print(plotCellChat(cellchat.a))
dev.off()

cellchat.b <- create.cellchat.obj(seur.b, CellChatDB.use)
saveRDS(cellchat.b, file = "Analysis1_before_cellchat.rds")
png("Cellchat_before.png", width=3200, height=1800, res = 300)
print(plotCellChat(cellchat.b))
dev.off()

cellchat.c <- create.cellchat.obj(seur.c, CellChatDB.use)
saveRDS(cellchat.c, file = "Analysis1_control_cellchat.rds")
png("Cellchat_control.png", width=3200, height=1800, res = 300)
print(plotCellChat(cellchat.c))
dev.off()

cellchat.w <- create.cellchat.obj(seur.w, CellChatDB.use)
saveRDS(cellchat.w, file = "Analysis1_without_cellchat.rds")
png("Cellchat_without.png", width=3200, height=1800, res = 300)
print(plotCellChat(cellchat.w))
dev.off()

# ----Combine and save----
object.list <- list(after = cellchat.a,
                    before = cellchat.b,
                    control = cellchat.c,
                    without = cellchat.w)

cellchat.all <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(cellchat.all, file = "Analysis1_all_cellchat.rds")

capture.output(sessionInfo(), "Analysis1_cellchat_sessionInfo.txt")
