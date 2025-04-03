setwd("/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/cellchat_BD/trimester2")

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
seur <- readRDS("trimester2_subset.rds")

# ---- Rename flare groups (None to before, control to uninfected)
seur$cr.flare_renamed <- seur$cr.flare
seur$cr.flare_renamed[seur$cr.flare == 'None'] <- 'before'
seur$cr.flare_renamed[seur$cr.flare == 'control'] <- 'uninfected'

print(table(seur@meta.data[,c('cr.flare', 'cr.flare_renamed')]))

seur.a <- subset(seur, subset = cr.flare_renamed == "after")
seur.b <- subset(seur, subset = cr.flare_renamed == "before")
seur.c <- subset(seur, subset = cr.flare_renamed == "uninfected")
seur.w <- subset(seur, subset = cr.flare_renamed == "without")

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
saveRDS(cellchat.a, file = "trimester2_after_cellchat.rds")
png("Cellchat_after.png", width=3200, height=1800, res = 300)
print(plotCellChat(cellchat.a))
dev.off()

cellchat.b <- create.cellchat.obj(seur.b, CellChatDB.use)
saveRDS(cellchat.b, file = "trimester2_before_cellchat.rds")
png("Cellchat_before.png", width=3200, height=1800, res = 300)
print(plotCellChat(cellchat.b))
dev.off()

cellchat.c <- create.cellchat.obj(seur.c, CellChatDB.use)
saveRDS(cellchat.c, file = "trimester2_control_cellchat.rds")
png("Cellchat_control.png", width=3200, height=1800, res = 300)
print(plotCellChat(cellchat.c))
dev.off()

cellchat.w <- create.cellchat.obj(seur.w, CellChatDB.use)
saveRDS(cellchat.w, file = "trimester2_without_cellchat.rds")
png("Cellchat_without.png", width=3200, height=1800, res = 300)
print(plotCellChat(cellchat.w))
dev.off()

# ----Combine and save----
object.list <- list(after = cellchat.a,
                    before = cellchat.b,
                    control = cellchat.c,
                    without = cellchat.w)

cellchat.all <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(cellchat.all, file = "trimester2_all_cellchat.rds")

capture.output(sessionInfo(), "trimester2_cellchat_sessionInfo.txt")
