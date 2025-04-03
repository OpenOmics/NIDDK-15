libloc = "/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/R4.3.0_cellchat/"

# devtools::install_github('immunogenomics/presto', lib = libloc)
# devtools::install_github("jinworks/CellChat", lib = libloc)

library(CellChat, lib.loc = libloc)
library(Seurat, lib.loc = libloc)

# ---- Seurat Object import and subset
full <- readRDS('/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/merge_processing/harmony/timepoint_patient/patient_timepoint-cr.multi_azimuth_sctransform_harmony_fix_crmetadata.rds')

tri2 <- subset(full, subset = cr.timepoint == '2_2Tri')
tri3 <- subset(full, subset = cr.timepoint == '3_3Tri')

# ---- Rename B cells
b_cells <- read.csv("/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/merge_processing/seurat_merge/B_cells/Bcells_2Tri_FCRL1_metadata.csv", row.names = 1)

tri2 <- AddMetaData(tri2, b_cells)
tri2@meta.data <- tri2@meta.data %>%
  mutate(predicted.celltype.l1.split = coalesce(FCRL1, predicted.celltype.l1))


b_cells <- read.csv("/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/merge_processing/seurat_merge/B_cells/Bcells_3Tri_FCRL1_metadata.csv", row.names = 1)

tri3 <- AddMetaData(tri3, b_cells)
tri3@meta.data <- tri3@meta.data %>%
  mutate(predicted.celltype.l1.split = coalesce(FCRL1, predicted.celltype.l1))

saveRDS(tri2, 'trimester2/trimester2_subset.rds')
saveRDS(tri3, 'trimester3/trimester3_subset.rds')

