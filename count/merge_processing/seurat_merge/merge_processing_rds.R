library(Seurat)
library(dplyr)
library(ggplot2)

patient_cmo <- read.csv('/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/patient_cmo_timepoint.csv')
patient_metadata <- read.csv('/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/patient_metadata.csv')


seur_list <- readRDS('/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/merge_processing/seurat_patient_timepoint_cr.multi_azimuth_sctransform_list.rds')

for (x in 1:length(seur_list)) {seur_list[[x]]$cr.patient <- patient_cmo %>% filter(sample == seur_list[[x]]$cr.sample[1]) %>% pull(patient)}
for (x in 1:length(seur_list)) {seur_list[[x]]$cr.flare <- patient_metadata %>% filter(patient == patient_cmo %>% filter(sample == seur_list[[x]]$cr.sample[1]) %>% pull(patient)) %>% pull(flare)}
for (x in 1:length(seur_list)) {seur_list[[x]]$cr.group <- patient_metadata %>% filter(patient == patient_cmo %>% filter(sample == seur_list[[x]]$cr.sample[1]) %>% pull(patient)) %>% pull(group)}


#for (x in 1:length(seur_list)) {seur_list[[x]]$cr.patient <- names(which(table(seur_list[[x]]$patient) == max(table(seur_list[[x]]$patient))))}

# Merge normalized samples
seur <- merge(x = seur_list[[1]],
                       y = seur_list[2:length(seur_list)],
                       merge.data = TRUE)

DefaultAssay(seur) <- 'RNA'

seur <- NormalizeData(seur)
seur <- FindVariableFeatures(seur)
all.genes <- rownames(seur)
seur <- ScaleData(seur, features=all.genes)
seur <- SCTransform(seur)
seur <- RunPCA(seur)
seur <- FindNeighbors(seur, dims = 1:30)
seur <- FindClusters(seur, resolution = 0.8, algorithm=3, verbose = FALSE)
seur <- RunUMAP(seur, reduction = 'pca', dims = 1:30, assay = 'RNA')

saveRDS(seur, 'seurat_patient_timepoint_cr.multi_merge_sctransform.rds')

