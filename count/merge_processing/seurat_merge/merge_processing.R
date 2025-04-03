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


library(colorBlindness, lib.loc='/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/R4.2.0/')
library(dichromat)

cell_cols <- c(colorBlindness::SteppedSequential5Steps, dichromat::colorschemes$SteppedSequential.5[16:21])[1:length(names(table(seur$predicted.celltype.l2)))]
names(cell_cols) <- names(table(seur$predicted.celltype.l2))


#write.table(grep('R', grep('-', grep('^IFN', row.names(seur), value=TRUE), value=TRUE, invert=TRUE), value=TRUE, invert=TRUE), 'seur_grepped_interferons.csv', row.names=FALSE, col.names=FALSE, quote=FALSE)

interferons <- read.csv('seur_grepped_interferons.csv', header=FALSE)
modules <- list(interferons = interferons$V1)
seur <- AddModuleScore(seur, modules, assay='RNA', name='module')
colnames(seur@meta.data) <- gsub('module1', 'interferon', colnames(seur@meta.data))

group <- 'cr.timepoint'

DefaultAssay(seur) <- 'RNA'

for (patient in unique(seur$cr.patient)) {
    cols = ceiling(sqrt(dim(unique(seur[[group]]))[[1]]))
    rows = dim(unique(seur[[group]]))[[1]] / cols

    sub <- subset(seur, cells = rownames(seur@meta.data)[which(seur$cr.patient == patient)])
    png(paste0('UMAP_SCT_CellType_', patient, '_Timepoint-Split.png'), width=1000+1300*cols, height=(1300*rows)+300, res = 300)
    print(DimPlot(sub, group.by = 'predicted.celltype.l2', split.by = group, ncol= cols, cols=cell_cols, repel = TRUE) + ggtitle(patient))
    dev.off()

    for (gene in c("GPT", "SRY")) {
	    png(paste0('FeaturePlot_RNA_', patient, '_', gene, '_Timepoint-Split.png'), width=1700*cols, height=(1300*rows)+300, res = 300)
	    print(FeaturePlot(sub, feature=gene, split.by=group, ncol=cols) & theme(legend.position = "right") )
	    dev.off()
    }

    png(paste0('FeaturePlot_ModuleScore_', patient, '_Interferon_Timepoint-Split.png'), width=1700*cols, height=(1300*rows)+300, res = 300)
    print(FeaturePlot(sub, feature='interferon', split.by=group, ncol=cols) & theme(legend.position = "right") )
    dev.off()

}

