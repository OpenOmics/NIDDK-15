library(Seurat)
library(dplyr)
library(ggplot2)

patient_cmo <- read.csv('/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/patient_cmo_timepoint.csv')
patient_metadata <- read.csv('/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/patient_metadata.csv')


seur_full <- readRDS('seurat_patient_timepoint_cr.multi_merge_sctransform.rds')


library(colorBlindness, lib.loc='/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/R4.2.0/')
library(dichromat)

cell_cols <- c(colorBlindness::SteppedSequential5Steps, dichromat::colorschemes$SteppedSequential.5[16:21])[1:length(names(table(seur_full$predicted.celltype.l2)))]
names(cell_cols) <- names(table(seur_full$predicted.celltype.l2))


#write.table(grep('R', grep('-', grep('^IFN', row.names(seur), value=TRUE), value=TRUE, invert=TRUE), value=TRUE, invert=TRUE), 'seur_grepped_interferons.csv', row.names=FALSE, col.names=FALSE, quote=FALSE)

group <- 'cr.timepoint'
timepoints <- c("1_1Tri", "2_2Tri", "3_3Tri")
seur <- subset(seur_full, cells=which(seur_full$cr.timepoint %in% timepoints))

DefaultAssay(seur) <- 'RNA'

for (patient in unique(seur$cr.patient)) {
    #cols = ceiling(sqrt(dim(unique(seur[[group]]))[[1]]))
    #rows = dim(unique(seur[[group]]))[[1]] / cols
    cols = min(dim(unique(seur[[group]]))[[1]], 3)
    rows = 1

    sub <- subset(seur, cells = rownames(seur@meta.data)[which(seur$cr.patient == patient)])
    png(paste0('UMAP_SCT_CellType.l2_', patient, '_Timepoint1-3-Split.png'), width=1000+1300*cols, height=(1300*rows)+300, res = 300)
    print(DimPlot(sub, group.by = 'predicted.celltype.l2', split.by = group, ncol= cols, cols=cell_cols, repel = TRUE) + ggtitle(patient))
    dev.off()

    png(paste0('UMAP_SCT_CellType.l1_', patient, '_Timepoint1-3-Split.png'), width=1000+1300*cols, height=(1300*rows)+300, res = 300)
    print(DimPlot(sub, group.by = 'predicted.celltype.l1', split.by = group, ncol= cols, repel = TRUE) + ggtitle(patient))
    dev.off()

    for (gene in c("GPT", "SRY")) {
            png(paste0('FeaturePlot_RNA_', patient, '_', gene, '_Order_Timepoint1-3-Split.png'), width=1700*cols, height=(1300*rows)+300, res = 300)
            print(FeaturePlot(sub, feature=gene, order=TRUE, split.by=group, ncol=cols) & theme(legend.position = "right") )
            dev.off()
    }

}

sessionInfo()
