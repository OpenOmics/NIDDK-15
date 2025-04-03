library(Seurat)

seur <- readRDS('patient_timepoint-cr.multi_azimuth_sctransform_harmony_fix_crmetadata.rds')

Idents(seur) <- seur$predicted.celltype.l2

allmarkers <- FindAllMarkers(seur, assay='RNA', test.use='MAST', only.pos=TRUE)

saveRDS(allmarkers, 'FindAllMarkers_RNA_MAST_positiveOnly.rds')
