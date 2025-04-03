library(Seurat)

seur <- readRDS('/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/merge_processing/harmony/timepoint_patient//patient_timepoint-cr.multi_azimuth_sctransform_harmony_fix_crmetadata.rds')

Idents(seur) <- "cr.timepoint"

for (timepoint in unique(seur$cr.timepoint)) {
	if (!timepoint %in% c("1_1Tri", "2_2Tri")) {
	sub <- subset(seur, idents=timepoint)
	saveRDS(sub, paste0('subset_', timepoint, '.rds'))
	}
}

sessionInfo()
