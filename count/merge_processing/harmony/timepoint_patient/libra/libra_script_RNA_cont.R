library(Libra, lib='/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/R4.3.0/')
library(Seurat)

seur <- readRDS('../patient_timepoint-cr.multi_azimuth_sctransform_harmony_fix_crmetadata.rds')
DefaultAssay(seur) <- 'RNA'

seur$bulk <- 'bulk'

seur$cr.active.flare <- "None"
seur$cr.active.flare[seur$cr.patient %in% c('5', '6', '10')] <- 'control'
seur$cr.active.flare[seur$cr.patient %in% c('9', '11', '22', '23')] <- 'no_flare'
seur$cr.active.flare[seur$cr.patient %in% c('4', '7', '12', '18')] <- 'no_flare'
seur$cr.active.flare[seur$cr.patient == '4' & seur$cr.timepoint == '5_Post6M'] <- 'flare'
seur$cr.active.flare[seur$cr.patient == '12' & seur$cr.timepoint == '4_Post2M'] <- 'flare'
seur$cr.active.flare[seur$cr.patient == '14' & seur$cr.timepoint %in% c("1_1Tri", "3_3Tri")] <- 'no_flare'
seur$cr.active.flare[seur$cr.patient == '14' & seur$cr.timepoint %in% c("2_2Tri")] <- 'flare'
seur$cr.active.flare[seur$cr.patient == '3' & seur$cr.timepoint %in% c("3_3Tri", '5_Post6M')] <- 'no_flare'
seur$cr.active.flare[seur$cr.patient == '3' & seur$cr.timepoint %in% c("1_1Tri", "2_2Tri")] <- 'flare'
seur$cr.active.flare[seur$cr.patient == '20' & seur$cr.timepoint %in% c("3_3Tri", '5_Post6M')] <- 'no_flare'
seur$cr.active.flare[seur$cr.patient == '20' & seur$cr.timepoint %in% c("1_1Tri", "2_2Tri", "4_Post2M")] <- 'flare'

timepoint <- 'AllTimepoints'
label <- "cr.active.flare"
replicate <- "cr.patient"
cell_type <- "bulk"
comparisons <- c("control-flare", "control-no_flare", "no_flare-flare")

for (comparison in comparisons) {
  time.sub <- subset(seur, cr.active.flare %in% strsplit(comparison, '-')[[1]])
  full_results <- c()
  for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
    full_results[[method]] <- run_de(time.sub, label_col=label, replicate_col=replicate, cell_type_col=cell_type, de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
  }
  for (method in names(full_results)) {
      write.csv(full_results[[method]], paste0(timepoint, '_label-', paste(unique(time.sub[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
  }
}

timepoint <- c(1, 2, 3)
for (comparison in comparisons) {
  time.sub <- subset(seur, cr.active.flare %in% strsplit(comparison, '-')[[1]] & cr.timepoint %in% paste(sapply(c(1, 2, 3), function(x) paste0(x, '_', x, 'Tri')), sep='-'))
  full_results <- c()
  for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
    full_results[[method]] <- run_de(time.sub, label_col=label, replicate_col=replicate, cell_type_col=cell_type, de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
  }
  for (method in names(full_results)) {
      write.csv(full_results[[method]], paste0(paste(paste0(timepoint, 'Tri'), collapse='-'), '_label-', paste(unique(time.sub[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
  }
}
