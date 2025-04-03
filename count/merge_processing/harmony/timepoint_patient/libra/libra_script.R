library(Libra, lib='/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/R4.3.0/')
library(Seurat)

analysis1 <- readRDS('../analysis1/analysis1_subset.rds')
seur <- readRDS('../patient_timepoint-cr.multi_azimuth_sctransform_harmony_fix_crmetadata.rds')


seur$bulk <- 'bulk'
analysis1$bulk <- 'bulk'

full_results <- c()
for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
  full_results[[method]] <- run_de(analysis1, label_col='cr.group', replicate_col='cr.patient', cell_type_col='bulk', de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
}

timepoint <- "1_1Tri"
label <- "cr.group"
replicate <- "cr.patient"
cell_type <- "bulk"
for (method in names(full_results)) {
    write.csv(full_results[[method]], paste0(strsplit(timepoint, '_')[[1]][[2]], '_label-', paste(unique(analysis1[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
}


timepoint1$bulk <- 'bulk'
timepoint <- "1_1Tri"
label <- "cr.flare"
replicate <- "cr.patient"
cell_type <- "bulk"
comparisons <- c("control-after", "control-without", "after-without")
for (comparison in comparisons) {
  time.sub <- subset(timepoint1, cr.flare %in% strsplit(comparison, '-')[[1]])
  full_results <- c()
  for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
    full_results[[method]] <- run_de(time.sub, label_col='cr.flare', replicate_col='cr.patient', cell_type_col='bulk', de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
  }
  for (method in names(full_results)) {
      write.csv(full_results[[method]], paste0(strsplit(timepoint, '_')[[1]][[2]], '_label-', paste(unique(time.sub[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
  }
}

timepoint2 <- readRDS('subset_2_2Tri.rds')

#timepoint2$cr.flare[timepoint2$cr.patient %in% c(9, 11, 22, 23)] <- 'without'
#timepoint2$cr.flare[timepoint2$cr.patient %in% c(4, 7, 12, 18)] <- 'after'
timepoint2$cr.flare[timepoint2$cr.patient %in% c(14, 3, 20)] <- 'before'
#timepoint2$cr.flare[timepoint2$cr.patient %in% c(5, 6, 10)] <- 'control'

saveRDS(timepoint2, 'subset_2_2Tri.rds')


timepoint2$bulk <- 'bulk'
timepoint <- "2_2Tri"
label <- "cr.flare"
replicate <- "cr.patient"
cell_type <- "bulk"
comparisons <- c("control-before", "control-after", "control-without", "after-before", "after-without", "before-without")

for (comparison in comparisons) {
  time.sub <- subset(timepoint2, cr.flare %in% strsplit(comparison, '-')[[1]])
  full_results <- c()
  for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
    full_results[[method]] <- run_de(time.sub, label_col='cr.flare', replicate_col='cr.patient', cell_type_col='bulk', de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
  }
  for (method in names(full_results)) {
      write.csv(full_results[[method]], paste0(strsplit(timepoint, '_')[[1]][[2]], '_label-', paste(unique(time.sub[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
  }
}



timepoint3 <- readRDS('subset_3_3Tri.rds')

#timepoint2$cr.flare[timepoint2$cr.patient %in% c(9, 11, 22, 23)] <- 'without'
#timepoint2$cr.flare[timepoint2$cr.patient %in% c(4, 7, 12, 18)] <- 'after'
timepoint3$cr.flare[timepoint3$cr.patient %in% c(14, 3, 20)] <- 'before'
#timepoint2$cr.flare[timepoint2$cr.patient %in% c(5, 6, 10)] <- 'control'

saveRDS(timepoint3, 'subset_3_3Tri.rds')


timepoint3$bulk <- 'bulk'
timepoint <- "3_3Tri"
label <- "cr.flare"
replicate <- "cr.patient"
cell_type <- "bulk"
comparisons <- c("control-before", "control-after", "control-without", "after-before", "after-without", "before-without")

for (comparison in comparisons) {
  time.sub <- subset(timepoint3, cr.flare %in% strsplit(comparison, '-')[[1]])
  full_results <- c()
  for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
    full_results[[method]] <- run_de(time.sub, label_col='cr.flare', replicate_col='cr.patient', cell_type_col='bulk', de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
  }
  for (method in names(full_results)) {
      write.csv(full_results[[method]], paste0(strsplit(timepoint, '_')[[1]][[2]], '_label-', paste(unique(time.sub[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
  }
}
