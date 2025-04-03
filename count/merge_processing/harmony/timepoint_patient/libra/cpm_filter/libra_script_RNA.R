library(Libra, lib='/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/R4.3.0/')
library(Seurat)
library(edgeR)
library(stringr)

make_matrix <- function(seur, sample) {
        rowSums(GetAssayData(seur, slot='count', assay='RNA')[,seur$cr.sample == sample])
}

analysis1 <- readRDS('../../analysis1/analysis1_subset.rds')

temp <- sapply(str_sort(unique(analysis1$cr.sample), numeric=TRUE), function(x) make_matrix(analysis1, x))
analysis1 <- subset(analysis1, features=names(which(rowSums(cpm(temp) > 0.5) > dim(temp)[[2]] / 2)))

#seur$bulk <- 'bulk'
analysis1$bulk <- 'bulk'

DefaultAssay(analysis1) <- 'RNA'

full_results <- c()
#for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
for (method in c("DESeq2-LRT")) {
  full_results[[method]] <- run_de(analysis1, label_col='cr.group', replicate_col='cr.patient', cell_type_col='bulk', de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
}

timepoint <- "1_1Tri"
label <- "cr.group"
replicate <- "cr.patient"
cell_type <- "bulk"
for (method in names(full_results)) {
    write.csv(full_results[[method]], paste0(strsplit(timepoint, '_')[[1]][[2]], '_label-', paste(unique(analysis1[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
}


timepoint <- "1_1Tri"
label <- "cr.group"
replicate <- "cr.patient"
cell_type <- "predicted.celltype.l1"

full_results <- c()
#for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
for (method in c("DESeq2-LRT")) {
  full_results[[method]] <- run_de(analysis1, label_col=label, replicate_col=replicate, cell_type_col=cell_type, de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
}

for (method in names(full_results)) {
    write.csv(full_results[[method]], paste0(strsplit(timepoint, '_')[[1]][[2]], '_label-', paste(unique(analysis1[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
}


timepoint1 <- readRDS('../subset_1_1Tri.rds')
DefaultAssay(timepoint1) <- 'RNA'

print(timepoint1)
print(analysis1)

temp <- sapply(str_sort(unique(timepoint1$cr.sample), numeric=TRUE), function(x) make_matrix(timepoint1, x))
timepoint1 <- subset(timepoint1, features=names(which(rowSums(cpm(temp) > 0.5) > dim(temp)[[2]] / 2)))


timepoint1$bulk <- 'bulk'
timepoint <- "1_1Tri"
label <- "cr.flare"
replicate <- "cr.patient"
cell_type <- "bulk"
comparisons <- c("control-after", "control-without", "after-without")
for (comparison in comparisons) {
  time.sub <- subset(timepoint1, cr.flare %in% strsplit(comparison, '-')[[1]])
  full_results <- c()
#  for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
  for (method in c("DESeq2-LRT")) {
    full_results[[method]] <- run_de(time.sub, label_col='cr.flare', replicate_col='cr.patient', cell_type_col='bulk', de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
  }
  for (method in names(full_results)) {
      write.csv(full_results[[method]], paste0(strsplit(timepoint, '_')[[1]][[2]], '_label-', paste(unique(time.sub[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
  }
}


timepoint1$bulk <- 'bulk'
timepoint <- "1_1Tri"
label <- "cr.flare"
replicate <- "cr.patient"
cell_type <- "predicted.celltype.l1"
comparisons <- c("control-after", "control-without", "after-without")
for (comparison in comparisons) {
  time.sub <- subset(timepoint1, cr.flare %in% strsplit(comparison, '-')[[1]])
  full_results <- c()
#  for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
  for (method in c("DESeq2-LRT")) {
    full_results[[method]] <- run_de(time.sub, label_col=label, replicate_col=replicate, cell_type_col=cell_type,  de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
  }
  for (method in names(full_results)) {
      write.csv(full_results[[method]], paste0(strsplit(timepoint, '_')[[1]][[2]], '_label-', paste(unique(time.sub[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
  }
}


rm(analysis1)
rm(timepoint1)
gc()

timepoint2 <- readRDS('../subset_2_2Tri.rds')

temp <- sapply(str_sort(unique(timepoint2$cr.sample), numeric=TRUE), function(x) make_matrix(timepoint2, x))
timepoint2 <- subset(timepoint2, features=names(which(rowSums(cpm(temp) > 0.5) > dim(temp)[[2]] / 2)))

DefaultAssay(timepoint2) <- 'RNA'

timepoint2$bulk <- 'bulk'
timepoint <- "2_2Tri"
label <- "cr.flare"
replicate <- "cr.patient"
cell_type <- "bulk"
comparisons <- c("control-before", "control-after", "control-without", "after-before", "after-without", "before-without")

for (comparison in comparisons) {
  time.sub <- subset(timepoint2, cr.flare %in% strsplit(comparison, '-')[[1]])
  full_results <- c()
#  for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
  for (method in c("DESeq2-LRT")) {
    full_results[[method]] <- run_de(time.sub, label_col=label, replicate_col=replicate, cell_type_col=cell_type, de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
  }
  for (method in names(full_results)) {
      write.csv(full_results[[method]], paste0(strsplit(timepoint, '_')[[1]][[2]], '_label-', paste(unique(time.sub[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
  }
}


timepoint <- "2_2Tri"
label <- "cr.flare"
replicate <- "cr.patient"
cell_type <- "predicted.celltype.l1"
comparisons <- c("control-before", "control-after", "control-without", "after-before", "after-without", "before-without")

for (comparison in comparisons) {
  time.sub <- subset(timepoint2, cr.flare %in% strsplit(comparison, '-')[[1]])
  full_results <- c()
#  for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
  for (method in c("DESeq2-LRT")) {
    full_results[[method]] <- run_de(time.sub, label_col=label, replicate_col=replicate, cell_type_col=cell_type, de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
  }
  for (method in names(full_results)) {
      write.csv(full_results[[method]], paste0(strsplit(timepoint, '_')[[1]][[2]], '_label-', paste(unique(time.sub[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
  }
}
rm(timepoint2)
gc()

timepoint3 <- readRDS('../subset_3_3Tri.rds')

temp <- sapply(str_sort(unique(timepoint3$cr.sample), numeric=TRUE), function(x) make_matrix(timepoint3, x))
timepoint3 <- subset(timepoint3, features=names(which(rowSums(cpm(temp) > 0.5) > dim(temp)[[2]] / 2)))

DefaultAssay(timepoint3) <- 'RNA'

timepoint3$bulk <- 'bulk'
timepoint <- "3_3Tri"
label <- "cr.flare"
replicate <- "cr.patient"
cell_type <- "bulk"
comparisons <- c("control-before", "control-after", "control-without", "after-before", "after-without", "before-without")

for (comparison in comparisons) {
  time.sub <- subset(timepoint3, cr.flare %in% strsplit(comparison, '-')[[1]])
  full_results <- c()
#  for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
  for (method in c("DESeq2-LRT")) {
    full_results[[method]] <- run_de(time.sub, label_col='cr.flare', replicate_col='cr.patient', cell_type_col='bulk', de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
  }
  for (method in names(full_results)) {
      write.csv(full_results[[method]], paste0(strsplit(timepoint, '_')[[1]][[2]], '_label-', paste(unique(time.sub[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
  }
}


timepoint3$bulk <- 'bulk'
timepoint <- "3_3Tri"
label <- "cr.flare"
replicate <- "cr.patient"
cell_type <- "predicted.celltype.l1"
comparisons <- c("control-before", "control-after", "control-without", "after-before", "after-without", "before-without")

for (comparison in comparisons) {
  time.sub <- subset(timepoint3, cr.flare %in% strsplit(comparison, '-')[[1]])
  full_results <- c()
#  for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
  for (method in c("DESeq2-LRT")) {
    full_results[[method]] <- run_de(time.sub, label_col=label, replicate_col=replicate, cell_type_col=cell_type, de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
  }
  for (method in names(full_results)) {
      write.csv(full_results[[method]], paste0(strsplit(timepoint, '_')[[1]][[2]], '_label-', paste(unique(time.sub[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
  }
}

seur <- readRDS('../../patient_timepoint-cr.multi_azimuth_sctransform_harmony_fix_crmetadata.rds')

temp <- sapply(str_sort(unique(seur$cr.sample), numeric=TRUE), function(x) make_matrix(seur, x))
seur <- subset(seur, features=names(which(rowSums(cpm(temp) > 0.5) > dim(temp)[[2]] / 2)))

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
#  for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
  for (method in c("DESeq2-LRT")) {
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
#  for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
  for (method in c("DESeq2-LRT")) {
    full_results[[method]] <- run_de(time.sub, label_col=label, replicate_col=replicate, cell_type_col=cell_type, de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
  }
  for (method in names(full_results)) {
      write.csv(full_results[[method]], paste0(paste(paste0(timepoint, 'Tri'), collapse='-'), '_label-', paste(unique(time.sub[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
  }
}


timepoint <- 'AllTimepoints'
label <- "cr.active.flare"
replicate <- "cr.patient"
cell_type <- "predicted.celltype.l1"
comparisons <- c("control-flare", "control-no_flare", "no_flare-flare")

for (comparison in comparisons) {
  time.sub <- subset(seur, cr.active.flare %in% strsplit(comparison, '-')[[1]])
  full_results <- c()
#  for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
  for (method in c("DESeq2-LRT")) {
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
#  for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
  for (method in c("DESeq2-LRT")) {
    full_results[[method]] <- run_de(time.sub, label_col=label, replicate_col=replicate, cell_type_col=cell_type, de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
  }
  for (method in names(full_results)) {
      write.csv(full_results[[method]], paste0(paste(paste0(timepoint, 'Tri'), collapse='-'), '_label-', paste(unique(time.sub[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
  }
}

