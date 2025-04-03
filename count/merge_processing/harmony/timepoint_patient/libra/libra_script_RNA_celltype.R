library(Libra, lib='/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/R4.3.0/')
library(Seurat)

analysis1 <- readRDS('../analysis1/analysis1_subset.rds')
DefaultAssay(analysis1) <- 'RNA'


timepoint <- "1_1Tri"
label <- "cr.group"
replicate <- "cr.patient"
cell_type <- "predicted.celltype.l1"

full_results <- c()
for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
  full_results[[method]] <- run_de(analysis1, label_col=label, replicate_col=replicate, cell_type_col=cell_type, de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
}

for (method in names(full_results)) {
    write.csv(full_results[[method]], paste0(strsplit(timepoint, '_')[[1]][[2]], '_label-', paste(unique(analysis1[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
}

rm(analysis1)

timepoint1 <- readRDS('subset_1_1Tri.rds')
DefaultAssay(timepoint1) <- 'RNA'

timepoint <- "1_1Tri"
label <- "cr.flare"
replicate <- "cr.patient"
cell_type <- "predicted.celltype.l1"
comparisons <- c("control-after", "control-without", "after-without")
for (comparison in comparisons) {
  time.sub <- subset(timepoint1, cr.flare %in% strsplit(comparison, '-')[[1]])
  full_results <- c()
  for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
    full_results[[method]] <- run_de(time.sub, label_col='cr.flare', replicate_col='cr.patient', cell_type_col=cell_type, de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
  }
  for (method in names(full_results)) {
      write.csv(full_results[[method]], paste0(strsplit(timepoint, '_')[[1]][[2]], '_label-', paste(unique(time.sub[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
  }
}

rm(timepoint1)

timepoint2 <- readRDS('subset_2_2Tri.rds')
DefaultAssay(timepoint2) <- 'RNA'

timepoint <- "2_2Tri"
label <- "cr.flare"
replicate <- "cr.patient"
cell_type <- "predicted.celltype.l1"
comparisons <- c("control-before", "control-after", "control-without", "after-before", "after-without", "before-without")

for (comparison in comparisons) {
  time.sub <- subset(timepoint2, cr.flare %in% strsplit(comparison, '-')[[1]])
  full_results <- c()
  for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
    full_results[[method]] <- run_de(time.sub, label_col=label, replicate_col=replicate, cell_type_col=cell_type, de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
  }
  for (method in names(full_results)) {
      write.csv(full_results[[method]], paste0(strsplit(timepoint, '_')[[1]][[2]], '_label-', paste(unique(time.sub[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
  }
}

rm(timepoint2)

timepoint3 <- readRDS('subset_3_3Tri.rds')
DefaultAssay(timepoint3) <- 'RNA'
timepoint <- "3_3Tri"
label <- "cr.flare"
replicate <- "cr.patient"
cell_type <- "predicted.celltype.l1"
comparisons <- c("control-before", "control-after", "control-without", "after-before", "after-without", "before-without")

for (comparison in comparisons) {
  time.sub <- subset(timepoint3, cr.flare %in% strsplit(comparison, '-')[[1]])
  full_results <- c()
  for (method in c("DESeq2-LRT", "DESeq2-Wald", "limma-trend", "limma-voom")) {
    full_results[[method]] <- run_de(time.sub, label_col='cr.flare', replicate_col='cr.patient', cell_type_col=cell_type, de_method=strsplit(method, '-')[[1]][[1]], de_type=strsplit(method, '-')[[1]][[2]])
  }
  for (method in names(full_results)) {
      write.csv(full_results[[method]], paste0(strsplit(timepoint, '_')[[1]][[2]], '_label-', paste(unique(time.sub[[label]][[1]]), collapse='vs'), '_replicate-', replicate, '_celltype-', cell_type, '_method-', method, '_alldegs.csv'))
  }
}
