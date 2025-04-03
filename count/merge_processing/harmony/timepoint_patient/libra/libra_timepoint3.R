library(Libra, lib='/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/R4.3.0/')
library(Seurat)


timepoint3 <- readRDS('subset_3_3Tri.rds')

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

sessionInfo()
