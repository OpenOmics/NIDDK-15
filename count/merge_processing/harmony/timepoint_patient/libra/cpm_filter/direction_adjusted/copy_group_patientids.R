library(dplyr)

new_dir <- '20230817_delivery'
info <- read.table('pseudobulk_groups_patient_sample_info.tsv', sep='\t', header=TRUE)
for (time in unique(info$Timepoint)) {
  for (filename in grep('bulk', Sys.glob(paste0(gsub('._', '', time), '_label*')), value=TRUE)) {
    groups <- strsplit(strsplit(strsplit(filename, '-')[[1]][[2]], '_')[[1]][[1]], 'vs')[[1]]
    new_file <- gsub('_replicate', paste0('_patients-', paste(sapply(groups, function(x) info %>% filter(Timepoint == time & Group == x) %>% pull(Patients)), collapse='vs'), '_replicate'), filename)
    file.copy(filename, file.path(new_dir, new_file))
  }
}
