render_one <- function(timepoint, comparison, savedir) {
  if (!dir.exists(savedir)) {
    dir.create(savedir, recursive=TRUE)
  }
  rmarkdown::render(
    '/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/merge_processing/harmony/timepoint_patient/libra/cpm_filter/direction_adjusted/LibraDEGreport_BTM_Hallmark_KEGG_params.Rmd',
    output_file = file.path(savedir, paste0('NIDDK-15_Pseudobulk_DEG_Enrichment_', timepoint, '_', comparison, '_Report_20230817.html')),
    params=list( libradir='/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/merge_processing/harmony/timepoint_patient/libra/cpm_filter/direction_adjusted/', enrichdir='/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/merge_processing/harmony/timepoint_patient/libra/cpm_filter/direction_adjusted/ClusterProfiler', timepoint=timepoint, comparison=comparison, groupinfo='/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/merge_processing/harmony/timepoint_patient/libra/direction_adjusted/pseudobulk_groups_patient_sample_info.tsv'),
    envir = new.env()
  )
}

savedir <- '/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/merge_processing/harmony/timepoint_patient/libra/cpm_filter/direction_adjusted/20230817_delivery'

for (filename in Sys.glob('2Tri_*DESeq2-LRT*')) {
  timepoint <- strsplit(filename, '_')[[1]][[1]]
  comparison <- strsplit(strsplit(filename, '_')[[1]][[2]], '-')[[1]][[2]]
  render_one(timepoint, comparison, savedir)
}
for (filename in Sys.glob('1Tri_*DESeq2-LRT*')) {
  timepoint <- strsplit(filename, '_')[[1]][[1]]
  comparison <- strsplit(strsplit(filename, '_')[[1]][[2]], '-')[[1]][[2]]
  render_one(timepoint, comparison, savedir)
}
for (filename in Sys.glob('3Tri_*DESeq2-LRT*')) {
  timepoint <- strsplit(filename, '_')[[1]][[1]]
  comparison <- strsplit(strsplit(filename, '_')[[1]][[2]], '-')[[1]][[2]]
  render_one(timepoint, comparison, savedir)
}

