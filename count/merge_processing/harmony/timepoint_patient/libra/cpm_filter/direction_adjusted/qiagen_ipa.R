gene_info <- read.table(gzfile('/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/Sample_1/outs/filtered_feature_bc_matrix/features.tsv.gz', 'rt'), stringsAsFactors=FALSE)
gene_info$unique <- make.unique(gene_info$V2)

#temp <- read.csv('2Tri_label-aftervswithout_patients-4,7,12,18vs9,11,22,23_replicate-cr.patient_celltype-bulk_method-DESeq2-LRT_alldegs.csv')

for (filename in Sys.glob('*DESeq2-LRT*.csv')) {
  temp <- read.csv(filename)
  temp$ensembl <- gene_info[match(temp$gene, gene_info$unique),'V1']
  temp$gene <- paste0("'", temp$gene)
  
  write.table(temp[!is.na(temp$p_val_adj),c('ensembl', 'gene', 'avg_logFC', 'p_val', 'p_val_adj')], file.path('qiagen_ipa', paste0(basename(tools::file_path_sans_ext(filename)), '_qiagen.txt')), row.names=FALSE, quote=FALSE, sep='\t')
}

