gene_info <- read.table(gzfile('/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/Sample_1/outs/filtered_feature_bc_matrix/features.tsv.gz', 'rt'), stringsAsFactors=FALSE)
gene_info$unique <- make.unique(gene_info$V2)

setwd('20230817_delivery')

dir.create('apostrophe_name_trim')
dir.create('apostrophe_name')
for (filename in Sys.glob('*DESeq2-LRT*.csv')) {
  temp <- read.csv(filename)
#  temp$ensembl <- gene_info[match(temp$gene, gene_info$unique),'V1']
#  temp$gene <-  sapply(temp$gene, function(x) strsplit(x, '\\.')[[1]][[1]])
#  temp$rank_pval <- sign(temp$avg_logFC) * -log10(temp$p_val)
#  temp <- temp[order(temp$rank_pval, decreasing=T),]
  temp$gene <- paste0("'", temp$gene)
  write.csv(temp[intersect(grep('\\.', temp$gene, invert=TRUE), which(!is.na(temp$p_val))),], file.path('apostrophe_name_trim', basename(filename)), row.names=FALSE, quote=FALSE)
  write.csv(temp, file.path('apostrophe_name', basename(filename)), row.names=FALSE, quote=FALSE)
}

