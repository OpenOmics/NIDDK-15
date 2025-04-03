library(Seurat)

aggr <- read.csv("/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/AggregatedDatasets/outs/aggregation.csv")
seur <- readRDS('/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/merge_processing/seurat_merge/seurat_patient_timepoint_cr.multi_merge_sctransform.rds')

metadata <- seur@meta.data[,keep_cols]
umap <- seur@reductions$umap@cell.embeddings

#Rename rows based on aggregation order
rows <- rownames(umap)
rows <- sapply(rows, function(x) gsub('1', which(aggr$sample_id == gsub('Library', 'Sample_', strsplit(x, '_')[[1]][1])), strsplit(x, '_')[[1]][2]))

#Rename rows in data 
rownames(metadata) <- rows
rownames(umap) <- rows

# put data into correct format for loupe browser
umap = cbind(rownames(umap), umap)
colnames(umap) = c('Barcode', 'UMAP-1', 'UMAP-2')
metadata = cbind('barcode' = rownames(metadata), metadata)

outdir <- '/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/aggr_seurat_export/20230901/seurat_merge'
write.table(umap, file = file.path(outdir, paste0('all_samples_merge_umap.csv')), quote = F, sep = ',', row.names = F, col.names = T)
