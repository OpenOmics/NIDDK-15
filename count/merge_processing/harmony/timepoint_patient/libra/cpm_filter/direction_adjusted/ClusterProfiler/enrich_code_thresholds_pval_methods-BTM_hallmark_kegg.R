library(clusterProfiler)
library(org.Hs.eg.db)
library(Seurat)
library(dplyr)
library(ggplot2)
library(optparse)
library(msigdbr)

option_list <- list(
  make_option(c("-p", "--pval"), type='double', action='store', default=NA,
              help = "Adjusted p-value threshold"),
  #make_option(c("-p", "--pop"), type='character', action='store', default='EUR',
  #            help = "Population used for downstream analysis, phenotypes extracted only if contains provided population. Default value is EUR"),
  make_option(c("-l", "--logfold"), type='double', action='store', default=NA,
              help = "Log 2 fold change threshold"),
  make_option(c("-t", "--timepoint"), type='character', action='store', default=NA,
              help = "Timepoints to consider, comma separated if multiple entries"),
  make_option(c("-m", "--method"), type='character', action='store', default=NA,
              help = "Methods to consider, comma separated if multiple entries"),
  make_option(c("-d", "--degdir"), type='character', action='store', default=NA,
              help = "Directory containing DEG results")
)

opt <- parse_args(OptionParser(option_list=option_list))

gene_info <- read.table(gzfile('/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/Sample_1/outs/filtered_feature_bc_matrix/features.tsv.gz', 'rt'), stringsAsFactors=FALSE)
gene_info$unique <- make.unique(gene_info$V2)

opt$timepoint = strsplit(opt$timepoint, ',')[[1]]
orgdb = "org.Hs.eg.db"
organism = "hsa"

m_h <- msigdbr(species = "Homo sapiens", category = "H")
hs <- org.Hs.eg.db

#downsamp <- readRDS('/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/merge_processing/harmony/timepoint_patient/dotplot/downsampled_celltype1000.rds')
#DefaultAssay(downsamp) <- 'RNA'
#universe <- rownames(downsamp)
#universe_entrez <- select(hs, keys=gene_info[match(universe, gene_info$unique),'V1'], columns= c('ENTREZID', "SYMBOL", "ENSEMBL"), keytype="ENSEMBL") %>% pull(ENTREZID)
#universe_entrez <- bitr(gene_info[match(universe, gene_info$unique),'V1'], fromType="ENSEMBL", toType="ENTREZID", OrgDb=orgdb) %>% pull(ENTREZID)

counts <- list()

btm <- read.table('/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/BTM/BTM_plus_data.tsv', sep='\t', header=TRUE)
modules <- c()
for(i in 1:dim(btm)[1]) {
  modules <- rbind(modules, cbind(btm[i,'name'], strsplit(btm[i,'genes'], ',')[[1]]))
}
modules <- as.data.frame(modules)

dir <- paste0('log2FC', opt$logfold, '_pvaladj', opt$pval, '_minBothPct', opt$pct)
if (!dir.exists(dir)) {
  dir.create(dir)
}

for (timepoint in opt$timepoint) {
  for (filename in Sys.glob(file.path(opt$degdir, paste0(timepoint, '_*', '_method-', opt$method, '*_alldegs.csv')))) {
    name <- basename(tools::file_path_sans_ext(filename))
    temp <- read.csv(filename)
    temp$ensembl <- gene_info[match(temp$gene, gene_info$unique),'V1']
    universe <- temp$gene
    universe_entrez <- bitr(gene_info[match(universe, gene_info$unique),'V1'], fromType="ENSEMBL", toType="ENTREZID", OrgDb=orgdb) %>% pull(ENTREZID)
    up <- temp %>% filter(avg_logFC > opt$logfold) %>% filter(p_val_adj < opt$pval)
    down <- temp %>% filter(avg_logFC < -opt$logfold) %>% filter(p_val_adj < opt$pval)
    counts[paste0(name, '_up')] <- dim(up)[1]
    counts[paste0(name, '_down')] <- dim(down)[1]
    print(paste(c(name, dim(up)[1], dim(down)[1]), collapse=' '))
    up_results <- c()
    down_results <- c()
    if (dim(up)[1] > 0){
      try({
        results <- enricher(up$gene, universe=universe, TERM2GENE=modules)
        up_results <- results
        write.csv(results@result, file.path(dir, paste0(name, ".enricher.BTM", ".Up.csv")))
        
        if (dim(results)[1] > 0) {
          png(file.path(dir, paste0(name, ".enricher.BTM", ".Up.DotPlot.png")), res=300, height=max(800, 250+(180*min(dim(results)[1], 20))), width=2000)
          print(dotplot(results, showCategory = 20))
          dev.off()
        }
      })
    }
    if (dim(down)[1] > 0){
      try({
        results <- enricher(down$gene, universe=universe, TERM2GENE=modules)
        down_results <- results
        write.csv(results@result, file.path(dir, paste0(name, ".enricher.BTM", ".Down.csv")))
        
        if (dim(results)[1] > 0) {
          png(file.path(dir, paste0(name, ".enricher.BTM", ".Down.DotPlot.png")), res=300, height=max(800, 250+(180*min(dim(results)[1], 20))), width=2000)
          print(dotplot(results, showCategory = 20))
          dev.off()
        }
      })
    }
    try({
      full_results <- c()
      if (!is.null(up_results))
        full_results <- head(mutate(up_results@result %>% filter(p.adjust < 0.05), logp.adjust = -log(p.adjust, base=10)), 10)
      if (!is.null(down_results)) {
        temp <- mutate(down_results@result %>% filter(p.adjust < 0.05), logp.adjust = -log(p.adjust, base=10))
        temp$logp.adjust <- -temp$logp.adjust
        full_results <- rbind(full_results, head(temp, 10))
      }
      if (!is.null(full_results)) {
        if (dim(full_results)[1] > 0) {
          png(file.path(dir, paste0(name, ".enricher.BTM.", "BarPlot.png")), res=300, height=2400, width=2400)
          print(ggplot(full_results, aes(x=reorder(Description, logp.adjust), y=logp.adjust, fill=logp.adjust)) +
                  geom_bar(stat="identity", alpha=.6, width=.4, aes(fill = logp.adjust < 0)) +
                  coord_flip() +
                  xlab("") +
                  theme_bw() + scale_fill_manual(guide = FALSE, breaks = c(TRUE, FALSE), values=c("blue", "red")))
          dev.off()
        }}
    })
    
    up_results <- c()
    down_results <- c()
    if (dim(up)[1] > 0){
      try({
        results <- enricher(up$gene, universe=universe, TERM2GENE=m_h[,c('gs_name', 'gene_symbol')])
        up_results <- results
        write.csv(results@result, file.path(dir, paste0(name, ".enricher.MSIGDB-HALLMARK", ".Up.csv")))
        
        if (dim(results)[1] > 0) {
          png(file.path(dir, paste0(name, ".enricher.MSIGDB-HALLMARK", ".Up.DotPlot.png")), res=300, height=max(800, 250+(180*min(dim(results)[1], 20))), width=2000)
          print(dotplot(results, showCategory = 20))
          dev.off()
        }
      })
    }
    if (dim(down)[1] > 0){
      try({
        results <- enricher(down$gene, universe=universe, TERM2GENE=m_h[,c('gs_name', 'gene_symbol')])
        down_results <- results
        write.csv(results@result, file.path(dir, paste0(name, ".enricher.MSIGDB-HALLMARK", ".Down.csv")))
        
        if (dim(results)[1] > 0) {
          png(file.path(dir, paste0(name, ".enricher.MSIGDB-HALLMARK", ".Down.DotPlot.png")), res=300, height=max(800, 250+(180*min(dim(results)[1], 20))), width=2000)
          print(dotplot(results, showCategory = 20))
          dev.off()
        }
      })
    }
    try({
      full_results <- c()
      if (!is.null(up_results))
        full_results <- head(mutate(up_results@result %>% filter(p.adjust < 0.05), logp.adjust = -log(p.adjust, base=10)), 10)
      if (!is.null(down_results)) {
        temp <- mutate(down_results@result %>% filter(p.adjust < 0.05), logp.adjust = -log(p.adjust, base=10))
        temp$logp.adjust <- -temp$logp.adjust
        full_results <- rbind(full_results, head(temp, 10))
      }
      if (!is.null(full_results)) {
        if (dim(full_results)[1] > 0) {
          png(file.path(dir, paste0(name, ".enricher.MSIGDB-HALLMARK.", "BarPlot.png")), res=300, height=2400, width=2400)
          print(ggplot(full_results, aes(x=reorder(Description, logp.adjust), y=logp.adjust, fill=logp.adjust)) +
                  geom_bar(stat="identity", alpha=.6, width=.4, aes(fill = logp.adjust < 0)) +
                  coord_flip() +
                  xlab("") +
                  theme_bw() + scale_fill_manual(guide = FALSE, breaks = c(TRUE, FALSE), values=c("blue", "red")))
          dev.off()
        }}
    })
    
    
    up_results <- c()
    down_results <- c()
    if (dim(up)[1] > 0){
      try({
        
        #results <- enrichKEGG(select(hs, keys=up$ensembl, columns= c('ENTREZID', "SYMBOL", "ENSEMBL"), keytype="ENSEMBL") %>% pull(ENTREZID), universe=universe_entrez, organism='hsa')
	results <- enrichKEGG(bitr(up$ensembl, fromType="ENSEMBL", toType="ENTREZID", OrgDb=orgdb) %>% pull("ENTREZID"), universe=universe_entrez, organism='hsa')
        up_results <- setReadable(results, hs, 'ENTREZID')
        write.csv(up_results@result, file.path(dir, paste0(name, ".enrichKEGG", ".Up.csv")))
        
        if (dim(results)[1] > 0) {
          png(file.path(dir, paste0(name, ".enrichKEGG", ".Up.DotPlot.png")), res=300, height=max(800, 250+(180*min(dim(results)[1], 20))), width=2000)
          print(dotplot(results, showCategory = 20))
          dev.off()
        }
      })
    }
    if (dim(down)[1] > 0){
      try({
        #results <- enrichKEGG(select(hs, keys=down$ensembl, columns= c('ENTREZID', "SYMBOL", "ENSEMBL"), keytype="ENSEMBL") %>% pull(ENTREZID), universe=universe_entrez, organism='hsa')
	results <- enrichKEGG(bitr(down$ensembl, fromType="ENSEMBL", toType="ENTREZID", OrgDb=orgdb) %>% pull("ENTREZID"), universe=universe_entrez, organism='hsa')
        down_results <- setReadable(results, hs, 'ENTREZID')
        write.csv(down_results@result, file.path(dir, paste0(name, ".enrichKEGG", ".Down.csv")))
        
        if (dim(results)[1] > 0) {
          png(file.path(dir, paste0(name, ".enrichKEGG", ".Down.DotPlot.png")), res=300, height=max(800, 250+(180*min(dim(results)[1], 20))), width=2000)
          print(dotplot(results, showCategory = 20))
          dev.off()
        }
      })
    }
    try({
      full_results <- c()
      if (!is.null(up_results))
        full_results <- head(mutate(up_results@result %>% filter(p.adjust < 0.05), logp.adjust = -log(p.adjust, base=10)), 10)
      if (!is.null(down_results)) {
        temp <- mutate(down_results@result %>% filter(p.adjust < 0.05), logp.adjust = -log(p.adjust, base=10))
        temp$logp.adjust <- -temp$logp.adjust
        full_results <- rbind(full_results, head(temp, 10))
      }
      if (!is.null(full_results)) {
        if (dim(full_results)[1] > 0) {
          png(file.path(dir, paste0(name, ".enrichKEGG.", "BarPlot.png")), res=300, height=2400, width=2400)
          print(ggplot(full_results, aes(x=reorder(Description, logp.adjust), y=logp.adjust, fill=logp.adjust)) +
                  geom_bar(stat="identity", alpha=.6, width=.4, aes(fill = logp.adjust < 0)) +
                  coord_flip() +
                  xlab("") +
                  theme_bw() + scale_fill_manual(guide = FALSE, breaks = c(TRUE, FALSE), values=c("blue", "red")))
          dev.off()
        }}
    })

  }
}

print(counts)
counts <- as.data.frame(counts)
write.csv(counts, file.path(dir, 'counts.csv'))
sink(file.path(dir, 'session.log'))
sessionInfo()
sink()
