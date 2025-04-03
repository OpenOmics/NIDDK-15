library(Seurat)
library(stringr)
library(FactoMineR)
library(sp)
library(edgeR)
library(colorBlindness, lib.loc='/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/R4.2.0/')
library(ggplot2)
library(PCAtools, lib.loc='/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/R4.2.0/')


patient_cmo_timepoint <- read.table('/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/patient_cmo_timepoint.csv', sep=',', header=TRUE)
patient_metadata <- read.table('/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/patient_metadata.csv', sep=',', header=TRUE)
library_patient <- read.table('/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/library_patient_full.csv', sep=',', header=TRUE)

seur <- readRDS('all_samples_marked_active-flare.rds')

sample_count_matrices <- sapply(str_sort(unique(seur$cr.sample), numeric=TRUE),function(x) NULL)

make_matrix <- function(cell_type) {
	if (sum(seur$cr.sample == sample & seur$predicted.celltype.l1 == cell_type) > 1) {
		rowSums(GetAssayData(seur, slot='count', assay='RNA')[,seur$cr.sample == sample & seur$predicted.celltype.l1 == cell_type])
	} else {
		GetAssayData(seur, slot='count', assay='RNA')[,seur$cr.sample == sample & seur$predicted.celltype.l1 == cell_type]
	}
}

for (sample in names(sample_count_matrices)) {
	#sample_count_matrices[[sample]] <- sapply(str_sort(unique(seur$predicted.celltype.l1[seur$cr.sample == sample])), function(cell_type) rowSums(GetAssayData(seur, slot='count', assay='RNA')[,seur$cr.sample == sample & seur$predicted.celltype.l1 == cell_type]))
	sample_count_matrices[[sample]] <- sapply(str_sort(unique(seur$predicted.celltype.l1[seur$cr.sample == sample])), make_matrix)
	write.csv(sample_count_matrices[[sample]], paste0('RNA_counts_rowSum_sample-', sample, '_celltype-l1.csv'))
}

for (sample in names(sample_count_matrices)) {
	colnames(sample_count_matrices[[sample]]) <- paste0(sample, '-', colnames(sample_count_matrices[[sample]]))
	
}

full_matrix <- do.call(cbind, sample_count_matrices)
write.csv(full_matrix, 'RNA_counts_rowSum_sample-all_celltype-l1.csv')


cols <- colorBlindness::SteppedSequential5Steps[seq(3, 25, 5)]
names(cols) <- c("before", "after", "without", "control", "None")

get_sample_flare <- function(names) {
	sapply(strsplit(names, '-'), function(x) c(patient_metadata$flare[patient_metadata$patient == patient_cmo_timepoint$patient[patient_cmo_timepoint$sample == x[[1]]]], patient_cmo_timepoint$timepoint[patient_cmo_timepoint$sample == x[[1]]]))
}

strsplit(colnames(full_matrix), '-')
full_metadata <- t(as.data.frame(get_sample_flare(colnames(full_matrix))))

full_matrix <- read.csv('RNA_counts_rowSum_sample-all_celltype-l1.csv', row.names=1)

get_sample_flare <- function(names) {
        sapply(strsplit(names, '\\.'), function(x) c(patient_metadata$flare[patient_metadata$patient == patient_cmo_timepoint$patient[patient_cmo_timepoint$sample == x[[1]]]], patient_cmo_timepoint$timepoint[patient_cmo_timepoint$sample == x[[1]]]))
}
full_metadata <- t(as.data.frame(get_sample_flare(gsub('X', '', colnames(full_matrix)))))
colnames(full_matrix) <- gsub('X', '', colnames(full_matrix))
rownames(full_metadata) <- colnames(full_matrix)

full_metadata <- cbind(full_metadata, sapply(rownames(full_metadata), function(x) paste(strsplit(x, '\\.')[[1]][2:length(strsplit(x, '\\.')[[1]])], collapse=' ')))

full_metadata <- cbind(full_metadata, patient_cmo_timepoint$patient[match(sapply(rownames(full_metadata), function(x) strsplit(x, '\\.')[[1]][[1]]), patient_cmo_timepoint$sample)])

colnames(full_metadata) <- c("flare", "timepoint", "cell type", "patient")
full_metadata[full_metadata[,'flare'] == 'None' & full_metadata[,'timepoint'] != '1_1Tri', 'flare'] <- "before"


celltype_cols <- colorBlindness::SteppedSequential5Steps[seq(1, 25, floor(25/length(unique(full_metadata[,"cell type"]))))][1:length(unique(full_metadata[,"cell type"]))]
names(celltype_cols) <- unique(full_metadata[,"cell type"])

p <- pca(cpm(full_matrix[,grep('other$', colnames(full_matrix), perl=TRUE, invert=TRUE)], log=TRUE))
biplot(p, lab=grep('other$', colnames(full_matrix), perl=TRUE, invert=TRUE, value=TRUE))
#metadata <- t(as.data.frame(get_sample_flare(grep('other$', colnames(full_matrix), perl=TRUE, invert=TRUE, value=TRUE))))
#flare_status <- as.data.frame(get_sample_flare(grep('other$', colnames(full_matrix), perl=TRUE, invert=TRUE, value=TRUE)))
#colnames(flare_status) <- 'flare'
#rownames(flare_status) <- grep('other$', colnames(full_matrix), perl=TRUE, invert=TRUE, value=TRUE)
p <- pca(cpm(full_matrix[,grep('other$', colnames(full_matrix), perl=TRUE, invert=TRUE)], log=TRUE), metadata=flare_status)
biplot(p, lab=grep('other$', colnames(full_matrix), perl=TRUE, invert=TRUE, value=TRUE), colby='flare', colkey=cols)


flare_status <- as.data.frame(get_sample_flare(grep('other$', colnames(full_matrix), perl=TRUE, invert=TRUE, value=TRUE)))
colnames(flare_status) <- 'flare'
rownames(flare_status) <- grep('other$', colnames(full_matrix), perl=TRUE, invert=TRUE, value=TRUE)
p <- pca(cpm(full_matrix[,grep('other$', colnames(full_matrix), perl=TRUE, invert=TRUE)], log=TRUE), metadata=flare_status)
biplot(p, lab=grep('other$', colnames(full_matrix), perl=TRUE, invert=TRUE, value=TRUE), colby='flare', colkey=cols)

p <- pca(cpm(full_matrix[,grep('other$', names(which(full_metadata[,'timepoint'] == '2_2Tri')), perl=TRUE, invert=TRUE, value=TRUE)], log=TRUE), metadata=full_metadata[grep('other$', names(which(full_metadata[,'timepoint'] == '2_2Tri')), perl=TRUE, invert=TRUE, value=TRUE),])
biplot(p, colby='flare', colkey=cols, shape='cell type', legendPosition = "right", title='Timepoint 2')
biplot(p, colby='cell type', colkey=celltype_cols, shape='flare', legendPosition = "right", title='Timepoint 2')


generate_plot <- function(full_matrix, full_metadata, used_rows) {
  p <- pca(cpm(full_matrix[,used_rows], log=TRUE), metadata=full_metadata[used_rows,])
  return(p)
}


biplot(pca.paec, lab = NULL, colby = 'Strain', shape = 'Time',
       legendPosition = "right", title = "PAE cells",
       gridlines.major = F, gridlines.minor = F)
