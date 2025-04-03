library(Seurat)
library(stringr)

seur <- readRDS('patient_timepoint-cr.multi_azimuth_sctransform_harmony_fix_crmetadata.rds')

data <- seur@meta.data[,c('predicted.celltype.l2', 'cr.group', 'cr.patient', 'cr.timepoint')]
data$label <- paste(data$cr.patient, data$cr.timepoint, sep='-')

tmp <- dplyr::rename(dplyr::count(data, label, predicted.celltype.l2), Freq=n)
cnt_output <- tidyr::pivot_wider(tmp, names_from = label, values_from = Freq)
cnt_output[is.na(cnt_output)] <- 0


str_sort(colnames(cnt_output)[2:61], numeric=TRUE)

cnt_reorder <- cnt_output[,c('predicted.celltype.l2', str_sort(colnames(cnt_output)[2:61], numeric=TRUE))]
colnames(cnt_reorder) <- gsub('^', 'Patient ', gsub('-', ' ', colnames(cnt_reorder)))
colnames(cnt_reorder)[1] <- 'predicted.celltype.l2'

cnt_trim <- cnt_reorder[,grep('-', colnames(cnt_reorder))]
colnames(cnt_trim) <- gsub('^', 'Patient ', gsub('-', ' ', colnames(cnt_trim)))
cnt_trim <- as.data.frame(cnt_trim)
rownames(cnt_trim) <- cnt_reorder$predicted.celltype.l2
write.csv(cnt_trim, "All_Samples_CellType_Count_l2.csv")

write.csv(sweep(cnt_trim, 2, colSums(cnt_trim), `/`), "All_Samples_CellType_Percentage_l2.csv")

patients <- c(str_sort(unique(data$cr.patient[data$cr.group == 2]), numeric=TRUE), str_sort(unique(data$cr.patient[data$cr.group == 3]), numeric=TRUE))

patient_index <- unlist(sapply(patients, function(x) grep(paste0('Patient ', x), colnames(cnt_trim), value=TRUE)))

write.csv(cnt_trim[,patient_index], "All_Samples_Ordered_CellType_Count_l2.csv")

percentage <- sweep(cnt_trim, 2, colSums(cnt_trim), `/`)
write.csv(percentage[, patient_index], "All_Samples_Ordered_CellType_Percentage_l2.csv")


patients <- str_sort(unique(data$cr.patient[data$cr.group == 2]), numeric=TRUE)

patient_index <- unlist(sapply(patients, function(x) grep(paste0('Patient ', x), colnames(cnt_trim), value=TRUE)))

write.csv(cnt_trim[,patient_index], "Group2_Samples_CellType_Count_l2.csv")

percentage <- sweep(cnt_trim, 2, colSums(cnt_trim), `/`)
write.csv(percentage[, patient_index], "Group2_Samples_CellType_Percentage_l2.csv")


patients <- str_sort(unique(data$cr.patient[data$cr.group == 3]), numeric=TRUE)

patient_index <- unlist(sapply(patients, function(x) grep(paste0('Patient ', x), colnames(cnt_trim), value=TRUE)))

write.csv(cnt_trim[,patient_index], "Group3_Samples_CellType_Count_l2.csv")

percentage <- sweep(cnt_trim, 2, colSums(cnt_trim), `/`)
write.csv(percentage[, patient_index], "Group3_Samples_CellType_Percentage_l2.csv")

