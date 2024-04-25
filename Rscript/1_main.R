##TCGA文件处理（HPV CESC数据）
setwd("/public6/hww/enhancer/TCGA/UVM")
library("rjson")
json <- jsonlite::fromJSON("metadata.cart.2024-03-11.json") #metadata文件名
sample_id <- sapply(json$associated_entities, function(x) { x[, 1] })
file_sample <- data.frame(sample_id, file_name = json$file_name)
count_file <- list.files('gdc_download_20240311_181008.535563', pattern = '*.tsv', recursive = TRUE) #Counts文件夹名
count_file_name <- strsplit(count_file, split = '/')
count_file_name <- sapply(count_file_name, function(x) { x[2] })
#下面的修改基因数
matrix = data.frame(matrix(nrow = 60660, ncol = 0))
#下面的修改样本例数
for (i in 1:80) {
  path = paste0('gdc_download_20240311_181008.535563/', count_file[i]) #Counts文件夹名
  data <- read.delim(path, fill = TRUE, header = FALSE, row.names = 1)
  colnames(data) <- data[2,]
  data <- data[-c(1:6),]
  data <- data[3]
  #数据类型，选择其中之一 3：unstranded；4：stranded_first；5：stranded_second；6：tpm_unstranded；7：fpkm_unstranded；8：fpkm_uq_unstranded
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name == count_file_name[i])]
  matrix <- cbind(matrix, data)
}
write.csv(matrix, 'UVM_Counts_matrix.csv', row.names = TRUE)