drug <- function(currentdir,output_dir) {
  options(stringsAsFactors = F)
  library(oncoPredict)
  library(data.table)
  library(gtools)
  library(reshape2)
  library(ggpubr)
  setwd(output_dir)
  th = theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
  dir = paste(currentdir,'/DataFiles/Training_Data/',sep='')
  GDSC2_Expr = readRDS(file = paste(dir, 'GDSC2_Expr_RMA.rds',sep=''))
  GDSC2_Res = readRDS(file = paste(dir, "GDSC2_Res.rds",sep=''))
  GDSC2_Res <- exp(GDSC2_Res)
  ###导入测试集
  testExpr <- read.delim(paste(output_dir, "/diff_gene_exp.txt", sep = ''), row.names = 1)
  testExpr <- as.matrix(testExpr)
  ###运行得到样本IC50
  calcPhenotype(trainingExprData = GDSC2_Expr,
                trainingPtype = GDSC2_Res,
                testExprData = testExpr,
                batchCorrect = 'eb', # "eb" for ComBat
                powerTransformPhenotype = TRUE,
                removeLowVaryingGenes = 0.2,
                minNumSamples = 10,
                printOutput = TRUE,
                removeLowVaringGenesFrom = 'rawData')
}