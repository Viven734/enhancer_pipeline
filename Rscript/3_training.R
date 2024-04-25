#训练集和测试集
training_and_testing <- function(output_dir) {
  options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  if (!require("BiocManager")) {
    install.packages("BiocManager")
    library(BiocManager)
  }
  if (!require(Matrix)) {
    install.packages("Matrix")
    library(Matrix)
  }
  if (!require(caret)) {
    install.packages("caret")
    library(caret)
  }
  if (!require(survival)) {
    install.packages("survival")
    library(survival)
  }

  library(Matrix)
  library(caret)
  library(survival)
  #读取表达（creat包按列分割，因此患者为横坐标）和临床txt文档
  exp <- read.delim(paste(output_dir, "/exp.txt", sep = ''), row.names = 1) #读取表达数据
  cl <- read.delim(paste(output_dir, "/cl.txt", sep = ''), row.names = 1) #读取生存数据
  # dat = data.matrix(exp)
  # lable = Surv(cl$OS, cl$Censor)
  #1.训练集和测试集队列拆分(7:3)：
  set.seed(820)
  inA <- createDataPartition(cl$Censor, p = 0.7, times = 1, list = F)
  #训练集：
  train_exp <- exp[, inA]
  train_cl <- cl[inA,]
  #测试集：
  test_exp <- exp[, -inA]
  test_cl <- cl[-inA,]
  #切割后两个队列患者生死比例保持基本一致：
  prop.table(table(train_cl$status))
  prop.table(table(test_cl$status))
  #导出矩阵保存
  write.csv(train_exp, paste(output_dir, "/train_exp.csv", sep = ''), row.names = TRUE)
  write.csv(train_cl, paste(output_dir, "/train_cl.csv", sep = ''), row.names = TRUE)
  write.csv(test_exp, paste(output_dir, "/test_exp.csv", sep = ''), row.names = TRUE)
  write.csv(test_cl, paste(output_dir, "/test_cl.csv", sep = ''), row.names = TRUE)
  # write.table(train_exp, paste(output_dir, "/train_exp.txt", sep = ''), quote=FALSE,sep = "\t")
  # write.table(train_cl, paste(output_dir, "/train_cl.txt", sep = ''), quote=FALSE,sep = "\t")
  # write.table(test_exp, paste(output_dir, "/test_exp.txt", sep = ''), quote=FALSE,sep = "\t")
  # write.table(test_cl, paste(output_dir, "/test_cl.txt", sep = ''), quote=FALSE,sep = "\t")
}