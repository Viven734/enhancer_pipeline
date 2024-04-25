cox <- function(output_dir) {
  options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  if (!require(survival)) {
    install.packages("survival")
    library(survival)
  }
  library(survival)
  #用一个循环对矩阵进行批量单因素cox
  train_exp <- read.csv(paste(output_dir, "/train_exp.csv", sep = ''), row.names = 1)
  train_cl <- read.csv(paste(output_dir, "/train_cl.csv", sep = ''), row.names = 1)
  test_exp <- read.csv(paste(output_dir, "/test_exp.csv", sep = ''), row.names = 1)
  test_cl <- read.csv(paste(output_dir, "/test_cl.csv", sep = ''), row.names = 1)
  cox <- apply(
    train_exp, 1, function(x) {
      train_cl$gene <- as.numeric(x)
      cox_genes <- coxph(Surv(OS, Censor) ~ gene, data = train_cl)
      coef <- coef(cox_genes) #回归系数
      SE <- sqrt(diag(vcov(cox_genes))) #标准误
      HR <- exp(coef) #风险比
      cox_need <- cbind(HR = HR,
                        HR.95L = exp(coef - qnorm(.975, 0, 1) * SE),
                        HR.95H = exp(coef + qnorm(.975, 0, 1) * SE),
                        pvalue = 1 - pchisq((coef / SE)^2, 1))
      return(cox_need['gene',])
    }
  )
  unicox <- t(cox)
  head(unicox)
  #提取预后差异基因列表(p<0.001)：
  diff_unicox <- unicox[unicox[, 4] < 0.001,]
  dim(diff_unicox) #从7936个差异eRNA中筛选出了36个候选预后DEHGs
  table(diff_unicox[, 1] < 1) #？个风险因子，？个保护因子
  head(diff_unicox)
  write.csv(diff_unicox, paste(output_dir, "/diff_unicox.csv", sep = ''), row.name = TRUE) #导出单因素cox，p<0.001的差异eRNAs矩阵
  #保存所需数据并清空工作环境：
  setwd(output_dir)
  save(train_exp, train_cl, test_exp, test_cl, diff_unicox, file = '04_diff_unicox.Rdata')
  rm(list = ls())
}