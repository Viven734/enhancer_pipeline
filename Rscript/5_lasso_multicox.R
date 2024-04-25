lasso <- function(input_dir, output_dir, out_dir2, out_dir3) {
  ###Lasso回归：十折交叉检验筛选最佳λ并获取最终DEHGs
  #继续载入lasso回归所需R包：
  options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  if (!require(glmnet)) {
    install.packages("glmnet")
  }
  if (!require(Matrix)) {
    install.packages("Matrix")
  }
  if (!require(caret)) {
    install.packages("caret")
  }
  if (!require(survival)) {
    install.packages("survival")
  }
  if (!require(survminer)) {
    BiocManager::install("survminer")
  }
  library(survminer)
  library(glmnet)
  library(Matrix)
  library(caret)
  library(survival)
  #载入训练集队列及候选预后DEHGs：
  setwd(input_dir)
  load('04_diff_unicox.Rdata')
  #简化训练集队列命名：
  exp <- train_exp
  cl <- train_cl
  #构建候选预后DEHGs表达矩阵：
  exp <- exp[rownames(diff_unicox),]
  #使用lasso回归进一步收缩单因素cox筛选的25个候选DEHGs：
  x <- t(exp)
  y <- data.matrix(Surv(time = cl$OS, event = cl$Censor))
  x[1:6, 1:6]
  head(y)
  ###构建模型：
  fit <- glmnet(x, y, family = 'cox', type.measure = "deviance", nfolds = 10)
  pdf(file = paste(output_dir, "/train_unicox_p0.001_lasso-1.pdf", sep = ''), height = 15, width = 20)
  plot(fit, xvar = 'lambda', label = T) #候选DEHGs的lasso系数
  dev.off()
  #十折交叉检验筛选最佳lambda：
  set.seed(007)
  lasso_fit <- cv.glmnet(x, y, family = 'cox', type.measure = 'deviance', nfolds = 10)
  pdf(file = paste(output_dir, "/train_unicox_p0.001_lasso-2.pdf", sep = ''), height = 15, width = 20)
  plot(lasso_fit)
  dev.off()
  ##提取最佳λ值(这里选择min或lse对应lambda)：
  lambda.min <- lasso_fit$lambda.min
  lambda.min
  #拎出建模使用基因：
  model_lasso_min <- coef(lasso_fit, s = lambda.min)
  Active.Index <- which(model_lasso_min != 0)
  Active.model_lasso_min <- model_lasso_min[Active.Index]
  Active.Index
  Active.model_lasso_min
  row.names(model_lasso_min)[Active.Index] #筛选出15个
  # 从exp中选出了15个变量，命名为beta
  beta <- row.names(model_lasso_min)[Active.Index]
  # 从cl中提取OS和Censor列
  cl_subset <- cl[, c("OS", "Censor")]
  # 从exp中提取beta选出的15个变量
  # 对exp进行转置
  exp <- t(exp)
  exp_subset <- exp[, beta, drop = FALSE]
  # 将cl_subset和exp_subset合并成一个新的数据框gene_min
  gene_min <- cbind(cl_subset, exp_subset)
  write.table(gene_min, paste(output_dir, "/train_lasso_min_eRNAs_unicox_p0.001.txt", sep = ''), quote = FALSE, row.names = TRUE, col.names = NA, sep = '\t')
  rt = read.table(paste(output_dir, "/train_lasso_min_eRNAs_unicox_p0.001.txt", sep = ''), header = T, sep = "\t", check.names = F, row.names = 1)
  cox <- coxph(Surv(OS, Censor) ~ ., data = rt)
  cox = step(cox, direction = "both") #找最优的eRNA个数
  p <- ggforest(cox, data = rt) #森林图，survminer包功能
  ggsave(paste(output_dir, "/riskScore_train_eRNA_forest.pdf", sep = ''), p, width = 15, height = 10)
  riskScore = predict(cox, type = "risk", newdata = rt)
  risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
  write.table(cbind(id = rownames(cbind(rt[, 1:2], riskScore, risk)), cbind(rt[, 1:2], riskScore, risk)), paste(output_dir, "/riskScore_train_eRNA.txt", sep = ''), row.names = TRUE, quote = FALSE, sep = '\t')
  ############################################################################################测试集风险评分
  #先输入test_exp和test_cl矩阵，此处通过Import data导入
  exp2 <- t(test_exp) # 先对exp进行转置，再从exp中提取beta选出的15个变量
  exp_subset2 <- exp2[, beta, drop = FALSE]
  # 将cl_subset和exp_subset合并成一个新的数据框gene_min
  rt2 <- cbind(test_cl, exp_subset2)
  write.table(rt2, paste(out_dir2, "/test_multicox.txt", sep = ''), row.names = TRUE, quote = FALSE, sep = '\t')
  riskScore = predict(cox, type = "risk", newdata = rt2)
  risk2 = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
  write.table(cbind(id = rownames(cbind(rt2[, 1:2], riskScore, risk2)), cbind(rt2[, 1:2], riskScore, risk2)), paste(out_dir2, "/riskScore_test_eRNA.txt", sep = ''), row.names = TRUE, quote = FALSE, sep = '\t')

  #输入整体集合exp,cl矩阵，此处通过Import data导入
  exp <- read.delim("exp.txt", row.names = 1) #读取表达数据
  cl <- read.delim("cl.txt", row.names = 1) #读取生存数据
  # 先对exp进行转置，再从exp中提取beta选出的15个变量
  exp3 <- t(exp)
  exp_subset3 <- exp3[, beta, drop = FALSE]
  # 将cl_subset和exp_subset合并成一个新的数据框gene_min
  rt3 <- cbind(cl, exp_subset3)
  write.table(rt3, paste(out_dir3, "/combined_multicox.txt", sep = ''), row.names = TRUE, quote = FALSE, sep = '\t')
  riskScore = predict(cox, type = "risk", newdata = rt3)
  risk3 = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
  write.table(cbind(id = rownames(cbind(rt3[, 1:2], riskScore, risk3)), cbind(rt3[, 1:2], riskScore, risk3)), paste(out_dir3, "/riskScore_combined_eRNA.txt", sep = ''), row.names = TRUE, quote = FALSE, sep = '\t')
}