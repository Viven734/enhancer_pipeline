KM <- function(input_dir, output_dir) {
  library("ggplot2")
  library("ggpubr")
  library("survminer")
  library("survival")
  setwd(input_dir)
  data = read.table("riskScore_test_eRNA.txt", header = T, sep = "\t", check.names = F, row.names = 1)
  # Fit survival curves
  fit <- survfit(Surv(OS, Censor) ~ risk2, data = data)
  # Plot informative survival curves
  p <- ggsurvplot(fit, data = data,
             # 标题
             title = "Testing set",
             # 增加p值显示
             pval = TRUE, pval.method = F,
             # 更改legend标题
             legend.title = "Risksocre",
             legend.labs = c("High", "Low"),
             # 更改legend标签
             palette = c("firebrick2", "dodgerblue3"),
             risk.table = TRUE,
             cumevents = TRUE,
             tables.height = 0.15,
             tables.theme = theme_cleantable(),
             tables.y.text = FALSE,
             xlab = "Survival(Days)",
             break.time.by = 2000)
  # 调整颜色
  #更改X轴标签
  # 修改X轴标尺刻度大小
  #ggexport(p,paste(output_dir, "/Training_set.pdf", sep = ''))
  ggsave(paste(output_dir, "/Testing_set.pdf", sep = ''), plot = p$plot)

}