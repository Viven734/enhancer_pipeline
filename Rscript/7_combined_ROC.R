ROC <- function(input_dir, output_dir) {
  if (!require(survivalROC)) {
    BiocManager::install("survivalROC")
  }
  setwd(input_dir)
  library(survivalROC)
  rt = read.table("riskScore_combined_eRNA.txt", header = T, sep = "\t", check.names = F, row.names = 1)
  pdf(file = paste(output_dir, "/Combined_ROC.pdf", sep = ''))
  par(oma = c(0.5, 1, 0, 1), font.lab = 1.5, font.axis = 1.5)
  roc = survivalROC(Stime = rt$OS, status = rt$Censor, marker = rt$riskScore,
                    predict.time = 2000, method = "KM")
  plot(roc$FP, roc$TP, type = "l", xlim = c(0, 1), ylim = c(0, 1), col = 'red',
       xlab = "False positive rate", ylab = "True positive rate",
       main = paste("ROC curve(", "AUC = ", round(roc$AUC, 3), ")"),
       lwd = 2, cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.2, font = 1.2)
  abline(0, 1)
  dev.off()
}