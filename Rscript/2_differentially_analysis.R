##差异分析
differentially_analysis <- function(output_dir) {
  options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  if (!require("BiocManager")) {
    install.packages("BiocManager")
    library(BiocManager)
  }
  if (!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
  }
  if (!require(limma)) {
    BiocManager::install("limma")
    library(limma)
  }
  if (!require(pheatmap)) {
    install.packages("pheatmap")
    library(pheatmap)
  }
  if (!require(ggsci)) {
    install.packages("ggsci")
    library(ggsci)
  }
  if (!require(dplyr)) {
    install.packages("dplyr")
    library(dplyr)
  }
  if (!require(ggrepel)) {
    install.packages("ggrepel")
    library(ggrepel)
  }
  #读取差异表达样本，正常样本放前面，癌症样本放后面。
  data = read.table(paste(output_dir, "/all_sample_eRNA.txt", sep = ''), sep = "\t", header = T, check.names = F)
  #过滤在至少在75%的样本中都有表达的基因,rowSums(rawcount>0) 是对矩阵取行的和
  data2 <- rowSums(data > 0) >= floor(0.75 * ncol(data))
  table(data2)
  filter_count <- data[data2,]
  filter_count[1:4, 1:4]
  dim(filter_count)
  write.table(filter_count, paste(output_dir, "/filter_count_eRNA.txt", sep = ''), quote=FALSE,sep = "\t",row.names = FALSE)


  ###清除数据，重新导入filter_count命名为data，后面的代码就不用改了
  data <- read.delim(paste(output_dir, "/filter_count_eRNA.txt", sep = ''), row.names = 1)
  #控制组的数量，这里是正常样本的数量11个，因为前面已经把正常组提到了最前面几列（一定要根据自己的数据集情况操作），所以直接提就行了
  afcon = 11
  conData = data[, as.vector(colnames(data)[1:afcon])]
  aftreat = afcon + 1
  treatData = data[, as.vector(colnames(data)[aftreat:ncol(data)])]
  rt = cbind(conData, treatData)
  conNum = ncol(conData)
  treatNum = ncol(treatData)
  #limma差异标准流程
  Type = c(rep("con", conNum), rep("treat", treatNum))
  design <- model.matrix(~0 + factor(Type))
  colnames(design) <- c("con", "treat")
  fit <- lmFit(rt, design)
  cont.matrix <- makeContrasts(treat - con, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  Diff = topTable(fit2, adjust = 'fdr', number = length(rownames(data)))
  #保存所有基因的差异结果
  DIFFOUT = rbind(id = colnames(Diff), Diff)
  write.table(DIFFOUT,paste(output_dir, "/DIFF_all.txt", sep = ''), sep = "\t", quote = F, col.names = F)
  ####热图展示差异最大的前50个基因
  Diff = Diff[order(as.numeric(as.vector(Diff$logFC))),]
  diffGene = as.vector(rownames(Diff))
  diffLength = length(diffGene)
  afGene = c()
  if (diffLength > (100)) {
    afGene = diffGene[c(1:50, (diffLength - 50 + 1):diffLength)]
  }else {
    afGene = diffGene
  }
  afExp = rt[afGene,]
  #分组标签
  Type = c(rep("N", conNum), rep("T", treatNum))
  names(Type) = colnames(rt)
  Type = as.data.frame(Type)
  #分组标签的注释颜色
  anncolor = list(Type = c(T = pal_npg()(1), N = pal_npg()(2)[2]))
  #热图
  pdf(file = paste(output_dir, "/DIFF_heatmap.pdf", sep = ''), height = 7, width = 8)
  pheatmap(afExp,
           #热图数据
           annotation = Type, #分组
           color = colorRampPalette(c(pal_npg()(2)[2], "white", pal_npg()(1)))(50), #热图颜色
           cluster_cols = F, #不添加列聚类树
           show_colnames = F,
           scale = "row",
           fontsize = 8,
           fontsize_row = 6,
           fontsize_col = 8,
           annotation_colors = anncolor #展示列名
  )
  dev.off()
  ####火山图差异标准设置
  adjP = 0.05
  aflogFC = 1
  Significant = ifelse((Diff$P.Value < adjP & abs(Diff$logFC) > aflogFC), ifelse(Diff$logFC > aflogFC, "Up", "Down"), "Not")
  #开始绘制
  p = ggplot(Diff, aes(logFC, -log10(P.Value))) +
    geom_point(aes(col = Significant), size = 3) +
    scale_color_manual(values = c(pal_npg()(2)[2], "#838B8B", pal_npg()(1))) +
    labs(title = " ") +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) +
    geom_hline(aes(yintercept = -log10(adjP)), colour = "gray", linetype = "twodash", linewidth = 1) +
    geom_vline(aes(xintercept = aflogFC), colour = "gray", linetype = "twodash", linewidth = 1) +
    geom_vline(aes(xintercept = -aflogFC), colour = "gray", linetype = "twodash", linewidth = 1)
  #添加标记，按照
  point.Pvalue = 0.0001
  point.logFc = 3.8
  #继续绘制
  Diff$symbol = rownames(Diff)
  p = p + theme_bw()
  for_label <- Diff %>%
    filter(abs(logFC) > point.logFc & P.Value < point.Pvalue)
  p +
    geom_point(size = 2, shape = 1, data = for_label) +
    ggrepel::geom_label_repel(
      aes(label = symbol),
      data = for_label,
      color = "black",
      label.size = 0.1
    )
  ggsave(paste(output_dir, "/DIFF_vol.pdf", sep = ''), p, width = 5.5, height = 5)
}