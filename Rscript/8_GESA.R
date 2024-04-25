GESA <- function(expressdir, output_dir) {
  options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  ##差异分析
  if (!require(org.Hs.eg.db)) {
    BiocManager::install("org.Hs.eg.db")
  }
  if (!require(clusterProfiler)) {
    BiocManager::install("clusterProfiler")
  }
  if (!require(enrichplot)) {
    BiocManager::install("enrichplot")
  }
  if (!require(patchwork)) {
    BiocManager::install("patchwork")
  }
  library(ggplot2)
  library(limma)
  library(pheatmap)
  library(ggsci)
  library(dplyr)
  #读取差异表达样本，正常样本放前面，癌症样本放后面。
  data = read.table(expressdir, sep = "\t", header = T, check.names = F)
  #过滤在至少在75%的样本中都有表达的基因,rowSums(rawcount>0) 是对矩阵取行的和
  data2 <- rowSums(data > 0) >= floor(0.75 * ncol(data))
  table(data2)
  filter_count <- data[data2,]
  filter_count[1:4, 1:4]
  dim(filter_count)
  write.table(filter_count, paste(output_dir, "/filter_count_gene_exp.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = '\t')

  ###清除数据，重新导入filter_count命名为data，后面的代码就不用改了
  data <- read.delim(paste(output_dir, "/filter_count_gene_exp.txt", sep = ''), row.names = 1)
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
  write.table(DIFFOUT, paste(output_dir, "/diff_gene.txt", sep = ''), sep = "\t", quote = F, col.names = F)
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
  pdf(file = paste(output_dir, "/diff_heatmap.pdf", sep = ''), height = 7, width = 8)
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
  ggsave(paste(output_dir, "/diff_vol.pdf", sep = ''), p, width = 5.5, height = 5)


  ####################################1.载入测试数据
  lapply(c('clusterProfiler', 'enrichplot', 'patchwork'),
         function(x) { library(x, character.only = T) })
  # Please go to https://yulab-smu.github.io/clusterProfiler-book/ for the full vignette.
  data(geneList, package = "DOSE")
  class(geneList)
  kk2 <- gseKEGG(geneList = geneList,
                 organism = 'hsa',
                 nPerm = 10000,
                 minGSSize = 10,
                 maxGSSize = 200,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "none")
  lapply(c('org.Hs.eg.db', 'stringr', 'dplyr'),
         function(x) { library(x, character.only = T) })
  deg = read.table(paste(output_dir, "/diff_gene.txt", sep = ''), header = T, sep = "\t", check.names = F, row.names = 1)
  logFC_t = 1.5
  deg$g = ifelse(deg$P.Value > 0.05, 'stable',
                 ifelse(deg$logFC > logFC_t, 'UP',
                        ifelse(deg$logFC < -logFC_t, 'DOWN', 'stable'))
  )
  deg$symbol = rownames(deg)
  df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
             toType = c("ENTREZID"),
             OrgDb = org.Hs.eg.db)
  DEG = deg
  DEG = merge(DEG, df, by.y = 'SYMBOL', by.x = 'symbol')
  data_all_sort <- DEG %>%
    arrange(desc(logFC))
  geneList = data_all_sort$logFC #把foldchange按照从大到小提取出来
  names(geneList) <- data_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID

  class(kk2)
  colnames(kk2@result)
  kegg_result <- as.data.frame(kk2)
  rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
  gseaplot2(kk2, geneSetID = rownames(kk2@result)[head(order(kk2@result$enrichmentScore))])
  gseaplot2(kk2, geneSetID = rownames(kk2@result)[tail(order(kk2@result$enrichmentScore))])
  p <- gseaplot2(kk2,
            title = "Focal adhesion", #设置title
            "hsa04510", #绘制hsa04510通路的结果
            color = "red", #线条颜色
            base_size = 20, #基础字体的大小
            subplots = 1:2, #展示上2部分
            pvalue_table = T) # 显示p值
  ggsave(paste(output_dir, "/gseaplot.pdf", sep = ''), p, width = 8, height = 5)
  #####山脊图可视化
  p1 <- ridgeplot(kk2, 10)
  ggsave(paste(output_dir, "/ridgeplot.pdf", sep = ''), p1, width = 8, height = 5)
  #####气泡图可视化
  p2 <- dotplot(kk2)
  ggsave(paste(output_dir, "/dotplot.pdf", sep = ''), p2, width = 8, height = 5)
}