#library(DistVis)
library(ggplot2)
library(stringr)
library(dplyr)

DistBoxplot <- function(dist, meta, dist_name = "Distance", outdir = ".", result_name){
  if(! dir.exists(outdir)) {
    dir.create(outdir)
  }
  print(identical(rownames(dist), colnames(dist)))
  print(identical(rownames(dist), rownames(meta)))
  
  groups <- c("ConfidentH_VS_ConfidentH", "RelativeH_VS_RelativeH", "Caries_VS_Caries", 
              "RelativeH_VS_Caries", "RelativeH_VS_ConfidentH", "ConfidentH_VS_Caries")
  dist$x <- rownames(dist)
  dist <- melt(dist)
  dist$group <- paste(meta[dist$x, "Future_Status_Tooth"], "VS", meta[dist$variable, "Future_Status_Tooth"], sep = "_")
  dist <- subset(dist, group %in% groups)
  
  # p<-ggplot(dist, aes(x = value, y = group)) + 
  #   geom_boxplot(outlier.shape = NA) + 
  #   geom_jitter(aes(color=group), position=position_jitterdodge(jitter.width= 0.2,dodge.width = 0.8),size=1,alpha=0.4) +
  #   scale_colour_manual(values = brewer.pal(8, "Pastel2")) +
  #   geom_signif(comparisons = list(c("ConfidentH_VS_Caries", "RelativeH_VS_Caries")), 
  #               map_signif_level = function(p) {if(p < 0.01) {p = "**"} else if(p < 0.05) {p = "*"} else {p = "NS"}},
  #               test = "wilcox.test", textsize = 4, step_increase = 0.1,
  #               test.args = list(exact = FALSE, correct = FALSE, conf.int = TRUE, 
  #                                conf.level = 0.95)) + # 添加wilcoxon test结果，并使不同分组的检验间隔0.01
  #   xlab(dist_name) + 
  #   ylab("Group") + 
  #   theme_minimal() +
  #   theme_bw() + 
  #   theme(axis.line = element_line(color = "black"),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.background = element_blank(),
  #         panel.border = element_blank(),
  #         strip.text = element_text(size = 14),
  #         strip.background = element_rect(colour = "white"),
  #         legend.title = element_text(size = 14), 
  #         legend.text = element_text(size = 14),
  #         axis.title.x = element_text(size = 14),
  #         axis.title.y = element_text(size = 14),
  #         axis.text.x = element_blank(),
  #         axis.text.y = element_text(size = 14),
  #         legend.position = "none")
  
  p <- ggplot(dist, aes(x = value, y = group, color = group)) +
    geom_boxplot(outlier.shape = NA) +
    scale_colour_manual(values = brewer.pal(8, "Pastel2")) +
    #geom_jitter(aes(color = group), shape = 1, width = 0.1, alpha = 0.5) +
    geom_signif(comparisons = list(c("ConfidentH_VS_Caries", "RelativeH_VS_Caries")), 
                map_signif_level = function(p) {if(p < 0.01) {p = "**"} else if(p < 0.05) {p = "*"} else {p = "NS"}},
                test = "wilcox.test", textsize = 4, step_increase = 0.1,
                test.args = list(exact = FALSE, correct = FALSE, conf.int = TRUE, 
                                 conf.level = 0.95)) + # 添加wilcoxon test结果，并使不同分组的检验间隔0.01
    xlab(dist_name) + 
    ylab("Group") + 
    theme_minimal() +
    theme_bw() + 
    theme(axis.text = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text = element_text(size = 14),
          strip.background = element_rect(colour = "white"),
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 14),
          plot.title = element_text(size = 14), 
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.position = "none")
  
  ggsave(filename=paste(outdir, result_name, sep = "/"),plot=p, width=6, height=4)
  
  # pdf(paste(outdir, result_name, sep = "/"), width = 8, height = 4)
  # par(mar = c(5, 15, 4, 2) + 0.1)  # 默认为c(5, 4, 4, 2) + 0.1，增加底部的空间
  # boxplot(value ~ group, data = dist, ylim = c(0, max(dist$value) + 1), horizontal = TRUE, xlab = "", ylab = "", las = 1)
  # # 添加x轴和标题，设置轴标题距离为2
  # title(xlab=dist_name, mgp=c(2, 1, 0))
  # # 添加y轴和标题，设置轴标题距离为8
  # title(ylab="Group", mgp=c(12, 1, 0))
  
  # kruskal_test_result <-kruskal.test(data = dist, value ~ group)
  # print(kruskal_test_result)
  # text(max(dist$value)+0.5, length(unique(dist$group)) * 0.9, paste0("Kruskal-Wallis\np-value = ", round(kruskal_test_result$p.value, 5)), col = "red", srt = 270, adj = 0)
  # 
  # # 对两两分组执行 Wilcoxon 秩和检验，并标注 p 值小于 0.05 的比较
  # groups <- unique(dist$group)
  # combinations <- combn(1:length(groups), 2)  # 生成所有两两组合
  # for (i in 1:ncol(combinations)) {
  #   group_a <- subset(dist, group == groups[combinations[1, i]])$value
  #   group_b <- subset(dist, group == groups[combinations[2, i]])$value
  #   
  #   test <- wilcox.test(group_a, group_b)
  #   
  #   # 如果 p 值小于 0.05，则在相应的位置上添加一个星号
  #   if (test$p.value < 0.05) {
  #     x_coord <- max(c(max(group_a), max(group_b))) + runif(1)
  #     y_coords <- combinations[, i]
  #     segments(x_coord, y_coords[1], x_coord, y_coords[2], col = "gray") # 添加一条线
  #     text(x_coord + 0.01, mean(y_coords), "*", col = "red")
  #   }
  # }
  # dev.off()
}

outpath <- "../../Results/Figure_S4/"

dist <- read.table("../../data/637/taxonomy/637_all_asv_phylo_rpca/phylo_rpca.distance_matrix/distance-matrix.tsv", sep = "\t", header = T, row.names = 1)
colnames(dist) <- str_sub(colnames(dist), 2)
identical(rownames(dist), colnames(dist))
meta <- read.table("../../data/637/taxonomy/637_taxonomy_metadata.tsv", sep = "\t", header = T, row.names = 1)
identical(rownames(dist), rownames(meta))
DistBoxplot(dist, meta, dist_name = "phylo-rPCA distance", outdir = outpath, result_name ="637_taxonomy_phylo-rPCA dist.boxplot.pdf")

dist <- read.table("../../data/637/function/637_all_ko_rpca/rpca.distance_matrix/distance-matrix.tsv")
identical(rownames(dist), colnames(dist))
meta <- read.table("../../data/637/function/637_function_metadata.tsv", sep = "\t", header = T, row.names = 1)
identical(rownames(dist), rownames(meta))
DistBoxplot(dist, meta, dist_name = "rPCA distance", outdir = outpath,  result_name = "637_function_rPCA dist.boxplot.pdf")

dist <- read.table("../../data/1867/taxonomy/1867_all_asv_phylo_rpca/phylo_rpca.distance_matrix/distance-matrix.tsv", sep = "\t", header = T, row.names = 1)
colnames(dist) <- str_sub(colnames(dist), 2)
identical(rownames(dist), colnames(dist))
meta <- read.table("../../data/1867/taxonomy/1867_taxonomy_metadata.tsv", sep = "\t", header = T, row.names = 1)
identical(rownames(dist), rownames(meta))
DistBoxplot(dist, meta, dist_name = "phylo-rPCA distance", outdir = outpath, result_name = "1867_taxonomy_phylo-rPCA dist.boxplot.pdf")

dist <- read.table("../../data/1867/function/1867_all_ko_rpca/rpca.distance_matrix/distance-matrix.tsv")
identical(rownames(dist), colnames(dist))
meta <- read.table("../../data/1867/function/1867_function_metadata.tsv", sep = "\t", header = T, row.names = 1)
identical(rownames(dist), rownames(meta))
DistBoxplot(dist, meta, dist_name = "rPCA distance", outdir = outpath, result_name = "1867_function_rPCA dist.boxplot.pdf")

