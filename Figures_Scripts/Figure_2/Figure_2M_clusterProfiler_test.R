if(! require("BiocManager", quietly = T)) {
  install.package("BiocManager")
}
BiocManager::install("clusterProfiler")

library(clusterProfiler)
gene_list <- read.table("sig_diff_KO.txt", sep = "\t", header = F)
gene_list <- gsub("^K0+", "", gene_list$V1)
ek <- enrichKEGG(gene_list)

ek <- mutate(ek, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

library(ggplot2)
library(forcats)
ggplot(ek, showCategory = 10,
       aes(richFactor, fct_reorder(Description, richFactor))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        trans = "log10",
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  scale_size_continuous(range=c(2, 10)) +
  labs(x = "Rich Factor", y = "Enriched pathway") # 设置标签

