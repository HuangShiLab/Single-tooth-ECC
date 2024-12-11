## install and load necessary libraries for data analyses
#-------------------------------
p <- c("reshape2","ggplot2","pheatmap","combinat","randomForest", "gridExtra", "crossRanger", "RColorBrewer",
       "vegan", "biomformat", "doMC", "cowplot", "pROC", "ggsignif", "ggpubr", "dplyr", "viridis", "patchwork")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
#-------------------------------

imp_score <- read.table("../../Results/Figure_6/MiC/taxonomy_Position_crossRF_out/all_Taxon_OTU.RF_imps.tsv", 
                        sep = "\t", header = T)
imp_score <- imp_score[order(imp_score$imp_rank_mean), ]
top_32_features <- imp_score[1:32, ]

Taxonomy <- read.table("../../data/1867/taxonomy/silva_taxonomy_collapse/level8/taxonomy.tsv", 
                       sep = "\t", header = T, row.names = 1)

getLastElement <- function(s) {
  elements <- strsplit(s, ";")[[1]]
  return(tail(elements, 1))
}

Taxonomy$Taxon <- with(Taxonomy, paste(Taxon, ASV_ID, sep = "_"))
Taxonomy$Taxon <- sapply(Taxonomy$Taxon, getLastElement)
Taxonomy$Taxon <- trimws(Taxonomy$Taxon)

imp_rank <- read.table("../../Results/Figure_6/MiC/taxonomy_Position_crossRF_out/all_Taxon_OTU.RF_imps.tsv",
                       sep = "\t", header = T)#, row.names = 1, check.names = F)
imp_rank <- subset(imp_rank, n_pos > 7 & imp_rank_min < 32)
imp_rank <- imp_rank[order(-imp_rank$imp_rank_mean), , drop = F]
imp_feature <- imp_rank$feature
matched_index <- match(imp_rank$feature, rownames(Taxonomy))
imp_rank$Taxon <- Taxonomy[matched_index, "Taxon"]
imp_tax <- imp_rank$Taxon

imp_s <- read.table("../../Results/Figure_6/MiC/taxonomy_Position_crossRF_out/all_Taxon_OTU.wilcox_res.tsv",
                    sep = "\t", header = T)
imp_s$index <- with(imp_s, paste(feature, dataset, sep = "_"))

glogfc <- read.table("../../Results/Figure_8/generalized_log2_fold_change/generalized_logfc_by_datasets.txt",
                     sep = "\t", header = T, row.names = 1)
glogfc <- subset(glogfc, rownames(glogfc) %in% imp_feature)
glogfc <- as.matrix(glogfc)
glogfc <- melt(glogfc)

matched_index <- match(glogfc$Var1, rownames(Taxonomy))
glogfc$Taxon <- Taxonomy[matched_index, "Taxon"]
glogfc$index <- with(glogfc, paste(Var1, Var2, sep = "_"))
matched_index <- match(glogfc$index, imp_s$index)
glogfc$imp_s <- imp_s[matched_index, "rf_imps_rank"]
glogfc$Taxon <- factor(glogfc$Taxon, levels = imp_tax, ordered = T)
glogfc$Var2 <- factor(glogfc$Var2,
                           levels = c("T55", "T54", "T5161", "T64", "T65", "T75", "T74", "T84", "T85"),
                           ordered = T)

HHRH_glogfc <- read.table("../../Results/Figure_8/generalized_log2_fold_change/HHRH_generalized_logfc_by_datasets.txt",
                          sep = "\t", header = T, row.names = 1)
HHRH_glogfc <- subset(HHRH_glogfc, rownames(HHRH_glogfc) %in% imp_feature)
HHRH_glogfc <- as.matrix(HHRH_glogfc)
HHRH_glogfc <- melt(HHRH_glogfc)
matched_index <- match(HHRH_glogfc$Var1, rownames(Taxonomy))
HHRH_glogfc$Taxon <- Taxonomy[matched_index, "Taxon"]
HHRH_glogfc$index <- with(HHRH_glogfc, paste(Var1, Var2, sep = "_"))
matched_index <- match(HHRH_glogfc$index, imp_s$index)
HHRH_glogfc$imp_s <- imp_s[matched_index, "rf_imps_rank"]
HHRH_glogfc$Taxon <- factor(HHRH_glogfc$Taxon, levels = imp_tax, ordered = T)
HHRH_glogfc$Var2 <- factor(HHRH_glogfc$Var2,
                       levels = c("T55", "T54", "T5161", "T64", "T65", "T75", "T74", "T84", "T85"),
                       ordered = T)

glogfc_plot <- ggplot(glogfc, aes(Var2, Taxon, fill = value)) +
  geom_tile() +
  geom_text(aes(label = imp_s), color = "white", size = 3) +
  scale_fill_viridis()+
  labs(fill = "Caries\nVS\nConfidentH\ngeneralized\nlog2(fold change)") +
  xlab("Tooth position") + 
  ylab("ASV") + 
  theme_bw() + 
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(0.9, "cm"))
ggsave(filename="../../Results/Figure_8/generalized_log2_fold_change/HHCC_top32_glogfc.pdf", 
       plot=glogfc_plot, width=7, height=16, limitsize = F)


HHRH_glogfc_plot <- ggplot(HHRH_glogfc, aes(Var2, Taxon, fill = value)) +
  geom_tile() +
  geom_text(aes(label = imp_s), color = "white", size = 3) +
  scale_fill_viridis()+
  labs(fill = "RelativeH\nVS\nConfidentH\ngeneralized\nlog2(fold change)") +
  xlab("Tooth position") + 
  ylab("ASV") + 
  theme_bw() + 
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(0.9, "cm"))
ggsave(filename="../../Results/Figure_8/generalized_log2_fold_change/HHCC_top32_HHRH_glogfc.pdf", 
       plot=HHRH_glogfc_plot, width=7, height=16, limitsize = F)


HHRH_glogfc_plot <- ggplot(HHRH_glogfc, aes(Var2, Taxon, fill = value)) +
  geom_tile() +
  # geom_text(aes(label = imp_s), color = "white", size = 3) +
  scale_fill_viridis()+
  #scale_fill_gradient2(low = "#330000", high = "#FF3300", mid = "gray") +
  #scale_fill_gradient(low="#EFEFFF", high="blue") +  
  labs(fill = "RelativeH\nVS\nConfidentH\ngeneralized\nlog2(fold change)") +
  xlab("Tooth position") + 
  ylab("ASV") + 
  theme_bw() + 
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(0.9, "cm"))

combined_plot <- glogfc_plot + HHRH_glogfc_plot + plot_layout(ncol = 2)
ggsave(filename="../../Results/Figure_8/generalized_log2_fold_change/generalized_log2fold_change.pdf", 
       plot=combined_plot, width=11, height=16, limitsize = F)

