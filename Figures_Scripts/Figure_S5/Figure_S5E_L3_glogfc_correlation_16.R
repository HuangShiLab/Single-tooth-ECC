## install and load necessary libraries for data analyses
#-------------------------------
p <- c("reshape2","ggplot2","pheatmap","combinat","randomForest", "gridExtra", "crossRanger", "RColorBrewer",
       "vegan", "biomformat", "doMC", "cowplot", "pROC", "ggsignif", "ggpubr", "colorspace", "stringr", "viridis",
       "patchwork", "dplyr", "ggsankey", "ggalluvial")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
source("../data_trimming_util.R")
#-------------------------------
## input args
#-------------------------------
rf.opts<-list(outdir=NA, ntree=5000, verbose=FALSE, nfolds=10)
#-------------------------------
datafile <- "../../data/1867/function/KO/1867_functions_category.l3.Abd.txt"
sample_metadata <- "../../data/1867/function//1867_function_metadata.tsv"
feature_metadata <- "../../data/1867/function//ML_features_table/1867_function_feature_taxon_ko.txt"
prefix_name<-"L3_pathway"
s_category<-"Position2"
c_category<-"Future_Status_Tooth"
AddTaxonomy=TRUE
p_cutoff=0.05
q_cutoff=0.05
p.adj.method="BH"
h=10
outpath <- paste("../../Results/Figure_8/MiC/", prefix_name,"_Position_crossRF_out/", sep="")
dir.create(outpath, recursive = T)

clinical_features <- c("sum_s_dt", "sum_ns_dt", "sum_dmfs", "spatial_dist_weighted_mean_dmfs")
mic_dist_features <- c("spatial_dist_weighted_mean_md", "mean_md", "mean_md_H_T51", "mean_md_H_T52",
                       "mean_md_H_T53", "mean_md_H_T54", "mean_md_H_T55", "mean_md_H_T61", "mean_md_H_T62",
                       "mean_md_H_T63", "mean_md_H_T64", "mean_md_H_T65", "mean_md_H_T71", "mean_md_H_T72",
                       "mean_md_H_T73", "mean_md_H_T74", "mean_md_H_T75", "mean_md_H_T81", "mean_md_H_T82",
                       "mean_md_H_T83", "mean_md_H_T84", "mean_md_H_T85", "mean_md_H_T5161")
if(grepl("biom$", datafile)){
  biom <- read_biom(datafile)
  df <- data.frame(as.matrix(biom_data(biom)), check.names = FALSE)
  df <- df[, which((!colnames(df) %in% clinical_features) & (!colnames(df) %in% mic_dist_features))]
}else{
  df<-read.table(datafile, header=T, row.names=1, sep="\t", quote="", comment.char = "")
  df <- df[, which((!colnames(df) %in% clinical_features) & (!colnames(df) %in% mic_dist_features))]
}
df<-df[order(rownames(df)), ]
#df<-sweep(df, 1, rowSums(df), "/")
#-------------------------------
# Feature metadata input
#-------------------------------
if(!is.na(feature_metadata)){
  #fmetadata<-read.table(feature_metadata,header=T, sep="\t", fill = TRUE, comment.char = "")
  fmetadata<-read.table(feature_metadata,header=T, sep='\t', quote = "",
                        row.names = NULL,
                        stringsAsFactors = FALSE, comment.char = "")
  fmetadata <- subset(fmetadata, fmetadata[, 1] %in% colnames(df))
}

add_ann<-function(tab, fmetadata, tab_id_col=1, fmetadata_id_col=1){
  fmetadata_matched<-fmetadata[which(fmetadata[, fmetadata_id_col] %in% tab[, tab_id_col]),]
  out<-merge(tab, fmetadata_matched, by.x=tab_id_col, by.y=fmetadata_id_col)
  out
}
#-------------------------------
# Sample Metadata input
#-------------------------------
allmetadata<-read.table(sample_metadata,header=T,sep="\t",row.names=1, quote="", comment.char="")
if(length(allmetadata)==1){metadata<-data.frame(allmetadata[order(rownames(allmetadata)),])
all_group<-colnames(metadata)<-colnames(allmetadata)
}else{
  metadata<-allmetadata[order(rownames(allmetadata)),]
  all_group<-colnames(metadata)
  all_group_f<-colnames(metadata)[sapply(metadata,class)=="factor"]
  all_group_n<-colnames(metadata)[sapply(metadata,class)!="factor"]
}
# heatmap of contingency tables
cols <- colorRampPalette(brewer.pal(3, "Blues"))
tab_ori<-table(metadata[, c(s_category, c_category)])
pheatmap(tab_ori, cluster_rows = F, cluster_cols = F, display_numbers = TRUE,  number_format = "%.0f", color = cols(100),
         filename = paste(outpath, "contingency_table", ".",c_category, "_by_", s_category,".pdf",sep=""), width=3, height=4)

print(identical(rownames(df), rownames(metadata)))

imp_r <- read.table("../../Results/Figure_6/MiC/function_l3_Position_crossRF_out/all_Taxon_OTU.wilcox_res.tsv",
                    sep = "\t", header = T)

data <- imp_r[, c("feature", "dataset", "generalized_logfc")]
df <- dcast(data, feature ~ dataset)
rownames(df) <- df$feature
df <- df[, -1]
cor_dist <- 1 - cor(df)

pdf(paste0(outpath, "/1867_rarefied_all_l3_function_gfc_cor_dist_heatmap.pdf"),width=4.5,height=4)
pheatmap(as.matrix(cor_dist),
         clustering_method = "average",
         border_color = NA,
         treeheight_row = 50,
         treeheight_col = 50,
         angle_col = 45)
dev.off()


imp_r <- read.table("../../Results/Figure_6/MiC/function_l3_Position_crossRF_out//all_Taxon_OTU.RF_imps.tsv",
                    sep = "\t", header = T)
most_imp_f <- subset(imp_r, imp_rank_min <= 16)$feature
df <- subset(df, rownames(df) %in% most_imp_f)
cor_dist <- 1 - cor(df)

pdf(paste0(outpath, "/1867_16_l3_function_gfc_cor_dist_heatmap.pdf"),width=4.5,height=4)
pheatmap(as.matrix(cor_dist),
         clustering_method = "average",
         border_color = NA,
         treeheight_row = 50,
         treeheight_col = 50,
         angle_col = 45)
dev.off()


data_HHCC_gfc <- read.table("l3_data_HHCC_res.txt", sep = "\t", header = T)
data <- data_HHCC_gfc[, c("feature", "dataset", "generalized_logfc")]
df <- dcast(data, feature ~ dataset)
rownames(df) <- df$feature
HHCC_gfc_df <- df[, -1]

data_HHRH_gfc <- read.table("L3_data_HHRH_res.txt", sep = "\t", header = T)
data <- data_HHRH_gfc[, c("feature", "dataset", "generalized_logfc")]
df <- dcast(data, feature ~ dataset)
rownames(df) <- df$feature
HHRH_gfc_df <- df[, -1]

identical(rownames(HHCC_gfc_df), rownames(HHRH_gfc_df))
identical(colnames(HHCC_gfc_df), colnames(HHRH_gfc_df))

HHCC_HHRH_cor <- matrix(NA, nrow = ncol(HHCC_gfc_df), ncol = 1)
rownames(HHCC_HHRH_cor) <- colnames(HHCC_gfc_df)
colnames(HHCC_HHRH_cor) <- "cor_dist"
for(i in 1:ncol(HHCC_gfc_df)) {
  HHCC_HHRH_cor[i, 1] <- cor(HHCC_gfc_df[, i], HHRH_gfc_df[, i])
}

data <- data.frame(feature = data_HHCC_gfc$feature,
                   dataset = data_HHCC_gfc$dataset,
                   HHCC_gfc = data_HHCC_gfc$generalized_logfc,
                   HHRH_gfc = data_HHRH_gfc$generalized_logfc)

gg <- ggplot(data, aes(x = HHCC_gfc, y = HHRH_gfc)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # 在 x = 0 处添加虚线
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +  # 在 y = 0 处添加虚线
  #scale_color_brewer(palette = "Set3") +
  facet_wrap(~dataset) +
  geom_smooth(method = lm, linetype = 1, se = F, span = 1, color = "darkgray") + # 趋势线
  stat_cor(method = "pearson", color = "black") +
  xlab("gFC (Caries V.S. ConfidentH)") + # 设置X轴名称
  ylab("gFC (RelativeH V.S. ConfidentH)") + # 设置Y轴名称
  #scale_color_discrete(name = "Group") + # 设置图例名称
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave(paste0(outpath, "/HHCC_VS_HHRH_gfc_pos_scattor_plot.pdf"), gg, height = 8, width = 8)

data <- subset(data, feature %in% most_imp_f)

gg <- ggplot(data, aes(x = HHCC_gfc, y = HHRH_gfc)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # 在 x = 0 处添加虚线
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +  # 在 y = 0 处添加虚线
  #scale_color_brewer(palette = "Set3") +
  facet_wrap(~feature) +
  geom_smooth(method = lm, linetype = 1, se = F, span = 1, color = "darkgray") + # 趋势线
  stat_cor(method = "pearson", color = "black") +
  xlab("gFC (Caries V.S. ConfidentH)") + # 设置X轴名称
  ylab("gFC (RelativeH V.S. ConfidentH)") + # 设置Y轴名称
  #scale_color_discrete(name = "Group") + # 设置图例名称
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave(paste0(outpath, "/HHCC_VS_HHRH_gfc_16_l3_function_scattor_plot.pdf"), gg, height = 40, width = 40)

