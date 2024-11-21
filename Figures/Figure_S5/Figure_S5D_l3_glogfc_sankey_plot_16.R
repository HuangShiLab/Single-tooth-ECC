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
  df<-read.table(datafile, header=T, row.names=1, sep="\t", quote="", comment.char = "", check.names = F)
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
                        row.names = NULL, check.names = F,
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
print(identical(colnames(df), fmetadata[, 1]))

data_HHCC_gfc <- read.table("../../Results/Figure_6/MiC/function_l3_Position_crossRF_out/all_Taxon_OTU.wilcox_res.tsv", sep = "\t", header = T)

data_HHCC_gfc <- subset(data_HHCC_gfc, rf_imps_rank <= 16)
most_imp_f <- unique(data_HHCC_gfc$feature)

data_HHCC_mean_gfc <- data_HHCC_gfc %>%
  group_by(feature) %>%
  summarise(mean_gfc = mean(generalized_logfc))
matched_index <- match(data_HHCC_gfc$feature, data_HHCC_mean_gfc$feature)
data_HHCC_gfc$mean_gfc <- unlist(data_HHCC_mean_gfc[matched_index, "mean_gfc"])
data_HHCC_gfc <- data_HHCC_gfc %>%
  mutate(status = case_when(
    mean_gfc > 0 ~ "Caries_enriched",     # 如果 mean_gfc > 0，设置为 Caries
    mean_gfc == 0 ~ "Neutral",   # 如果 mean_gfc == 0，设置为 Neutral
    mean_gfc < 0 ~ "ConfidentH_enriched"      # 如果 mean_gfc < 0，设置为 ConfidentH
  ))
data_HHCC_gfc$status <- with(data_HHCC_gfc, paste(dataset, status, sep = "_"))

meta_confidentH <- subset(metadata, Future_Status_Tooth == "ConfidentH")
df_confidentH <- subset(df, rownames(df) %in% rownames(meta_confidentH))
identical(rownames(meta_confidentH), rownames(df_confidentH))
df_confidentH <- df_confidentH[, which(colSums(df_confidentH) != 0)]
pos_imp_r <- rf_clf.comps(df_confidentH, factor(meta_confidentH$NicheIM), "I")
pos_gfc <- pos_imp_r$feature_imps_list[[1]]
write.table(pos_gfc, "l3_pos_feature_imps_result_16.txt", quote = F, sep = "\t", row.names = F)
pos_gfc <- subset(pos_gfc, feature %in% most_imp_f)
pos_gfc <- pos_gfc %>%
  mutate(pos = case_when(
    generalized_logfc > 0 ~ "ConfidentH_Incisor_enriched",     # 如果 mean_gfc > 0，设置为 Caries
    generalized_logfc == 0 ~ "Neutral",   # 如果 mean_gfc == 0，设置为 Neutral
    generalized_logfc < 0 ~ "ConfidentH_Molar_enriched"      # 如果 mean_gfc < 0，设置为 ConfidentH
  ))

matched_index <- match(data_HHCC_gfc$feature, pos_gfc$feature)
data_HHCC_gfc$pos <- pos_gfc[matched_index, "pos"]

data <- data_HHCC_gfc[, c("status", "pos")]

# 为每行添加一个计数列
data$count <- 1

sum_count <- data %>%
  group_by(status, pos) %>%
  summarise(sum_count = sum(count), .groups = 'drop')  # 使用 .groups='drop' 来避免分组信息

# unique_labels <- data %>%
#   distinct(status, pos, .keep_all = TRUE)  # 保留每个组的第一次出现

# 绘制三层的 Sankey 图
g <- ggplot(data = data,
       aes(axis1 = status, axis2 = pos, y = count)) +
  geom_alluvium(aes(fill = status)) +   # 绘制流动线，并根据 status 进行颜色填充
  geom_stratum() +                      # 绘制层（节点）
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  # 为节点添加标签
  scale_x_discrete(limits = c("Status", "Position"), expand = c(0.15, 0.15)) +  # 设置X轴显示的层次
  theme_minimal() +
  labs(x = "Category", y = "Count") +
  theme(legend.position = "none")

g

ggsave(paste0(outpath, "1867_rarefied_sankey_plot_16.pdf"), g, height = 12, width = 11, limitsize = F)

summarise <- table(data[, c("status", "pos")])
write.table(summarise, paste0(outpath, "/1867_rarefied_sankey_plot_GTDB_summarise_16.txt"), quote = F, sep = "\t", row.names = T, col.names = NA)

