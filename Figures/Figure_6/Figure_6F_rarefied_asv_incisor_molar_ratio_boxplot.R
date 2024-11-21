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
datafile <- "../../data/1867/taxonomy/silva_taxonomy_collapse/level8/rarefied-feature-table.biom"
sample_metadata <- "../../data/1867/taxonomy/1867_taxonomy_metadata.tsv"
feature_metadata <- "../../data/1867/taxonomy/GTDB_taxonomy_collaspse/taxonomy.tsv"
prefix_name<-"taxonomy"
s_category<-"Position2"
c_category<-"Future_Status_Tooth"
AddTaxonomy=TRUE
p_cutoff=0.05
q_cutoff=0.05
p.adj.method="BH"
h=10
outpath <- paste("../../Results/Figure_8/MiC/rarefied_",prefix_name,"_Position_crossRF_out/", sep="")
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
print(identical(colnames(df), fmetadata[, 1]))

log.mat<-function(mat, base=2){
  if(any(mat == 0)) {
    if(sum(mat == 0) == length(unlist(mat))) {
      mat <- mat + 1e-5
    }
    else{
      v <- as.vector(mat)
      minval <- min(v[v > 0])/2
      mat <- mat + minval
    }
  }
  out<-log(mat, base)
  return(out)
}

df <- as.matrix(df)
data <- melt(df)
colnames(data) <- c("sample", "feature", "abundance")

pos_imp_r <- read.table("~/Projects/Teng/github/ECC/Results/Figure_8/MiC/1867_H2H_rarefied_clr_differential_result.txt",
                        sep = "\t", header = T)

matched_index <- match(data$feature, pos_imp_r$X)
data$incisor_vs_molar_gfc <- pos_imp_r[matched_index, "generalized_logfc"]
data <- data %>%
  mutate(feature_group = case_when(
    incisor_vs_molar_gfc > 0 ~ "Incisor_enriched",     
    incisor_vs_molar_gfc < 0 ~ "Molar_enriched",   
    incisor_vs_molar_gfc == 0 ~ "Neutral"    
  ))

  
data_all <- data %>% 
  group_by(sample, feature_group) %>%
  summarise(sum_abd = sum(abundance))

data_all <- subset(data_all, feature_group != "Neutral")
df_all <- dcast(data_all, sample ~ feature_group)
df_all$ratio <- with(df_all, Incisor_enriched / Molar_enriched)
log_ration <- log.mat(df_all[, c("Incisor_enriched", "Molar_enriched")], base = 10)
df_all$log_ratio <- with(log_ration, Incisor_enriched - Molar_enriched)

matched_index <- match(df_all$sample, rownames(metadata))
df_all$Position <- metadata[matched_index, "Position2"]
df_all$Group <- metadata[matched_index, "Future_Status_Tooth"]
df_all <- subset(df_all, Group != "C_H")

df_all$Group <- factor(df_all$Group, 
                       levels = c("ConfidentH", "RelativeH", "Caries"),
                       ordered = T)
g <- ggplot(df_all, aes(x = Group, y = log_ratio, color = Group)) +
  geom_boxplot(outlier.shape = NA) +
  scale_colour_manual(values = brewer.pal(8, "Pastel2")) +
  geom_jitter(aes(color = Group), shape = 1, width = 0.2) +
  facet_wrap(~Position, scales = "free_x", nrow = 1) +
  geom_signif(comparisons = list(c("ConfidentH", "RelativeH"), c("ConfidentH", "Caries")), 
              test = "wilcox.test", textsize = 2, step_increase = 0.1, map_signif_level = TRUE,
              test.args = list(exact = FALSE, correct = FALSE, conf.int = TRUE, 
                               conf.level = 0.95), color = "black") + # 添加wilcoxon test结果，并使不同分组的检验间隔0.01
  xlab("Group") + 
  labs(y = expression(log(frac("Incisor-enriched bacteria", "Molar-enriched bacteria")))) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

g

ggsave("../../Results/Figure_8/MiC/1867_rarefied_asv_incisor_molar_ratio_boxplot.pdf", g, height = 4, width = 10)

matched_index <- match(data$sample, rownames(metadata))
data$dataset <- metadata[matched_index, "Position2"]
data$index <- with(data, paste(feature, dataset, sep = "_"))
HHCC_imp_r <- read.table("../../Results/Figure_6/MiC/rarefied_taxonomy_Position_crossRF_out/HHCC_all_Taxon_OTU.wilcox_res.tsv", sep = "\t", header = T)
HHCC_imp_r <- subset(HHCC_imp_r, rf_imps_rank < 32)
HHCC_imp_r$index <- with(HHCC_imp_r, paste(feature, dataset, sep = "_"))
data <- subset(data, index %in% HHCC_imp_r$index)

data_32 <- data %>% 
  group_by(sample, feature_group) %>%
  summarise(sum_abd = sum(abundance))

data_32 <- subset(data_32, feature_group != "Neutral")
df_32 <- dcast(data_32, sample ~ feature_group)
df_32$ratio <- with(df_32, Incisor_enriched / Molar_enriched)
log_ration <- log.mat(df_32[, c("Incisor_enriched", "Molar_enriched")], base = 10)
df_32$log_ratio <- with(log_ration, Incisor_enriched - Molar_enriched)

matched_index <- match(df_32$sample, rownames(metadata))
df_32$Position <- metadata[matched_index, "Position2"]
df_32$Group <- metadata[matched_index, "Future_Status_Tooth"]
df_32 <- subset(df_32, Group != "C_H")

df_32$Group <- factor(df_32$Group, 
                       levels = c("ConfidentH", "RelativeH", "Caries"),
                       ordered = T)
g <- ggplot(df_32, aes(x = Group, y = log_ratio, color = Group)) +
  geom_boxplot(outlier.shape = NA) +
  scale_colour_manual(values = brewer.pal(8, "Pastel2")) +
  geom_jitter(aes(color = Group), shape = 1, width = 0.2) +
  facet_wrap(~Position, scales = "free_x", nrow = 1) +
  geom_signif(comparisons = list(c("ConfidentH", "RelativeH"), c("ConfidentH", "Caries")), 
              test = "wilcox.test", textsize = 2, step_increase = 0.1, map_signif_level = TRUE,
              test.args = list(exact = FALSE, correct = FALSE, conf.int = TRUE, 
                               conf.level = 0.95), color = "black") + # 添加wilcoxon test结果，并使不同分组的检验间隔0.01
  xlab("Group") + 
  labs(y = expression(log(frac("Incisor-enriched bacteria", "Molar-enriched bacteria")))) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

g

ggsave("../../Results/Figure_8/MiC/1867_rarefied_32_asv_incisor_molar_ratio_boxplot.pdf", g, height = 4, width = 10)

