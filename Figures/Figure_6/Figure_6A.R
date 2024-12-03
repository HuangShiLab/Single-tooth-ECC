## install and load necessary libraries for data analyses
#-------------------------------
p <- c("reshape2","ggplot2","pheatmap","combinat","randomForest", "gridExtra", "crossRanger", "RColorBrewer",
       "vegan", "biomformat", "doMC", "cowplot", "pROC", "ggsignif", "ggpubr", "colorspace")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
#-------------------------------
## input args
#-------------------------------
rf.opts<-list(outdir=NA, ntree=5000, verbose=FALSE, nfolds=10)
#source("ranger_util_20190520.R")
source("../Utilities/data_trimming_util.R")
#-------------------------------
datafile <- "../../Data/1867/taxonomy/1867_taxonomic_rarefied_feature_abd_table.biom"
sample_metadata <- "../../Data/1867/1867_metadata.txt"
feature_metadata <- "../../Data/1867/taxonomy/1867_taxonomic_feature_taxon_GTDB.txt"
prefix_name<-"taxonomy"
s_category<-"Position2"
c_category<-"Future_Status_Tooth"
AddTaxonomy=TRUE
p_cutoff=0.05
q_cutoff=0.05
p.adj.method="BH"
h=10
outpath <- paste0("../../Results/Figure_6/")
if(!dir.exists(outpath)) {
  dir.create(outpath, recursive = T)
}
#-------------------------------
# Biom table input
#-------------------------------
clinical_features <- c("sum_s_dt", "sum_ns_dt", "sum_dmfs", "spatial_dist_weighted_mean_dmfs")
mic_dist_features <- c("spatial_dist_weighted_mean_md", "mean_md", "mean_md_H_T51", "mean_md_H_T52",
                       "mean_md_H_T53", "mean_md_H_T54", "mean_md_H_T55", "mean_md_H_T61", "mean_md_H_T62",
                       "mean_md_H_T63", "mean_md_H_T64", "mean_md_H_T65", "mean_md_H_T71", "mean_md_H_T72",
                       "mean_md_H_T73", "mean_md_H_T74", "mean_md_H_T75", "mean_md_H_T81", "mean_md_H_T82",
                       "mean_md_H_T83", "mean_md_H_T84", "mean_md_H_T85", "mean_md_H_T5161")
if(grepl("biom$", datafile)){
  df <- read_hdf5_biom(datafile)
  df <- biom(df)
  df <- biom_data(df)
  df <- t(as.matrix(df))
  df <- data.frame(df, check.names = F)
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
#'-------------------------------
#' Train data: filtering
#'-------------------------------
data_list<-filter_samples_by_sample_ids_in_metadata(df, dm = F, metadata)
#'-------------------------------
#' Train data: filter out samples with null values in target_field
#'-------------------------------
data_list<-filter_samples_by_NA_in_target_field_of_metadata(data_list$data, data_list$metadata, target_field = s_category)
#'-------------------------------
#' Train data: filter out samples in particular groups
#'-------------------------------
data_HHCC_list<-filter_samples_by_groups_in_target_field_of_metadata(data_list$data, data_list$metadata,
                                                                     target_field = c_category, negate=FALSE,
                                                                     groups = c("ConfidentH", "Caries"))
#-------------------------------
# rf_clf.by_datasets
#-------------------------------
## "rf_clf.by_datasets" runs standard random forests with oob estimation for classification of
## c_category in each the sub-datasets splited by the s_category.
## The output includes a summary of rf models in the sub datasets
## and all important statistics for each of features.

res_file<-paste("../../Results/RF_results/MiC/rarefied_",prefix_name,"_Position_crossRF_out/", prefix_name, "_rf_clf.by_datasets_res.RData", sep="")
if(file.exists(res_file)){
  load(res_file)
}else{
  rf_clf_res<-rf_clf.by_datasets(data_HHCC_list$data, data_HHCC_list$metadata, s_category, c_category,
                                 positive_class="Caries", clr_transform = TRUE, rf_imp_pvalues = FALSE,
                                 p.adj.method=p.adj.method, q_cutoff=q_cutoff, nfolds=5, verbose=FALSE, ntree=500)
  save(rf_clf_res, file=res_file)
}

rf_clf_res$feature_imps_list<-lapply(rf_clf_res$feature_imps_list, function(x) {x[x$mean_all==0, "rf_imps"]<-NA; return(x)})

## RF classification performance VS number of features used
# Prediction performances at increasing number of features obtained by
# retraining the random forest regressor on the top-ranking features identified
# with a first random forest model training in a cross-validation setting
top_n_perf_list<-list()
for(n in 1:length(rf_clf_res$rf_model_list)){
  rf_imps<-rf_clf_res$feature_imps_list[[n]][, "rf_imps"]
  rf_imps_rank<-rank(-rf_imps, na.last = "keep")
  x<-rf_clf_res$x_list[[n]]
  max_n<-ncol(x)
  n_features<-c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
  n_features<-n_features[n_features < max_n]
  nrow<-length(n_features)+1
  top_n_perf<-matrix(NA, ncol=2, nrow=nrow)
  colnames(top_n_perf)<-c("n_features", "AUROC")
  rownames(top_n_perf)<-top_n_perf[,1]<-c(n_features, max_n)
  for(i in 1:length(n_features)){
    idx<-which(rf_imps_rank<=n_features[i])
    x_n<-x[, idx]
    y_n<-rf_clf_res$y_list[[n]]
    top_n_rf<-rf.out.of.bag(x_n, y_n, ntree=5000)
    rf_AUROC<-get.auroc(top_n_rf$probabilities[, "Caries"], y_n, positive_class="Caries")
    top_n_perf[i, 1]<-n_features[i]
    top_n_perf[i, 2]<-rf_AUROC
  }
  top_n_perf[nrow, ]<-c(max_n, rf_clf_res$rf_AUROC[n])
  top_n_perf_list[[n]]<-top_n_perf
}
names(top_n_perf_list)<-rf_clf_res$datasets
top_n_perf_list<-lapply(1:length(top_n_perf_list),
                        function(x) data.frame(Dataset=rep(names(top_n_perf_list)[x], nrow(top_n_perf_list[[x]])), top_n_perf_list[[x]]))

top_n_perf_comb<-do.call(rbind, top_n_perf_list)
top_n_perf_comb$n_features<-as.numeric(as.character(top_n_perf_comb$n_features))
top_n_perf_comb_m<-melt(top_n_perf_comb, id.vars = c("n_features", "Dataset"))
breaks<-top_n_perf_comb_m$n_features
p<-ggplot(subset(top_n_perf_comb_m, variable=="AUROC"), aes(x=n_features, y=value)) +
  xlab("# of features used")+
  ylab("AUROC")+
  scale_x_continuous(trans = "log",breaks=breaks)+
  ylim(c(0.5, 1))+
  geom_point(aes(color = Dataset)) +
  geom_line(aes(color = Dataset)) +
  scale_color_manual(values = darken(brewer.pal(9, "Pastel1"), amount = 0.2)) +
  #facet_wrap(~Dataset) +
  theme_bw()+
  theme(axis.line = element_line(color="black"),
        axis.title = element_text(size=18),
        strip.background = element_rect(colour = "white"),
        panel.border = element_blank(),
        axis.text.x=element_text(angle=90, hjust=1))
ggsave(filename=paste(outpath,"Figure_6A_AUROC__top_rankings.scatterplot.pdf",sep=""), plot=p, width=5, height=4)
