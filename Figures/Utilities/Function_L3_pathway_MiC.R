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
datafile <- "../../Data/1867/function/1867_functional_l3_abd_table.biom"
sample_metadata <- "../../Data/1867/1867_metadata.txt"
feature_metadata <- NA
prefix_name<-"Function_L3"
s_category<-"Position2"
c_category<-"Future_Status_Tooth"
AddTaxonomy=TRUE
p_cutoff=0.05
q_cutoff=0.05
p.adj.method="BH"
h=10
outpath <- paste("../../Results/RF_results/MiC/",prefix_name,"_Position_crossRF_out/", sep="")
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

res_file<-paste(outpath, prefix_name, "_rf_clf.by_datasets_res.RData", sep="")
if(file.exists(res_file)){
  load(res_file)
}else{
  rf_clf_res<-rf_clf.by_datasets(data_HHCC_list$data, data_HHCC_list$metadata, s_category, c_category,
                                 positive_class="Caries", clr_transform = TRUE, rf_imp_pvalues = FALSE,
                                 p.adj.method=p.adj.method, q_cutoff=q_cutoff, nfolds=5, verbose=FALSE, ntree=500)
  save(rf_clf_res, file=res_file)
}

rf_clf_res$feature_imps_list<-lapply(rf_clf_res$feature_imps_list, function(x) {x[x$mean_all==0, "rf_imps"]<-NA; return(x)})

rf_clf_res.summ<-plot_clf_res_list(rf_clf_res, p_cutoff=0.05, p.adj.method = p.adj.method, q_cutoff=0.05, outdir=outpath)

wilcox_res<-rf_clf_res.summ$feature_res
rf_models<-rf_clf_res$rf_model_list
# Add feature annotations using feature metadata
rbind.na<-function(l){
  max_len<-max(unlist(lapply(l, length)))
  c_l<-lapply(l, function(x) {c(x, rep(NA, max_len - length(x)))})
  do.call(rbind, c_l)
}
expand_Taxon<-function(df, Taxon){
  taxa_df <- rbind.na(strsplit(as.character(df[, Taxon]), ';'))
  colnames(taxa_df) <- c("phylum","class","order","family","genus","species") #"kingdom",
  data.frame(df, taxa_df)
}

if(AddTaxonomy==TRUE){
  wilcox_res<-add_ann(wilcox_res, fmetadata)
  #wilcox_res<-expand_Taxon(wilcox_res, "Taxon")
}else{
  wilcox_res$Taxon=wilcox_res$feature
}
## modify the rf imp scores
#wilcox_res[wilcox_res$mean_all==0, "rf_imps"]<-NA
library("dplyr")
wilcox_res<-wilcox_res %>% dplyr::group_by(dataset) %>% mutate(rf_imps_rank=rank(-rf_imps, na.last = "keep"))

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
  ylim(c(0, 1)) +
  geom_point() + geom_line() +facet_wrap(~Dataset) +
  theme_bw()+
  theme(axis.line = element_line(color="black"),
        axis.title = element_text(size=18),
        strip.background = element_rect(colour = "white"),
        panel.border = element_blank(),
        axis.text.x=element_text(angle=90, hjust=1))
ggsave(filename=paste(outpath,"AUROC__top_rankings.facets.scatterplot.pdf",sep=""), plot=p, width=7, height=6)
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
ggsave(filename=paste(outpath,"AUROC__top_rankings.scatterplot.pdf",sep=""), plot=p, width=5, height=4)
sink(paste(outpath,"top_n_feature_perf.xls",sep=""));write.table(top_n_perf_comb,quote=FALSE,sep="\t", row.names = F);sink()



## Non-specific features across datasets (at least present in two of datasets)
res_spcf<-id_non_spcf_markers(wilcox_res, positive_class="Caries", other_class="ConfidentH", p.adj.method, outdir=outpath)
#summary(res_spcf)
wilcox_res_spcf<-res_spcf$feature_res_spcf


## keeps only OTUs which were significant in at least one datasets.
keep_sig_markers<-function(wilcox_res, q_cutoff){
  wilcox_res_sig<-wilcox_res[which(wilcox_res[, "feature"] %in% unique(as.character(subset(wilcox_res, non.param.test_p.adj < q_cutoff)[, "feature"]))), ]
  wilcox_res_sig
}

wilcox_res_spcf_sig<-keep_sig_markers(wilcox_res_spcf, q_cutoff=q_cutoff)
sink(paste(outpath,"BetweenGroupTest_out_all.xls",sep=""));write.table(wilcox_res_spcf,quote=FALSE,sep="\t", row.names = F);sink(NULL)
sink(paste(outpath,"BetweenGroupTest_out_sig_",p.adj.method,"_q_cutoff_",q_cutoff,".xls",sep=""));write.table(wilcox_res_spcf_sig,quote=FALSE,sep="\t", row.names = F);sink(NULL)
#-------------------------------
# Plotting
#-------------------------------
## customised colors for Enr's 3 factors: XX-depleted, XX-enriched and Neutral
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
my3cols<-c("grey60", rev(gg_color_hue(2)))
#-------------------------------
# The performance of rf models
#-------------------------------
m_s<-split(data_HHCC_list$metadata, data_HHCC_list$metadata[, s_category])
rocobj_list<-list()
ciobj_list<-list()
for(i in 1:length(rf_models)){
  sub_group<-m_s[[i]][, c_category]
  #--------------------------------------------------
  #  ROC plot using "pROC" package
  #--------------------------------------------------
  pdf(paste(outpath, "rf_", names(m_s)[i],".pROC.ci.pdf",sep=""), width=4, height=4)
  rocobj_list[[i]]<-rocobj <- plot.roc(sub_group, rf_models[[i]]$probabilities[,2],main="", percent=TRUE,ci=TRUE) # print the AUROC (will contain the CI)
  ciobj_list[[i]]<-ciobj <- ci.se(rocobj,specificities=seq(0, 100, 5)) # over a select set of specificities
  plot(ciobj, type="shape", col="#1c61b6AA", print.auc=TRUE) # plot as a blue shape
  text(70,25, paste0(levels(sub_group), collapse = " VS "), pos=4)
  text(70,15, paste("AUROC = ",formatC(rocobj$auc,digits=2,format="f"),sep=""),pos=4)
  ci.lower<-formatC(rocobj$ci[1],digits=2,format="f")
  ci.upper<-formatC(rocobj$ci[3],digits=2,format="f")
  text(70,5, paste("95% CI: ",ci.lower,"-",ci.upper,sep=""),pos=4)
  dev.off()
}

#-------------------------------
# AUROC distribution
#-------------------------------
p<-ggplot(wilcox_res_spcf_sig, aes(AUROC, fill=dataset)) + geom_histogram(alpha=0.4) +
  #scale_fill_manual(values = my3cols)+
  facet_wrap(~Enr, nrow=1)
ggsave(filename=paste(outpath,"Markers_AUROC_distr_",p.adj.method, ".",c_category, "_by_", s_category,".pdf",sep=""),plot=p, width=8, height=3)
#-------------------------------
# Mean VS sd of markers in each of datasets
#-------------------------------
l_m<-wilcox_res_spcf_sig[, c("feature", "Taxon","dataset","mean_all", "sd_all", "IfSig", "IfSigEnr", "Enr")]
p<-ggplot(l_m, aes(x=log(mean_all), y=log(sd_all), color=Enr)) + scale_color_manual(values = my3cols) +
  geom_point(stat="identity", alpha=.5) + xlab("log10(mean)")+ ylab("log10(sd)")+
  geom_abline(slope = 1, intercept = 0)+
  coord_flip()+ # if want to filp coordinate
  theme_bw()+
  facet_wrap(~ dataset, nrow=1)
ggsave(filename=paste(outpath,"Markers_mean_VS_sd_",p.adj.method,"_",c_category, "_by_", s_category,".pdf",sep=""),plot=p, width=8, height=4)

#-------------------------------
# Statistics summary of features in all datasets
#-------------------------------
rf_clf_res.summ_all<-plot_clf_res_list(rf_clf_res, p_cutoff=1, p.adj.method = p.adj.method, q_cutoff=1, outdir=outpath)
all_wilcox_res<-rf_clf_res.summ_all$feature_res
if(AddTaxonomy==TRUE){
  all_wilcox_res<-add_ann(all_wilcox_res, fmetadata)
  #wilcox_res<-expand_Taxon(wilcox_res, "Taxon")
}else{
  all_wilcox_res$Taxon=all_wilcox_res$feature
}
#-------------------------------
# RF imps
#-------------------------------
if(any(grepl("package:plyr", search()))) detach("package:plyr") else message("plyr not loaded")
## plyr not loaded
library("dplyr")
all_wilcox_res$Taxon<-as.factor(all_wilcox_res$Taxon)
if(any(nchar(all_wilcox_res$feature)>20)){
  all_wilcox_res<-all_wilcox_res %>% dplyr::group_by(dataset) %>% mutate(asv_ids=paste("ASV_", seq(1, length(feature)),sep = ""))
  all_wilcox_res<-all_wilcox_res %>% dplyr::group_by(dataset) %>% mutate(rf_imps_rank=rank(-rf_imps, na.last = "keep"),
                                                                         feature_Taxon=paste(asv_ids, Taxon, sep="__"))
}else{
  all_wilcox_res<-all_wilcox_res %>% dplyr::group_by(dataset) %>% mutate(rf_imps_rank=rank(-rf_imps, na.last = "keep"),
                                                                         feature_Taxon=paste(feature, Taxon, sep="__"))
}

write.table(all_wilcox_res, paste0(outpath, "/HHCC_all_Taxon_OTU.wilcox_res.tsv"), sep = "\t", quote = F,
            row.names = T, col.names = NA)

all_Taxon_df<- all_wilcox_res %>%
  dplyr::group_by(Taxon, feature) %>%
  dplyr::summarise(n=length(rf_imps), n_NA=sum(is.na(rf_imps)), n_neg=sum(rf_imps<0, na.rm = T), n_pos=sum(rf_imps>0, na.rm = T),
                   imp_sum = sum(rf_imps, na.rm = T), imp_mean=mean(rf_imps, na.rm = T),
                   imp_rank_sum= sum(rf_imps_rank, na.rm = T), imp_rank_mean=mean(rf_imps_rank, na.rm = T),
                   imp_rank_min= min(rf_imps_rank, na.rm = T) #,
                   #imp_rank_min_Enr=ifelse(length(rf_imp_rank)>1 && !all(is.na(rf_imps_rank)), Enr[which.min(rf_imp_rank)], as.character(Enr))
                   #imp_rank_min_mean_logfc=mean_logfc[which.min(rf_imp_rank)]
  )

write.table(all_Taxon_df, paste(outpath, "all_Taxon_OTU.RF_imps.tsv", sep=""), sep="\t", quote=FALSE, row.names = FALSE)
cat("# of sum values: ", dim(all_Taxon_df)[1], '\n')
cat("# of Taxa: ", length(levels(as.factor(all_wilcox_res$Taxon))))
Taxon_asc_levels<-levels(as.factor(all_wilcox_res$Taxon))[with(all_Taxon_df,  order(imp_sum))]
Taxon_desc_levels<-levels(as.factor(all_wilcox_res$Taxon))[with(all_Taxon_df,  order(-imp_sum))]
# plotting imp rank of top-n features in multiple classifiers
top_N<-32
# select features that minimum rank < top_N and n_pos > 4
all_Taxon_df_top_n<-subset(all_Taxon_df, imp_rank_min< top_N & n_pos>7)
top_features <- as.character(all_Taxon_df_top_n$feature)
all_wilcox_res_top_N <- all_wilcox_res %>% filter(feature %in% top_features)
all_Taxon_df_top_n_acast<-data.frame(acast(all_wilcox_res_top_N, feature~dataset, value.var="rf_imps_rank"))

row_sum<-rowSums(all_Taxon_df_top_n_acast, na.rm = T)
NA_perc<-rowSums(is.na(all_Taxon_df_top_n_acast))/ncol(all_Taxon_df_top_n_acast)
NA_perc_idx<-which(NA_perc<=0.8); taxa_kept<-names(NA_perc)[NA_perc_idx]
row_peudosum4NA<-rowMeans(all_Taxon_df_top_n_acast, na.rm = T)*apply(all_Taxon_df_top_n_acast, 1, function(x) sum(is.na(x)))
row_sum_comb<-row_sum+row_peudosum4NA
ord<-order(row_sum_comb, decreasing = F)
all_wilcox_res_top_N$feature_Taxon<-factor(all_wilcox_res_top_N$feature_Taxon,
                                           levels=levels(factor(all_wilcox_res_top_N$feature_Taxon))[ord], ordered=TRUE)
# ggplot
p<-ggplot(all_wilcox_res_top_N, aes(dataset, feature_Taxon, fill = rf_imps_rank)) +
  geom_tile() + # adding "-" can switch the coloring
  geom_text(aes(label = round(rf_imps_rank, 1))) +
  #scale_y_discrete(label=f_ann_m$Taxon)+
  #scale_fill_continuous(type = "viridis", na.value = 'white')+
  #coord_flip()+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  scale_fill_gradient2(high ="#f46d43", mid="white", low= "steelblue", na.value = NA)

ggsave(filename=paste(outpath,"/all_Top_N_min_imp_rank_by_Taxon_VS_dataset_heatmap.ggplot.pdf",sep=""),plot=p, width=15, height=18)

#-------------------------------
# ALE plot
#-------------------------------
#library(ALEPlot)
#imp_rank_df<- wilcox_res %>%
#  dplyr::group_by(dataset) %>%  filter(rf_imps_rank<=32)
#summary(wilcox_res$rf_imps_rank)
#
#feature_imps_rank_list<-lapply(rf_clf_res$feature_imps_list, function(x)
#  with(x, data.frame(x, rf_imps_rank=rank(-rf_imps, na.last = "keep") ) ) )
#out_df_list<-list()
#for(n in 1:length(rf_clf_res$x_list)){
#  feature_imps_rank_n<-feature_imps_rank_list[[n]][,"rf_imps_rank"]
#  rf_model_n<-rf_clf_res$rf_model_list[[n]]
#  x_n<-rf_clf_res$x_list[[n]]
#  dataset_n<-rf_clf_res$datasets[[n]]
#  N=64
#  #pdf(paste(outpath,"dataset_", n,"_top_n_feature_ALEPlot.pdf",sep=""), width=40, height=20)
#  #par(mfrow = c(4,8))
#  ale_out_n<-list()
#  out_n<-matrix(NA, ncol=7, nrow=N)
#  colnames(out_n)<-c("Feature_rank", "Sparsity", "Prevalence", "N", "K", "Corr", "P_value")
#  yhat <- function(X.model, newdata) predict(X.model, newdata)
#  for(i in 1:N){
#    FeatureIdx<-which(feature_imps_rank_n==i)
#    FeatureID<-names(FeatureIdx)
#    out<-ALEPlot(x_n, rf_model_n, pred.fun = yhat, J=FeatureIdx, K=10000)
#    ale_out_n[[i]]<-data.frame(FeatureID, out)
#    out_corr<-cor.test(out[[2]], out[[3]], method="spearman")
#    out_n[i,1]<-feature_imps_rank_n[FeatureIdx]
#    out_n[i,2]<-sum(x_n[,FeatureIdx]==0)/nrow(x_n)
#    out_n[i,3]<-sum(x_n[,FeatureIdx]!=0)/nrow(x_n)
#    out_n[i,4]<-sum(x_n[,FeatureIdx]!=0)
#    out_n[i,5]<-out$K
#    out_n[i,6]<-out_corr$estimate
#    out_n[i,7]<-out_corr$p.value
#  }
#  #names(ale_list)<-names(feature_imps_rank_n)[which(feature_imps_rank_n<=N)]
#  ale_out_n<-do.call(rbind, ale_out_n)
#  p<-ggplot(ale_out_n, aes(x=x.values, y=f.values))+geom_point(alpha=0.5) + geom_line()+
#    geom_rug() + facet_wrap(~FeatureID, scales="free") + theme_bw()
#  ggsave(filename=paste(outpath,"dataset_",n,"_top_n_feature_ALEPlot.pdf",sep=""),plot=p, width=20, height=20)
#  #dev.off()
#  out_df_list[[n]]<-out_df_n<-data.frame(Dataset=rep(dataset_n, N), FeatureID=names(sort(feature_imps_rank_n[which(feature_imps_rank_n<=N)])), out_n)
#}
#
#
#ale_out_df<-reshape::merge_all(ale_out_df_list)
#ale_out_df
#sink(paste(outpath,"RF_ALEPlot_summ.xls",sep=""));write.table(ale_out_df,quote=FALSE,sep="\t", row.names = F);sink()

#-------------------------------
## fold change
#-------------------------------
#---pheatmap: Occurence rate heatmap of all  biomarkers
OccRate_dataset_df<-acast(wilcox_res_spcf_sig, Taxon~dataset, value.var="OccRate_all", fun.aggregate = mean)
pheatmap(OccRate_dataset_df, cluster_cols = F, filename = paste(outpath,"Taxa_OccRate_all_VS_dataset_",p.adj.method,"_",c_category, "_by_", s_category,".pdf",sep=""), width = 20)
# pheatmap(OccRate_dataset_df, cluster_cols = F, filename = paste(outpath,"Taxa_OccRate_all_VS_dataset_",p.adj.method,"_",c_category, "_by_", s_category,".pdf",sep=""))
## pheatmap: heatmap in mean log2 fold change of markers across datasets
logfc_dataset_df<-acast(wilcox_res_spcf_sig, Taxon~dataset, value.var="mean_logfc", fun.aggregate = mean)
head(logfc_dataset_df)
pheatmap(logfc_dataset_df, cluster_cols = F, filename = paste(outpath,"Taxa_mean_logfc_VS_dataset_",p.adj.method,"_",c_category, "_by_", s_category,".pdf",sep=""), width = 20)
# pheatmap(logfc_dataset_df, cluster_cols = F, filename = paste(outpath,"Taxa_mean_logfc_VS_dataset_",p.adj.method,"_",c_category, "_by_", s_category,".pdf",sep=""))
#---ggplot2: heatmap in mean log2 fold change of markers across datasets
feature_logfc_dataset_df<-acast(wilcox_res_spcf_sig, feature~dataset, value.var="mean_logfc")
ord <- hclust(dist(feature_logfc_dataset_df))$order
mean_logfc_m<-wilcox_res_spcf_sig[, c("dataset", "feature", "Taxon", "IfSig", "Enr", "mean_logfc")]
#-------------------------------barchart
p_logfc<-ggplot(mean_logfc_m, aes(x=Taxon, y=mean_logfc, fill=Enr)) + scale_fill_manual(values = my3cols) +
  geom_bar(stat="identity", alpha=.5) + ylab("mean log2(fold change)")+
  coord_flip()+ # if want to filp coordinate
  theme_bw()+
  facet_wrap(~ dataset, nrow=1)+
  theme(axis.line = element_line(color="black"),
        strip.background = element_rect(colour = "white"),
        panel.border = element_blank(),
        axis.text.x=element_text(angle=90, hjust=1))
ggsave(filename=paste(outpath,"Markers_logfc_",p.adj.method,"_",c_category, "_by_", s_category,".barplot.pdf",sep=""),plot=p_logfc, width=18, height=h)
#-------------------------------geom_tile
mean_logfc_m[which(mean_logfc_m$IfSig=="NotSig"), "mean_logfc"]<-NA
mean_logfc_m$feature<-factor(mean_logfc_m$feature,levels=rownames(feature_logfc_dataset_df)[ord], ordered=TRUE)
#-------------------------------
p<-ggplot(mean_logfc_m, aes(dataset, feature) ) +
  geom_tile(data=mean_logfc_m, na.rm = FALSE, aes(fill = -mean_logfc), colour = 'black') + # adding "-" can switch the coloring
  scale_y_discrete(label=mean_logfc_m$Taxon[ord])+
  #scale_fill_continuous(type = "viridis", na.value = 'salmon')+
  #coord_flip()+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  scale_fill_gradient2(high ="#f46d43", mid="white", low= "steelblue", na.value = NA)
ggsave(filename=paste(outpath,"Markers_mean_logfc_VS_dataset_",p.adj.method,"_",c_category, "_by_", s_category,".heatmap.ggplot.pdf",sep=""),plot=p, width=18, height=h)


#-------------------------------
# AUROC
#-------------------------------
l_m<-melt(wilcox_res_spcf_sig[, c("feature", "Taxon","dataset","AUROC", "IfSigEnr", "Enr")])
p_auc<-ggplot(l_m, aes(x=Taxon, y=value, color=Enr)) + ylim(0.5, 1) + ylab("AUROC")+ scale_color_manual(values = my3cols) +
  geom_abline(slope=0, intercept=0.5,  col = "black", lty=2) +
  geom_point(aes(colour=Enr), shape="diamond", size=4) +
  coord_flip()+ # if want to filp coordinate
  theme_bw()+
  facet_wrap(~ dataset, nrow=1)+
  theme(axis.line = element_line(color="black"),
        strip.background = element_rect(colour = "white"),
        panel.border = element_blank(),
        axis.text.x=element_text(angle=90, hjust=1))
ggsave(filename=paste(outpath,"Markers_AUROC_sig_",p.adj.method, ".", c_category, "_by_", s_category,".scatterplot.pdf",sep=""),plot=p_auc,width=18,height=h)
#-------------------------------
# -log10(q values)
#-------------------------------
l_m<-melt(wilcox_res_spcf_sig[, c("feature", "Taxon","dataset","non.param.test_p.adj", "IfSigEnr", "Enr")])
p_q<-ggplot(l_m, aes(x=Taxon, y=-log10(value), fill=Enr)) + scale_fill_manual(values = my3cols) +
  geom_bar(stat="identity", alpha=.5) + ylab("-log10(q-value) (Wilcoxon.test)")+
  coord_flip()+ # if want to filp coordinate
  theme_bw()+
  facet_wrap(~ dataset, nrow=1)+
  theme(axis.line = element_line(color="black"),
        strip.background = element_rect(colour = "white"),
        panel.border = element_blank(),
        axis.text.x=element_text(angle=90, hjust=1))
ggsave(filename=paste(outpath,"Markers_Wilcoxon.test_p.adj_",p.adj.method,"_",c_category, "_by_", s_category,".barplot.pdf",sep=""),plot=p_q, width=18, height=h)
#-------------------------------
# OccRate
#-------------------------------
OccRateCols<-colnames(wilcox_res_spcf_sig)[grep("^(?=.*OccRate)(?!.*all)", colnames(wilcox_res_spcf_sig), perl=TRUE)]
l_m<-melt(wilcox_res_spcf_sig[, c("feature", "Taxon", "dataset", OccRateCols, "IfSigEnr", "Enr")])
tmp<-do.call(rbind, strsplit(as.character(l_m[, "variable"]), "__")); colnames(tmp)<-c("OccRate", c_category)
l_m_OccRate<-data.frame(tmp, l_m)

p_occ<-ggplot(l_m_OccRate, aes(x=Taxon, y=value, fill=get(c_category))) +
  geom_bar(position=position_dodge(), stat="identity", alpha=.5) +
  ylab("Occurrence rate")+ #scale_fill_manual(values = my3cols[2:3]) +
  coord_flip()+ # if want to filp coordinate
  theme_bw()+
  facet_grid( ~ dataset, scales="free")+
  theme(axis.line = element_line(color="black"),
        strip.background = element_rect(colour = "white"),
        panel.border = element_blank(),
        axis.text.x=element_text(angle=90, hjust=1))
ggsave(filename=paste(outpath,"OccRate_sig_",p.adj.method, ".",c_category, "_by_", s_category,".pdf",sep=""), plot=p_occ, width=18, height=h)

#-------------------------------
# Combined plots: require(cowplot)
#-------------------------------
p_logfc_ <- p_logfc + theme(legend.position="none")
p_auc_ <- p_auc + theme(axis.text.y =element_blank(), axis.title.y = element_blank())+ theme(legend.position="none")
p_occ_ <- p_occ + theme(axis.text.y =element_blank(), axis.title.y = element_blank())+ theme(legend.position="none")
p_q <- p_q + theme(axis.text.y =element_blank(), axis.title.y = element_blank())
combined_plot=plot_grid(p_logfc_, p_auc_, p_occ_, p_q, nrow=1, rel_widths = c(5, 2, 2, 2.5), labels="AUTO")
combined_plot
ggsave(filename=paste(outpath,"Stats_combined_plot.pdf", sep=""), plot=combined_plot, width=25, height = 6)
#-------------------------------
# Mean +/- sd abundance
#-------------------------------
# Mean abundance
meanCols<-colnames(wilcox_res_spcf_sig)[grep("^(?=.*mean)(?!.*all)(?!.*logfc)", colnames(wilcox_res_spcf_sig), perl=TRUE)]
l_m<-melt(wilcox_res_spcf_sig[, c("feature", "Taxon", "dataset", meanCols, "IfSigEnr", "Enr")])
tmp<-do.call(rbind, strsplit(as.character(l_m[, "variable"]), "__")); colnames(tmp)<-c("MeanAbd", c_category)
l_m_MeanAbd<-data.frame(tmp, l_m)
# sd
sdCols<-colnames(wilcox_res_spcf_sig)[grep("^(?=.*sd)(?!.*all)(?!.*logfc)", colnames(wilcox_res_spcf_sig), perl=TRUE)]
l_m<-melt(wilcox_res_spcf_sig[, c("feature", "Taxon", "dataset", sdCols, "IfSigEnr", "Enr")])
tmp<-do.call(rbind, strsplit(as.character(l_m[, "variable"]), "__")); colnames(tmp)<-c("sd", c_category)
l_m_sd<-data.frame(tmp, l_m)
# Mean +/- sd abundance
l_m_MeanAbd_sd<-data.frame(l_m_MeanAbd[, -ncol(l_m_MeanAbd)+1 : -ncol(l_m_MeanAbd)], mean=l_m_MeanAbd$value, sd=l_m_sd$value)

p_mean_sd<-ggplot(l_m_MeanAbd_sd, aes(x=Taxon, y=mean, group=get(c_category), fill=Enr)) +
  geom_bar(position=position_dodge(), stat="identity", alpha=.5) + scale_fill_manual(values = my3cols) +
  ylab("Relative abundance")+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  coord_flip()+ # if want to filp coordinate
  theme_bw()+
  facet_grid( ~ dataset, scales="free")
ggsave(filename=paste(outpath,"MeanAbd_sig_",p.adj.method, ".",c_category, "_by_", s_category,".pdf",sep=""),plot=p_mean_sd, width=18, height=h)


#-------------------------------
# rf.cross.datasets
#-------------------------------
# "rf.cross.datasets" runs standard random forests with oob estimation for classification of
# c_category in each the sub-datasets splited by the s_category,
# and apply the model to all the other datasets. The output includes
# accuracy, auc and Kappa statistics.
# perf_summ<-rf.cross.datasets(df_k, metadata_k, s_category, c_category, nfolds=10, verbose=FALSE, ntree=5000)

crossRF_res<-rf_clf.cross_appl(rf_clf_res$rf_model_list, rf_clf_res$x_list, rf_clf_res$y_list, positive_class="Caries")
perf_summ<-crossRF_res$perf_summ
p<-ggplot(perf_summ, aes(x=Validation_type, y=AUROC, color=Test_data)) +
  geom_point(aes(shape=Validation_type), size=3, alpha=0.4) +
  facet_grid(~Train_data)+
  #geom_line(aes(x=Validation_type, y=AUROC, color=Train_data, group=Train_data))+
  theme_bw()+
  theme(axis.text.x= element_text(angle = 45, hjust=1))
ggsave(filename=paste(outpath,"AUROC_cross_application_",c_category, "_among_", s_category,".pdf",sep=""),plot=p, width=10, height=5)

sink(paste(outpath,"crossRF_perf_summ.xls",sep=""));write.table(perf_summ,quote=FALSE,sep="\t", row.names = F);sink()

#' The performance of cross-applications
self_validation=as.factor(perf_summ$Train_data==perf_summ$Test_data)
library(viridis)
p_AUROC<-ggplot(perf_summ, aes(x=as.factor(Test_data), y=as.factor(Train_data), z=AUROC)) +
  xlab("Test data") +
  ylab("Train data") +
  labs(fill = "AUROC\n(ConfidentH\nVS\nCaries)") +
  geom_tile(aes(fill = AUROC, color = self_validation, width=0.9, height=0.9), size=1) +
  scale_color_manual(values=c("white","grey80")) +
  geom_text(aes(label = round(AUROC, 2)), color = "white") +
  scale_fill_viridis(limit = c(0.5, 1))+
  theme_bw() + theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.line = element_blank(),
        axis.ticks = element_blank())
p_AUROC
ggsave(filename=paste(outpath,"AUROC_cross_appl_matrix_",c_category, "_among_", s_category,".heatmap.pdf",sep=""),plot=p_AUROC, width=6, height=4.6)


# Only CH and RH groups in the StatusToothPre kept for prediction
#'-------------------------------
#' Train data: filter out samples in particular groups
#'-------------------------------
data_CHRH_list<-filter_samples_by_groups_in_target_field_of_metadata(data_list$data, data_list$metadata,
                                                                     target_field = c_category, negate=FALSE,
                                                                     groups = c("C_H", "RelativeH"))
# data list
pred_x_list<-split(data_CHRH_list$data, data_CHRH_list$metadata[, s_category])
pred_y_list<-split(data_CHRH_list$metadata[, c_category], data_CHRH_list$metadata[, s_category])

rf_clf.pred<-function(rf_model_list, x_list, y_list, newx_list, newy_list, positive_class=NA)
{
  L<-length(rf_model_list)
  positive_class<-ifelse(is.na(positive_class), levels(factor(y_list[[1]]))[1], positive_class)
  try(if(!identical(L, length(x_list), length(y_list))) stop("The length of x list, y list and rf model list should be identical."))
  perf_summ<-data.frame(matrix(NA, ncol=17, nrow=L))
  colnames(perf_summ)<-c("Train_data", "Test_data", "Validation_type", "Accuracy", "AUROC", "Kappa",
                         "Sensitivity", "Specificity", "Pos_Pred_Value","Neg_Pred_Value", "Precision", "Recall",
                         "F1", "Prevalence", "Detection_Rate", "Detection_Prevalence", "Balanced_Accuracy")
  predicted<-matrix(list(), ncol=2, nrow=L)
  colnames(predicted)<-c("train_predicted", "test_predicted")
  for(i in 1:L){
    y<-y_list[[i]]
    x<-x_list[[i]]
    try(if(nlevels(y)==1) stop("Less than one level in the subgroup for classification"))
    #rf_model<-randomForest(x, y, ntree=5000, importance=T)
    #oob<-rf.out.of.bag(x, y, nfolds=nfolds, verbose=verbose, ntree=ntree)
    oob<-rf_model_list[[i]]
    #---
    #  RF Training accuracy
    #---
    cat("\nTraining dataset: ", names(x_list)[i] ,"\n\n")
    conf<-caret::confusionMatrix(data=oob$predicted, oob$y, positive=positive_class)
    acc<-conf$overall[1]
    kappa_oob<-conf$overall[2]
    cat("Accuracy in the self-validation: ", acc ,"\n")
    #---
    #  AUROC computation using "pROC" package
    #---
    auc<-get.auroc(oob$probabilities[, positive_class], oob$y, positive_class)
    cat("AUROC in the self-validation: ", auc ,"\n")
    perf_summ[i, 1:3]<-c(names(x_list)[i], names(x_list)[i], "self_validation")
    perf_summ[i, 4:17]<-c(acc, auc, kappa_oob, conf$byClass)
    #predictions
    cat("Predictions: confusion matrix \n")
    newx<-newx_list[[i]]
    newy<-newy_list[[i]]
    #pred_prob<-predict(oob$rf.model, x_list[[j]], type="prob") # for regular rf.out.of.bag (randomForest) function
    # add a 0 matrix with features that are not in the rf model
    if(class(oob)=="rf.cross.validation"){
      model_features <- oob$rf.model[[1]]$forest$independent.variable.names
    }else{
      model_features <- oob$rf.model$forest$independent.variable.names
    }
    undected_ids <- model_features[!model_features %in% colnames(newx)]
    zero_mat<-matrix(0, nrow=nrow(newx), ncol=length(undected_ids))
    colnames(zero_mat) <- undected_ids
    newx <- data.frame(newx, zero_mat, check.names = FALSE)
    if(class(oob)=="rf.cross.validation"){
      pred_prob <- get.predict.probability.from.forest(oob$rf.model[[1]], newx) # only use one of n models for prediction
      pred_newy<-factor(predict(oob$rf.model[[1]], newx, type="response")$predictions)
    }else{
      pred_prob <- get.predict.probability.from.forest(oob$rf.model, newx) # ranger only
      pred_newy<-factor(predict(oob$rf.model, newx, type="response")$predictions)
    }
    pred_prob<-pred_prob[, order(colnames(pred_prob))] # to avoid unanticipated order of numeric levels of factor y
    # ranger only
    colnames(pred_prob)<- levels(oob$y)
    levels(pred_newy)<- levels(oob$y)
    print(table(newy,pred_newy))
    predicted[i, 1][[1]]<-data.frame(y=y, pred_y=oob$predicted, oob$probabilities); names(predicted[i, 1]) <-names(x_list)[i]
    predicted[i, 2][[1]]<-data.frame(y=newy, pred_y=pred_newy, pred_prob); names(predicted[i, 2]) <-names(x_list)[i]
  }
  res<-list()
  res$perf_summ<-perf_summ
  res$predicted<-predicted
  res
}


RFpred_res<-rf_clf.pred(rf_model_list=rf_clf_res$rf_model_list,
                        x_list=rf_clf_res$x_list,
                        y_list=rf_clf_res$y_list,
                        newx_list=pred_x_list,
                        newy_list=pred_y_list, positive_class="Caries")
# The predicted values and probability of prediction samples merged
pred_list<-lapply(RFpred_res$predicted[, "test_predicted"], as.data.frame)
names(pred_list)<-NULL
#pred_list<-lapply(1:length(pred_list), function(x) data.frame(dataset=rep(names(pred_list)[x], nrow(pred_list[[x]])), pred_list[[x]]))
pred_comb<-do.call(rbind, pred_list)
pred_comb<-data.frame(dataset=rep("prediction", nrow(pred_comb)), pred_comb)
pred_comb<-merge(pred_comb, data_CHRH_list$metadata, by=c(0, 0))
# The predicted values and probability of train samples merged
train_list<-lapply(RFpred_res$predicted[, "train_predicted"], as.data.frame)
names(train_list)<-NULL
train_comb<-do.call(rbind, train_list)
train_comb<-data.frame(dataset=rep("train", nrow(train_comb)), train_comb)
train_comb<-merge(train_comb, data_list$metadata, by=c(0, 0))
all_comb<-rbind(train_comb, pred_comb)
all_comb$y=factor(all_comb$y, levels=c("ConfidentH", "C_H", "RelativeH", "Caries"), ordered = TRUE)
sink(paste(outpath,"RF_pred_summ.xls",sep=""));write.table(all_comb,quote=FALSE,sep="\t", row.names = F);sink()


# Plotting
# boxplot of CC probability for all status
p <- ggplot(all_comb, aes(x=y, y=Caries)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=y), position=position_jitterdodge(jitter.width= 0.2,dodge.width = 0.8),size=1,alpha=0.4) +
  scale_colour_manual(values = viridis(4), name = "Future_status_of_tooth") +
  # geom_signif(comparisons = list(c("ConfidentH", "RelativeH"), c("ConfidentH", "C_H"),
  #                                c("RelativeH", "C_H")),
  #             map_signif_level = function(p) {if(p < 0.01) {p = "**"} else if(p < 0.05) {p = "*"} else {p = "NS"}},
  #             test = "wilcox.test", textsize = 4, step_increase = 0.1,
  #             test.args = list(exact = FALSE, correct = FALSE, conf.int = TRUE,
  #                              conf.level = 0.95)) + # 添加wilcoxon test结果，并使不同分组的检验间隔0.01
  ylab("MiC")+ xlab("Future_status_of_tooth")+
  ylim(c(0, 1))+
  #geom_hline(yintercept=0.5, linetype="dashed")+
  theme_bw()+
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        strip.text = element_text(size = 15),
        strip.background = element_rect(colour = "white"),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 15),
        legend.position = "none")
p
ggsave(filename=paste(outpath,"Pred_in_",c_category,".boxplot.pdf",sep=""),plot=p, width=3, height=4)

# facted by datasets (i.e. positions)
p<-ggplot(all_comb, aes(x=y, y=Caries)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=y), position=position_jitterdodge(jitter.width= 0.2,dodge.width = 0.8),size=1,alpha=0.4) +
  scale_colour_manual(values = viridis(4), name = "Future_status_of_tooth") +
  ylab("MiC")+ xlab("Future_status_of_tooth")+
  ylim(c(0, 1))+
  geom_hline(yintercept=0.5, linetype="dashed")+
  facet_wrap(~Position2, nrow = 1)+
  theme_bw() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        strip.text = element_text(size = 14),
        strip.background = element_rect(colour = "white"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        legend.position = "right")
p
ggsave(filename=paste(outpath,"Pred_in_",c_category, "_among_", s_category,".facets.boxplot.pdf",sep=""),plot=p, width=10, height=4)

# facted by Timepoint
p<-ggplot(all_comb, aes(x=y, y=Caries)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=y), position=position_jitterdodge(jitter.width= 0.2,dodge.width = 0.8),size=1,alpha=0.4) +
  scale_colour_manual(values = viridis(4), name = "Future_status_of_tooth") +
  ylab("MiC")+ xlab("Future_status_of_tooth")+
  ylim(c(0, 1))+
  #geom_hline(yintercept=0.5, linetype="dashed")+
  facet_wrap(~Timepoint, nrow = 1)+
  theme_bw() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        strip.text = element_text(size = 14),
        strip.background = element_rect(colour = "white"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        legend.position = "right")
p
ggsave(filename=paste(outpath,"Pred_in_",c_category, "_among_Timepoint.facets.boxplot.pdf",sep=""),plot=p, width=10, height=4)

# facted by StatusHostChange
p<-ggplot(all_comb, aes(x=y, y=Caries)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=y), position=position_jitterdodge(jitter.width= 0.2,dodge.width = 0.8),size=1,alpha=0.4) +
  scale_colour_manual(values = viridis(4), name = "Future_status_of_tooth") +
  ylab("MiC")+ xlab("Future_status_of_tooth")+
  ylim(c(0, 1))+
  geom_hline(yintercept=0.5, linetype="dashed")+
  facet_wrap(~StatusHostChange, nrow = 1)+
  theme_bw() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        strip.text = element_text(size = 14),
        strip.background = element_rect(colour = "white"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        legend.position = "right")
p
ggsave(filename=paste(outpath,"Pred_in_",c_category, "_among_StatusHostChange.facets.boxplot.pdf",sep=""),plot=p, width=10, height=4)

# facted by StatusToothChange
p<-ggplot(all_comb, aes(x=y, y=Caries)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=y), position=position_jitterdodge(jitter.width= 0.2,dodge.width = 0.8),size=1,alpha=0.4) +
  scale_colour_manual(values = viridis(4), name = "Future_status_of_tooth") +
  ylab("MiC")+ xlab("Future_status_of_tooth")+
  ylim(c(0, 1))+
  geom_hline(yintercept=0.5, linetype="dashed")+
  facet_wrap(~StatusToothChange, nrow = 1)+
  theme_bw() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        strip.text = element_text(size = 14),
        strip.background = element_rect(colour = "white"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        legend.position = "right")
p
ggsave(filename=paste(outpath,"Pred_in_",c_category, "_among_StatusToothChange.facets.boxplot.pdf",sep=""),plot=p, width=10, height=4)

# facted by StatusHostChange + Timepoint
p<-ggplot(all_comb, aes(x=y, y=Caries)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=y), position=position_jitterdodge(jitter.width= 0.2,dodge.width = 0.8),size=1,alpha=0.4) +
  scale_colour_manual(values = viridis(4), name = "Future_status_of_tooth") +
  ylab("MiC")+ xlab("Future_status_of_tooth")+
  ylim(c(0, 1))+
  geom_hline(yintercept=0.5, linetype="dashed")+
  facet_grid(Timepoint~StatusHostChange)+
  theme_bw() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        strip.text = element_text(size = 14),
        strip.background = element_rect(colour = "white"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        legend.position = "bottom")
p
ggsave(filename=paste(outpath,"Pred_in_",c_category, "_among_StatusHostChange_Timepoint.facets.boxplot.pdf",sep=""),plot=p, width=9, height=8)


# faceted by StatusToothChange and Timepoint
p<-ggplot(all_comb, aes(x=y, y=Caries)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=y), position=position_jitterdodge(jitter.width= 0.2,dodge.width = 0.8),size=1,alpha=0.4) +
  scale_colour_manual(values = viridis(4), name = "Future_status_of_tooth") +
  ylab("MiC")+ xlab("Future_status_of_tooth")+
  ylim(c(0, 1))+
  geom_hline(yintercept=0.5, linetype="dashed")+
  facet_grid(Timepoint~StatusToothChange)+
  theme_bw() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        strip.text = element_text(size = 14),
        strip.background = element_rect(colour = "white"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        legend.position = "bottom")
p
ggsave(filename=paste(outpath,"Pred_in_",c_category, "_among_StatusToothChange_Timepoint.facets.boxplot.pdf",sep=""),plot=p, width=9, height=8)
