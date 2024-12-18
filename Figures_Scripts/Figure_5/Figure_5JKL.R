## install and load necessary libraries for data analyses
#-------------------------------
p <- c("reshape2","ggplot2","pheatmap","combinat","randomForest", "gridExtra", "crossRanger", "RColorBrewer",
       "vegan", "biomformat", "doMC", "cowplot", "pROC", "ggsignif", "ggpubr", "colorspace", "viridis")
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
datafile <- "../../Data/1867/function/1867_functional_feature_abd_table.biom"
sample_metadata <- "../../Data/1867/1867_metadata.txt"
feature_metadata <- "../../Data/1867/function/1867_functional_annotation.txt"
prefix_name<-"Function"
s_category<-"Position2"
c_category<-"Future_Status_Tooth"
AddTaxonomy=TRUE
p_cutoff=0.05
q_cutoff=0.05
p.adj.method="BH"
h=10
outpath <- "../../Results/Figure_3/"
if(! dir.exists(outpath)) {
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

res_file<-paste("../../Results/RF_results/sMiC/",prefix_name,"_Position_crossRF_out/", prefix_name, "_rf_clf.by_datasets_res.RData", sep="")
if(file.exists(res_file)){
  load(res_file)
}else{
  rf_clf_res<-rf_clf.by_datasets(data_HHCC_list$data, data_HHCC_list$metadata, s_category, c_category,
                                 positive_class="Caries", clr_transform = TRUE, rf_imp_pvalues = FALSE,
                                 p.adj.method=p.adj.method, q_cutoff=q_cutoff, nfolds=5, verbose=FALSE, ntree=500)
  save(rf_clf_res, file=res_file)
}

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
ggsave(filename=paste(outpath,"Figure_5J_AUROC_cross_appl_matrix_",c_category, "_among_", s_category,".heatmap.pdf",sep=""),plot=p_AUROC, width=6, height=4.6)


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
ggsave(filename=paste(outpath,"Figure_5K_Pred_in_",c_category, "_among_", s_category,".facets.boxplot.pdf",sep=""),plot=p, width=10, height=4)

cal_auc_by_position <- function(x, y) {
  rocobj <- roc(x, y, smooth = F)
  auc <- auc(rocobj)[1]
  return (auc)
}

AUROC <- function(df, label) {
  positions <- unique(df$Position2)
  auc <- matrix(0, nrow = length(positions), ncol = 1)
  auc <- as.data.frame(auc)
  rownames(auc) <- positions
  colnames(auc) <- "AUROC"
  for(pos in positions) {
    data <- df[which(df$Position2 == pos), ]
    x <- data$y
    y <- data$Caries
    # rocobj <- roc(x, y, smooth = F)
    # auc[pos, 1] <- auc(rocobj)[1])
    auc[pos, 1] <- cal_auc_by_position(x, y)
  }
  
  auc$position <- rownames(auc)
  auc$feature <- paste("AUROC", label, sep = "\n")
  p_AUROC<-ggplot(auc, aes(x=as.factor(position), y=as.factor(feature))) +
    xlab("Tooth Position") +
    geom_tile(aes(fill = AUROC, width=0.9, height=0.9), size=6) +
    labs(names = "AUROC") +
    scale_color_manual(values=c("white","grey80"))+
    geom_text(aes(label = round(AUROC, 2)), size = 10, color = "white") +
    scale_fill_viridis(limit = c(0.5, 1))+
    theme_bw() + theme_classic() +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          strip.text = element_text(size = 14),
          strip.background = element_rect(colour = "white"),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          plot.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.position = "right")
  return (p_AUROC)
}

df = read.table(paste0(outpath, "/RF_pred_summ.xls"), sep = "\t", header = T, row.names = 1)
HHRH_df <- subset(df, y %in% c("ConfidentH", "RelativeH"))
HHRH_df$y = factor(HHRH_df$y, levels = c("ConfidentH", "RelativeH"), ordered = T)
plot <- AUROC(HHRH_df, "(ConfidentH\nVS\nRelativeH)")
ggsave(filename=paste0(outpath, "Figure_5L_HHRH_AUROC_by_position.heatmap.pdf"),plot=plot, width=12, height=2)

