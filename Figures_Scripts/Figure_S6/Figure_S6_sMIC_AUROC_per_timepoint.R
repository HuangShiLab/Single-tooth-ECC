# ---
# title: "ECC single-tooth prediction"
# author: "ShiHuang"
# date: "4/17/2019"
# output: html_document
# ---
#-------------------------------
## install and load necessary libraries for data analyses
#-------------------------------
p <- c("reshape2","ggplot2","pheatmap","combinat","randomForest", "gridExtra", "crossRanger", "RColorBrewer",
       "vegan", "biomformat", "doMC", "cowplot", "pROC", "ggsignif", "aplot", "ggpubr", "viridis", "dplyr")
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
source("../../../data_trimming_util.R")
#-------------------------------
datafile <- "../../../../data/1867/taxonomy/ML_features_table/1867_taxonomy_feature_table.biom"
sample_metadata <- "../../../../data/1867/taxonomy/1867_taxonomy_metadata.tsv"
feature_metadata <- "../../../../data/1867/taxonomy/ML_features_table/1867_taxonomy_feature_taxon.txt"
prefix_name<-"taxonomy"
s_category<-"Position2" 
c_category<-"Status_ToothPre"
AddTaxonomy=TRUE
p_cutoff=0.05
q_cutoff=0.05
p.adj.method="BH"
h=10
outpath <- paste(prefix_name,"_Position_crossRF_out/", sep="")
dir.create(outpath)

cal_auc_by_position <- function(x, y) {
  rocobj <- roc(x, y, smooth = F)
  auc <- auc(rocobj)[1]
  return (auc)
}

# AUROC计算
# AUROC计算
AUROC <- function(df, outfile, label) {
  rocobj <- roc(df$y, df$Caries, smooth = F)       # 曲线是否光滑，当光滑时，无法计算置信区间
  
  # 计算AUROC值
  auc<-auc(rocobj)[1]
  # AUROC的置信区间
  auc_low<-ci(rocobj,of="auc")[1]
  auc_high<-ci(rocobj,of="auc")[3]
  print(auc)
  print(auc_low)
  print(auc_high)
  
  # 计算置信区间
  ciobj <- ci.se(rocobj,specificities=seq(0, 1, 0.01))
  data_ci<-ciobj[1:101,1:3]
  data_ci<-as.data.frame(data_ci)
  x=as.numeric(rownames(data_ci))
  data_ci<-data.frame(x,data_ci)
  
  timepoints <- unique(df$Timepoint)
  for(time in timepoints) {
    data <- df[which(df$Timepoint == time), ]
    x <- data$y
    y <- data$Caries
    if(length(unique(x)) >= 2) {
      rocobj <- roc(x, y, smooth = F)
      # 计算置信区间
      ciobj <- ci.se(rocobj,specificities=seq(0, 1, 0.01))
      data_ci<-ciobj[1:101,1:3]
      data_ci<-as.data.frame(data_ci)
      x=as.numeric(rownames(data_ci))
      data_ci<-data.frame(x,data_ci)
      
      plot <- ggroc(rocobj,
                    color="red",
                    size=1,
                    legacy.axes = F # FALSE时 横坐标为1-0 specificity；TRUE时 横坐标为0-1 1-specificity
      )+
        theme_classic()+
        geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1),        # 绘制对角线
                     colour='grey', 
                     linetype = 'dotdash'
        ) +
        geom_ribbon(data = data_ci,                                # 绘制置信区间
                    aes(x=x,ymin=X2.5.,ymax=X97.5.),               # 当legacy.axes=TRUE时， 把x=x改为x=1-x
                    fill = 'lightblue',
                    alpha=0.5) +
        geom_text(aes(x = 0.4, y = 0.2, label = paste0("AUROC = ", round(auc, 2), "\n", label, "\n", time, " timepoint")), size = 5)
      #geom_point(aes(x = cutOffPoint[[2]],y = cutOffPoint[[3]]))+ # 绘制临界点/阈值
      #geom_text(aes(x = cutOffPoint[[2]],y = cutOffPoint[[3]],label=cutOffPointText),vjust=-1) # 添加临界点/阈值文字标签
      
      ggsave(filename=paste0(outfile, "_", time, ".pdf"), plot=plot, width=3.5, height=3.5)
    }
  }
  
  positions <- unique(df$Position2)
  auc <- matrix(0, nrow = length(positions), ncol = length(unique(df$Timepoint)))
  auc <- as.data.frame(auc)
  rownames(auc) <- positions
  colnames(auc) <- unique(df$Timepoint)
  for(pos in positions) {
    data <- df[which(df$Position2 == pos), ]
    timepoints <- unique(data$Timepoint)
    for(time in timepoints) {
      data2 <- data[which(data$Timepoint == time), ]
      x <- data2$y
      y <- data2$Caries
      if(length(unique(x)) < 2) {
        auc[pos, time] <- NA
      } else {
        auc[pos, time] <- cal_auc_by_position(x, y)
      }
    }
  }
  print(auc)
  write.table(auc, paste0(outfile, "_by_position.xls"), sep = "\t", quote = F, row.names = T, col.names = NA)
  
  auc$Position <- rownames(auc)
  auc <- melt(auc)
  colnames(auc) <- c("Position", "Timepoint", "AUROC")
  
  auc$Position <- factor(auc$Position, levels = c("T55", "T54", "T5161", "T64", "T65", "T75", "T74", "T84", "T85"), order = T)
  auc$Timepoint <- factor(auc$Timepoint, levels = c("T1", "T2", "T3", "T4", "T5"), order = T)
  
  auc <- auc %>%
         arrange(Position, Timepoint)

  return (auc)
}

# 读取ROC数据文件
df = read.table(paste0("../", outpath, "RF_pred_summ.xls"), sep = "\t", header = T, row.names = 1)
HHRH_df <- subset(df, y %in% c("ConfidentH", "RelativeH"))
outfile <- paste0(outpath, "/taxonomy_HHRH_RF_pred_summ_ROC")
HHRH_auc <- AUROC(HHRH_df, outfile, "ConfidentH V.S. RelativeH")
plot <- ggplot(HHRH_auc, aes(Position, Timepoint, fill = AUROC)) +
  geom_tile() +
  scale_color_manual(values=c("white","grey80")) +
  scale_fill_viridis(limit = c(0.5, 1))+ 
  geom_text(aes(label = round(AUROC, 2)), color = "white") +
  geom_text(data = HHRH_auc[is.na(HHRH_auc$AUROC), ], aes(label = "NA"), color = "white") +
  #scale_fill_gradient(low="#EFEFFF", high="blue") + 
  xlab("Tooth position") + 
  ylab("Timepoint") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
plot
ggsave(filename=paste0(outfile, "_by_position.heatmap.pdf"),plot=plot, width=5, height=3)

HHCC_df <- subset(df, y %in% c("ConfidentH", "Caries"))
outfile <- paste0(outpath, "/taxonomy_HHCC_RF_pred_summ_ROC")
HHCC_auc <- AUROC(HHCC_df, outfile, "ConfidentH V.S. Caries")
plot <- ggplot(HHCC_auc, aes(Position, Timepoint, fill = AUROC)) +
  geom_tile() +
  scale_color_manual(values=c("white","grey80")) +
  scale_fill_viridis(limit = c(0.5, 1))+ 
  geom_text(aes(label = round(AUROC, 2)), color = "white") +
  #geom_text(data = HHCC_auc[is.na(HHCC_auc$AUROC), ], aes(label = "NA"), color = "white") +
  #scale_fill_gradient(low="#EFEFFF", high="blue") + 
  xlab("Tooth position") + 
  ylab("Timepoint") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
plot
ggsave(filename=paste0(outfile, "_by_position.heatmap.pdf"),plot=plot, width=5, height=3)
