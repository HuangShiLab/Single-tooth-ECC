## install and load necessary libraries for data analyses
#-------------------------------
p <- c("reshape2","ggplot2","pheatmap","combinat","randomForest", "gridExtra", "crossRanger", "RColorBrewer",
       "vegan", "biomformat", "doMC", "cowplot", "pROC")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
#-------------------------------
# input args
#-------------------------------
source("../data_trimming_util.R")
#-------------------------------
sample_metadata <-"../../data/1867/taxonomy/1867_taxonomy_metadata.tsv"
s_category<-"Position2" #c("Group", "Timepoint")
c_category<-"Future_Status_Tooth"
outpath <- "../../Results/Figure_S1/"
if (!dir.exists(outpath)) {
  dir.create(outpath)
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
#colnames(metadata)
#with(metadata, table(HostGroup, HostID))
metadata$Future_Status_Tooth <- factor(metadata$Future_Status_Tooth, 
                                   levels = c("ConfidentH", "C_H", "RelativeH", "Caries"),
                                   order = T)
p <- ggplot(metadata, aes(x=Timepoint, y=HostID)) + 
     geom_point(aes(size=dmfs_Tooth, color=Future_Status_Tooth)) + 
     scale_color_manual(values = brewer.pal(4, "Pastel2")[c(1, 3, 2, 4)]) +
     scale_size(range = c(1, 6), breaks = c(0, 2, 4, 6, 8, 10)) +
     facet_grid(HostGroup~Position2, scales = "free", space = "free") +
     # labs(color = "Future status (Tooth)") +
     theme_bw() + 
     theme(panel.border = element_rect(color = "white", fill = NA),
           strip.background.x=element_rect(color = NA,  fill="grey90"), 
           strip.background.y=element_rect(color = NA,  fill="grey90"))
p
ggsave(filename=paste(outpath,"FigS1_1867_sample_collection.pdf",sep="/"),plot=p, width=12, height=10)


#-------------------------------
# Sample Metadata input in the 637 cohort
#-------------------------------
vld_sample_metadata <- "../../data/637/taxonomy/637_taxonomy_metadata.tsv"
vld_allmetadata<-read.table(vld_sample_metadata,header=T,sep="\t",row.names=1, quote="", comment.char="")
if(length(vld_allmetadata)==1){vld_metadata<-data.frame(vld_allmetadata[order(rownames(vld_allmetadata)),])
all_vld_group<-colnames(vld_metadata)<-colnames(vld_allmetadata)
}else{
  vld_metadata<-vld_allmetadata[order(rownames(vld_allmetadata)),]
}

#colnames(vld_metadata)
#with(vld_metadata, table(HostGroup, HostID))
vld_metadata$Future_Status_Tooth <- factor(vld_metadata$Future_Status_Tooth, 
                                   levels = c("ConfidentH", "C_H", "Caries"),
                                   order = T)
p <- ggplot(vld_metadata, aes(x=Timepoint, y=HostID)) + 
  geom_point(aes(size=dmfs_Tooth, color=Future_Status_Tooth)) + 
  scale_color_manual(values = brewer.pal(4, "Pastel2")[c(1, 3, 4)]) +
  scale_size(range = c(1, 2), breaks = c(0, 1))+
  facet_grid(HostGroup~Position2, scales = "free", space = "free") + 
  # labs(color = "Future status (Tooth)") +
  theme_bw() + 
  theme(panel.border = element_rect(color = "white", fill = NA), 
        strip.background.x=element_rect(color = NA,  fill="grey90"), 
        strip.background.y=element_rect(color = NA,  fill="grey90"))
p
ggsave(filename=paste(outpath,"FigS1_637_sample_collection.pdf",sep="/"),plot=p, width=12, height=4)




