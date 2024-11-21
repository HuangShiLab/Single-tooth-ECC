p <- c("reshape2","ggplot2","pheatmap","combinat","randomForest", "gridExtra", "crossRanger", "RColorBrewer",
       "vegan", "biomformat", "doMC", "cowplot", "pROC", "ggsignif", "aplot", "ggpubr", "viridis")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

source("../data_trimming_util.R")
outdir <- "../../Results/Figure_5/"

surrounding_teeth = list(
  "T55" = c("T54", "T85", "T84"),
  "T54" = c("T55", "T53", "T85", "T84", "T83"),
  "T53" = c("T54", "T52", "T84", "T83", "T82"),
  "T52" = c("T53", "T51", "T83", "T82", "T81"),
  "T51" = c("T52", "T61", "T82", "T81", "T71"),
  "T61" = c("T51", "T62", "T81", "T71", "T72"),
  "T5161" = c("T82", "T81", "T71", "72"),
  "T62" = c("T61", "T63", "T71", "T72", "T73"),
  "T63" = c("T62", "T64", "T72", "T73", "T74"),
  "T64" = c("T63", "T65", "T73", "T74", "T75"),
  "T65" = c("T64", "T75", "T74"),
  "T85" = c("T54", "T55", "T84"),
  "T84" = c("T55", "T53", "T85", "T54", "T83"),
  "T83" = c("T54", "T52", "T84", "T53", "T82"),
  "T82" = c("T53", "T51", "T83", "T52", "T81"),
  "T81" = c("T52", "T61", "T82", "T51", "T71"),
  "T71" = c("T51", "T62", "T81", "T61", "T72"),
  "T72" = c("T61", "T63", "T71", "T62", "T73"),
  "T73" = c("T62", "T64", "T72", "T63", "T74"),
  "T74" = c("T63", "T65", "T73", "T64", "T75"),
  "T75" = c("T64", "T65", "T74")
)

metadata_file="../../data/1867/taxonomy/1867_taxonomy_metadata.tsv"
metadata<-read.table(metadata_file,header=T,sep="\t",row.names=1, quote="", comment.char="")

dist <- dist(metadata[, c("X", "Y", "Z")])
dist <- as.matrix(dist)

all_samples <- rownames(metadata)
healthy_teeth_positions <- unique(subset(metadata, HostGroup == "H2H")$Position2)
clinical_new_features <- c("sum_s_dt", "sum_ns_dt", "sum_dmfs", "spatial_dist_weighted_mean_dmfs")

clinical_features <- as.data.frame(matrix(0, nrow=length(all_samples), ncol=length(clinical_new_features), byrow = TRUE, 
                                          dimnames = list(all_samples, clinical_new_features)))

for (sample in all_samples) {
  t <- metadata[sample, "Timepoint"]
  h <- metadata[sample, "HostID"]
  p <- metadata[sample, "Position2"]
  clinical_features[sample, "sum_s_dt"] <- sum(with(metadata, Timepoint == t & 
                                                      HostID == h & 
                                                      Position2 %in% surrounding_teeth[[p]] & 
                                                      Status_Tooth == "C"))
  clinical_features[sample, "sum_ns_dt"] <- sum(with(metadata, Timepoint == t & 
                                                       HostID == h & 
                                                       ! (Position2 %in% c(p, surrounding_teeth[[p]])) & 
                                                       Status_Tooth == "C"))
}

identical(all_samples, rownames(metadata))
for (s1 in all_samples) {
  samples <- all_samples[which(metadata$Timepoint == metadata[s1, "Timepoint"] & metadata$HostID == metadata[s1, "HostID"])]
  clinical_features[s1, "sum_dmfs"] = 0
  clinical_features[s1, "spatial_dist_weighted_mean_dmfs"] = 0
  sum_of_weights = 0
  for (s2 in samples) {
    if(s1 == s2) next
    clinical_features[s1, "sum_dmfs"] <- clinical_features[s1, "sum_dmfs"] + metadata[s2, "dmfs_Tooth"]
    clinical_features[s1, "spatial_dist_weighted_mean_dmfs"] <- 
      clinical_features[s1, "spatial_dist_weighted_mean_dmfs"] + metadata[s2, "dmfs_Tooth"] / dist[s1, s2]
    sum_of_weights <- sum_of_weights + 1.0 / dist[s1, s2]
  }
  if(sum_of_weights != 0) {
    clinical_features[s1, "spatial_dist_weighted_mean_dmfs"] <- clinical_features[s1, "spatial_dist_weighted_mean_dmfs"] / sum_of_weights
  }
  else {
    clinical_features[s1, "spatial_dist_weighted_mean_dmfs"] <- 0
  }
}

clinical_feature <- t(clinical_features)
write.table(clinical_features, paste0(outdir, "/1867_clinical_features.txt"),
            sep = "\t", quote = F, row.names = T, col.names = NA)

identical(rownames(clinical_features), rownames(metadata))
df <- data.frame(clinical_features, Future_Status_Tooth = metadata$Future_Status_Tooth)
df <- melt(df)
df$Future_Status_Tooth <- factor(df$Future_Status_Tooth, levels = c("ConfidentH", "C_H", "RelativeH", "Caries"), ordered = TRUE)

df1 <- subset(df, variable != "spatial_dist_weighted_mean_dmfs")
df2 <- subset(df, variable == "spatial_dist_weighted_mean_dmfs")

#"#440154FF" "#31688EFF" "#35B779FF" "#FDE725FF"

p <- ggplot(df1, aes(x = Future_Status_Tooth, y = value, color = Future_Status_Tooth)) + 
  geom_boxplot(outlier.shape = NA) + 
  #geom_jitter(aes(color=Future_Status_Tooth), position=position_jitterdodge(jitter.width= 0.2,dodge.width = 0.8),size=1,alpha=0.4) +
  scale_colour_manual(values = viridis(4)) +
  ylab("value")+ xlab("Actual status of tooth")+
  ylim(0, max(df1$value)+6) + 
  geom_signif(comparisons = list(c("ConfidentH", "RelativeH"),
                                 c("ConfidentH", "Caries")),
              map_signif_level = function(p) {if(p < 0.01) {p = "**"} else if(p < 0.05) {p = "*"} else {p = "NS"}},
              test = "wilcox.test", textsize = 5, step_increase = 0.2,
              test.args = list(exact = FALSE, correct = FALSE, conf.int = TRUE, conf.level = 0.95)) +
  facet_wrap(~variable, scales="free_x")+
  #theme_bw()+
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(), 
        panel.grid = element_blank())
p
ggsave(filename=paste0(outdir, "/Fig5A1_1867_clinical_features1.pdf"),plot=p, width=10, height=3)

p <- ggplot(df2, aes(x = Future_Status_Tooth, y = value, color = Future_Status_Tooth)) + 
  geom_boxplot(outlier.shape = NA) + 
  #geom_jitter(aes(color=Future_Status_Tooth), position=position_jitterdodge(jitter.width= 0.2,dodge.width = 0.8),size=1,alpha=0.4) +
  scale_colour_manual(values = viridis(4)) +
  ylab("value")+ xlab("Actual status of tooth")+
  ylim(0, max(df2$value)+2) + 
  geom_signif(comparisons = list(c("ConfidentH", "RelativeH"),
                                 c("ConfidentH", "Caries")),
              map_signif_level = function(p) {if(p < 0.01) {p = "**"} else if(p < 0.05) {p = "*"} else {p = "NS"}},
              test = "wilcox.test", textsize = 5, step_increase = 0.2,
              test.args = list(exact = FALSE, correct = FALSE, conf.int = TRUE, conf.level = 0.95)) +
  facet_wrap(~variable, scales="free_x")+
  #theme_bw()+
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(), 
        panel.grid = element_blank())        
p
ggsave(filename=paste0(outdir, "/Fig5A2_1867_clinical_features2.pdf"),plot=p, width=4, height=3)



# clinical features of functional profiling
metadata_file="../../data/1867/function/1867_function_metadata.tsv"
metadata<-read.table(metadata_file,header=T,sep="\t",row.names=1, quote="", comment.char="")

dist <- dist(metadata[, c("X", "Y", "Z")])
dist <- as.matrix(dist)

all_samples <- rownames(metadata)
healthy_teeth_positions <- unique(subset(metadata, HostGroup == "H2H")$Position2)
clinical_new_features <- c("sum_s_dt", "sum_ns_dt", "sum_dmfs", "spatial_dist_weighted_mean_dmfs")

clinical_features <- as.data.frame(matrix(0, nrow=length(all_samples), ncol=length(clinical_new_features), byrow = TRUE, 
                                          dimnames = list(all_samples, clinical_new_features)))

for (sample in all_samples) {
  t <- metadata[sample, "Timepoint"]
  h <- metadata[sample, "HostID"]
  p <- metadata[sample, "Position2"]
  clinical_features[sample, "sum_s_dt"] <- sum(with(metadata, Timepoint == t & 
                                                      HostID == h & 
                                                      Position2 %in% surrounding_teeth[[p]] & 
                                                      Status_Tooth == "C"))
  clinical_features[sample, "sum_ns_dt"] <- sum(with(metadata, Timepoint == t & 
                                                       HostID == h & 
                                                       ! (Position2 %in% c(p, surrounding_teeth[[p]])) & 
                                                       Status_Tooth == "C"))
}

identical(all_samples, rownames(metadata))
for (s1 in all_samples) {
  samples <- all_samples[which(metadata$Timepoint == metadata[s1, "Timepoint"] & metadata$HostID == metadata[s1, "HostID"])]
  clinical_features[s1, "sum_dmfs"] = 0
  clinical_features[s1, "spatial_dist_weighted_mean_dmfs"] = 0
  sum_of_weights = 0
  for (s2 in samples) {
    if(s1 == s2) next
    clinical_features[s1, "sum_dmfs"] <- clinical_features[s1, "sum_dmfs"] + metadata[s2, "dmfs_Tooth"]
    clinical_features[s1, "spatial_dist_weighted_mean_dmfs"] <- 
      clinical_features[s1, "spatial_dist_weighted_mean_dmfs"] + metadata[s2, "dmfs_Tooth"] / dist[s1, s2]
    sum_of_weights <- sum_of_weights + 1.0 / dist[s1, s2]
  }
  if(sum_of_weights != 0) {
    clinical_features[s1, "spatial_dist_weighted_mean_dmfs"] <- clinical_features[s1, "spatial_dist_weighted_mean_dmfs"] / sum_of_weights
  }
  else {
    clinical_features[s1, "spatial_dist_weighted_mean_dmfs"] <- 0
  }
}

clinical_feature <- t(clinical_features)
write.table(clinical_features, paste0(outdir, "/1867_clinical_features_function.txt"),
            sep = "\t", quote = F, row.names = T, col.names = NA)

identical(rownames(clinical_features), rownames(metadata))

