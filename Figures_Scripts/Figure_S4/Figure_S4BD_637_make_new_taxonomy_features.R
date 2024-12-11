library(ggplot2)
library(ggpubr)
library(viridis)
source("../data_trimming_util.R")

outdir <- "../../Results/Figure_S5/"

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

metadata_file="../../data/637/taxonomy/637_taxonomy_metadata.tsv"
metadata<-read.table(metadata_file,header=T,sep="\t",row.names=1, quote="", comment.char="")

dist <- dist(metadata[, c("X", "Y", "Z")])
dist <- as.matrix(dist)

all_samples <- rownames(metadata)
healthy_teeth_positions <- unique(subset(metadata, HostGroup == "H2H")$Position2)

taxonomy_dist <- read.table("../../data/637/taxonomy/637_all_asv_phylo_rpca/phylo_rpca.distance_matrix/distance-matrix.tsv",
                            sep = "\t", header = T, row.names = 1, check.names = F)
identical(rownames(taxonomy_dist), colnames(taxonomy_dist))

mic_dist_new_features <- c("mean_md", "spatial_dist_weighted_mean_md", paste("mean_md_H", healthy_teeth_positions, sep = "_"))
mic_dist_features <- as.data.frame(matrix(0, nrow=length(all_samples), ncol=length(mic_dist_new_features), byrow = TRUE, 
                                          dimnames = list(all_samples, mic_dist_new_features)))
for(sample in all_samples) {
  other_samples <- all_samples[which(metadata$Timepoint == metadata[sample, "Timepoint"] &
                                       metadata$HostID == metadata[sample, "HostID"] &
                                       all_samples != sample)]
  if(length(other_samples) > 0) {
    mic_dist_features[sample, "mean_md"] <- mean(unlist(taxonomy_dist[sample, other_samples]))
    mic_dist_features[sample, "spatial_dist_weighted_mean_md"] <-
      mean(dist[sample, other_samples] * unlist(taxonomy_dist[sample, other_samples])) 
  }
}

library(stringr)
taxonomy_dist <- subset(taxonomy_dist, metadata[rownames(taxonomy_dist), "HostGroup"] == "H2H")
mean_all_dist <- apply(taxonomy_dist, 2, function(x) {aggregate(x, by=list(type=metadata[rownames(taxonomy_dist), "Position2"]), mean)})
mean_all_dist <- as.data.frame(mean_all_dist)
rownames(mean_all_dist) <- mean_all_dist[, 1]
mean_all_dist <- mean_all_dist[, -seq(1,ncol(mean_all_dist),2)]
colnames(mean_all_dist) <- str_sub(colnames(mean_all_dist), 2, nchar(colnames(mean_all_dist))-2)
mean_all_dist <- t(mean_all_dist)
identical(rownames(mic_dist_features), rownames(mean_all_dist))
mic_dist_features[, paste("mean_md_H", colnames(mean_all_dist), sep = "_")] <- mean_all_dist
write.table(mic_dist_features, paste0(outdir, "/637_all_taxonomy_dist_features.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)

identical(rownames(mic_dist_features), rownames(metadata))
df <- data.frame(mic_dist_features, Future_Status_Tooth = metadata$Future_Status_Tooth)
df <- melt(df)
df$Future_Status_Tooth <- factor(df$Future_Status_Tooth, levels = c("ConfidentH", "C_H", "Caries"), ordered = TRUE)

df1 <- subset(df, variable %in% c("mean_md", "spatial_dist_weighted_mean_md"))
p <- ggplot(df1, aes(x = Future_Status_Tooth, y = value)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=Future_Status_Tooth), position=position_jitterdodge(jitter.width= 0.2,dodge.width = 0.8),size=1,alpha=0.4) +
  scale_colour_manual(values = viridis(4)[c(1, 2, 4)]) +
  ylab("value")+ xlab("Actual status of tooth")+
  ylim(0, max(df1$value)+1) + 
  geom_signif(comparisons = list(c("ConfidentH", "C_H"),
                                 c("ConfidentH", "Caries")),
              map_signif_level = function(p) {if(p < 0.01) {p = "**"} else if(p < 0.05) {p = "*"} else {p = "NS"}},
              test = "wilcox.test", textsize = 3, step_increase = 0.3,
              test.args = list(exact = FALSE, correct = FALSE, conf.int = TRUE, conf.level = 0.95)) +
  facet_wrap(~variable, scales="free_x")+
  theme_bw()+
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
        panel.background = element_blank())

p
ggsave(filename= paste0(outdir, "/FigS5B_637_all_taxonomy_dist_features1.pdf"),plot=p, width=7, height=3)

df2 <- subset(df, ! variable %in% c("mean_md", "spatial_dist_weighted_mean_md"))
p <- ggplot(df2, aes(x = Future_Status_Tooth, y = value)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=Future_Status_Tooth), position=position_jitterdodge(jitter.width= 0.2,dodge.width = 0.8),size=1,alpha=0.4) +
  scale_colour_manual(values = viridis(4)[c(1, 2, 4)]) +
  ylab("value")+ xlab("Actual status of tooth")+
  ylim(0, max(df2$value)+1) + 
  geom_signif(comparisons = list(c("ConfidentH", "C_H"),
                                 c("ConfidentH", "Caries")),
              map_signif_level = function(p) {if(p < 0.01) {p = "**"} else if(p < 0.05) {p = "*"} else {p = "NS"}},
              test = "wilcox.test", textsize = 3, step_increase = 0.3,
              test.args = list(exact = FALSE, correct = FALSE, conf.int = TRUE, conf.level = 0.95)) +
  facet_wrap(~variable, scales="free_x")+
  theme_bw()+
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
        panel.background = element_blank())
p
ggsave(filename= paste0(outdir, "/FigS5D_637_all_taxonomy_dist_features2.pdf"),plot=p, width=10, height=10)

