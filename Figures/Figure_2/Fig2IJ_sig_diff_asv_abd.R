## install and load necessary libraries for data analyses
#-------------------------------
p <- c("reshape2","ggplot2","pheatmap","combinat","randomForest", "gridExtra", "crossRanger", "RColorBrewer",
       "vegan", "biomformat", "doMC", "cowplot", "pROC", "ggsignif", "ggpubr", "dplyr", "viridis")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
#-------------------------------
getLastElement <- function(s) {
  elements <- strsplit(s, ";")[[1]]
  return(tail(elements, 1))
}

outdir <- "../../Results/Figure_3/"

Taxonomy <- read.table("../../data/1867/taxonomy/silva_taxonomy_collapse/level8/taxonomy.tsv", 
                       sep = "\t", header = T, row.names = 1)

Taxonomy$Taxon <- with(Taxonomy, paste(Taxon, ASV_ID, sep = "_"))
Taxonomy$Taxon <- sapply(Taxonomy$Taxon, getLastElement)
Taxonomy$Taxon <- trimws(Taxonomy$Taxon)

meta_637 <- read.table("../../data/637/taxonomy/637_taxonomy_metadata.tsv", sep = "\t", header = T, row.names = 1)
mapping <- unique(meta_637[, c("Tooth_num2", "Position2")])
colnames(mapping) <- c("Tooth_num2", "position")
rownames(mapping) <- mapping$position

sig_diff_asv_1867 <- read.table("../../data/1867/taxonomy/sig_differental_results/level_8/1867_H2H_clr_sig_feature.txt",
                                sep = "\t", header = T, row.names = 1)

asv_abd_637 <- read_hdf5_biom("../../data/637/taxonomy/silva_taxonomy_collapse/level8/feature-table.biom")
asv_abd_637 <- biom(asv_abd_637)
asv_abd_637 <- data.frame(t(as.matrix(biom_data(asv_abd_637))), check.names = FALSE)
asv_abd_637 <- asv_abd_637 / rowSums(asv_abd_637)

sig_diff_asv_table_637 <- asv_abd_637[, which(colnames(asv_abd_637) %in% colnames(sig_diff_asv_1867))]
sig_diff_asv_table_637 <- sig_diff_asv_table_637[order(rownames(sig_diff_asv_table_637)), ]
identical(rownames(sig_diff_asv_table_637), rownames(meta_637))
colnames(sig_diff_asv_table_637) <- Taxonomy[colnames(sig_diff_asv_table_637), "Taxon"]
sig_diff_asv_table_637$samples <- rownames(sig_diff_asv_table_637)

sig_diff_asv_table_637 <- melt(sig_diff_asv_table_637)
sig_diff_asv_table_637$Position <- meta_637[sig_diff_asv_table_637$samples, "Position2"]
sig_diff_asv_table_637 <- sig_diff_asv_table_637 %>% 
  group_by(Position, variable) %>% 
  summarise(value = mean(value))

ranked_data <- sig_diff_asv_table_637 %>%
  group_by(variable) %>%
  mutate(rank = rank(value)) %>%
  ungroup()

df_637 <- ranked_data[, c("Position", "variable", "rank")]
colnames(df_637) <- c("Position", "feature", "Rank")

ordered_features <- subset(df_637, Position == "T51")
ordered_features$feature <- as.character(ordered_features$feature)
ordered_features <- ordered_features %>% arrange(Rank, feature)
ordered_features <- ordered_features$feature

write.table(ordered_features, paste0(outdir, "/ordered_spatial_related_features.txt"),
            sep = "\t", quote = F, row.names = F)

df_637$Tooth_num <- mapping[df_637$Position, "Tooth_num2"]
df_637$feature <- factor(df_637$feature, levels = ordered_features, ordered = T)
df_637 <- df_637 %>% arrange(Tooth_num, feature)
df_637$Position <- factor(df_637$Position, levels = unique(df_637$Position), ordered = T)
  
plot <- ggplot(df_637, aes(Position, feature, fill = Rank)) +
  geom_tile() +
  scale_fill_viridis()+
  xlab("Tooth position") + 
  ylab("Log-scaled Abundance of ASV") + 
  theme_bw() + 
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave(filename=paste0(outdir, "/Fig3C_637_spatial_sig_diff_asv_abd.pdf"), plot=plot, width=10, height=3.5)



meta_1867 <- read.table("../../data/1867/taxonomy/1867_taxonomy_metadata.tsv", sep = "\t", header = T, row.names = 1)
sig_diff_asv_1867 <- read.table("../../data/1867/taxonomy/sig_differental_results/level_8/1867_H2H_clr_sig_feature.txt",
                                sep = "\t", header = T, row.names = 1)

asv_abd_1867 <- read_hdf5_biom("../../data/1867/taxonomy/silva_taxonomy_collapse/level8/feature-table.biom")
asv_abd_1867 <- biom(asv_abd_1867)
asv_abd_1867 <- data.frame(t(as.matrix(biom_data(asv_abd_1867))), check.names = FALSE)
asv_abd_1867 <- asv_abd_1867 / rowSums(asv_abd_1867)

sig_diff_asv_table_1867 <- asv_abd_1867[, which(colnames(asv_abd_1867) %in% colnames(sig_diff_asv_1867))]
sig_diff_asv_table_1867 <- sig_diff_asv_table_1867[order(rownames(sig_diff_asv_table_1867)), ]
identical(rownames(sig_diff_asv_table_1867), rownames(meta_1867))
colnames(sig_diff_asv_table_1867) <- Taxonomy[colnames(sig_diff_asv_table_1867), "Taxon"]
sig_diff_asv_table_1867$samples <- rownames(sig_diff_asv_table_1867)

sig_diff_asv_table_1867 <- melt(sig_diff_asv_table_1867)
sig_diff_asv_table_1867$Position <- meta_1867[sig_diff_asv_table_1867$samples, "Position2"]
sig_diff_asv_table_1867 <- sig_diff_asv_table_1867 %>% 
  group_by(Position, variable) %>% 
  summarise(value = mean(value))

ranked_data <- sig_diff_asv_table_1867 %>%
  group_by(variable) %>%
  mutate(rank = rank(value)) %>%
  ungroup()

df_1867 <- ranked_data[, c("Position", "variable", "rank")]
df_1867 <- dcast(df_1867, variable ~ Position)
df_1867$T51 <- df_1867$T5161
df_1867$T61 <- df_1867$T5161
df_1867 <- df_1867[, which(colnames(df_1867) != "T5161")]

other_positions <- mapping[which(! mapping$position %in% colnames(df_1867)), "position"]
other_positions_table <- as.data.frame(matrix(NA, nrow = nrow(df_1867), ncol = length(other_positions),
                                              dimnames = list(rownames(df_1867), other_positions)))
df_1867 <- cbind(df_1867, other_positions_table)
df_1867 <- melt(df_1867, id.vars="variable")
colnames(df_1867) <- c("feature", "Position", "Rank")

df_1867$Position <- as.character(df_1867$Position)
df_1867$Tooth_num <- mapping[df_1867$Position, "Tooth_num2"]
df_1867$feature <- factor(df_1867$feature, levels = ordered_features, ordered = T)
df_1867 <- df_1867 %>% arrange(Tooth_num, feature)
df_1867$Position <- factor(df_1867$Position, levels = unique(df_1867$Position), ordered = T)

plot <- ggplot(df_1867, aes(Position, feature, fill = Rank)) +
  geom_tile() +
  scale_fill_viridis()+
  xlab("Tooth position") + 
  ylab("ASV") + 
  theme_bw() + 
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave(filename=paste0(outdir, "/Fig3D_1867_spatial_sig_diff_asv_abd.pdf"), plot=plot, width=10, height=3.5)
