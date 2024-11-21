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

meta_637 <- read.table("../../data/637/function/637_function_metadata.tsv", sep = "\t", header = T, row.names = 1)
mapping <- unique(meta_637[, c("Tooth_num2", "Position2")])
colnames(mapping) <- c("Tooth_num2", "position")
rownames(mapping) <- mapping$position

sig_diff_ko_table_637 <- read.table("../../data/637/function/sig_differental_results/functions_category.l3.Abd",
                                sep = "\t", header = T, row.names = 1)

sig_diff_ko_table_637 <- sig_diff_ko_table_637[order(rownames(sig_diff_ko_table_637)), ]
identical(rownames(sig_diff_ko_table_637), rownames(meta_637))
sig_diff_ko_table_637$samples <- rownames(sig_diff_ko_table_637)

sig_diff_ko_table_637 <- melt(sig_diff_ko_table_637)
sig_diff_ko_table_637$Position <- meta_637[sig_diff_ko_table_637$samples, "Position2"]
sig_diff_ko_table_637$variable <- as.character(sig_diff_ko_table_637$variable)
sig_diff_ko_table_637$Position <- as.character(sig_diff_ko_table_637$Position)
sig_diff_ko_table_637 <- sig_diff_ko_table_637 %>% 
  dplyr::group_by(Position, variable) %>% 
  dplyr::summarise(value = mean(value))

ranked_data <- sig_diff_ko_table_637 %>%
  dplyr::group_by(variable) %>%
  dplyr::mutate(rank = rank(value)) %>%
  dplyr::ungroup()

df_637 <- ranked_data[, c("Position", "variable", "rank")]
colnames(df_637) <- c("Position", "feature", "Rank")

ordered_features <- subset(df_637, Position == "T51")
ordered_features$feature <- as.character(ordered_features$feature)
ordered_features <- ordered_features %>% arrange(Rank, feature)
ordered_features <- ordered_features$feature

df_637$Tooth_num <- mapping[df_637$Position, "Tooth_num2"]
df_637$feature <- factor(df_637$feature, levels = ordered_features, ordered = T)
df_637 <- df_637 %>% arrange(Tooth_num, feature)
df_637$Position <- factor(df_637$Position, levels = unique(df_637$Position), ordered = T)
  
plot <- ggplot(df_637, aes(Position, feature, fill = Rank)) +
  geom_tile() +
  scale_fill_viridis()+
  xlab("Tooth position") + 
  ylab("Level 3 KEGG pathway") + 
  theme_bw() + 
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave(filename=paste0(outdir, "/Fig3I_637_L3_spatial_sig_diff_ko_abd.pdf"), plot=plot, width=10, height=15)



meta_1867 <- read.table("../../data/1867/function/1867_function_metadata.tsv", sep = "\t", header = T, row.names = 1)
sig_diff_ko_table_1867 <- read.table("../../data/1867/function/sig_differental_results/functions_category.l3.Abd",
                                sep = "\t", header = T, row.names = 1)
sig_diff_ko_table_1867 <- sig_diff_ko_table_1867[order(rownames(sig_diff_ko_table_1867)), ]

identical(rownames(sig_diff_ko_table_1867), rownames(meta_1867))
sig_diff_ko_table_1867$samples <- rownames(sig_diff_ko_table_1867)

sig_diff_ko_table_1867 <- melt(sig_diff_ko_table_1867)
sig_diff_ko_table_1867$Position <- meta_1867[sig_diff_ko_table_1867$samples, "Position2"]
sig_diff_ko_table_1867 <- sig_diff_ko_table_1867 %>% 
  dplyr::group_by(Position, variable) %>% 
  dplyr::summarise(value = mean(value))

ranked_data <- sig_diff_ko_table_1867 %>%
  dplyr::group_by(variable) %>%
  dplyr::mutate(rank = rank(value)) %>%
  dplyr::ungroup()

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
  ylab("Level 3 KEGG pathway") + 
  theme_bw() + 
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave(filename=paste0(outdir, "/Fig3J_1867_L3_spatial_sig_diff_ko_abd.pdf"), plot=plot, width=10, height=15)

