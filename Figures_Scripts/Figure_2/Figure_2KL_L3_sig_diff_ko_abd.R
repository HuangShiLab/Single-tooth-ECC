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

outdir <- "."

meta_637 <- read.table("~/Projects/Teng/github/ECC/data/637/function/637_function_metadata.tsv", sep = "\t", header = T, row.names = 1)
mapping <- unique(meta_637[, c("Tooth_num2", "Position2")])
colnames(mapping) <- c("Tooth_num2", "position")
rownames(mapping) <- mapping$position

ko_abd <- read.table("~/Projects/Teng/github/ECC/data/637/function/KO/637_function_kos_abd.tsv",
                     sep = "\t", header = T, row.names = 1, comment.char = "")

identical(colnames(ko_abd), rownames(meta_637))

ko_richiness_results <- read.table("KOs.txt", sep = "\t", header = T, row.names = 1) 

# 初始化一个空列表来存储KOs
kos_list <- list()

# 遍历每行数据
for (i in 1:nrow(ko_richiness_results)) {
  # 分割KOs列中的字符串，这里假设KO之间是通过'/'分隔的
  kos <- strsplit(ko_richiness_results$KOs[i], "/")[[1]]
  
  # 将分割后的KOs添加到列表中，使用Description作为名称
  kos_list[[ko_richiness_results$Description[i]]] <- kos
}

kos_list <- stack(kos_list)

write.table(kos_list, "KO_ko.txt", sep = "\t", quote = F, row.names = F)

df <- merge(ko_abd, kos_list, by.x = 0, by.y = 1, all.x = F, all.y = T)

df <- df[, -1]
summarized_df <- df %>%
  group_by(ind) %>%
  summarize(across(everything(), sum))

write.table(summarized_df, "637_sig_diff_enriched_kos_abd.txt", sep = "\t", quote = F, row.names = F)

df <- melt(summarized_df)
df <- merge(df, meta_637, by.x = 2, by.y = 0, all = T)

summarized_df <- df %>%
  group_by(ind, Position2) %>%
  summarize(abd = mean(value))

ranked_data <- summarized_df %>%
  dplyr::group_by(ind) %>%
  dplyr::mutate(rank = rank(abd)) %>%
  dplyr::ungroup()

ordered_inds <- subset(ranked_data, Position2 == "T51")
ordered_inds$ind <- as.character(ordered_inds$ind)
ordered_inds <- ordered_inds %>% arrange(rank, ind)
ordered_inds <- ordered_inds$ind

ranked_data$ind <- factor(ranked_data$ind,
                          levels = ordered_inds,
                          ordered = T)
ranked_data$Position2 <- factor(ranked_data$Position2, 
                                levels = c("T55", "T54", "T53", "T52", "T51",
                                           "T61", "T62", "T63", "T64", "T65", 
                                           "T75", "T74", "T73", "T72", "T71",
                                           "T81", "T82", "T83", "T84", "T85"),
                                ordered = T)

plot1 <- ggplot(ranked_data, aes(Position2, ind, fill = rank)) +
  geom_tile() +
  scale_fill_viridis()+
  xlab("Tooth position") + 
  ylab("Enriched pathway") + 
  theme_bw() + 
  theme(strip.text = element_text(size = 15),
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom")

plot
ggsave(filename="Fig3I_637_L3_enriched_spatial_sig_diff_ko_abd.pdf", plot=plot1, width=12, height=10)


meta_1867 <- read.table("~/Projects/Teng/github/ECC/data/1867/function/1867_function_metadata.tsv", sep = "\t", header = T, row.names = 1)
mapping <- unique(meta_637[, c("Tooth_num2", "Position2")])
colnames(mapping) <- c("Tooth_num2", "position")
rownames(mapping) <- mapping$position

ko_abd <- read.table("~/Projects/Teng/github/ECC/data/1867/function/KO/1867_function_kos_abd.tsv",
                     sep = "\t", header = T, row.names = 1, comment.char = "")

identical(colnames(ko_abd), rownames(meta_1867))

df <- merge(ko_abd, kos_list, by.x = 0, by.y = 1, all.x = F, all.y = T)

df <- df[, -1]
summarized_df <- df %>%
  group_by(ind) %>%
  summarize(across(everything(), sum))

write.table(summarized_df, "1867_sig_diff_enriched_kos_abd.txt", sep = "\t", quote = F, row.names = F)

df <- melt(summarized_df)
df <- merge(df, meta_1867, by.x = 2, by.y = 0, all = T)

summarized_df <- df %>%
  group_by(ind, Position2) %>%
  summarize(abd = mean(value))

ranked_data <- summarized_df %>%
  dplyr::group_by(ind) %>%
  dplyr::mutate(rank = rank(abd)) %>%
  dplyr::ungroup()

ranked_data <- dcast(ranked_data, ind ~ Position2)

other_positions <- unique(meta_637[which(! meta_637$Position2 %in% colnames(ranked_data)), "Position2"])
other_positions_table <- as.data.frame(matrix(NA, nrow = nrow(ranked_data), ncol = length(other_positions),
                                              dimnames = list(rownames(ranked_data), other_positions)))
ranked_data <- cbind(ranked_data, other_positions_table)
ranked_data$T51 <- ranked_data$T5161
ranked_data$T61 <- ranked_data$T5161
ranked_data <- ranked_data[, which(colnames(ranked_data) != "T5161")]

ranked_data <- melt(ranked_data, id.vars="ind")
colnames(ranked_data) <- c("ind", "Position2", "rank")

ranked_data$ind <- factor(ranked_data$ind,
                          levels = ordered_inds,
                          ordered = T)
ranked_data$Position2 <- factor(ranked_data$Position2, 
                                levels = c("T55", "T54", "T53", "T52", "T51",
                                           "T61", "T62", "T63", "T64", "T65", 
                                           "T75", "T74", "T73", "T72", "T71",
                                           "T81", "T82", "T83", "T84", "T85"),
                                ordered = T)
# ranked_data$Position2 <- factor(ranked_data$Position2, 
#                                   levels = c("T55", "T54", "T5151", "T64", "T65", "T75", "T74", "T84", "T85"),
#                                   ordered = T)

plot2 <- ggplot(ranked_data, aes(Position2, ind, fill = rank)) +
  geom_tile() +
  scale_fill_viridis()+
  xlab("Tooth position") + 
  ylab("Enriched pathway") + 
  theme_bw() + 
  theme(strip.text = element_text(size = 15),
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom")

plot
ggsave(filename="Fig3J_1867_L3_enriched_spatial_sig_diff_ko_abd.pdf", plot=plot2, width=12, height=10)


library(ggplot2)
library(forcats)
ko_richiness_results <- mutate(ko_richiness_results, 
        richFactor = as.numeric(sub("/\\d+", "", GeneRatio)) / as.numeric(sub("/\\d+", "", BgRatio)))

ko_richiness_results$Description <- 
  factor(ko_richiness_results$Description,
         levels = ordered_inds,
         ordered = T)
ko_richiness_results <- ko_richiness_results[order(ko_richiness_results$Description), ]
ko_richiness_results$order <- 1:nrow(ko_richiness_results)

plot3 <- ggplot(ko_richiness_results, showCategory = 32,
  aes(richFactor, fct_reorder(Description, order))) +
  #aes(richFactor, fct_reorder(Description, richFactor))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        trans = "log10",
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  scale_size_continuous(range=c(2, 10)) +
  xlab("Rich Factor") +
  ylab("Enriched pathway") +
  theme_bw() + 
  theme(strip.text = element_text(size = 15),
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom")

ggsave(filename="Fig3H_KO_L3_enrichness_result.pdf", plot=plot3, width=8, height=10)


plot <- ggarrange(plot_grid(plot1, plot2, plot3, ncol = 3, nrow = 1, align = "hv"),
                  common.legend = F, legend = "bottom")
ggsave(filename="Fig3IJH_IM_sig_diff_KO_enrichness_result.pdf", plot=plot, width=32, height=10)

