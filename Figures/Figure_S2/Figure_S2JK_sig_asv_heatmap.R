rm(list = ls())

library("circlize")
library("tidyverse")
library("ComplexHeatmap")
library("gridBase")
library("reshape2")
library("biomformat")
library("dplyr")
library("plyr")
outpath <- "Fig2"
if (!dir.exists(outpath)) {
  dir.create(outpath)
} 

sig_table_1867 <- read.table("../../data/1867/taxonomy/sig_differental_results/level_8/1867_H2H_clr_sig_feature.txt", sep = "\t", header = T, row.names = 1)
fmetadata_1867 <- read.table("../../data/1867/taxonomy/silva_taxonomy_collapse/level8/taxonomy.tsv", sep = "\t", header = T, row.names = 1, quote = "", comment.char = "")
sig_species <- subset(fmetadata_1867, rownames(fmetadata_1867) %in% colnames(sig_table_1867))$Taxon
#head(sig_table)
meta_637 <- read.table("../../data/637/taxonomy/637_taxonomy_metadata.tsv", sep = "\t", header = T, row.names = 1)
mapping <- unique(meta_637[, c("Tooth_num2", "Position2")])
colnames(mapping) = c("Tooth_num2", "position")

biom_table <- read_biom("../../data/1867/taxonomy/silva_taxonomy_collapse/level8/feature-table.biom")
otu <- biom_data(biom_table)
otu <- data.matrix(otu)
table <- t(otu)
rowsums <- rowSums(table)
table <- table / rowsums
table <- table[order(rownames(table)), colnames(sig_table_1867)]

meta <- read.table("../../data/1867/taxonomy/1867_taxonomy_metadata.tsv", sep = "\t", header = T, row.names = 1)
meta <- meta[order(rownames(meta)), ]
identical(rownames(table), rownames(meta))

asv <- read.table("../../data/1867/taxonomy/silva_taxonomy_collapse/level8/taxonomy.tsv", sep = "\t", header = F, row.names = 1, quote = "")
#head(asv)
asv[1] <- str_split_fixed(unlist(asv[1]), "; ", 6)[, 6]
#print(colnames(table) %in% rownames(asv))

h_table <- subset(table, meta[rownames(table), "HostGroup"] == "H2H")
h_meta <- subset(meta, HostGroup == "H2H")
identical(rownames(h_table), rownames(h_meta))

h_data <- aggregate(h_table, by = list(h_meta[rownames(h_table), "Position2"]), FUN = mean)
rownames(h_data) <- h_data[, 1]
h_data <- h_data[, -1]
h_mat_dataR_1867 <- t(h_data)
h_mat_dataR_1867 <- merge(h_mat_dataR_1867, asv[1], by = 0, all.x = T, all.y = F)
rownames(h_mat_dataR_1867) <- h_mat_dataR_1867[, 1]
h_mat_dataR_1867 <- h_mat_dataR_1867[, -1]
h_mat_dataR_1867 <- aggregate(h_mat_dataR_1867[1:ncol(h_mat_dataR_1867)-1], by = list(h_mat_dataR_1867$V2), FUN = sum)
rownames(h_mat_dataR_1867) <- h_mat_dataR_1867[, 1]
h_mat_dataR_1867 <- h_mat_dataR_1867[, -1]
h_mat_dataR_1867$T51 <- h_mat_dataR_1867$T5161
h_mat_dataR_1867$T61 <- h_mat_dataR_1867$T5161
h_mat_dataR_1867 <- h_mat_dataR_1867[, -which(colnames(h_mat_dataR_1867) == "T5161")]

t_h_mat_dataR_1867 <- t(h_mat_dataR_1867)
t_h_mat_dataR_1867 <- apply(t_h_mat_dataR_1867, 2, rank)
t_h_mat_dataR_1867 <- floor(t_h_mat_dataR_1867)
h_mat_dataR_1867 <- t(t_h_mat_dataR_1867)

meta_637 <- read.table("../../data/637/taxonomy/637_taxonomy_metadata.tsv", sep = "\t", header = T, row.names = 1)
positions <- unique(meta_637$Position2)
other_positions <- positions[which(! positions %in% colnames(h_mat_dataR_1867))]
other_positions_table <- matrix(data = NA, nrow = nrow(h_mat_dataR_1867), ncol = length(other_positions),
                                dimnames = list(rownames(h_mat_dataR_1867), other_positions))
h_mat_dataR_1867 <- cbind(h_mat_dataR_1867, other_positions_table)

data <- melt(h_mat_dataR_1867)
colnames(data) <- c("asv", "position", "rank")

dataR_1867 <- merge(data, mapping, by = "position", all.x = T, all.y = F)

dataR_1867 <- dataR_1867 %>%
  arrange(Tooth_num2, rank)
dataR_1867$position <- factor(dataR_1867$position, levels=unique(dataR_1867$position), order = T)
dataR_1867$asv <- factor(dataR_1867$asv, levels=unique(dataR_1867$asv), order = T)
plot <- ggplot(dataR_1867, aes(as.factor(position), as.factor(asv), fill = rank)) +
  geom_tile() +
  scale_fill_gradient(low = "#EFEFFF", high="blue", na.value = "lightgray", breaks = seq(0, 10, by = 2)) + 
  guides(override.aes = guide_colorbar(list(fill = "lightgray"))) +
  xlab("Tooth position") + 
  ylab("ASV") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))

plot
ggsave(filename=paste(outpath, "/Fig2D_1867_sig_differential_asv_heatmap_2.pdf", sep=""), plot=plot, width=10, height=3)


biom_table <- read_biom("../../data/637/taxonomy/silva_taxonomy_collapse/level8/feature-table.biom")
otu <- biom_data(biom_table)
otu <- data.matrix(otu)
table <- t(otu)
rowsums <- rowSums(table)
table <- table / rowsums
fmetadata_637 <- read.table("../../data/637/taxonomy/silva_taxonomy_collapse/level8/taxonomy.tsv", sep = "\t", header = T, row.names = 1, quote = "", comment.char = "")
colnames(table) <- fmetadata_637[colnames(table), "Taxon"]
table <- melt(table)
sum_data <- table %>%
  group_by(Var1, Var2) %>%
  reframe(sum_value = sum(value))
table <- acast(sum_data, Var1 ~ Var2)
table <- table[order(rownames(table)), unique(sig_species)]

meta <- read.table("../../data/637/taxonomy/637_taxonomy_metadata.tsv", sep = "\t", header = T, row.names = 1)
meta <- meta[order(rownames(meta)), ]
identical(rownames(table), rownames(meta))

colnames(table) <- str_split_fixed(unlist(colnames(table)), "; ", 6)[, 6]
#print(colnames(table) %in% rownames(asv))

h_table <- subset(table, meta[rownames(table), "HostGroup"] == "H2H")
h_meta <- subset(meta, HostGroup == "H2H")
identical(rownames(h_table), rownames(h_meta))

h_data <- aggregate(h_table, by = list(h_meta[rownames(h_table), "Position2"]), FUN = mean)
rownames(h_data) <- h_data[, 1]
h_data <- h_data[, -1]
h_mat_dataR_637 <- t(h_data)

t_h_mat_dataR_637 <- t(h_mat_dataR_637)
t_h_mat_dataR_637 <- apply(t_h_mat_dataR_637, 2, rank)
t_h_mat_dataR_637 <- floor(t_h_mat_dataR_637)
h_mat_dataR_637 <- t(t_h_mat_dataR_637)

data <- melt(h_mat_dataR_637)
colnames(data) <- c("asv", "position", "rank")

dataR_637 <- merge(data, mapping, by = "position", all.x = T, all.y = F)
dataR_637 <- dataR_637 %>%
  arrange(Tooth_num2, -rank)

dataR_637$position <- factor(dataR_637$position, levels=unique(dataR_637$position), order = T)
dataR_637$asv <- factor(dataR_637$asv, levels=unique(dataR_1867$asv), order = T)
#print(dataR_637)
plot <- ggplot(dataR_637, aes(as.factor(position), as.factor(asv), fill = rank)) +
  geom_tile() +
  scale_fill_gradient(low="#EFEFFF", high="blue") + 
  xlab("Tooth position") + 
  ylab("ASV") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))

plot
ggsave(filename=paste(outpath, "/Fig2D_637_sig_differential_asv_heatmap_2.pdf", sep=""), plot=plot, width=10, height=3)


