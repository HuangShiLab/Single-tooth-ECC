## install and load necessary libraries for data analyses
#-------------------------------
p <- c("reshape2","ggplot2","pheatmap","combinat","randomForest", "gridExtra", "crossRanger", "RColorBrewer",
       "vegan", "biomformat", "doMC", "cowplot", "pROC", "ggsignif", "ggpubr", "plyr", "viridis")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

shannon_diversity_position_group_boxplot <- function(abd_table, metadata, group, mapping = data.frame(), outpath = ".") {
  diversity_index <- diversity(abd_table, index = "shannon")
  diversity_index <- data.frame(shannon = diversity_index)
  write.table(diversity_index, paste0(outpath, "/function_", group, "_diversity.txt"), sep = "\t", quote = F, row.names = T, col.names = F)
  
  identical(rownames(diversity_index), rownames(metadata))
  
  diversity_index <- merge(metadata[, c("Tooth_num2", "Position2", "NicheICM", "HostGroup", "UL")], diversity_index, by="row.names", all=FALSE)
  if("T5161" %in% diversity_index$Position2) {
    x <- subset(diversity_index, Position2 == "T5161")
    x$Position2 <- "T51"
    diversity_index <- rbind(diversity_index, x)
    x$Position2 <- "T61"
    diversity_index <- rbind(diversity_index, x)
    diversity_index <- subset(diversity_index, Position2 != "T5161")
    diversity_index <- rbind.fill(diversity_index, mapping)
  }
  
  diversity_index <- diversity_index[order(diversity_index$Tooth_num2), ]
  diversity_index$Position2 <- factor(diversity_index$Position2, levels=unique(diversity_index$Position2), order = T)
  
  plot <- ggplot(diversity_index, aes(x = Position2, y = shannon)) +
    geom_boxplot(outlier.shape = NA) +
    scale_colour_manual(values = brewer.pal(8, "Pastel2")) +
    geom_jitter(aes(color=NicheICM), shape = 1, width = 0.2) +
    geom_smooth(mapping = aes(x = Position2, y = shannon, group = 1), data = subset(diversity_index, UL == "Upper_teeth"), method = "loess", se=TRUE) + 
    facet_wrap(~HostGroup, scales = "free_x") +
    xlab("Position") + 
    ylab("Shannon diversity") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
  
  plot
  
  ggsave(filename=paste0(outpath, "/function_", group, "_position_HostGroup.pdf"), plot=plot, width=8, height=3)
}

metadata <- read.table("../../data/637/function/637_function_metadata.tsv", sep = "\t", header = T, row.names = 1)
mapping <- unique(metadata[, c("Tooth_num2", "Position2", "NicheICM", "UL")])
mapping$HostGroup <- "H2H"
x <- mapping
x$HostGroup <- "H2C / C2C"
mapping <- rbind(mapping, x)
rownames(mapping) <- 1:nrow(mapping)

outpath <- "../../Results/Figure_S2/"

abd_table <- read_hdf5_biom("../../data/1867/function/KO/1867_function_kos_abd.biom")
abd_table <- biom(abd_table)
abd_table <- data.frame(t(as.matrix(biom_data(abd_table))), check.names = FALSE)
abd_table <- abd_table[order(rownames(abd_table)), ]
metadata <- read.table("../../data/1867/function/1867_function_metadata.tsv", sep = "\t", header = T, row.names = 1)
metadata[which(metadata$HostGroup %in% c("H2C", "C2C")), "HostGroup"] = "H2C / C2C"
identical(rownames(abd_table), rownames(metadata))
shannon_diversity_position_group_boxplot(abd_table, metadata, "1867", mapping, outpath = outpath)


abd_table <- read_hdf5_biom("../../data/637/function/KO/637_function_kos_abd.biom")
abd_table <- biom(abd_table)
abd_table <- data.frame(t(as.matrix(biom_data(abd_table))), check.names = FALSE)
abd_table <- abd_table[order(rownames(abd_table)), ]
metadata <- read.table("../../data/637/function/637_function_metadata.tsv", sep = "\t", header = T, row.names = 1)
identical(rownames(abd_table), rownames(metadata))
shannon_diversity_position_group_boxplot(abd_table, metadata, "637", mapping, outpath = outpath)

