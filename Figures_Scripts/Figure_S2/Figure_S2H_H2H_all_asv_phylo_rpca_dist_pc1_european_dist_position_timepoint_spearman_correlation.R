# ---
# title: "alpha diversity analysis"
# author: "ShiHuang"
# date: "12/29/2019"
# output: html_document
# ---
#-------------------------------
# install and load necessary libraries for data analyses
#-------------------------------
p <- c("ggplot2", "RColorBrewer", "viridis", "palmerpenguins", "tidyverse", "colorspace", "plyr")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
source("../data_trimming_util.R")

outpath <- "Fig2"
if (!dir.exists(outpath)) {
  dir.create(outpath)
} 

heatmap <- function(data, file_name) {
  data$x <- rownames(data)
  dataR <- melt(data)

  plot <- ggplot(dataR, aes(as.factor(variable), as.factor(x), fill = value)) +
    geom_tile() +
    scale_fill_gradient(low="white", high="blue", limits = c(0, 0.8)) + 
    xlab("Tooth position") + 
    ylab("Tooth position") + 
    theme_bw() + 
    geom_text(label = round(dataR$value, 3)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
  
  plot
  ggsave(filename=file_name, plot=plot, width=5, height=4)
  
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y, method = "spearman")) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}

same_host_timepoint_delta_pc1 <- function(df, metadata, group) {
  identical(rownames(df), rownames(metadata))
  euro_dist <- as.matrix(dist(df$PC1))
  rownames(euro_dist) <- rownames(df)
  colnames(euro_dist) <- rownames(df)
  
  df <- melt(euro_dist)
  df[, c("h1", "t1", "p1", "g1")] <- metadata[df$Var1, c("HostID", "Timepoint", "Position2", "HostGroup")]
  df[, c("h2", "t2", "p2", "g2")] <- metadata[df$Var2, c("HostID", "Timepoint", "Position2", "HostGroup")]
  
  df <- subset(df, h1 == h2 & t1 == t2)
  df <- df[order(df$t1), ]

  r <- df |>
    group_by(t1) |>
    nest()
  
  df <- subset(df, h1 == h2 & t1 == t2)
  df <- df[order(df$t1), ]
  
  df$g <- with(df, paste(h1, p1, h2, p2, sep = "_"))
  r <- df |>
    group_by(t1) |>
    nest()
  
  groups <- unique(r$data[[1]]$g)
  for(i in 2:length(r$t1)) {
    groups <- intersect(groups, unique(r$data[[i]]$g))
  }
  
  data <- matrix(NA, nrow = length(groups), ncol = length(r$t1))
  data <- as.data.frame(data)
  rownames(data) <- groups[order(groups)]
  colnames(data) <- r$t1
  for(i in 1:length(r$t1)) {
    r$data[[i]] <- subset(r$data[[i]], g %in% groups)
    r$data[[i]] <- r$data[[i]][order(r$data[[i]]$g), ]
    print(identical(r$data[[i]]$g, rownames(data)))
    data[r$t1[i]] <- r$data[[i]]$value
  }
  
  pdf(paste0(outpath, "/Fig2C_", group, "_H2H_all_asv_phylo_rpca_dist_pc1_european_dist_spearman_correlation_matrix.pdf"))
  pairs(data,
        upper.panel = panel.cor,    # Correlation panel
        lower.panel = panel.smooth) # Smoothed regression lines
  dev.off()
  
  
}



scatter_plot <- function(df, metadata, group) {
  identical(rownames(df), rownames(metadata))
  euro_dist <- as.matrix(dist(df$PC1))
  rownames(euro_dist) <- rownames(df)
  colnames(euro_dist) <- rownames(df)
  
  df <- melt(euro_dist)
  df[, c("h1", "t1", "p1", "g1")] <- metadata[df$Var1, c("HostID", "Timepoint", "Position2", "HostGroup")]
  df[, c("h2", "t2", "p2", "g2")] <- metadata[df$Var2, c("HostID", "Timepoint", "Position2", "HostGroup")]
  
  df <- subset(df, h1 == h2 & t1 == t2)
  df <- df[order(df$t1), ]
  
  r <- df |>
    group_by(t1) |>
    nest()
  
  df$g <- with(df, paste(h1, p1, h2, p2, sep = "_"))
  r <- df |>
    group_by(t1) |>
    nest()
  
  groups <- unique(r$data[[1]]$g)
  for(i in 2:length(r$t1)) {
    groups <- intersect(groups, unique(r$data[[i]]$g))
  }
  
  data <- matrix(NA, nrow = length(groups), ncol = length(r$t1))
  data <- as.data.frame(data)
  rownames(data) <- groups[order(groups)]
  colnames(data) <- r$t1
  for(i in 1:length(r$t1)) {
    r$data[[i]] <- subset(r$data[[i]], g %in% groups)
    r$data[[i]] <- r$data[[i]][order(r$data[[i]]$g), ]
    print(identical(r$data[[i]]$g, rownames(data)))
    data[r$t1[i]] <- r$data[[i]]$value
  }
  
  T1 <- unlist(data[r$t1[[1]]])
  T2 <- unlist(data[r$t1[[2]]])
  cor_res <- cor(T1, T2, method = "spearman")
  pdf(paste0(outpath, "/Fig2B_", group, "_H2H_all_asv_phylo_rpca_dist_pc1_european_dist_spearman_correlation_matrix.pdf"))
  #plot(T1, T2, pch = 19) +
  plot(T1, T2, type = "n") + 
  points(T1, T2, cex = .5, col = rgb(0, 0, 0, 70, maxColorValue=255)) +
  abline(lm(T2 ~ T1), col = "red", lwd = 3) + # Regression line
  text(0.9 * max(T1), 0.9 * max(T2), paste("Correlation:", round(cor_res, 2))) # Pearson correlation
  dev.off()
}
  
  
metadata_file="../../data/637/taxonomy/637_taxonomy_metadata.tsv"
options(scipen = 4)
data_file="../../data/637/taxonomy/637_all_asv_phylo_rpca/phylo_rpca.biplot/PCs.txt"
df<-read.table(data_file, header=T, row.names=1, quote="", comment.char = "")
metadata<-read.table(metadata_file,header=T,sep="\t",row.names=1, quote="", comment.char="")
metadata<-metadata[which(rownames(metadata) %in% rownames(df)), ]
identical(rownames(df), rownames(metadata))
metadata<-subset(metadata, HostGroup == "H2H")
df<-subset(df, rownames(df) %in% rownames(metadata))
identical(rownames(df), rownames(metadata))
scatter_plot(df, metadata, 637)
#same_host_timepoint_delta_pc1(df, metadata, 637)

metadata_file="../../data/1867/taxonomy/1867_taxonomy_metadata.tsv"
options(scipen = 4)
data_file="../../data/1867/taxonomy/1867_all_asv_phylo_rpca/phylo_rpca.biplot/PCs.txt"
df<-read.table(data_file, header=T, row.names=1, quote="", comment.char = "")
df<-df[order(rownames(df)), ]
metadata<-read.table(metadata_file,header=T,sep="\t",row.names=1, quote="", comment.char="")
metadata<-metadata[which(row.names(metadata) %in% row.names(df)), ]
identical(rownames(df), rownames(metadata))
metadata<-subset(metadata, HostGroup == "H2H")
df<-subset(df, rownames(df) %in% rownames(metadata))
identical(rownames(df), rownames(metadata))
same_host_timepoint_delta_pc1(df, metadata, 1867)


