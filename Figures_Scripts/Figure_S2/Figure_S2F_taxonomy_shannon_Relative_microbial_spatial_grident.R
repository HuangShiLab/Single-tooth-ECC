# install and load necessary libraries for data analyses
#-------------------------------
p <- c("ggplot2", "RColorBrewer", "viridis", "palmerpenguins", "tidyverse", "colorspace", "plyr", "ggpubr", "biomformat")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
source("../data_trimming_util.R")

outpath <- "../../Results/Figure_S2/"
if (!dir.exists(outpath)) {
  dir.create(outpath)
} 

relative_microbial_Gradient_boxplot <- function(df, metadata) {
  df<-df[order(rownames(df)), ]
  metadata<-metadata[order(rownames(metadata)),]
  identical(rownames(df), rownames(metadata))
  
  metadata <- subset(metadata, UL == "Upper_teeth")
  df <- subset(df, rownames(df) %in% rownames(metadata))
  identical(rownames(df), rownames(metadata))
  
  diversity_index <- diversity(df, index = "shannon")
  df <- data.frame(shannon = diversity_index)

  df <- merge(metadata[, c("Tooth_num2", "Position2", "NicheICM", "HostGroup", "UL")], df, by="row.names", all=FALSE)
  if("T5161" %in% df$Position2) {
    x <- subset(df, Position2 == "T5161")
    x$Position2 <- "T51"
    df <- rbind(df, x)
    x$Position2 <- "T61"
    df <- rbind(df, x)
    df <- subset(df, Position2 != "T5161")
  }
  
  df <- df[order(df$Tooth_num2), ]
  table(df$Position2)
  df$Position2 <- factor(df$Position2, levels=unique(df$Position2), order = T)

  
  H2H_df <- subset(df, HostGroup == "H2H")
  # 将因子水平手动映射到数值
  H2H_df$Position_number <- as.numeric(H2H_df$Position2)
  
  # loess拟合
  fit <- loess(shannon ~ Position_number, data = H2H_df)
  
  # 获取拟合曲线的数据
  x_vals <- unique(H2H_df$Position_number)
  fit_data <- data.frame(Position_number = x_vals)
  fit_data$predicted_shannon <- predict(fit, newdata = fit_data)
  
  # 将拟合结果与原始x值匹配
  fit_data$position <- factor(fit_data$Position_number, labels = levels(H2H_df$Position2))
  
  # 显示拟合曲线的y值
  # print(fit_data)
  
  match_indices <- match(df$Position2, fit_data$position)
  df$predicted_shannon <- fit_data[match_indices, "predicted_shannon"]
  df$relative_microbial_Gradient <- with(df, shannon - predicted_shannon)
  
  match_indices <- match(df$Row.names, rownames(metadata))
  df$Status_Tooth <- metadata[match_indices, "Status_Tooth"]

  print(head(df))

  plot1 <- ggplot(df, aes(x = Status_Tooth, y = relative_microbial_Gradient)) +
    ylim(min(df$relative_microbial_Gradient) * 1.1, max(df$relative_microbial_Gradient) * 1.1) + 
    geom_boxplot(outlier.shape = NA) +
    scale_color_manual(values = c("C" = "#943735", "H" = "#41663A")) +
    geom_jitter(aes(color=Status_Tooth), shape = 1, width = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    # geom_signif(data = subset(df, HostGroup %in% c("H2C", "C2C")), 
    #                           comparisons = list(c("C", "H")),
    #             map_signif_level = function(p) {if(p < 0.01) {p = "**"} else if(p < 0.05) {p = "*"} else {p = "NS"}},
    #             test = "wilcox.test", textsize = 5, step_increase = 0.1,
    #             test.args = list(exact = FALSE, correct = FALSE, conf.int = TRUE, conf.level = 0.95)) +
    facet_grid(. ~ HostGroup, scales = "free", space = "free") +
    xlab("Tooth Status") + 
    ylab("Relative Microbial\nSpatial Gradient (shannon)") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "gray", fill = NA),
          strip.text = element_text(size = 14),
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.position = "bottom")
  print(plot1)

  if("H2C" %in% df$HostGroup) {
    df$x <- with(df, interaction(HostGroup, Status_Tooth))
    df$x <- factor(df$x, levels = c("C2C.C", "C2C.H", "H2C.C", "H2C.H", "H2H.H"), ordered = T)
    plot2 <- ggplot(df, aes(x = x, y = relative_microbial_Gradient)) +
      ylim(min(df$relative_microbial_Gradient) * 1.1, max(df$relative_microbial_Gradient) * 1.8) + 
      geom_boxplot(outlier.shape = NA) +
      scale_color_manual(values = brewer.pal(8, "Pastel1")) +
      geom_jitter(aes(color=Status_Tooth), shape = 1, width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      geom_signif(comparisons = list(c("C2C.C", "C2C.H"), c("H2C.C", "H2C.H"),
                                     c("H2C.C", "H2H.H"), c("C2C.C", "H2H.H")),
                  map_signif_level = function(p) {if(p < 0.01) {p = "**"} else if(p < 0.05) {p = "*"} else {p = "NS"}},
                  test = "wilcox.test", textsize = 5, step_increase = 0.1,
                  test.args = list(exact = FALSE, correct = FALSE, conf.int = TRUE, conf.level = 0.95)) +
      xlab("Tooth Status") + 
      ylab("Relative Microbial\nSpatial Gradient (shannon)") + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            strip.text = element_text(size = 14),
            legend.title = element_text(size = 14), 
            legend.text = element_text(size = 14),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            legend.position = "bottom")
  }
  else {
    df$x <- with(df, interaction(HostGroup, Status_Tooth))
    df$x <- factor(df$x, levels = c("C2C.C", "C2C.H", "H2H.H"), ordered = T)
    plot2 <- ggplot(df, aes(x = interaction(HostGroup, Status_Tooth), y = relative_microbial_Gradient)) +
      ylim(min(df$relative_microbial_Gradient) * 1.1, max(df$relative_microbial_Gradient) * 1.5) + 
      geom_boxplot(outlier.shape = NA) +
      scale_color_manual(values = brewer.pal(8, "Pastel1")) +
      geom_jitter(aes(color=Status_Tooth), shape = 1, width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      geom_signif(comparisons = list(c("C2C.C", "C2C.H"), c("C2C.C", "H2H.H")),
                  map_signif_level = function(p) {if(p < 0.01) {p = "**"} else if(p < 0.05) {p = "*"} else {p = "NS"}},
                  test = "wilcox.test", textsize = 5, step_increase = 0.1,
                  test.args = list(exact = FALSE, correct = FALSE, conf.int = TRUE, conf.level = 0.95)) +
      xlab("Tooth Status") + 
      ylab("Relative Microbial\nSpatial Gradient (shannon)") + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            strip.text = element_text(size = 14),
            legend.title = element_text(size = 14), 
            legend.text = element_text(size = 14),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            legend.position = "bottom") 
  }
  print(plot2)
  
  if("H2C" %in% metadata$HostGroup) {
    plot3 <- ggplot(df, aes(x = HostGroup, y = relative_microbial_Gradient)) +
      ylim(min(df$relative_microbial_Gradient) * 1.1, max(df$relative_microbial_Gradient) * 1.5) + 
      geom_boxplot(outlier.shape = NA) +
      scale_color_manual(values = c("C" = "#943735", "H" = "#41663A")) +
      geom_jitter(aes(color=Status_Tooth), shape = 1, width = 0.1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      geom_signif(comparisons = list(c("H2C", "H2H"), c("C2C", "H2H")),
                  map_signif_level = function(p) {if(p < 0.01) {p = "**"} else if(p < 0.05) {p = "*"} else {p = "NS"}},
                  test = "wilcox.test", textsize = 5, step_increase = 0.1,
                  test.args = list(exact = FALSE, correct = FALSE, conf.int = TRUE, conf.level = 0.95)) +
      xlab("Tooth Status") + 
      ylab("Relative Microbial\nSpatial Gradient (shannon)") + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(color = "gray", fill = NA),
            strip.text = element_text(size = 14),
            legend.title = element_text(size = 14), 
            legend.text = element_text(size = 14),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            legend.position = "bottom")
  } else {
    plot3 <- ggplot(df, aes(x = HostGroup, y = relative_microbial_Gradient)) +
      ylim(min(df$relative_microbial_Gradient) * 1.1, max(df$relative_microbial_Gradient) * 1.5) + 
      geom_boxplot(outlier.shape = NA) +
        scale_color_manual(values = c("C" = "#943735", "H" = "#41663A")) +
      geom_jitter(aes(color=Status_Tooth), shape = 1, width = 0.1) +
      geom_signif(comparisons = list(c("C2C", "H2H")),
                  map_signif_level = function(p) {if(p < 0.01) {p = "**"} else if(p < 0.05) {p = "*"} else {p = "NS"}},
                  test = "wilcox.test", textsize = 5, step_increase = 0.1,
                  test.args = list(exact = FALSE, correct = FALSE, conf.int = TRUE, conf.level = 0.95)) +
      xlab("Tooth Status") + 
      ylab("Relative Microbial\nSpatial Gradient (shannon)") + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(color = "gray", fill = NA),
            strip.text = element_text(size = 14),
            legend.title = element_text(size = 14), 
            legend.text = element_text(size = 14),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            legend.position = "bottom")
  }
  
  print(plot3)
  
  return(list(plot1 = plot1, plot2 = plot2, plot3 = plot3))
}

metadata_file="../../data/637/taxonomy/637_taxonomy_metadata.tsv"
options(scipen = 4)
data_file="../../data/637/taxonomy/silva_taxonomy_collapse/level8/feature-table.biom"
df <- read_biom(data_file)
df <- as.data.frame(t(as.matrix(biom_data(df))))
df <- df / rowSums(df)
#df<-read.table(data_file, header=T, row.names=1, quote="", comment.char = "")
metadata<-read.table(metadata_file,header=T,sep="\t",row.names=1, quote="", comment.char="")
metadata<-metadata[which(row.names(metadata) %in% row.names(df)), ]
r <- relative_microbial_Gradient_boxplot(df, metadata)
plot1 <- r$plot1 
ggsave(filename=paste(outpath, "/FigS2_637_taxonomy_shonnon_relative_microbial_Gradient.pdf", sep=""), plot=plot1, width=4, height=4)
plot2 <- r$plot2
ggsave(filename=paste(outpath, "/FigS2_637_taxonomy_shonnon_relative_microbial_Gradient2.pdf", sep=""), plot=plot2, width=4, height=4)
plot3 <- r$plot3
ggsave(filename=paste(outpath, "/FigS2_637_taxonomy_shonnon_relative_microbial_Gradient_HostGroup.pdf", sep=""), plot=plot3, width=4, height=4)


metadata_file="../../data/1867/taxonomy/1867_taxonomy_metadata.tsv"
options(scipen = 4)
data_file="../../data/1867/taxonomy/silva_taxonomy_collapse/level8/feature-table.biom"
df <- read_biom(data_file)
df <- as.data.frame(t(as.matrix(biom_data(df))))
df <- df / rowSums(df)
#df<-read.table(data_file, header=T, row.names=1, quote="", comment.char = "")
metadata<-read.table(metadata_file,header=T,sep="\t",row.names=1, quote="", comment.char="")
metadata<-metadata[which(row.names(metadata) %in% row.names(df)), ]
df <- df[order(rownames(df)), ]
metadata <- metadata[order(rownames(metadata)), ]
#metadata[which(metadata$HostGroup %in% c("H2C", "C2C")), "HostGroup"] = "H2C / C2C"
r <- relative_microbial_Gradient_boxplot(df, metadata)
plot1 <- r$plot1
ggsave(filename=paste(outpath, "/FigS2_1867_taxonomy_shonnon_relative_microbial_Gradient.pdf", sep=""), plot=plot1, width=5, height=4)
plot2 <- r$plot2
ggsave(filename=paste(outpath, "/FigS2_1867_taxonomy_shonnon_relative_microbial_Gradient2.pdf", sep=""), plot=plot2, width=7, height=4)
plot3 <- r$plot3
ggsave(filename=paste(outpath, "/FigS2_1867_taxonomy_shonnon_relative_microbial_Gradient_HostGroup.pdf", sep=""), plot=plot3, width=4, height=4)


