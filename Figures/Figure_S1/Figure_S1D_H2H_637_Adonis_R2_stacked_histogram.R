p <- c("reshape2","ggplot2", "dplyr", "viridis", "colorspace", "RColorBrewer", "cowplot", "grid", "gridExtra", "gtable",
       "ggpubr")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

all_group <- c("Timepoint", "Position_Tooth", "HostID", "Age", "NicheICM")

outpath <- "../../../Results/Figure_2/"
stat_summ_taxo_637_all_group <- read.table(paste0(outpath, 
                                        "637_taxonomy_dist_H2H.PerMANOVA_out/All_metadata", 
                                        length(all_group), ".Beta_diversity_summ.xls"),
                                 sep='\t', header = T)

stat_summ_func_637_all_group <- read.table(paste0(outpath, 
                                        "637_function_dist_H2H.PerMANOVA_out/All_metadata", 
                                        length(all_group), ".Beta_diversity_summ.xls"),
                                 sep='\t', header = T)

identical(stat_summ_taxo_637_all_group$feature, stat_summ_func_637_all_group$feature)
adonis_R2_637 <- data.frame(feature = stat_summ_taxo_637_all_group$feature, 
                            Function = stat_summ_func_637_all_group$Adonis.R2,
                            Taxonomy = stat_summ_taxo_637_all_group$Adonis.R2)
summ <- adonis_R2_637$Taxonomy + adonis_R2_637$Function
adonis_R2_637 <- adonis_R2_637[order(summ), ]
ordered_features <- unique(adonis_R2_637$feature)
adonis_R2_637$feature <- factor(adonis_R2_637$feature, 
                                levels=ordered_features, 
                                order = T)
melt_adonis_R2_637 <- melt(adonis_R2_637)

stacked_histogram <- ggplot(melt_adonis_R2_637, aes(feature, value, fill = variable)) +
  scale_fill_manual(values = brewer.pal(8, "Pastel2")) + 
  geom_bar(stat="identity", position="stack", color="lightgray", width = 0.7) +
  coord_flip() +
  scale_x_discrete(position = "bottom") +
  scale_y_continuous(position = "right") +
  labs(x = "Group", y = bquote('Adonis '~R^2), fill = "Distance") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 30, color = "black"),
        axis.title.x = element_text(size = 30, color = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 30, color = "black"),
        legend.text = element_text(size = 30, color = "black"),
        legend.key.size = unit(1, "cm"),
        panel.background = element_rect(fill="white"),
        legend.position = "bottom")
stacked_histogram

stat_summ_taxo_637_block_hostid <- read.table(paste0(outpath, 
                                                     "637_taxonomy_dist_H2H_block_hostid.PerMANOVA_out/All_metadata", 
                                                     length(all_group), ".Beta_diversity_summ.xls"),
                                              sep='\t', header = T)
stat_summ_func_637_block_hostid <- read.table(paste0(outpath, 
                                                     "637_function_dist_H2H_block_hostid.PerMANOVA_out/All_metadata", 
                                                     length(all_group), ".Beta_diversity_summ.xls"),
                                              sep='\t', header = T)

stat_summ_taxo_637_block_hostid_timepoint <- read.table(paste0(outpath, 
                                                     "637_taxonomy_dist_H2H_block_hostid_timepoint.PerMANOVA_out/All_metadata", 
                                                     length(all_group), ".Beta_diversity_summ.xls"),
                                              sep='\t', header = T)
stat_summ_func_637_block_hostid_timepoint <- read.table(paste0(outpath, 
                                                     "637_function_dist_H2H_block_hostid_timepoint.PerMANOVA_out/All_metadata", 
                                                     length(all_group), ".Beta_diversity_summ.xls"),
                                              sep='\t', header = T)

identical(stat_summ_taxo_637_all_group$feature, stat_summ_func_637_all_group$feature)
identical(stat_summ_taxo_637_all_group$feature, stat_summ_taxo_637_block_hostid$feature)
identical(stat_summ_taxo_637_all_group$feature, stat_summ_func_637_block_hostid$feature)
identical(stat_summ_taxo_637_all_group$feature, stat_summ_taxo_637_block_hostid_timepoint$feature)
identical(stat_summ_taxo_637_all_group$feature, stat_summ_func_637_block_hostid_timepoint$feature)

p_value <- data.frame(feature = stat_summ_taxo_637_all_group$feature,
                      Taxonomy1 = stat_summ_taxo_637_all_group$Adonis.P,
                      Function1 = stat_summ_func_637_all_group$Adonis.P,
                      Taxonomy2 = stat_summ_taxo_637_block_hostid$Adonis.P,
                      Function2 = stat_summ_func_637_block_hostid$Adonis.P,
                      Taxonomy3 = stat_summ_taxo_637_block_hostid_timepoint$Adonis.P,
                      Function3 = stat_summ_func_637_block_hostid_timepoint$Adonis.P)

melt_p_value <- melt(p_value)
melt_p_value$feature <- factor(melt_p_value$feature,
                               levels = ordered_features,
                               ordered = T)

melt_p_value[which(melt_p_value$value > 0.05), "value"] <- "NS"
melt_p_value[which(melt_p_value$value < 0.05), "value"] <- ""
melt_p_value[which(melt_p_value$value < 0.01), "value"] <- ""
host_features <- c("HostID", "HostGroup", "Status_Host", "dmfs_Host")
melt_p_value[with(melt_p_value, variable %in% c("Taxonomy2", "Function2") & 
                    feature %in% host_features), "value"] <- "NA"
timepoint_features <- c("Timepoint", "Age")
melt_p_value[with(melt_p_value, variable %in% c("Taxonomy3", "Function3") & 
                    feature %in% c(host_features, timepoint_features)), "value"] <- "NA"

heatmap <- ggplot(melt_p_value, aes(x = variable, y = feature)) + 
  geom_tile(aes(fill = value, width = 0.9, height = 0.9)) +
  scale_fill_manual(values = c("lightgreen", "lightgray", "lightyellow"),
                    labels = c("p < 0.01", "NA", "NS"))+
  scale_x_discrete(labels=c("Without\nblocking\nTaxonomy", 
                            "Without\nblocking\nFunction", 
                            "HostID\nblocked\nTaxonomy", 
                            "HostID\nblocked\nFunction", 
                            "HostID -\nTimepoint\nblocked\nTaxonomy",
                            "HostID -\nTimepoint\nblocked\nFunction"),
                   position = "top") +
  geom_text(aes(label = value), size = 10, color = "black") +
  theme(axis.text.x=element_text(angle=90, hjust=1), panel.border = element_blank())+
  labs(fill = "Adonis.P") +
  ylab("Group")+
  xlab("Adonis p-value")+
  theme_bw() + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"),
        axis.title.y = element_text(size = 30, color = "black"),
        legend.title = element_text(size = 30, color = "black"),
        legend.text = element_text(size = 30, color = "black"),
        legend.key.size = unit(1, "cm"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")
heatmap

pdf(paste0(outpath, "/637_H2H_Adonis_results.pdf"), width = 22, height = 9)
ggarrange(plot_grid(heatmap, stacked_histogram, align = "hv"), 
          common.legend = TRUE, legend = "right")
dev.off()


