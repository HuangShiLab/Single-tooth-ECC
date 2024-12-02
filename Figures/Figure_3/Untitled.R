clinical_features <- c("sum_s_dt", "sum_ns_dt", "sum_dmfs", "spatial_dist_weighted_mean_dmfs")
mic_dist_features <- c("spatial_dist_weighted_mean_md", "mean_md", "mean_md_H_T51", "mean_md_H_T52",
                       "mean_md_H_T53", "mean_md_H_T54", "mean_md_H_T55", "mean_md_H_T61", "mean_md_H_T62",
                       "mean_md_H_T63", "mean_md_H_T64", "mean_md_H_T65", "mean_md_H_T71", "mean_md_H_T72",
                       "mean_md_H_T73", "mean_md_H_T74", "mean_md_H_T75", "mean_md_H_T81", "mean_md_H_T82",
                       "mean_md_H_T83", "mean_md_H_T84", "mean_md_H_T85", "mean_md_H_T5161")

df <- read_hdf5_biom("../../Data/1867/taxonomy/1867_taxonomic_rarefied_ASV_abd_table.biom")
df <- biom(df)
df <- biom_data(df)
df <- t(as.matrix(df))
df <- data.frame(df, check.names = F)
asv <- df[, which((!colnames(df) %in% clinical_features) & (!colnames(df) %in% mic_dist_features))]

df <- read_hdf5_biom("../../Data/1867/taxonomy/1867_taxonomic_feature_abd_table.biom")
df <- biom(df)
df <- biom_data(df)
df <- t(as.matrix(df))
df <- data.frame(df, check.names = F)
df <- df[, which(colnames(df) %in% c(clinical_features, mic_dist_features))]

identical(rownames(asv), rownames(df))

new_df <- cbind(asv, df)
new_df <- new_df[, order(colnames(new_df))]
new_df <- t(new_df)
write.table(new_df, "1867_taxonomic_rarefied_feature_abd_table.tsv", sep = "\t", quote = F, row.names = T, col.names = NA)
