#--------------------------------------------------
p <- c("reshape2", "ade4", "vegan")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE, repos = "http://cran.us.r-project.org")
  suppressWarnings(suppressMessages(invisible(require(p, character.only = TRUE))))
}
invisible(lapply(p, usePackage))

#' @title check_metadata
#' @param metadata A dataframe with > two columns corresponds to samples (rownames) in the biological data.
#' @param more_missing_values A optional string(s) can be added to define the missing values.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' x0 <- data.frame(rbind(t(rmultinom(7, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(8, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289)))))
#' a<-factor(c(rep("A", 29), NA, rep("B", 29), NA))
#' b<-factor(c(rep("A", 27), NA, "Not applicable", "Missing:not collected", rep("B", 28), NA, NA))
#' c<-factor(c(rep("A", 20), rep("B", 20), rep("C", 20)))
#' d<-c(sample(1:18), NA, "Not applicable", sample(5:44))
#' e<-c(rnorm(40), NA, "Not applicable", rnorm(18, 4))
#' e0<-c(NA, "Not applicable", rnorm(58, 4))
#' f<-rep("C", 60)
#' g<-rep(4, 60)
#' h<-c(rep("C", 59), "B")
#' metadata<-data.frame(a, b, c, d, e, e0, f, g, h)
#' check_metadata(metadata)
#' @author Shi Huang
#' @export
check_metadata <- function(metadata, more_missing_values=NULL, unique_rate_thres=0.2){
  if(is.null(more_missing_values)){
    missing_values<-c("not applicable", "Not applicable", "Missing:not collected",
                      "Not provided", "missing: not provided", "unknown",
                      "not provided", "not collected", "NA", NA, "")
  }

  check_integers <- function(vector, all=TRUE){
    checker <- grepl("^[0-9]+$", as.character(vector), perl = T)
    total_len <- length(vector)
    n_integers <- sum(checker)
    n_non_integers <- total_len - n_integers
    #cat("Number of integers: ", n_integers, "\n")
    #cat("Number of non-integers: ", n_non_integers, "\n")
    if(all) if_all_integers <- all(checker) ## if any TRUE in the checker
    out <- list(if_all_integers=if_all_integers,
                n_integers=n_integers,
                n_non_integers=n_non_integers)
    return(out)
  }

  check_characters <- function(vector){
    checker <- grepl("[A-Za-z-$_]+", as.character(vector), perl = T) # [A-Za-z[:punct:]]+
    total_len <- length(vector)
    n_characters <- sum(checker)
    n_non_characters <- total_len - n_characters
    #cat("Number of characters: ", n_characters, "\n")
    #cat("Number of non-characters: ", n_non_characters, "\n")
    if_character_existed <- any(checker) ## if any TRUE in the checker

    out <- list(if_character_existed=if_character_existed,
                n_characters=n_characters,
                n_non_characters=n_non_characters)
    return(out)
  }

  check_numeric <- function(vector){
    checker <- suppressWarnings(as.numeric(as.character(vector)) %% 1 != 0)
    total_len <- length(vector)
    n_numeric <- sum(checker[!is.na(checker)]) # if any NA values it means a string appear among the numeric values.
    n_non_numeric <- total_len - n_numeric
    #cat("Number of numeric: ", n_numeric, "\n")
    #cat("Number of non-numeric: ", n_non_numeric, "\n")
    if_numeric_existed <- any(checker[!is.na(checker)])
    ## if any TRUE in the checker
    out <- list(if_numeric_existed=if_numeric_existed,
                n_numeric=n_numeric,
                n_non_numeric=n_non_numeric)
    return(out)
  }

  n_unique_values <- sapply(metadata, function(x) nlevels(factor(x)))
  unique_values <- sapply(metadata, function(x) paste0(levels(factor(x)), collapse="|"))
  n_missing_values <- sapply(metadata, function(x) sum(x %in% missing_values))
  n_real_values <- sapply(metadata, function(x) sum(!x %in% missing_values))
  completeness <- n_real_values/nrow(metadata)
  unique_rate <- ifelse(n_real_values==0, 0, n_unique_values/n_real_values)
  #n_unique_real_values <- sapply(metadata, function(x) nlevels(as.factor(x[!x %in% missing_values])))
  if_character_existed <- sapply(metadata, function(x) check_characters(x[!x %in% missing_values])[[1]])
  n_characters <- sapply(metadata, function(x) check_characters(x[!x %in% missing_values])[[2]])
  if_all_integers <- sapply(metadata, function(x) check_integers(x[!x %in% missing_values])[[1]])
  n_integers <- sapply(metadata, function(x) check_integers(x[!x %in% missing_values])[[2]])
  if_numeric_existed <- sapply(metadata, function(x) check_numeric(x[!x %in% missing_values])[[1]])
  n_numeric <- sapply(metadata, function(x) check_numeric(x[!x %in% missing_values])[[2]])
  all_values_identical <- apply(metadata, 2, function(x) length(unique(x))==1)
  numeric_var <- if_numeric_existed | (if_all_integers & unique_rate >= unique_rate_thres)
  categorical_var <- if_character_existed | (if_all_integers & unique_rate < unique_rate_thres)
  metadata_summ<-data.frame(metadata=colnames(metadata),
                            all_values_identical,
                            total_n=nrow(metadata),
                            n_missing_values,
                            n_real_values,
                            completeness,
                            n_unique_values,
                            unique_values,
                            unique_rate,
                            if_character_existed,
                            n_characters,
                            if_all_integers,
                            n_integers,
                            if_numeric_existed,
                            n_numeric,
                            numeric_var,
                            categorical_var)

  metadata_summ
}

normalize_NA_in_metadata<-function(md){
  apply(md, 1, function(x) {
   idx <- which(x=="not provided" | x=="Not provided" | x=="Not Provided"
          | x=="not applicable" | x=="Not applicable"
          | x=="Missing:not collected"
          | x=="NA" | x=="na" | x=="Na"
          | x=="none" | x=="None" | x=="NONE")
   x[idx]<-NA
  })
  md
}

discard_uninfo_columns_in_metadata<-function(md){
  noninfo_idx<-which(apply(md, 2, function(x) length(unique(x))==1))
  md<-md[-noninfo_idx]
  md
}

trim_metadata <- function(md, completeness_threshold=0.5, filter_cols_by_type="numeric"){
  if(!is.element(filter_cols_by_type, c("numeric", "categorical", NA)))
    stop("Only 'numeric', 'categorical', or NA are allowed for 'filter_cols_by_type'!")
  md<-normalize_NA_in_metadata(md)
  md_summ<-check_metadata(metadata = md)
  filtered_md_summ <- md_summ %>% 
    filter(n_real_values!=n_unique_values) %>% 
    filter(completeness>completeness_threshold) %>% 
    filter(n_unique_values>1)
  if(is.na(filter_cols_by_type)){
    cat("No column was filtered by data type!\n")
  }else if(filter_cols_by_type=="numeric"){
    filtered_md_summ <- filtered_md_summ %>% filter(numeric_var==TRUE)
  }else{
    filtered_md_summ <- filtered_md_summ %>% filter(categorical_var==TRUE)
  }
  vars<-filtered_md_summ[, "metadata"]
  return(md[, vars])
}

chop_seq_to_x_nt<-function(df, start=1, nt=100){
  cat("The length of feature sequences all equal to ", nt, " before :", all(nchar(colnames(df))==nt), "\n")
  if(!all(nchar(colnames(df))==nt)){
    if(start==1){end <- nt}else{end <- start + nt -1 }
    colnames(df)<-substr(colnames(df), start, end)
  }
  cat("The length of feature sequences all equal to ", nt, " after :", all(nchar(colnames(df))==nt), "\n")
  df
}




filter_features_allzero<-function(data, samples=TRUE, features=TRUE){
  if(samples & features){
    result<-data[which(apply(data, 1, sum)!=0), ]
    result<-data[, which(apply(result, 2, sum)!=0)]
  }else if(samples & !features){
    result<-data[which(apply(data, 1, sum)!=0), ]
  }else if(!samples & features){
    result<-data[, which(apply(data, 2, sum)!=0)]
  }else{
    stop("Nothing has been done!")
  }
  result
}

filter_features_by_prev <- function(data, prev=0.001){
  data<-data[, which(colSums(data!=0) > prev * nrow(data))]
  data
}

filter_features_by_abundance <- function(data, mean_abd_cutoff=0.001){
  hist(colMeans(data))
  data<-data[, which(colMeans(data) > mean_abd_cutoff)]
  data
}

remove_noisy_feature_by_count_cutoff <- function(data, count_cutoff=1000){
  if(class(data)=="data.frame") data<-data.matrix(data)
  cat("The total number of zeros in the table:", sum(data==0), "\n")
  data[data<=count_cutoff]<-0
  cat("The total number of zeros in the filtered table:", sum(data==0))
  data<-data.frame(data)
  data
}

filter_features_with_NA <- function(data, min_prev=0){
  idx <- which(colSums(is.na(data)) <= min_prev * nrow(data))
  data<-data[, idx]
  data
}

filter_samples_with_NA <- function(data, min_prev=0){
  idx <- which(rowSums(is.na(data)) <= min_prev * ncol(data))
  data<-data[idx, ]
  data
}


filter_samples_by_NA_in_y <- function(data, y){
  y_k<-y[which(!is.na(y))]
  data_k<-data[which(!is.na(y)) ,]
  result<-list()
  result$data_k<-data_k
  result$y_k<-y_k
  result
}

#' @title filter_samples_by_groups_in_target_field_of_metadata
#' @param data A data.frame.
#' @param metadata A data.frame including multiple metadata variables corresponds to samples in the data.
#' @param target_field A character string indicating a target variable in the metadata.
#' @param groups A character string(s) indicating one or multiple groups in the target metadata variable specified.
#' @param negate A bool value indicates if samples in the specified groups should be removed (TRUE) or kept (FALSE).
#' @examples
#' set.seed(123)
#' data <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' z<-factor(c(rep("E", 20), rep("C", 20), rep("D", 20)))
#' y<-factor(c(rep("A", 30), rep("B", 30)))
#' y0<-factor(c(rep("A", 5), rep("B", 55)))
#' metadata <- data.frame(z, y, y0)
#' filter_samples_by_groups_in_target_field_of_metadata(data, metadata, target_field="z", groups=c("D", "A"), negate=TRUE)
#' @export
filter_samples_by_groups_in_target_field_of_metadata <- function(data, metadata, target_field, groups, negate=FALSE, ids_col=NA){
  if(is.na(ids_col)){
    SampleIDs<-rownames(metadata)
  }else{
    SampleIDs<-metadata[, ids_col]
  }
  target_field_levels <- levels(factor(metadata[, target_field]))
  if(!identical(rownames(data), SampleIDs)) stop("The sample IDs should be idenical in feature table and metadata!")
  if_groups_exist<-target_field_levels %in% groups
  if(all(if_groups_exist==FALSE)) stop("All specified groups do not exist in the target field!")
  if(any(if_groups_exist==FALSE)) cat("Only specified groups:",  target_field_levels[if_groups_exist], "existed in the target field!\n")
  groups<-target_field_levels[if_groups_exist]
  if(negate){
    idx<-which(!metadata[, target_field] %in% groups)
  }else{
    idx<-which(metadata[, target_field] %in% groups)
  }
  metadata_k<-metadata[idx, ]
  metadata_k[, target_field] <- factor(metadata_k[, target_field])
  data_k<-data[idx, ]
  cat("The number of kept samples (only samples related to ", groups," in ",target_field,"): ", nrow(metadata_k) ,"\n")
  result<-list()
  result$data<-data_k
  result$metadata<-metadata_k
  result
}

filter_samples_by_NA_in_target_field_of_metadata <- function(data, metadata, target_field, ids_col=NA){
  if(is.na(ids_col)){
    SampleIDs<-rownames(metadata)
  }else{
    SampleIDs<-metadata[, ids_col]
  }
  if(!identical(rownames(data), SampleIDs)) stop("The sample IDs should be idenical in feature table and metadata!")
  NAN_values<-c("not applicable", "Not applicable", "Missing:not collected",
                "Not provided", "missing: not provided",
                "not provided", "not collected", "NA", NA, "")
  idx<-which(!metadata[, target_field] %in% NAN_values)
  metadata_k<-metadata[idx, ]
  data_k<-data[idx, ]
  cat("The number of kept samples (after filtering out samples with NA values in ",target_field,"): ", nrow(metadata_k) ,"\n")
  result<-list()
  result$data<-data_k
  result$metadata<-metadata_k
  result
}

split_dm_by_metadata<-function(dm, metadata, split_factor){
  #if(is.null(rownames(dm))#
  if(!is.element(split_factor, colnames(metadata)))
    stop("The split_factor you specified should be one of the column names of input metadata.")
  if(class(dm)=="dist") dm <- data.matrix(dm)
  f<-metadata[, split_factor]
  sub_dm_name_list<-split(1:nrow(dm), f, drop=TRUE)
  sub_dm_list<-
    lapply(sub_dm_name_list, function(x){
      list(dm[x, x], metadata[x, ])
    })

  sub_dm_list

}

filter_dm_by_NA_in_target_field_of_metadata <- function(dm, metadata, target_field, ids_col=NA){
  if(is.na(ids_col)){
    SampleIDs<-rownames(metadata)
  }else{
    SampleIDs<-metadata[, ids_col]
  }
  if(!identical(rownames(dm), SampleIDs)) stop("The sample IDs should be idenical in feature table and metadata!")
  NAN_values<-c("not applicable", "Not applicable", "Missing:not collected",
                "Not provided", "missing: not provided", "unknown", "Unknown",
                "not provided", "not collected", "NA", NA, "")
  idx<-which(!metadata[, target_field] %in% NAN_values)
  metadata_k<-metadata[idx, ]
  dm_k<-dm[idx, idx]
  cat("The number of kept samples (after filtering out samples with NA values in ",target_field,"): ", nrow(metadata_k) ,"\n")
  result<-list()
  result$dm<-dm_k
  result$metadata<-metadata_k
  result
}

filter_samples_by_sample_ids_in_metadata <- function(data, metadata, ids_col=NA){
  if(is.na(ids_col)){
    shared_ids<-intersect(rownames(data), rownames(metadata))
    metadata_idx<-which(rownames(metadata) %in% shared_ids)
  }else{
    shared_ids<-intersect(rownames(data), metadata[, ids_col])
    metadata_idx<-which(metadata[, ids_col] %in% shared_ids)
  }
  data_matched<-data[shared_ids, ]
  data_matched<-data_matched[order(rownames(data_matched)),]
  cat("The number of samples in feature table (after filtering out samples with no metadata): ",
      nrow(data_matched) ,"\n")
  metadata_matched<-metadata[metadata_idx, ]

  cat("The number of samples metadata (after filtering out samples with no metadata): ",
      nrow(metadata_matched) ,"\n")

  if(is.na(ids_col)){
    metadata_matched<-metadata_matched[order(rownames(metadata_matched)),]
    cat("The sample IDs are idenical in feature table and metadata: ",
        identical(rownames(data_matched), rownames(metadata_matched)), "\n")
  }else{
    metadata_matched<-metadata_matched[order(as.character(metadata_matched[, ids_col])),]
    cat("The sample IDs are idenical in feature table and metadata: ",
        identical(rownames(data_matched), as.character(metadata_matched[, ids_col])), "\n")
  }


  result<-list()
  result$data<-data_matched
  result$metadata<-metadata_matched
  return(result)
}

filter_samples_by_sample_ids_in_metadata <- function(data, dm=FALSE, metadata, ids_col=NA){
  data<-data[order(rownames(data)), order(colnames(data))]
  metadata<-metadata[order(rownames(metadata)),]
  if(is.na(ids_col)){
    shared_ids<-intersect(rownames(data), rownames(metadata))
    metadata_idx<-which(rownames(metadata) %in% shared_ids)
  }else{
    shared_ids<-intersect(rownames(data), metadata[, ids_col])
    metadata_idx<-which(metadata[, ids_col] %in% shared_ids)
  }
    if(dm==TRUE){
     data_matched<-data[shared_ids, shared_ids]
     data_matched<-data_matched[order(rownames(data_matched)),order(colnames(data_matched))]
     cat("The number of samples in distance matrix (after filtering out samples with no metadata): ",
      nrow(data_matched) ,"\n")
    }else{
     data_matched<-data[shared_ids, ]
     data_matched<-data_matched[order(rownames(data_matched)),]
     cat("The number of samples in feature table (after filtering out samples with no metadata): ",
      nrow(data_matched) ,"\n")
    }

  metadata_matched<-metadata[metadata_idx, ]
  
  cat("The number of samples metadata (after filtering out samples with no metadata): ",
      nrow(metadata_matched) ,"\n")
  
  if(is.na(ids_col)){
    metadata_matched<-metadata_matched[order(rownames(metadata_matched)),]
    cat("The sample IDs are idenical in feature table and metadata: ",
        identical(rownames(data_matched), rownames(metadata_matched)), "\n")
  }else{
    metadata_matched<-metadata_matched[order(as.character(metadata_matched[, ids_col])),]
    cat("The sample IDs are idenical in feature table and metadata: ",
        identical(rownames(data_matched), as.character(metadata_matched[, ids_col])), "\n")
  }

  result<-list()
  result$data<-data_matched
  result$metadata<-metadata_matched
  return(result)
}

filter_samples_by_seq_depth<-function(data, metadata, cutoff=1000){
  seq_dep_idx<-which(rowSums(data) > cutoff)
  cat("The number of kept samples with more than ",cutoff," reads: ", length(seq_dep_idx), "\n")
  if(length(seq_dep_idx)>0){
    metadata_k<-metadata[seq_dep_idx, ]
    data_k<-data[seq_dep_idx, ]
  }else{
    stop("The read count of all samples is less than sequencing depth threshold!")
  }
  cat("The sample IDs are idenical in feature table and metadata: ", identical(rownames(data_k), rownames(metadata_k)), "\n")
  #metadata_k$Seq_depth<-rowSums(df_k)
  #p<-ggplot(metadata_k, aes(Seq_depth)) + geom_histogram() + xlim(1000, 70000) + theme_bw()
  #p
  result<-list()
  result$data<-data_k
  result$metadata<-metadata_k
  return(result)
}

convert_y_to_numeric<-function(y, reg=TRUE){
  if(reg & !is.numeric(y)){
    y=as.numeric(as.character(y))
  }else if(!reg){
    y=factor(y)
  }else{
    y=y
  }
  y
}

keep_shared_features<-function(train_x, test_x){
  common_idx<-intersect(colnames(train_x), colnames(test_x))
  train_x_shared<-train_x[, common_idx]
  test_x_shared<-test_x[, common_idx]
  cat("Number of features kept:", length(common_idx), "\n")
  cat("The proportion of commonly shared features in train and test dataset respectively: \n")
  cat("Train data: ", length(common_idx)/ncol(train_x), "\n")
  cat("Test data: ", length(common_idx)/ncol(test_x), "\n")
  result<-list()
  result$train_x_shared<-train_x_shared
  result$test_x_shared<-test_x_shared
  result
}


#' feature metadata
add_ann<-function(tab, fmetadata, tab_id_col=1, fmetadata_id_col=1){
  fmetadata[, fmetadata_id_col]<-as.character(fmetadata[, fmetadata_id_col])
  tab[, tab_id_col]<-as.character(tab[, tab_id_col])
  matched_idx<-which(fmetadata[, fmetadata_id_col] %in% tab[, tab_id_col])
  uniq_features_len<-length(unique(tab[, tab_id_col]))
  if(uniq_features_len>length(matched_idx)){
    warning("# of features has no matching IDs in the taxonomy file: ", uniq_features_len-length(matched_idx), "\n")
  }
  fmetadata_matched<-fmetadata[matched_idx,]
  out<-merge(tab, fmetadata_matched, by.x=tab_id_col, by.y=fmetadata_id_col)
  out
}

rbind.na<-function(l){
  max_len<-max(unlist(lapply(l, length)))
  c_l<-lapply(l, function(x) {c(x, rep(NA, max_len - length(x)))})
  do.call(rbind, c_l)
}
expand_Taxon<-function(df, Taxon){
  taxa_df <- rbind.na(strsplit(as.character(df[, Taxon]), '; '))
  colnames(taxa_df) <- c("kingdom","phylum","class","order","family","genus","species") #"kingdom",
  data.frame(df, taxa_df)
}
