# Ligand Receptor Analysis  by BulkSignaR : Bulk Data,Geomx Full Transcriptomic signal

# parameters for BulkSignalR

# paired_only = FALSE # TODO I haven't tried this : for paired samples only set to TRUE
# TODO rmv - in the master script
#qval_threshold = 0.01 # filter significant LR pairs

combined_Data = TRUE # To run for all data as well. you should keep this TRUE if you want to generate the signature score heatmap

# unnormalized expression data (exprs) or normalized data, best harmony_q3_norm
# unnormalised data will be normalised by default with quartile norm UQ (upper quartile - 0.75)
data_type <- "harmony_batch_corr_q3_norm"  

# variables to merge the final csv with
meta_names <- c(aoi_id, roi_id, aoi_segment_var, sample_name, main_experimental_condition, 
                grouping_var_col_ids)

# you can specify one of the null model from this list (see learnParameters() documentation): 
# c("automatic", "mixedNormal", "normal", "kernelEmpirical", "empirical",  "stable")
# or 'automatic' to get the best fitting null model
null_model = 'automatic' 

norm_is_log <- ifelse(data_type %in% c('exprs', 'q3_norm', 'deseq2_norm'), FALSE, TRUE)
normalize_needed <- ifelse(data_type == 'exprs', TRUE, FALSE)

# gene count below this thr is considered unexpressed
# for unnormalised (exprs) may be set eg to 10 - best to check on expr histogram for selected data_type
# for normalised/batch corrected etc best set to 0 to not disrupt the distribution
min_expr_counts <- ifelse(data_type == 'exprs', 10, 0)

# load data ---------------------------------------------------------------

geomx_obj <<- readRDS(geomx_norm_batch_eff_rm_path) # batch effect corrected Geomx Object

count_geomx <- as.data.frame(geomx_obj@assayData[[data_type]])
meta_data_all <- sData(geomx_obj)[, unique(meta_names)]

# make variable with all categories from grouping_var_col_ids
meta_data_all$comparison_group <- apply(meta_data_all, 1, function(row){
  group <- sapply(grouping_var_col_ids, function(var){
    paste(row[var])
  })
  group <- paste(group, collapse = '_')
  return(group)
})

# TODO change it to filter through Paired_status column
# if (paired_only == TRUE){
#   
#   paired_samples = filter_paired_data(geomx_obj, main_experimental_condition, paired_id)
#   meta_data_all = meta_data_all %>% filter(!!sym(sample_name) %in% paired_samples)
#   meta_data_all[,aoi_id] = gsub('-', '.', meta_data_all[,aoi_id])
#   col_ids = colnames(count_geomx)  %in% meta_data_all[,aoi_id]
#   count_geomx = count_geomx[,col_ids]
#   
# }

# Run BulkSignalR predictions ---------------------------------------------

BulkSignalR_Output <- list()

# run for all data
if (combined_Data == TRUE){
  print("combined data")

  # TODO add norm type to output name
  BulkSignalR_Output[['combined']] = BulkSignalR_LR_prediction(count_geomx, 
                                                              meta_data_all, 
                                                              normalize_needed, 
                                                              norm_is_log,
                                                              min_expr_counts, 
                                                              null_model,
                                                              qval_threshold, 
                                                              group = 'combined', 
                                                              file.path(output_dir, 'lr_interactions', 'bulk_signalr'))
} 

# separately for groups
for(group in unique(meta_data_all$comparison_group)) {
  print(group)
  
  #filter to group 
  meta_data_group <- meta_data_all[meta_data_all$comparison_group == group, ]
  count_geomx_group <- count_geomx[, meta_data_group$dcc_filename]
  
  BulkSignalR_Output[[group]] = BulkSignalR_LR_prediction(count_geomx_group, 
                                                              meta_data_group, 
                                                              normalize_needed, 
                                                              norm_is_log,
                                                              min_expr_counts, 
                                                              null_model,
                                                              qval_threshold, 
                                                              group = group, 
                                                              file.path(output_dir, 'lr_interactions', 'bulk_signalr'))

      
}

saveRDS(BulkSignalR_Output, geomx_BulkSignalR_path)  

# combine all LR predictions to a single dataframe list -------------------

lr_df_all_groups <- lapply(names(BulkSignalR_Output), function(group){
  print(group)
  # loop through all bsrinf objects with different reducing options (see vignette)
  lr_df_bsrinf <- lapply(grep('bsrinf', names(BulkSignalR_Output[[group]]), value = T), function(bsrinf_name){
    print(bsrinf_name)
    
    bsrinf <- BulkSignalR_Output[[group]][[bsrinf_name]]
    reduction_type <- gsub('bsrinf_', '', bsrinf_name)
    
    # extracting and filtering LR dataframes
    LRinter.df <- LRinter(bsrinf) %>%
      filter(qval <= qval_threshold) %>%
      arrange(desc(qval)) %>%
      mutate(qval = ifelse(qval == 0, 1e-70, qval),
             neg_log10_p_adj = -log10(qval),
             reduction = reduction_type,
             group = group) 
    
    rownames(LRinter.df) <- NULL
    
    return(LRinter.df)
  })
  return(lr_df_bsrinf)
})


lr_df_combined <- do.call(Map, c(f = rbind, lr_df_all_groups))
names(lr_df_combined) <- grep('bsrinf', names(BulkSignalR_Output[[1]]), value = T)

# save BRS output and LR dfs
saveRDS(lr_df_combined, file = lr_output_path)

