# Ligand Receptor Analysis  by BulkSignaR : Bulk Data,Geomx Full Transcriptomic signal

# parameters for BulkSignalR


# you can specify one of the null model from this list (see learnParameters() documentation): 
# c("automatic", "mixedNormal", "normal", "kernelEmpirical", "empirical",  "stable")
# or 'automatic' to get the best fitting null model
null_model = 'automatic' 
paired_only = FALSE # TODO I haven't tried this : for paired samples only set to TRUE
qval_threshold = 0.01 # filter significant LR pairs

# TODO make it useful + add else to the loop if FALSE
combined_Data = TRUE # To run for combined data as well. you should keep this TRUE if you want to generate the signature score heatmap

# unnormalized expression data (exprs) or normalized data, best harmony_q3_norm
# unnormalised data will be normalised y default with UQ (upper quartile - 0.75)
data_type <- "harmony_batch_corr_q3_norm"  

# variables to merge the final csv with
meta_names <- c(aoi_id, roi_id, aoi_segment_var, sample_name, main_experimental_condition, 
                grouping_var_col_ids)

norm_is_log <- ifelse(data_type %in% c('exprs', 'q3_norm', 'deseq2_norm'), FALSE, TRUE)
normalize_needed <- ifelse(data_type == 'exprs', TRUE, FALSE)

# gene count below this thr is considered unexpressed
# for unnormalised (exprs) may be set eg to 10 - best to check on expr histogram for selected data_type
# for normalised/batch corrected etc best set to 0 to not disrupt the distribution
min_expr_counts <- ifelse(data_type == 'exprs', 10, 0)

# load data ---------------------------------------------------------------

geomx_obj <<- readRDS(geomx_norm_batch_eff_rm_path) # batch effect corrected Geomx Object

count_geomx <- data.frame(geomx_obj@assayData[[data_type]])
meta_data_all <- sData(geomx_obj)[, unique(meta_names)]

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

BulkSignaR_Output <- list()

# run for all data
if (combined_Data == TRUE){
  print("combined data")

  # TODO add norm type to output name
  BulkSignaR_Output[[output_name]] = BulkSignalR_LR_prediction(count_geomx, 
                                                              meta_data_all, 
                                                              normalize_needed, 
                                                              norm_is_log,
                                                              min_expr_counts, 
                                                              null_model,
                                                              qval_threshold, 
                                                              group = 'combined', 
                                                              file.path(output_dir, 'lr_interactions', 'bulk_signalr'))
} 


for (group in comparison) {
  
  print(group)
  
  # TODO simplify by passing metadata colnames as 1 argument      
  count_geomx_filtered_list <- Filter_for_BulkSignaR_LR_prediction(geomx_obj, 
                                                                aoi_id, 
                                                                sample_name, 
                                                                aoi_segment_var, 
                                                                main_experimental_condition, 
                                                                grouping_var_col_ids,
                                                                count_geomx,
                                                                group,
                                                                paired_only,
                                                                paired_id = NULL)
        
  BulkSignaR_Output[[group]] = BulkSignaR_LR_prediction(count_geomx_filtered_list,
                                                          normalize_needed,
                                                          normalize_method,
                                                          UQ_pc, 
                                                        file.path(output_dir, 'lr_interactions', 'bulk_signalr'),
                                                          qval_threshold,
                                                          group,
                                                          null_model)
      
      
}
  
# combining all the LR predictions from both groups to a single dataframe

lr_df_list = list()

for (group in comparison){
  
  df = BulkSignaR_Output[[group]]$LRinter_pairs_best_pws
  df$group = group
  lr_df_list[[group]] = df
  
}

df_combined = do.call(rbind, lr_df_list)
df_combined$lr_interaction = paste0(df_combined$L,"-",df_combined$R)
df_combined$qval[df_combined$qval == 0] <- 1e-70
df_combined$neg_log10_p_adj = -log(df_combined$qval)
rownames(df_combined) <- NULL


BulkSignaR_Output_list = list(
  
  BulkSignaR_Output = BulkSignaR_Output,
  count_geomx_filtered = count_geomx_filtered_list$count_geomx,
  meta_data_filtered = count_geomx_filtered_list$meta_data,
  unfiltered_LR_df_for_plotting = df_combined
  
)
  

saveRDS(BulkSignaR_Output_list, file = geomx_BulkSignalR_path)
  





