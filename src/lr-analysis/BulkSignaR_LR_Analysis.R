# Ligand Receptor Analysis  by BulkSignaR : Bulk Data,Geomx Full Transcriptomic signal

# parameters for BulkSignalR


# you can specify one of the null model from this list: c("mixedNormal", "normal", "kernelEmpirical","empirical", "stable")
# or leave null to get the best fitting null model
null_model = NULL 
paired_only = FALSE # TODO I haven't tried this : for paired samples only set to TRUE
UQ_pc = 0.75 # Upper quantile percentage
qval_threshold = 0.01 # filter significant LR pairs

# TODO make it useful + add else to the loop if FALSE
combined_Data = TRUE # To run for combined data as well. you should keep this TRUE if you want to generate the signature score heatmap

data_type <- "exprs"  # unnormalized expression data or normalized data, eg : "q3_norm" - if q3_norm normalize_needed = FALSE
normalize_needed = TRUE # FALSE if the provided data are normalized. Default "UQ", 
normalize_method = "UQ" # c("UQ",TC") or user defined(eg: q3_norm) if the provided data are normalized. ('UQ' for upper quartile or 'TC' for total count. UQ.pc = 0.75.


# load data
geomx_obj <<- readRDS(geomx_norm_batch_eff_rm_path) # batch effect corrected Geomx Object

count_geomx <- data.frame(geomx_obj@assayData[[data_type]])

# TODO remove?
# plot_dir = file.path(output_dir, output_dir, 'lr_interactions', 'bulk_signalr','plots_and_csv_files')
# dir.create(plot_dir , recursive = T, showWarnings = F)

# Run BulkSignalR predictions


BulkSignaR_Output <- list()

# run for all data
if (combined_Data == TRUE){
  print("combined data")
  
  #output_name = paste(comparison, collapse = "-combined-")
  output_name = "combined"
  
  meta_data_all = sData(geomx_obj)
  meta_data_all = meta_data_all %>% select(!!sym(aoi_id), !!sym(sample_name), !!sym(aoi_segment_var), !!sym(main_experimental_condition), all_of(grouping_var_col_ids)) 
  
  
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
  
  count_geomx_list = list(count_geomx = count_geomx,
                          meta_data = meta_data_all)
  
  BulkSignaR_Output[[output_name]] = BulkSignaR_LR_prediction(count_geomx_list,
                                                              normalize_needed, 
                                                              normalize_method, 
                                                              UQ_pc, 
                                                              file.path(output_dir, 'lr_interactions', 'bulk_signalr'), 
                                                              qval_threshold,
                                                              group = NULL,
                                                              null_model
                                                              ) 
  
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
  





