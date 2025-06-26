# Ligand Receptor Analysis  by BulkSignaR : Bulk Data,Geomx Full Transcriptomic signal

# parameters for BulkSignalR

# TODO : NULL model check !!

null.model = NULL # c("automatic", "mixedNormal", "normal", "kernelEmpirical","empirical", "stable")
normalize_needed = TRUE # Default "UQ", FALSE if the provided data are normalized
paired_only = FALSE # for paired samples only TRUE
normalize_method = "UQ" # c("UQ",TC") or user defined if the provided data are normalized. ('UQ' for upper quartile or 'TC' for total count. UQ.pc = 0.75.

UQ_pc = 0.75
qval_threshold = 0.01 # filter significant LR pairs

combined_Data = TRUE # To run for combined data as well

# load data



count_geomx = data.frame(geomx_obj@assayData$exprs) # count data
#count_geomx  = data.frame(geomx_obj@assayData$q3_norm) # q3 normalized data


plot_dir = file.path(output_dir, BulkSignalR_folder_name,'plots_and_csv_files')
dir.create(plot_dir , recursive = T, showWarnings = F)



# Run BulkSignalR predictions


BulkSignaR_Output <- list()

if (combined_Data == TRUE){
  print("combined data")
  
  #output_name = paste(comparison, collapse = "-combined-")
  output_name = "combined"
  
  meta_data_all = sData(geomx_obj)
  meta_data_all = meta_data_all %>% select(!!sym(aoi_id), !!sym(sample_name), !!sym(aoi_segment_var), !!sym(main_experimental_condition), all_of(grouping_var_col_ids)) 
  
  count_geomx_list = list(count_geomx = count_geomx,
                          meta_data = meta_data_all)
  
  BulkSignaR_Output[[output_name]] = BulkSignaR_LR_prediction(count_geomx_list,
                                                              normalize_needed, 
                                                              normalize_method, 
                                                              UQ_pc, 
                                                              plot_dir, 
                                                              qval_threshold,
                                                              group = NULL
                                                              ) 
  
}


  
for (group in comparison) {
  
  print(group)
        
  count_geomx_filtered_list <- Filter_for_BulkSignaR_LR_prediction(geomx_obj, 
                                                                aoi_id, 
                                                                sample_name, 
                                                                aoi_segment_var, 
                                                                main_experimental_condition, 
                                                                grouping_var_col_ids, 
                                                                paired_only = FALSE, 
                                                                count_geomx,
                                                                group)
        
  BulkSignaR_Output[[group]] = BulkSignaR_LR_prediction(count_geomx_filtered_list,
                                                          normalize_needed,
                                                          normalize_method,
                                                          UQ_pc, 
                                                          plot_dir,
                                                          qval_threshold,
                                                          group)
      
      
}
  
BulkSignaR_Output_list = list(
  
  BulkSignaR_Output = BulkSignaR_Output,
  count_geomx_filtered = count_geomx_filtered_list$count_geomx,
  meta_data_filtered = count_geomx_filtered_list$meta_data
  
)
  

saveRDS(BulkSignaR_Output_list, file = geomx_BulkSignalR_path)
  





