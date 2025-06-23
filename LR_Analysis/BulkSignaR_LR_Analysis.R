# Ligand Receptor Analysis  by BulkSignaR : Bulk Data,Geomx Full Transcriptomic signal

# parameters for BulkSignalR

# TODO : NULL model check !!

null.model = NULL # c("automatic", "mixedNormal", "normal", "kernelEmpirical","empirical", "stable")
normalize_needed = TRUE # Default "UQ", FALSE if the provided data are normalized
paired_only = FALSE # for paired samples only TRUE
normalize_method = "UQ" # c("UQ",TC") or user defined if the provided data are normalized. ('UQ' for upper quartile or 'TC' for total count. UQ.pc = 0.75.

UQ_pc = 0.75
qval_threshold = 0.01 # filter significant LR pairs

# load data


geomx_obj <- readRDS(paste0(data_dir,'/geomx_qc_norm_batch_eff_rm.RDS'))
count_geomx = data.frame(geomx_obj@assayData$exprs) # count data
#count_geomx  = data.frame(geomx_obj@assayData$q3_norm) # q3 normalized data


# create directory for output files

output_folder_name  = "/BulkSignalR_objects"

if (!dir.exists(paste0(output_dir,output_folder_name))) {
  new_dir = paste0(output_dir,output_folder_name)
  dir.create(new_dir,recursive = TRUE)
  message("Directory created")
  output_dir = new_dir
} else {
  output_dir = paste0(output_dir,output_folder_name)
}


parts <- c(nact_status, segment, annotation)
parts_non_NA <- parts[!sapply(parts, is.na)]


# Run BulkSignalR predictions


if (length(parts_non_NA) == 0) {
  
  nact_status = NULL
  segment = NULL
  annotation = NULL
  
  # Filteration
  
  count_geomx_filtered <- Filter_for_BulkSignaR_LR_prediction(geomx_obj, 
                                                              nact_status, 
                                                              segment, 
                                                              annotation, 
                                                              paired_only, 
                                                              count_geomx)
  
  BulkSignaR_Output = BulkSignaR_LR_prediction(count_geomx_filtered$count_geomx,
                                               count_geomx_filtered$meta_data,
                                               normalize_needed, 
                                               normalize_method, 
                                               UQ_pc, 
                                               output_dir , 
                                               nact_status, 
                                               segment, 
                                               annotation, 
                                               qval_threshold
                                               )
  
  saveRDS(BulkSignaR_Output, file = paste0(output_dir,'/BulkSignalR_combined_output.RDS'))
  
  
} else {
  
  BulkSignaR_Output <- list()
  
  for (seg in segment) {
    
    if (is.na(seg)){seg = NULL}
    
    for (nact in nact_status) {
      
      if (is.na(nact)){nact = NULL}
      
      for (ann in annotation) {
        
        if (is.na(ann)){ann = NULL}
        
        key <- paste0(seg, "_", nact, "_", ann)
        
        count_geomx_filtered <- Filter_for_BulkSignaR_LR_prediction(geomx_obj, 
                                                                    nact, 
                                                                    seg, 
                                                                    ann, 
                                                                    paired_only, 
                                                                    count_geomx)
        
        BulkSignaR_Output[[key]] = BulkSignaR_LR_prediction(count_geomx_filtered$count_geomx,
                                                            count_geomx_filtered$meta_data,
                                                            normalize_needed, 
                                                            normalize_method, 
                                                            UQ_pc, 
                                                            output_dir , 
                                                            nact_status = nact, 
                                                            segment = seg, 
                                                            annotation = ann, 
                                                            qval_threshold
                                                            )
        
        file_name = paste(parts_non_NA, collapse = "_")
        saveRDS(BulkSignaR_Output, file = paste0(output_dir,'/BulkSignalR_',file_name,'_output.RDS'))
        
        
      }
    }
  }
  
}



