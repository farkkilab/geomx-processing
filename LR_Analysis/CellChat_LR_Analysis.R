# Ligand Receptor Analysis  by CellChat : Geomx Deconvoluted data

# TODO If All are NA (nact_status ,segment , annotation) error  : Run for the whole dataset 


# parameters for Cellchat


cell_frac_cutoff = 0.01 # 0.01 or 0.005
min_cells = 10
# row_variance_cutoff = NULL  # I did not filter based on row variance



# load data


geomx_obj = readRDS(paste0(data_dir,"/geomx_qc_norm_batch_eff_rm.RDS"))
bprism_res = readRDS(paste0(data_dir,"/bp_res_mid_lvl_ct.RDS"))
cell_fractions_df = read.csv(paste0(data_dir,"/bp_res_mid_lvl_ct_ct_fraction.csv"))


# create directory for output files

output_folder_name  = "/CellChat_objects"

if (!dir.exists(paste0(output_dir,output_folder_name))) {
  new_dir = paste0(output_dir,output_folder_name)
  dir.create(new_dir,recursive = TRUE)
  message("Directory created")
  output_dir = new_dir
} else {
  output_dir = paste0(output_dir,output_folder_name)
}




##  preprocessing and generating log normalized expr matrix


#  Extracting expression data of the desired cell types and combining 

if (is.null(cell_types)) { # If the cell types are not defined take all the cell types in the prism object
  ct_names <- colnames(bprism_res@posterior.theta_f@theta.cv)
  
}else{
  ct_names <- cell_types
}


temp <- extract_and_combine_expression_data(cell_types, bprism_res, ct_names)

deconv_ct_list_int_all = temp$deconv_ct_list_int_all
metadt_all = temp$metadt_all

###############################     normalization           ################

# TODO check the warning

# 1: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors

# Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = t(deconv_ct_list_int_all),
                              colData = metadt_all,
                              design = formula(~ Segment + NACT_status))

# TODO examine eg if add Annotation_cell or NACT status?

dds <- estimateSizeFactors(dds, type = "poscounts")

#dds <- estimateSizeFactors(dds, type = "iterate")

#dds <- estimateSizeFactors(dds) # did not work

deseq2_norm_counts <- counts(dds, normalized=TRUE)
expr_deseq2_norm_log <- log2(deseq2_norm_counts + 1) # log transformation



################   filtering based on cell fraction #####################



formatted_ct_names <- gsub(" ", ".", ct_names) #because the cell fraction dataframe colnames are different: dot instead of space

filtered_sample_list <- setNames(lapply(formatted_ct_names, function(ct_name) {
  
  # Filter out cells less than the cell_frac_cutoff
  if(!is.null(cell_frac_cutoff)){
    
    cells_greater_than_cutoff <- cell_fractions_df[cell_fractions_df[, ct_name] >= cell_frac_cutoff,]$dcc_filename
    
    cells_greater_than_cutoff <- gsub('\\.dcc', paste0('_', ct_name), cells_greater_than_cutoff)
    
    return(cells_greater_than_cutoff)
    
  }
  
  return(NULL)
  
}), ct_names) 

# Get a unique list of filtered samples
filtered_samples <- unlist(filtered_sample_list)

# Subset the expression matrix
expr_deseq2_norm_log <- expr_deseq2_norm_log[, colnames(expr_deseq2_norm_log) %in% filtered_samples]




#############  calculating the cellchat probabiities

cellchat_results <- list()

for (seg in segment) {
  
  if (is.na(seg)){seg = NULL}
  
  for (nact in nact_status) {
    
    if (is.na(nact)){nact = NULL}
    
    for (ann in annotation) {
      
      if (is.na(ann)){ann = NULL}
      
      key <- paste0(seg, "_", nact, "_", ann)
      cellchat_results[[key]] <- cellchat_predict_prob(
        metadt_all,
        expr_deseq2_norm_log,
        nact_status = nact,
        segment = seg,
        annotation = ann,
        output_dir = output_dir
      )

    }
  }
}


# save final LR predictions

parts <- c(nact_status, segment, annotation)
parts_non_NA <- parts[!sapply(parts, is.na)]
file_name = paste(parts_non_NA, collapse = "_")

cellchat_output = list(
  expr_deseq2_norm_log = expr_deseq2_norm_log,
  metadt_all = metadt_all,
  cellchat_results = cellchat_results
  )

saveRDS(cellchat_output, file = paste0(output_dir,'/CellChat_',file_name,'_output.RDS'))



