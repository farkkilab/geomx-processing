# Ligand Receptor Analysis  by CellChat : Geomx Deconvoluted data

# TODO If All are NA (nact_status ,segment , annotation) error  : Run for the whole dataset 


# parameters for Cellchat


cell_frac_cutoff = 0.01 # 0.01 or 0.005
min_cells = 10 # Number of minimum cells in each cell group


#  Extracting expression data and meta data of the desired cell types and combining 

if (is.null(cell_types)) { # If the cell types are not defined take all the cell types in the prism object
  ct_names <- colnames(bprism_res@posterior.theta_f@theta.cv)
  
  
}else{
  ct_names <- cell_types
}

# expr_and_meta_list

combined_expression_data <- extract_and_combine_expression_data(bprism_res, ct_names, geomx_obj,aoi_id, sample_name, aoi_segment_var, main_experimental_condition, grouping_var_col_ids)
deconv_ct_list_int_all <- combined_expression_data$deconv_ct_list_int_all
metadt_all <- combined_expression_data$metadt_all



# save files

file_name = paste(ct_names, collapse = "_")
saveRDS(combined_expression_data, file = file.path(output_dir,paste0(file_name,'_expr_and_meta_list.RDS')))



# normalization           

# TODO check the warning
# 1: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
# Create DESeq2Dataset object
design_formula <- as.formula(
  paste("~", aoi_segment_var, "+", main_experimental_condition)
)

dds <- DESeqDataSetFromMatrix(countData = t(deconv_ct_list_int_all),
                              colData = metadt_all,
                              design = design_formula) # TODO check the warning message

# TODO examine eg if add Annotation_cell or NACT status?

dds <- estimateSizeFactors(dds, type = "poscounts")
#dds <- estimateSizeFactors(dds, type = "iterate")
#dds <- estimateSizeFactors(dds) # did not work
deseq2_norm_counts <- counts(dds, normalized=TRUE)
expr_deseq2_norm_log <- log2(deseq2_norm_counts + 1) # log transformation



# filtering based on cell fraction 

expr_deseq2_norm_log_cf_filtered = filter_based_on_cell_fraction(ct_names, cell_frac_cutoff, cell_fractions_df, expr_deseq2_norm_log)

  
#  calculating the cellchat probabilties

cellchat_results <- list()

for (group in comparison) {
  
  print(group)
  cellchat_results[[group]] <- cellchat_predict_prob(
  metadt_all,
  expr_deseq2_norm_log_cf_filtered,
  sample_name = sample_name,
  grouping_var_col_ids = grouping_var_col_ids,
  group = group,
  output_dir =  output_dir
  )
  }


# combining all the LR predictions from both groups to a single dataframe


lr_df_list = list()

for (group in comparison){
  
  df = cellchat_results[[group]]$df.net
  df$group = group
  lr_df_list[[group]] = df
  
}

df_combined = do.call(rbind, lr_df_list)
df_combined$sender_receiver = paste(df_combined$source, df_combined$target, sep = " -> ")
rownames(df_combined) <- NULL

# save final LR predictions


cellchat_output = list(
  expr_deseq2_norm_log = expr_deseq2_norm_log,
  metadt_all = metadt_all,
  cellchat_results = cellchat_results,
  unfiltered_LR_df_for_plotting = df_combined
  )

saveRDS(cellchat_output, file = geomx_CellChat_path)



