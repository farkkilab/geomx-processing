# Ligand Receptor Analysis  by CellChat : Geomx Deconvoluted data

# Cellchat expects Normalized data (e.g., library-size normalization and then log-transformed with a pseudocount of 1)
# negative values are not accepted, so harmony batch correction is not suitable
# normalisation will be done using normalizeData() provided by CellChat
# first full pseudo scRNAseq dataset has to be normalised, and later subsetted to the cell types of interest
# https://github.com/sqjin/CellChat/issues/300
# https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html


# define parameters -------------------------------------------------------

# common params
scrna_anno <<- 'mid_lvl_ct_updated' #either 'cell_type' / 'mid_lvl_ct' / 'mid_lvl_ct_updated' / 'low_lvl_ct'
sample_name <- 'Sample'
grouping_var_col_ids <- c("Segment") # define the meta data column names of the groups that needed to be compared separately eg: c("Segment","NACT_status")

# parameters for CellChat and MultiNicheNet : Single cell approaches
# names of cells to fin
#cell_types_selected = c("Tcells_CD8","Macrophages_Monocytes") # set to NULL to get all the cell types : ct_of_interest
cell_types_selected <- NULL

# parameters for Cellchat
cell_frac_cutoff = 0.005 # 0.01 or 0.005
min_cells = 10 # Number of minimum cells in each cell group

# variables to merge the final csv with
meta_names <- c(aoi_id, roi_id, aoi_segment_var, sample_name, main_experimental_condition, 
                grouping_var_col_ids, main_batch_var, secondary_batch_var)


# path to raw BP results and cell fractions
bp_res_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', 
                         paste0('bp_res_', scrna_anno, '.RDS'))

bp_ct_frac_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                       paste0('bp_res_', scrna_anno, '_ct_fraction.csv'))

# path to the prepared normalized pseudo scRNaseq dataset from all bpres
bp_pseudosc_path <- file.path(output_dir, 'lr_interactions', 
                              paste0('bp_res_pseudosc_', scrna_anno, 'ct_frac_', cell_frac_cutoff, '_norm.csv'))

# output paths
outct <- ifelse(!is.null(cell_types_selected), paste(cell_types_selected, collapse = '_'), 'all')

# path to cellchat output file
geomx_CellChat_path <- file.path(output_dir, 'lr_interactions', 'cell_chat', 
                                 paste0('CellChat_output_',paste(grouping_var_col_ids, collapse = '_'),
                                        '_', outct, '.RDS'))

# path to output lr dataframe
cc_lr_df_path <- file.path(output_dir, 'lr_interactions', 'cell_chat', 
                           paste0('CellChat_df_',paste(grouping_var_col_ids, collapse = '_'),
                                  '_', outct, '_lr.csv'))

cc_path_df_path <- file.path(output_dir, 'lr_interactions', 'cell_chat', 
                             paste0('CellChat_df_', paste(grouping_var_col_ids, collapse = '_'),
                                    '_', outct, '_pathway.csv'))

# load geomx metadata -----------------------------------------------------

geomx_obj <<- readRDS(geomx_norm_batch_eff_rm_path) # batch effect corrected Geomx Object
meta_data_all <- sData(geomx_obj)[, unique(meta_names)]

# TODO code repetition from BSR Analysis
# make variable with all categories from grouping_var_col_ids
meta_data_all$comparison_group <- apply(meta_data_all, 1, function(row){
  group <- sapply(grouping_var_col_ids, function(var){
    paste(row[var])
  })
  group <- paste(group, collapse = '_')
  return(group)
})

rm(geomx_obj)
gc()

# create normalised pseudo scRNAseq dataset -------------------------------

# create a combined 'artificial pseudo-bulk scRNAseq' dataset with all ct specific counts 
if(file.exists(bp_pseudosc_path)){
  bprism_res_norm <- as.matrix(fread(bp_pseudosc_path), rownames=1)
} else{
  bprism_res_norm <- create_norm_pseudosc_from_deconv(bp_res_path, bp_ct_frac_path, scrna_anno, cell_frac_cutoff, bp_pseudosc_path)
}

# make metadata -----------------------------------------------------------

# make metadata with dcc_cell type
dcc_ct <- data.frame('dcc_filename' = gsub('_.*', '', colnames(bprism_res_norm)), 
                     'dcc_ct' = colnames(bprism_res_norm))
meta_data_ct <- left_join(dcc_ct, meta_data_all)
meta_data_ct$ct_label <- gsub('^[^_]*', '', meta_data_ct$dcc_ct)
meta_data_ct$ct_label <- gsub('^_', '', meta_data_ct$ct_label)
rownames(meta_data_ct) <- meta_data_ct$dcc_ct
meta_data_ct$samples <- meta_data_ct[[sample_name]]

# filter to ct of interest ------------------------------------------------

if(!is.null(cell_types_selected)) { 
  # filter to interesting cells
  meta_data_sel <- meta_data_ct[grepl(paste0(cell_types_selected, collapse = '|'), meta_data_ct$dcc_ct), ]
  bprism_res_sel <- bprism_res_norm[, meta_data_sel$dcc_ct]
} else{
  # If the cell types are not defined take all the cell types in the bprism object
  meta_data_sel <- meta_data_ct
  bprism_res_sel <- bprism_res_norm
}

rm(bprism_res_norm)
gc()

bprism_res_sel <- bprism_res_sel[, !grepl('Mast_cells', colnames(bprism_res_sel))]
meta_data_sel <- meta_data_sel[meta_data_sel$ct_label != 'Mast_cells', ]

# calculate cellchat probabilities ----------------------------------------

#  calculating the cellchat probabilties

cellchat_results <- list()

for(group in unique(meta_data_sel$comparison_group)) {
  print(group)
  
  # filtering to group
  meta_data_sel_group <- meta_data_sel[meta_data_sel$comparison_group == group, ]
  bprism_res_sel_group <- bprism_res_sel[, meta_data_sel_group$dcc_ct]
  
  # return error if not enough cells for comparison
  # TODO if too little cells in a group - remove this ct + omit it while plotting
  if(any(as.vector(table(meta_data_sel_group$ct_label)) < min_cells)){
    print(table(meta_data_sel_group$ct_label))
    next(paste0('nr of cells in comparison group smaller than min_cells'))
  }

  cellchat_results[[group]] <- cellchat_predict_prob(
    meta_data_sel_group,
    bprism_res_sel_group,
    thresh_fc = 0.1, 
    thresh_p = 0.05, 
    min_cells = min_cells
  )
}

saveRDS(cellchat_results, geomx_CellChat_path)


# combine all LR predictions to a single dataframe list -------------------

lr_df_all_groups <- lapply(names(cellchat_results), function(group){
  print(group)
  
  cellchat_obj <- cellchat_results[[group]]
  
  # calculate communication network
  df.net <- subsetCommunication(cellchat_obj, slot.name = "net") # L-R lvl
  df.path <- subsetCommunication(cellchat_obj, slot.name = "netP") # pathway lvl
  
  df.net$group <- group
  df.path$group <- group
  
  return(list(net = df.net, path = df.path))
})

lr_df_combined <- do.call(Map, c(f = rbind, lr_df_all_groups))

# save LR and pathway dfs
fwrite(lr_df_combined[[1]], file = cc_lr_df_path)
fwrite(lr_df_combined[[2]], file = cc_path_df_path)
