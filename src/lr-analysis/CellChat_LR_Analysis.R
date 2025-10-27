# Ligand Receptor Analysis  by CellChat : Geomx Deconvoluted data

# TODO If All are NA (nact_status ,segment , annotation) error  : Run for the whole dataset 


# define parameters -------------------------------------------------------

scrna_anno <<- 'mid_lvl_ct_updated' #either 'cell_type' / 'mid_lvl_ct' / 'mid_lvl_ct_updated' / 'low_lvl_ct'

# common parameters
grouping_var_col_ids <- c("Segment") # define the meta data column names of the groups that needed to be compared separately eg: c("Segment","NACT_status")

# TODO BSR - handled with code - check for NN and CC
# define the groups from  "grouping_var_col_ids" that needed to be compared eg: c("pre_stroma","pre_tumor") order matters. 
# Can compare only two groups at a time
# comparison <- c("stroma","tumor") 

# parameters for CellChat and MultiNicheNet : Single cell approaches
# names of cells to fin
cell_types_selected = c("Tcells_CD8","Macrophages_Monocytes") # set to NULL to get all the cell types : ct_of_interest

# parameters for Cellchat
cell_frac_cutoff = 0.005 # 0.01 or 0.005
min_cells = 10 # Number of minimum cells in each cell group

# normalised + in log form (best after batch effect correction - then always in log form)
# data_type <- "harmony_batch_corr_deseq2_vst"  

# variables to merge the final csv with
meta_names <- c(aoi_id, roi_id, aoi_segment_var, sample_name, main_experimental_condition, 
                grouping_var_col_ids, main_batch_var, secondary_batch_var)

# TODOthis to rm
# bp_res_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', 
#                          paste0('bp_res_', scrna_anno, '_expr_mtx_cleaned_deseq2_vst_harmony_corr.RDS'))

# path to raw BP results and cell fractions
bp_res_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', 
                         paste0('bp_res_', scrna_anno, '.RDS'))

bp_ct_frac_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                       paste0('bp_res_', scrna_anno, '_ct_fraction.csv'))


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

# load and filter raw bayesprism results ----------------------------------

bprism_res <<- readRDS(bp_res_path) # normalised and batch corrected bp results
ct_names <- colnames(get.fraction (bp=bprism_res, which.theta="final", state.or.type="type"))

# make metadata with dcc_cell type 
dcc_ct <- data.frame('dcc_filename' = gsub('_.*', '', colnames(bprism_res_norm)), 'dcc_ct' = colnames(bprism_res_norm))
meta_data_ct <- left_join(dcc_ct, meta_data_all)

cell_fractions_df <<- read.csv(bp_ct_frac_path) # cell fractions from BayesPrism

# create a combined 'artificial pseudo-bulk scRNAseq' ---------------------
#TODO make a separate function for cleaning + normalising + batch corr
# create a combined 'artificial pseudo-bulk scRNAseq' dataset with all ct specific counts 

# filtering out aois with very low cell fraction and combining into 1 mtx
bprism_res_filtered <- lapply(ct_names, function(ct_name){
  
  bprism_ct <- t(BayesPrism::get.exp(bp=bprism_res,
                                   state.or.type="type",
                                   cell.name=ct_name))
  
  ct_fraq <- cell_fractions_df[!is.na(cell_fractions_df[[ct_name]]), c('dcc_filename', ct_name)]
  aoi_with_ct <- ct_fraq$dcc_filename[ct_fraq[[ct_name]] >= cell_frac_cutoff]
  
  bprism_ct_filt <- bprism_ct[, colnames(bprism_ct) %in% aoi_with_ct]
  colnames(bprism_ct_filt) <- paste0(colnames(bprism_ct_filt), '_', ct_name)
  
  return(bprism_ct_filt)
})

cbind.fill <- function(df_list){
  nm <- lapply(df_list, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

bprism_res_filtered <- cbind.fill(bprism_res_filtered)

dim(bprism_res_filtered)

rm(bprism_res)
gc()

# normalisation 
# from Cellchat vignette: (e.g., library-size normalization and then log-transformed with a pseudocount of 1)
# for our purpose: deseq2+vst

# remove genes with only 0 counts
bprism_res_filtered <- bprism_res_filtered[rowSums(bprism_res_filtered) != 0, ]

#add pseudocount 1 to avoid vst error with log geo means
# https://help.galaxyproject.org/t/error-with-deseq2-every-gene-contains-at-least-one-zero/564/2
bprism_res_filtered <- bprism_res_filtered + 1

# do vst normalisation
bprism_res_norm <- varianceStabilizingTransformation(round(bprism_res_filtered))

fwrite(bprism_res_norm, file.path(output_dir, 'lr_interactions', 
                             paste0('bp_res_pseudosc_', scrna_anno, 'ct_frac_', cell_frac_cutoff, '_deseq2_vst.csv')))


# do batch effect correction with harmony
bprism_res_norm_batch_rm <- t(HarmonyMatrix(bprism_res_norm, 
                                   meta_data = meta_data_ct,
                                   vars_use = c(primary_batch_var, secondary_batch_var)))

fwrite(bprism_res_norm_batch_rm, file.path(output_dir, 'lr_interactions', 
                                  paste0('bp_res_pseudosc_', scrna_anno, 'ct_frac_', cell_frac_cutoff, '_deseq2_vst_harmony.csv')))

# plot distributions
plot_expr_distribution(bprism_res_filtered, paste('pseudo scRNAseq from deconv ct raw'), 
                       file.path(output_dir, 'lr_interactions', 
                                 paste0('expr_hist_bp_res_pseudosc_', scrna_anno, 'ct_frac_', cell_frac_cutoff, '_raw.png')), is_log = F)

plot_expr_distribution(bprism_res_norm, paste('pseudo scRNAseq from deconv ct deseq2 vst norm'), 
                       file.path(output_dir, 'lr_interactions', 
                                 paste0('expr_hist_bp_res_pseudosc_', scrna_anno, 'ct_frac_', cell_frac_cutoff, '_deseq2_vst.png')), is_log = T)

plot_expr_distribution(bprism_res_norm_batch_rm, paste('pseudo scRNAseq from deconv ct deseq2 vst norm harmony batch corr'), 
                       file.path(output_dir, 'lr_interactions', 
                                 paste0('expr_hist_bp_res_pseudosc_', scrna_anno, 'ct_frac_', cell_frac_cutoff, '_deseq2_vst_harmony.png')), is_log = T)

rm(bprism_res_filtered)
rm(bprism_res_norm)
gc()



# filter to ct of interest ------------------------------------------------


if(!is.null(cell_types_selected)) { 
  # filter to interesting cells
  meta_data_sel <- meta_data_ct[grepl(paste0(cell_types_selected, collapse = '|'), meta_data_ct$dcc_ct), ]
  bprism_res_sel <- bprism_res_norm_batch_rm[, meta_data_sel$dcc_ct]
} else{
  # If the cell types are not defined take all the cell types in the bprism object
  meta_data_sel <- meta_data_ct
  bprism_res_sel <- bprism_res_norm_batch_rm
}


# calculate cellchat probabilities ----------------------------------------

#  calculating the cellchat probabilties

group <- 'stroma'

cellchat_results <- list()

for(group in unique(meta_data_all$comparison_group)) {
  
  print(group)
  meta_data_group <- meta_data_all[meta_data_all$comparison_group == group, ]
  
  # filtering to group
  bprism_res_filt_group <- lapply(ct_names, function(ct_name){
    bprism_ct <- bprism_res_filtered[[ct_name]]
    bprism_ct_group <- bprism_ct[, colnames(bprism_ct) %in% meta_data_group[[aoi_id]]]
    return(bprism_ct_group)
  })
  
  names(bprism_res_filt_group) <- ct_names
  
  # TODO cell chat expects 1 df with all ct specific expr together as single cell !!
  # TODO check distribution of such df and decide if renormalisation + batch corr is needed..
  
  cellchat_results[[group]] <- cellchat_predict_prob(
    meta_data_all,
    expr_deseq2_norm_log_cf_filtered,
    sample_name = sample_name,
    grouping_var_col_ids = grouping_var_col_ids,
    group = group,
    output_dir =  output_dir
  )
}


#########################################3
#########################################

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



####################################################################
####################################################################
#  Extracting expression data and meta data of the desired cell types and combining 

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

  
