# README: script for normalising qc-ed raw data and making umap projections

# define variables --------------------------------------------------------

umap_vars <- c(aoi_segment_var, main_roi_label, main_experimental_condition, sample_name, 
               main_batch_var, batch_var, other_vars_bio, other_vars_tech)

exp_design <- as.formula(paste('~', aoi_segment_var, '+', main_experimental_condition))

outliers_to_rm <- c('DSP-1001660037684-F-A12.dcc') # dcc filename of outliers to manually remove

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'qc', 'umap_tsne', 'all'), showWarnings = T, recursive = T)

# load qc geomx data ------------------------------------------------------

geomx_obj <- readRDS(geomx_qc_path)

# remove selected outliers after manual inspection (eg anormally low counts nr)
if(length(outliers_to_rm) > 0){
  geomx_obj <- geomx_obj[, which(!(colnames(geomx_obj) %in% outliers_to_rm))]
}

# Q3 normalisation --------------------------------------------------------

# this plot only makes sense for Q3 norm since it explores q3 value against NegGeoMean 
plot_q3_stats(geomx_obj, aoi_segment_var, file.path(output_dir,'qc', 'q3_stats.png'))

expr_q3_norm <- do_normalisation(geomx_obj@assayData$exprs, norm_type = 'q3')

# DESeq2 normalisation ----------------------------------------------------

expr_deseq2_norm <- do_normalisation(geomx_obj@assayData$exprs, norm_type = 'deseq2', meta_data = sData(geomx_obj), 
                       aoi_segment_var, main_experimental_condition)

# do variance stabilising transformation - for PCA and other downstream analysis
# output is in log-like-scale !
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#dimensional-reduction
expr_deseq2_vst_norm <- do_normalisation(geomx_obj@assayData$exprs, norm_type = 'deseq2_vst', meta_data = sData(geomx_obj), 
                                         aoi_segment_var, main_experimental_condition)

# scaling is better for PCA
expr_deseq2_vst_norm_scaled <- scale(expr_deseq2_vst_norm)

# add norm matrices to geomx obj ------------------------------------------

# hacking GeoMx class object 
newassay <- new.env(parent=geomx_obj@assayData)
newassay$exprs <- geomx_obj@assayData$exprs
newassay$q3_norm <- expr_q3_norm
newassay$deseq2_norm <- expr_deseq2_norm
newassay$deseq2_vst_norm <- expr_deseq2_vst_norm
newassay$deseq2_vst_norm_scaled <- expr_deseq2_vst_norm_scaled

geomx_obj@assayData <- newassay

# plot effects of normalisation -------------------------------------------

for(norm_type in c('exprs', 'q3_norm', 'deseq2_norm', 'deseq2_vst_norm', 'deseq2_vst_norm_scaled')) {
  plt_title <- ifelse(norm_type == 'exprs', 'raw_counts', norm_type)
  islog <- ifelse(norm_type %in% c('deseq2_vst_norm', 'deseq2_vst_norm_scaled'), T, F)

  # plot_norm_effect(assayDataElement(geomx_obj[,1:10], elt = norm_type),
  #                  plt_title, file.path(output_dir, 'qc', paste0('norm_', plt_title, '.png')),
  #                  is_log = islog)


  # plots with xlim = 0.99 percentile to rmv long tail
  plot_expr_distribution(geomx_obj@assayData[[norm_type]], plt_title, 
                         file.path(output_dir, 'qc', paste0('expr_hist_', plt_title, '.png')),
                                   is_log = islog)
}

# make UMAP and t-SNE -----------------------------------------------------

# divide for segment and do dimentionality reduction for all
seg_types <- unique(sData(geomx_obj)[, aoi_segment_var])

geomx_obj_seg_list <- lapply(seg_types, function(seg){
  dir.create(file.path(output_dir,'qc', 'umap_tsne', seg), showWarnings = T, recursive = T)
  
  geomx_obj_seg <- geomx_obj[, geomx_obj@phenoData@data[[aoi_segment_var]] == seg]
  
  return(geomx_obj_seg)
})

names(geomx_obj_seg_list) <- seg_types
geomx_list <- c(all = geomx_obj, geomx_obj_seg_list)


geomx_list_dim_red <- lapply(1:length(geomx_list), function(n){
  geomx <- geomx_list[[n]]

  norm_type <- 'deseq2_vst_norm_scaled' # just one to check befor batch effect rmv  
  # run UMAP and tSNE on norm counts
  geomx <- make_umap_tsne(geomx, norm_type, assay_is_log = T)
    
  # generate umap and tsne plots and color by variables
  for(method in c('UMAP', 'tSNE')){
    for(color_var in umap_vars){
      print(color_var)
      plot_umap_tsne(pData(geomx), method_type = method, 
                      norm_type = norm_type, color_var = color_var,
                      output_name = file.path(output_dir,'qc', 'umap_tsne', names(geomx_list)[n], 
                                              paste0(method, '_', norm_type, '_', color_var, '.pdf')))
    }
  }
  return(geomx)
})

# update objects - only keep whole geomx object
geomx_obj <- geomx_list_dim_red[[1]]

rm(geomx_list)
rm(geomx_list_dim_red)

# save geomx as RDS
saveRDS(geomx_obj, file = geomx_norm_path)
