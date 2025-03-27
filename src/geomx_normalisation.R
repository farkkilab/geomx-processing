# README: script for normalising qc-ed raw data and making umap projections

# define variables --------------------------------------------------------

umap_vars <- c(aoi_segment_var, main_roi_label, main_experimental_condition, sample_name, 
               main_batch_var, batch_var, other_vars_bio, other_vars_tech)

exp_design <- as.formula(paste('~', aoi_segment_var, '+', main_experimental_condition))

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'umap_tsne', 'all'), showWarnings = T, recursive = T)

# load qc geomx data ------------------------------------------------------

geomx_obj <- readRDS(geomx_qc_path)

# Q3 normalisation --------------------------------------------------------

# this plot only makes sense for Q3 norm since it explores q3 value against NegGeoMean 
plot_q3_stats(geomx_obj, main_roi_label, file.path(output_dir,'qc', 'q3_stats.png'))

geomx_obj <- normalize(geomx_obj ,
                       norm_method = "quant", 
                       desiredQuantile = .75,
                       toElt = "q3_norm")

# DESeq2 normalisation ----------------------------------------------------

#change to integers
expr_int <- apply(geomx_obj@assayData$exprs, c(1, 2), function(x) {(as.integer(x))})

## Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = expr_int,
                              colData = sData(geomx_obj),
                              design = exp_design)

# normalise
dds <- estimateSizeFactors(dds)
deseq2_norm_counts <- counts(dds, normalized=TRUE)
dimnames(deseq2_norm_counts) = dimnames(geomx_obj@assayData$exprs)

# sizeFactors(dds)[1:10] # have a look at size factors

# do variance stabilising transformation - for PCA and other downstream analysis
# output is in log-like-scale !
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#dimensional-reduction
deseq2_vst <- vst(dds, blind = FALSE)
deseq2_vst_counts <- assay(deseq2_vst)
dimnames(deseq2_vst_counts) = dimnames(geomx_obj@assayData$exprs)

# scaling is better for PCA
deseq2_vst_scaled <- scale(deseq2_vst_counts)

# add norm matrices to geomx obj ------------------------------------------

# hacking GeoMx class object 
newassay <- new.env(parent=geomx_obj@assayData)
newassay$exprs <- geomx_obj@assayData$exprs
newassay$q3_norm <- geomx_obj@assayData$q3_norm
newassay$deseq2_norm <- deseq2_norm_counts
newassay$deseq2_vst <- deseq2_vst_counts
newassay$deseq2_vst_scaled <- deseq2_vst_scaled

geomx_obj@assayData <- newassay

# plot effects of normalisation -------------------------------------------

for(norm_type in c('exprs', 'q3_norm', 'deseq2_norm', 'deseq2_vst', 'deseq2_vst_scaled')){
  plt_title <- ifelse(norm_type == 'exprs', 'raw_counts', norm_type)
  islog <- ifelse(norm_type %in% c('deseq2_vst', 'deseq2_vst_scaled'), T, F)

  plot_norm_effect(assayDataElement(geomx_obj[,1:10], elt = norm_type),
                   plt_title, file.path(output_dir, 'qc', paste0('norm_', plt_title, '.png')),
                   is_log = islog)
  
  
  # plots with xlim = 0.99 percentile to rmv long tail
  plot_expr_distribution(geomx_obj@assayData[[norm_type]], plt_title, 
                         file.path(output_dir, 'qc', paste0('expr_hist_', plt_title, '.png')),
                                   is_log = islog)
}

# make UMAP and t-SNE -----------------------------------------------------

# TODO for segment, UMAP/tSNE should be done on mtx filtered for segment
# divide for segment and do dimentionality reduction for all
seg_types <- unique(sData(geomx_obj)[, aoi_segment_var])

geomx_obj_seg_list <- lapply(seg_types, function(seg){
  dir.create(file.path(output_dir, 'umap_tsne', seg), showWarnings = T, recursive = T)
  geomx_obj_seg <- geomx_obj[, geomx_obj@phenoData@data[[aoi_segment_var]] == seg]
  
  return(geomx_obj_seg)
})

names(geomx_obj_seg_list) <- seg_types
geomx_list <- c(all = geomx_obj, geomx_obj_seg_list)


norm <- 'deseq2_vst_scaled'

geomx_list_dim_red <- lapply(1:length(geomx_list)[1], function(n){
  geomx <- geomx_list[[n]]
  
  # run UMAP and tSNE on deseq2 vst counts
  geomx <- make_umap_tsne(geomx, norm, assay_is_log = T)

  # generate umap and tsne plots and color by variables
  for(method in c('UMAP', 'tSNE')){
    for(color_var in umap_vars){
      print(color_var)
      plot_umap_tsne(pData(geomx), method_type = method, 
                     norm_type = norm, color_var = color_var,
                     output_name = file.path(output_dir, 'umap_tsne', names(geomx_list)[n], 
                                             paste0(method, '_', norm, '_', color_var, '.pdf')))
    }
  }
  
  return(geomx)
})

# update objects
geomx_obj <- geomx_list_dim_red[[1]]

rm(geomx_list)
rm(geomx_list_dim_red)

# save geomx as RDS
saveRDS(geomx_obj, file = geomx_norm_path)
