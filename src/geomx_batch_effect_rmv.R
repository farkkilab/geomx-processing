# nice explanation
# https://www.biostars.org/p/366403/

# define variables --------------------------------------------------------

# main experimental conditions
exp_design <- formula(~ Segment + NACT_status)

# biological covariates which effect should be ignored by limma 
# if NULL no cov are added to limma rmv batch eff
# TODO check if this is beneficial 
cov_design <- formula(~ Patient + Site) 

# all variables to check for variance
tech_vars <- c('Slide_Name', 'batch_nr') # 'batch_nr_sample_collection'
bio_vars <- c('Patient', 'Site', 'Sample', 'Annotation_cell', 'NACT_status', 'Segment')
batch_vars <- c(tech_vars, bio_vars)

aoi_segment_var <- 'Segment' # not included in batch_vars bcs it will take most of the variance

dir.create(file.path(output_dir, 'batch_correction'), showWarnings = T, recursive = T)

# main cause of the batch effect, from 1st PVCA plot
main_batch_var <- 'batch_nr'
secondary_batch_var <- NULL # it has to be INDEPENDENT from the main_batch_var

# normalisation used for batch effect correction calculation
norm_type <- 'deseq2_vst' # best to use vst data, eventually deseq2_norm
norm_is_log <- TRUE # vst is already in the log-like scale, deseq2_norm not

# PVCA threshold
pct_threshold <- 0.6 

dir.create(file.path(output_dir, 'batch_correction'), recursive = T, showWarnings = F)

# load geomx object and create expression set -----------------------------

geomx_obj <- readRDS(geomx_norm_path)

# change vars into factors
for(colname in batch_vars){
  pData(geomx_obj)[[paste0(colname, '_factor')]] <- as.factor(pData(geomx_obj)[[colname]])
}
batch_factors_names <- paste0(batch_vars, '_factor')

# make expression sets for PVCA
phenoData <- new("AnnotatedDataFrame", data=geomx_obj@phenoData@data, 
                 varMetadata=geomx_obj@phenoData@varMetadata)
featureData <- new("AnnotatedDataFrame", data=geomx_obj@featureData@data, 
                 varMetadata=geomx_obj@featureData@varMetadata)

# !! by default deseq2_norm is used to check for initial batch effect by pvca
exprset_deseq2_norm <- ExpressionSet(assayData=geomx_obj@assayData$deseq2_norm, 
                              phenoData = phenoData,
                              featureData = featureData)

if(!norm_is_log){
  # make log2 transformed normalised counts if norm_type not in log scale
  expr_norm_log <- log2(geomx_obj@assayData[[norm_type]] + 1)
} else{
  expr_norm_log <- geomx_obj@assayData[[norm_type]]
}

# check initial batch effect with PVCA ------------------------------------

pvcaObj_ini <- pvcaBatchAssess(exprset_deseq2_norm, batch_factors_names, pct_threshold) 

plot_pvca(pvcaObj_ini, 'before_correction_deseq2_norm', file.path(output_dir, 'batch_correction'))

# remove batch effect with limma ------------------------------------------

design <- model.matrix(exp_design, data = sData(geomx_obj))
batch <-  sData(geomx_obj)[[main_batch_var]]

batch2 <- ifelse((!is.null(secondary_batch_var)), 
                 sData(geomx_obj)[[secondary_batch_var]], NULL)

cov <- ifelse((!is.null(cov_design)), 
              model.matrix(cov_design, data = sData(geomx_obj)), NULL)


# if there are 1 main variable responsible for batcheffect
limma_res <- limma::removeBatchEffect(expr_norm_log, batch = batch, batch2 = batch2,
                                      covariates = cov, design = design)

##########################
# https://portals.broadinstitute.org/harmony/articles/quickstart.html
# remove batch effect with harmony
meta_dt <- pData(geomx_obj)[, c(batch_vars, aoi_segment_var, 'dcc_filename')]

harmony_res <- t(HarmonyMatrix(expr_norm_log, 
                               meta_data = meta_dt,
                               vars_use = c(main_batch_var, secondary_batch_var)))

# additional scaling and PCA before harmony -------------------------------
# https://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/doc/detailedWalkthrough.html
# TODO not so clear how to use it later - so far we'll stay with just norm-log data
# may be important in the future for clustering etc
# expr_norm_log_scaled <- scale(expr_norm_log)
# expr_norm_log_scaled_pca <- prcomp(t(expr_norm_log_scaled))
# expr_norm_log_scaled_pca_top20 <- t(expr_norm_log_scaled_pca$x)[1:20, ]
# 
# 
# harmony_res_scaledpca <- t(HarmonyMatrix(expr_norm_log_scaled_pca_top20, 
#                                          meta_data = meta_dt,
#                                          vars_use = main_batch_var))



# check PVCA after batch effect removal -----------------------------------

exprset_after_limma <- ExpressionSet(assayData=limma_res, 
                                         phenoData = phenoData,
                                         featureData = featureData)

exprset_after_harmony <- ExpressionSet(assayData=harmony_res, 
                                     phenoData = phenoData,
                                     featureData = featureData)


pvcaObj_limma <- pvcaBatchAssess(exprset_after_limma, batch_factors_names, pct_threshold) 
pvcaObj_harmony <- pvcaBatchAssess(exprset_after_harmony, batch_factors_names, pct_threshold) 

plot_pvca(pvcaObj_limma, paste0('after_correction_limma_', cov_name,  main_batch_var, secondary_batch_var), 
          file.path(output_dir, 'batch_correction'))

plot_pvca(pvcaObj_harmony, paste0('after_correction_harmony_', main_batch_var, secondary_batch_var), 
          file.path(output_dir, 'batch_correction'))

# add limma and harmony res to umap object --------------------------------

geomx_obj@assayData$limma_batch_corr <- limma_res
geomx_obj@assayData$harmony_batch_corr <- harmony_res

# make UMAP and visualise batch-corrected results -------------------------

geomx_obj <- make_umap_tsne(geomx_obj, 'limma_batch_corr', assay_is_log = T)
geomx_obj <- make_umap_tsne(geomx_obj, 'harmony_batch_corr', assay_is_log = T)

# generate umap and tsne plots and color by variables
for(corr_type in c('limma_batch_corr', 'harmony_batch_corr')){
  for(method in c('UMAP', 'tSNE')){
    for(color_var in batch_vars){
      print(color_var)
      plot_umap_tsne(pData(geomx_obj), method_type = method, 
                     norm_type = corr_type, color_var = color_var,
                     output_name = file.path(output_dir, 'batch_correction', 
                                             paste0(method, '_', corr_type, '_',  color_var, '.pdf')))
    }
  }
}

# save geomx_obj with batch eff correction --------------------------------

saveRDS(geomx_obj, file = geomx_norm_batch_eff_rm_path)
