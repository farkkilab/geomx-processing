# README: script for removing batch effect using limma and/or harmony
# nice explanation
# https://www.biostars.org/p/366403/

# define variables --------------------------------------------------------

# main experimental conditions
exp_design <- as.formula(paste('~', aoi_segment_var, '+', main_experimental_condition))
exp_design <- as.formula(paste('~', 'tls_status'))

# all variables to check for variance
batch_vars <- c(primary_batch_var, secondary_batch_var, other_vars_tech, 
                aoi_segment_var, main_roi_label, sample_name, main_experimental_condition, 
                other_vars_bio)

batch_vars <- c(primary_batch_var, secondary_batch_var, other_vars_tech, 
                sample_name,  "Segment_geomx", "Patient", 'tls_status')

# normalisation used for batch effect correction calculation
norm_type <- 'q3_norm' # best to use deseq2_vst data, eventually deseq2_norm


# whether or not compute pvca - it takes awful amount of time 
# and is needed only 1nce in a given batch
calculate_pvca <- TRUE
# PVCA threshold
pct_threshold <- 0.6 

# pre-umap filtering params (if no filtering set to NULL)
top_var <- 2000 # filter to top variable genes
top_pca <- 50 # do PCA and filter to top components

# biological covariates which effect should be ignored by limma 
# if NULL no cov are added to limma rmv batch eff
# TODO check if this is beneficial 
# cov_design <- formula(~ Patient + Site) 
cov_design <- NULL

# make dirs and set additional vars ---------------------------------------

dir.create(file.path(output_dir, 'batch_correction'), showWarnings = T, recursive = T)

norm_is_log <- ifelse(norm_type %in% c('exprs', 'q3_norm', 'deseq2_norm'), FALSE, TRUE)
covname <- ifelse(is.null(cov_design), 'no', gsub(' ', '', as.character(cov_design)[2]))

# load geomx object and create expression set -----------------------------

geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)

# change vars into factors and remove vars with only 1 unique value
batch_vars_filt <- c()
for(colname in batch_vars){
  if(length(unique(pData(geomx_obj)[[colname]])) > 1){
    pData(geomx_obj)[[paste0(colname, '_factor')]] <- as.factor(pData(geomx_obj)[[colname]])
    batch_vars_filt <- c(batch_vars_filt, colname)
  }
}
batch_factors_names <- paste0(batch_vars_filt, '_factor')

# make expression sets for PVCA
phenoData <- new("AnnotatedDataFrame", data=geomx_obj@phenoData@data, 
                 varMetadata=geomx_obj@phenoData@varMetadata)
featureData <- new("AnnotatedDataFrame", data=geomx_obj@featureData@data, 
                 varMetadata=geomx_obj@featureData@varMetadata)

exprset_norm <- ExpressionSet(assayData=geomx_obj@assayData[[norm_type]], 
                              phenoData = phenoData,
                              featureData = featureData)

if(!norm_is_log){
  # make log2 transformed normalised counts if norm_type not in log scale
  expr_norm_log <- log2(geomx_obj@assayData[[norm_type]] + 1)
} else{
  expr_norm_log <- geomx_obj@assayData[[norm_type]]
}

# check initial batch effect with PVCA ------------------------------------

if(calculate_pvca){
  pvcaObj_ini <- pvcaBatchAssess(exprset_norm, batch_factors_names, pct_threshold) 
  
  plot_pvca(pvcaObj_ini, paste0('before_correction_', norm_type), file.path(output_dir, 'batch_correction'))
}


# remove batch effect with limma ------------------------------------------

design <- model.matrix(exp_design, data = sData(geomx_obj))
prim_batch <-  sData(geomx_obj)[[primary_batch_var]]

if(!is.null(secondary_batch_var)){
  second_batch <- sData(geomx_obj)[[secondary_batch_var]]
}else{second_batch <- NULL} 

if(!is.null(cov_design)){
  cov <- model.matrix(cov_design, data = sData(geomx_obj))
  }else{cov <- NULL} 

limma_res <- limma::removeBatchEffect(expr_norm_log, batch = prim_batch, batch2 = second_batch,
                                      covariates = cov, design = design)

##########################
# https://portals.broadinstitute.org/harmony/articles/quickstart.html
# remove batch effect with harmony
meta_dt <- pData(geomx_obj)[, unique(c(batch_vars_filt, aoi_segment_var, aoi_id))]

harmony_res <- t(HarmonyMatrix(expr_norm_log, 
                               meta_data = meta_dt,
                               vars_use = c(primary_batch_var, secondary_batch_var)))

# additional scaling and PCA before harmony -------------------------------
# https://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/doc/detailedWalkthrough.html
# TODO for now we want to keep all genes, not PCA
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


if(calculate_pvca){
  pvcaObj_limma <- pvcaBatchAssess(exprset_after_limma, batch_factors_names, pct_threshold) 
  pvcaObj_harmony <- pvcaBatchAssess(exprset_after_harmony, batch_factors_names, pct_threshold) 
  
  plot_pvca(pvcaObj_limma, paste0('after_correction_limma_', primary_batch_var, secondary_batch_var,
                                  '_cov_', covname), file.path(output_dir, 'batch_correction'))
  
  plot_pvca(pvcaObj_harmony, paste0('after_correction_harmony_', primary_batch_var, secondary_batch_var), 
            file.path(output_dir, 'batch_correction'))
}

# add limma and harmony res to umap object --------------------------------

geomx_obj@assayData[[paste0('limma_batch_corr_', norm_type)]] <- limma_res
geomx_obj@assayData[[paste0('harmony_batch_corr_', norm_type)]] <- harmony_res

# plot counts distribution ------------------------------------------------

plot_expr_distribution(geomx_obj@assayData[[paste0('limma_batch_corr_', norm_type)]], paste0('limma_batch_corr_', norm_type), 
                       file.path(output_dir, 'batch_correction', 
                                 paste0('expr_hist_limma_batch_corr_', 
                                        primary_batch_var, secondary_batch_var,
                                        '_cov_', covname, '.png')), is_log = T)


plot_expr_distribution(geomx_obj@assayData[[paste0('harmony_batch_corr_', norm_type)]], paste0('harmony_batch_corr_', norm_type), 
                       file.path(output_dir, 'batch_correction',
                                 paste0('expr_hist_harmony_batch_corr_', 
                                        primary_batch_var, secondary_batch_var,
                                        '_cov_', covname, '.png')), is_log = T)

# make UMAP and visualise batch-corrected results -------------------------
# TODO simplify code (as in sanity_check)
# divide for segment and do dimentionality reduction for all

# make separate geomx obj for all + each segment
seg_types <- unique(sData(geomx_obj)[, aoi_segment_var])

geomx_obj_seg_list <- lapply(seg_types, function(seg){
  dir.create(file.path(output_dir, 'batch_correction', seg), showWarnings = T, recursive = T)
  geomx_obj_seg <- geomx_obj[, geomx_obj@phenoData@data[[aoi_segment_var]] == seg]
  
  return(geomx_obj_seg)
})

names(geomx_obj_seg_list) <- seg_types
geomx_list <- c(all = geomx_obj, geomx_obj_seg_list)
dir.create(file.path(output_dir, 'batch_correction', 'all'), showWarnings = T, recursive = T)

# iterate through all objects 
geomx_list_dim_red <- lapply(1:length(geomx_list), function(n){
  geomx <- geomx_list[[n]]
  
  # run UMAP and tSNE on limma and harmony batch effect correction
  geomx <- make_umap_tsne(geomx, paste0('limma_batch_corr_', norm_type), assay_is_log = T, top_var = top_var, top_PCA = top_pca)
  geomx <- make_umap_tsne(geomx, paste0('harmony_batch_corr_', norm_type), assay_is_log = T, top_var = top_var, top_PCA = top_pca)
  
  # generate umap and tsne plots and color by variables
  for(corr_type in c(paste0('limma_batch_corr_', norm_type), paste0('harmony_batch_corr_', norm_type))){
    for(method in c('UMAP', 'tSNE')){
      for(color_var in batch_vars_filt){
        print(color_var)
        
        plot_umap_tsne(pData(geomx), method_type = method, 
                       norm_type = corr_type, color_var = color_var,
                       output_name = file.path(output_dir, 'batch_correction', names(geomx_list)[n], 
                                               paste0(method, '_', corr_type, '_', color_var, 
                                                      '_topvargenes_', ifelse(is.null(top_var), 'NULL', as.character(top_var)),
                                                      '_toppca_', ifelse(is.null(top_var), 'NULL', as.character(top_pca)), '.pdf')))
      }
    }
  }
  
  return(geomx)
})

# update objects - only keep whole geomx object
geomx_obj <- geomx_list_dim_red[[1]]

rm(geomx_list)
rm(geomx_list_dim_red)



# save geomx_obj with batch eff correction --------------------------------

saveRDS(geomx_obj, file = geomx_norm_batch_eff_rm_path)

# save logs
writeLines(c('batch effect rmv logs:',
             '; normalisation type : ', norm_type,
             '; primary batch effect variable : ', primary_batch_var,
             '; secondary batch effect variable : ', secondary_batch_var,
             '; limma experimental design : ', as.character(exp_design)[2],
             '; limma covariate : ', as.character(cov_design)[2]), 
           file.path(output_dir, 'batch_correction', 'batch_correction_logs.txt'))
