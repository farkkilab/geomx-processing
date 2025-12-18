# README: script for removing batch effect using limma and/or harmony
# nice explanation
# https://www.biostars.org/p/366403/

# define variables --------------------------------------------------------

# main experimental conditions
exp_design <- as.formula(paste('~', aoi_segment_var, '+', main_experimental_condition))
#exp_design <- as.formula(paste('~', 'tls_status'))

# all variables to check for variance
batch_vars <- unique(c(primary_batch_var, secondary_batch_var, other_vars_tech, 
                aoi_segment_var, main_roi_label, sample_name, main_experimental_condition, 
                other_vars_bio))

# batch_vars <- c(primary_batch_var, secondary_batch_var, other_vars_tech, 
#                 sample_name,  "Segment_geomx", "Patient", 'tls_status')

# normalisation used for batch effect correction calculation best: 'deseq2_vst_norm' / 'q3_norm'
norm_type <- 'deseq2_vst_norm' 
batch_rm_type <- 'harmony' # either 'limma' or 'harmony'

# whether or not compute pvca - it takes awful amount of time
# and is needed only 1nce in a given batch
calculate_pvca <- FALSE
# PVCA threshold
pct_threshold <- 0.6 

# pre-umap filtering params (if no filtering set to NULL)
top_var <- 2000 # filter to top variable genes
top_pca <- 50 # do PCA and filter to top components

# make dirs and set additional vars ---------------------------------------

dir.create(file.path(output_dir, 'batch_correction'), showWarnings = T, recursive = T)

norm_is_log <- ifelse(norm_type %in% c('exprs', 'q3_norm', 'deseq2_norm'), FALSE, TRUE)

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

if(!norm_is_log){
  # make log2 transformed normalised counts if norm_type not in log scale
  expr_norm_log <- log2(geomx_obj@assayData[[norm_type]] + 1)
} else{
  expr_norm_log <- geomx_obj@assayData[[norm_type]]
}

# make phenoData and featureData in correct class for PVCA and limma 
phenoData <- new("AnnotatedDataFrame", data=geomx_obj@phenoData@data, 
                 varMetadata=geomx_obj@phenoData@varMetadata)
featureData <- new("AnnotatedDataFrame", data=geomx_obj@featureData@data, 
                   varMetadata=geomx_obj@featureData@varMetadata)

# check initial batch effect with PVCA ------------------------------------

if(calculate_pvca){
  phenoData <- as.data.frame(geomx_obj@phenoData@data)
  rownames(phenoData) <- phenoData[[aoi_id]]
  
  # make expression sets for PVCA
  exprset_norm = ExpressionSet(assayData=geomx_obj@assayData[[norm_type]],
                                phenoData = AnnotatedDataFrame(phenoData),
                                featureData = AnnotatedDataFrame(as.data.frame(geomx_obj@featureData@data)))
  
  pvcaObj_ini <- pvcaBatchAssess(exprset_norm, batch_factors_names, pct_threshold) 
  
  plot_pvca(pvcaObj_ini, paste0('before_correction_', norm_type), file.path(output_dir, 'batch_correction'))
}


# remove batch effect with limma or harmony -------------------------------

if(batch_rm_type == 'limma'){
  # remove batch effect with limma
  design <- model.matrix(exp_design, data = sData(geomx_obj))
  prim_batch <-  sData(geomx_obj)[[primary_batch_var]]
  
  if(!is.null(secondary_batch_var)){
    second_batch <- sData(geomx_obj)[[secondary_batch_var]]
  }else{second_batch <- NULL} 
  
  batch_rm_res <- limma::removeBatchEffect(expr_norm_log, batch = prim_batch, batch2 = second_batch,
                                        design = design)
  
} else if(batch_rm_type == 'harmony'){
  # remove batch effect with harmony
  # https://portals.broadinstitute.org/harmony/articles/quickstart.html
  meta_dt <- pData(geomx_obj)[, unique(c(batch_vars_filt, aoi_segment_var, aoi_id))]
  
  batch_rm_res <- t(HarmonyMatrix(expr_norm_log, 
                                 meta_data = meta_dt,
                                 vars_use = c(primary_batch_var, secondary_batch_var)))
} else{
  stop("batch_rm_type can be either limma or harmony")
}

# add batch effect rm res to umap object --------------------------------

geomx_obj@assayData[[paste0(batch_rm_type, '_', norm_type)]] <- batch_rm_res

# plot counts distribution ------------------------------------------------

plot_expr_distribution(geomx_obj@assayData[[paste0(batch_rm_type, '_', norm_type)]], paste0(batch_rm_type, '_', norm_type), 
                       file.path(output_dir, 'batch_correction', 
                                 paste0('expr_hist_', batch_rm_type, '_', norm_type, '.png')), is_log = T)

# check PVCA after batch effect removal -----------------------------------

if(calculate_pvca){
  exprset_after_batch_rm <- ExpressionSet(assayData=batch_rm_res,
                                       phenoData = AnnotatedDataFrame(phenoData),
                                       featureData = AnnotatedDataFrame(as.data.frame(geomx_obj@featureData@data)))

  pvcaObj_after_batch_rm <- pvcaBatchAssess(exprset_after_batch_rm, batch_factors_names, pct_threshold) 
  
  plot_pvca(pvcaObj_after_batch_rm, paste0('after_correction_', batch_rm_type, '_', norm_type), 
            file.path(output_dir, 'batch_correction'))
}


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
  
  # run UMAP and tSNE on batch effect correction
  geomx <- make_umap_tsne(geomx, paste0(batch_rm_type, '_', norm_type), assay_is_log = T, top_var = top_var, top_PCA = top_pca)

  # generate umap and tsne plots and color by variables
  for(method in c('UMAP', 'tSNE')){
    for(color_var in batch_vars_filt){
      print(color_var)
      
      plot_umap_tsne(pData(geomx), method_type = method, 
                     assay_name = paste0(batch_rm_type, '_', norm_type), color_var = color_var, shape_var = aoi_segment_var,
                     output_name = file.path(output_dir, 'batch_correction', names(geomx_list)[n], 
                                             paste0(method, '_', batch_rm_type, '_', norm_type, '_', color_var, 
                                                    '_topvargenes_', ifelse(is.null(top_var), 'NULL', as.character(top_var)),
                                                    '_toppca_', ifelse(is.null(top_var), 'NULL', as.character(top_pca)), '.pdf')))
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
             '; batch effect removal type : ', batch_rm_type,
             '; primary batch effect variable : ', primary_batch_var,
             '; secondary batch effect variable : ', secondary_batch_var,
             '; limma experimental design : ', as.character(exp_design)[2]), 
           file.path(output_dir, 'batch_correction', 'batch_correction_logs.txt'))
