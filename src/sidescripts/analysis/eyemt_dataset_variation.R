library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library(GeomxTools)
library(pvca)

# define paths ------------------------------------------------------------

proj_dir <<- '~/Documents/phd/st'
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808')

metadt_path <- file.path(output_dir, 'metadata_full_SENSITIVE.csv')

# expression files
geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS')
deconv_bp_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', 'bp_res_mid_lvl_ct_updated_expr_mtx_cleaned_deseq2_vst_harmony_corr.RDS')

signal <- 'Macrophages_Monocytes' # either 'all' or deconv ct

pct_threshold <- 0.6 # threshold for PVCA
top_var <- 2000 # filter to top variable genes
top_pca <- 50 # do PCA and filter to top components

source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))
out_dir <- file.path(output_dir, 'downstream', 'variation', signal)
dir.create(out_dir)

# define variables --------------------------------------------------------

# TODO  myeloids, lymphoids?
ct_of_interest <- c("Tcells_CD8","Tcells_CD4", "Macrophages_Monocytes", "DCs", "Fibroblasts_Mesothelial", "tumor")
ct_of_interest_immune <- c("Tcells_CD8","Tcells_CD4", "Macrophages_Monocytes", "DCs", "other_immune")

pt_id <- 'Patient'
sample_id <- 'Sample'
roi_id <- 'sample_roi'
aoi_id <- 'dcc_filename'
aoi_segment_var <- 'Segment'

vars_id <- c(pt_id, sample_id, roi_id, aoi_id)

# vars names from each categories, discr and cont separately
vars_survival_cont <- c('OS_days', 'PFS_days')

vars_pt_disc <- c('stage', 'primary_surgery_residual', 'primary_treatment_response',
                   'treatment_bevacizumab_1st_line', 'treatment_PARPi', 
                  'HRP_status', 'BRCA_status') # 'PARPi_line',
vars_pt_cont <- c('age_at_diagnosis', 'ovaHRDscar_score', 'TMB')

vars_sample_disc <- c('Site', 'NACT_status')

vars_roi_disc <- c('Segment_geomx', 'roi_cluster_label') #'tls_status'
vars_roi_cont <- c(paste0('ct_frac_sd_roi_', ct_of_interest), paste0('ct_immunefrac_sd_roi_', ct_of_interest_immune))

vars_aoi_disc <- c('Segment')
vars_aoi_cont <- paste0('ct_frac_sd_aoi_', ct_of_interest)

vars_technical_disc <- c('main_batch_nr')

vars_all_disc <- c(vars_pt_disc, vars_sample_disc, vars_roi_disc, vars_aoi_disc, vars_technical_disc)
vars_all_cont <- c(vars_survival_cont, vars_pt_cont, vars_roi_cont, vars_aoi_cont)
  
# load data ---------------------------------------------------------------

geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path) # needed for feature data

if(signal == 'all'){
  expr_mtx <- geomx_obj@assayData$harmony_batch_corr_q3_norm
} else{
  expr_mtx <- readRDS(deconv_bp_path)[[signal]]
}

metadt <- fread(metadt_path, select = c(vars_id, vars_all_disc, vars_all_cont)) %>%
  filter(dcc_filename %in% colnames(expr_mtx)) %>%
  as.data.frame()

metadt[metadt == ''] <- NA

# do PVCA -----------------------------------------------------------------

# create pheno and feature data
phenoData <- metadt[, c(vars_id, vars_all_disc)] %>% # only disrete vars can be used in pvca
  select(-sample_roi) %>% 
  mutate_all( as.factor) %>%
  column_to_rownames('dcc_filename')
phenoData <- phenoData[match(colnames(expr_mtx), rownames(phenoData)),]

featureData <- as.data.frame(geomx_obj@featureData@data) 
featureData <- featureData[match(rownames(expr_mtx), rownames(featureData)),]

# make expression sets for PVCA
exprset = ExpressionSet(assayData=expr_mtx,
                             phenoData = AnnotatedDataFrame(phenoData),
                             featureData = AnnotatedDataFrame(featureData))

pvcaObj <- pvcaBatchAssess(exprset, colnames(phenoData), pct_threshold) 

plot_pvca(pvcaObj, plot_name = signal, output_dir = out_dir)

saveRDS(pvcaObj, file.path(out_dir, paste0('pvca_obj_', signal, '.RDS')))


# do PCA + UMAP + tSNE coloured by vars -----------------------------------

# get deconv df and filter metadata
ct_name <- signal
deconv_ct <- expr_mtx
metadt_ct <- metadt[metadt$dcc_filename %in% colnames(deconv_ct),]

# iterate through all + different segments
seg_types <- c('all', unique(metadt[, aoi_segment_var]))

sapply(seg_types, function(seg){

  print(seg)
  dir.create(file.path(out_dir, paste0('deconv_umap_tsne_', seg)), showWarnings = T, recursive = T)

  if(seg != 'all'){
    deconv_seg <- deconv_ct[, metadt_ct$dcc_filename[metadt_ct[[aoi_segment_var]] == seg]]
    metadt_seg <- metadt_ct[metadt_ct$dcc_filename %in% colnames(deconv_seg),]
  } else{
    deconv_seg <- deconv_ct
    metadt_seg <- metadt_ct
  }

  print(dim(deconv_seg))

  # run UMAP and tSNE
  ###########################
  # get top N variable genes
  if(!is.null(top_var)){
    per_gene_variance <- apply(deconv_seg, 1, stats::var)
    top_var_genes <- names(sort(per_gene_variance, decreasing = T)[1:top_var])

    deconv_seg <- deconv_seg[rownames(deconv_seg) %in% top_var_genes, ]
  }

  # do PCA
  pca_obj <- pca(deconv_seg, scale = T)
  pca_res <- t(-1*pca_obj$rotated) # reverse the signs of eigen vectors
  pca_loads <- -1*pca_obj$loadings
  #pca_vars <- pca_obj$variance
  metadt_seg[, c("PCA1","PCA2")] <- t(pca_res)[, c(1,2)]

  if(!is.null(top_pca)){
    deconv_seg <- pca_res[1:top_pca, ]
  }

  # make umap
  custom_umap <- umap::umap.defaults
  custom_umap$random_state <- 42
  umap_out <- umap(t(deconv_seg), config = custom_umap)
  metadt_seg[, c("UMAP1","UMAP2")] <- umap_out$layout[, c(1,2)]

  # make tsne
  set.seed(42)
  tsne_out <- Rtsne(t(deconv_seg), perplexity = ncol(deconv_seg)*.15)
  metadt_seg[, c("tSNE1","tSNE2")] <- tsne_out$Y[, c(1,2)]

  for(method in c('UMAP', 'tSNE', 'PCA')){

    # for discrete labels
    for(color_var in c(vars_all_disc 'roi_cluster_label')){
      print(color_var)

      sub <- ifelse(method == 'PCA', paste0('% of variance explained: PC1= ', as.character(round(pca_obj$variance[1], 2)),
                                            ' PC2= ', as.character(round(pca_obj$variance[2], 2))), '')

      ggplot(metadt_seg,
             aes(x = get(paste0(method, '1')),
                 y = get(paste0(method, '2')),
                 color = get(color_var), shape = get(aoi_segment_var))) +
        geom_point(size = 3) +
        xlab(paste0(method, '1')) +
        ylab(paste0(method, '2')) +
        theme(plot.subtitle = sub) +
        scale_color_discrete(name = color_var) +
        scale_shape_discrete(name = aoi_segment_var) +
        theme_bw()

      ggsave(file.path(out_dir, paste0('deconv_umap_tsne_', seg),
                       paste0(ct_name, '_', method, '_topvargenes_', ifelse(is.null(top_var), 'NULL', as.character(top_var)),
                              '_toppca_', ifelse(is.null(top_pca), 'NULL', as.character(top_pca)),
                              '_', color_var, '.png')),
             width = 2000, height = 1500, unit='px', device='png')
    }

    # for continuous variables
    for(color_var in vars_all_cont){
      print(color_var)

      sub <- ifelse(method == 'PCA', paste0('% of variance explained: PC1= ', as.character(round(pca_obj$variance[1], 2)),
                                            ' PC2= ', as.character(round(pca_obj$variance[2], 2))), '')

      ggplot(metadt_seg,
             aes(x = get(paste0(method, '1')),
                 y = get(paste0(method, '2')),
                 color = get(color_var), shape = get(aoi_segment_var))) +
        geom_point(size = 3) +
        xlab(paste0(method, '1')) +
        ylab(paste0(method, '2')) +
        theme(plot.subtitle = sub) +
        scale_color_continuous(name = color_var) +
        scale_shape_discrete(name = aoi_segment_var) +
        theme_bw()

      ggsave(file.path(out_dir, paste0('deconv_umap_tsne_', seg),
                       paste0(ct_name, '_', method, '_topvargenes_', ifelse(is.null(top_var), 'NULL', as.character(top_var)),
                              '_toppca_', ifelse(is.null(top_pca), 'NULL', as.character(top_pca)),
                              '_', color_var, '.png')),
             width = 2000, height = 1500, unit='px', device='png')
    }
  }
  # save PCA res df
  fwrite(pca_loads[, 1:10], file.path(out_dir, paste0('deconv_umap_tsne_', seg),
                                      paste0('pca_loads_', ct_name, '.csv')), row.names = F)
})

