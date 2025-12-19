# README: script for performing deconvolution using SpatialDecon (output: cell_types fractions)
# and BayesPrism (output: cell types fractions and cell-type-specific transcriptional profiles
# which undergo vst normalisation and batch effect correction)

# recommended usage for scRNAseq reference dataset is raw counts
# although not log transofrmation of both sc and bulk is also ok

# TODO improve adjust synonym genes function..

# define variables --------------------------------------------------------

# variables for reference scRNAseq dataset
ct_nr_thr <- 50 # min recommended is 20 - rmv cell states lower than thr in scrnaseq
# tumor ct label in scrna_anno
tumor_ct_name <- 'Epithelial cells' #vaharautio
#tumor_ct_name <- 'tumor' #hautaniemi
# patient col name in reference seq
pt_colname <- 'publication_patient_code_final' #vaharautio
#pt_colname <- 'patient_id' #hautaniemi
# single cell id colname in scrnaseq ref
cell_id <- 'cell_name' #vaharautio
#cell_id <- 'cell' #hautaiemi

adjust_synonym_gene_names <- T # whether or not to adjust synonymical gene names between scRNAsea and GeoMX
# that help rescue typically around 300 genes with synonym names, but sometimes Ensembl not work
raw_counts_layer <- "counts" # raw counts slot (layer) name in reference scRNAseq dataset

# which norm should be used for bayesprism results
deconv_norm_type <- 'q3_norm' # c('q3_norm', 'log_norm', 'deseq2', 'deseq2_vst', 'libsize_log')
batch_rm_type <- 'harmony' # c('harmony', 'limma')
#cell_frac_cutoff = 0.005 # aois with lower ct frec will be removed for given ct for bp results processing

sd_norm_type <- 'q3_norm' # suggested normalisation type for SpatialDecon (CANNOT BE IN LOG FORM)

# variables to merge the final csv with
meta_names <- unique(c(aoi_id, roi_id, aoi_segment_var, sample_name, main_experimental_condition, 
                main_roi_label, other_vars_bio))

# main experimental conditions for limma batch eff rmv
exp_design <- as.formula(paste('~', aoi_segment_var, '+', main_experimental_condition))

# make dirs and set additional vars ---------------------------------------

dir.create(file.path(output_dir, 'deconvolution'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'deconvolution', 'spatial_decon', scrna_anno), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'deconvolution', 'bayes_prism', scrna_anno), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'deconvolution', 'bayes_prism', scrna_anno, 'hist'), showWarnings = T, recursive = T)

scrna_ref_cleaned_path <- file.path(output_dir, 'deconvolution', gsub('.RDS', '_cleaned_for_deconv.RDS', basename(scrna_ref_path)))
bp_res_path <- file.path(output_dir,'deconvolution', 'bayes_prism', paste0('bp_res_', scrna_anno, '.RDS'))
bp_ct_frac_path <- file.path(output_dir,'deconvolution', 'bayes_prism', paste0('bp_res_', scrna_anno, '_ct_fraction.csv')) 

norm_is_log <- ifelse(deconv_norm_type %in%  c('q3_norm', 'deseq2_norm'), FALSE, TRUE)

# prepare scrnaseq reference dataset --------------------------------------

if(!file.exists(scrna_ref_cleaned_path)){
  
  adjust_scrna_ref(scrna_ref_path, scrna_ref_cleaned_path, raw_counts_layer, 
                   tumor_ct_name, ct_nr_thr, pt_colname,
                   adjust_synonym_gene_names, geomx_norm_batch_eff_rm_path)
  gc()
} else{
  print(paste0("cleaned reference scRNAseq dataset made from ", scrna_ref_path,
               " already exists under the path: ", scrna_ref_cleaned_path))
}

# deconvolution by bayesprism ---------------------------------------------
# https://github.com/Danko-Lab/BayesPrism/blob/main/tutorial_deconvolution.html

# load cleaned scrna and geomx
geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)
# make metadata for batch effect correction
meta_dt <- sData(geomx_obj)[, c(meta_names, primary_batch_var, secondary_batch_var)]
scrna_ref_obj <- readRDS(scrna_ref_cleaned_path)
scrna_ref_raw_counts_mtx <- GetAssayData(object = scrna_ref_obj[["RNA"]], layer = raw_counts_layer)

print('cell types and cell states to use:')
print(table(scrna_ref_obj@meta.data[[scrna_anno]], scrna_ref_obj@meta.data$cell_state))

# make a prism object
prism_obj <- new.prism(
  reference=t(scrna_ref_raw_counts_mtx), 
  mixture=t(geomx_obj@assayData$exprs),
  input.type="count.matrix", 
  cell.type.labels = scrna_ref_obj@meta.data[[scrna_anno]], 
  cell.state.labels = scrna_ref_obj@meta.data$cell_state,
  key="tumor",
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

# save some memory
rm(scrna_ref_obj)
rm(scrna_ref_raw_counts_mtx)
gc()

# run bayesprism
bprism_res <- run.prism(prism = prism_obj, n.cores = detectCores()-2)

# save res
saveRDS(bprism_res, file = bp_res_path)


# extract and save ct fractions  ------------------------------------------

cell_frac_cv <- as.data.frame(bprism_res@posterior.theta_f@theta.cv)
# mask ct_frac results if cv > 0.2-0.5 (0.1 thr for bulk, 0.5 for Visium, GeoMx should be in the middle)
# histogram from batch 1 suggests 0.2 as thr

ct_frac <- get.fraction (bp=bprism_res,
                         which.theta="final",
                         state.or.type="type")

ct_names <- colnames(ct_frac)

# mask  unreliable results
ct_frac[cell_frac_cv > 0.2] <- NA

ct_frac <- rownames_to_column(as.data.frame(ct_frac), aoi_id)
ct_frac <- left_join(ct_frac, sData(geomx_obj)[, meta_names],
                     by = aoi_id)

fwrite(ct_frac, bp_ct_frac_path)


# normalise deconvolution expr mtx  ---------------------------------------

deconv_ct_norm_list <- lapply(ct_names, function(ct_name){
  print(ct_name)
  cell_frac_cv <- as.data.frame(bprism_res@posterior.theta_f@theta.cv)
  # mask ct_frac results if cv > 0.2-0.5 (0.1 thr for bulk, 0.5 for Visium, GeoMx should be in the middle)
  # histogram from batch 1 suggests 0.2 as thr
  cell_to_rm <- rownames(cell_frac_cv)[cell_frac_cv[[ct_name]] > 0.2]
  
  deconv_ct <- t(BayesPrism::get.exp(bp=bprism_res,
                                     state.or.type="type",
                                     cell.name=ct_name))
  
  # remove unreliable cells predictions
  deconv_ct_cleaned <- deconv_ct[,!(colnames(deconv_ct) %in% cell_to_rm)]
  # remove genes with only 0 counts
  deconv_ct_cleaned <- deconv_ct_cleaned[rowSums(deconv_ct_cleaned) != 0,]
  # remove aois with only 0 counts
  deconv_ct_cleaned <- deconv_ct_cleaned[, colSums(deconv_ct_cleaned) != 0]
  
  plot_expr_distribution(deconv_ct, paste0(ct_name, '_raw'),
                         file.path(output_dir, 'deconvolution', 'bayes_prism',
                                   scrna_anno, 'hist', paste0('expr_hist_', ct_name, '_raw.png')), is_log = F)

  # do normalisation
  if(deconv_norm_type == 'deseq2_vst'){
    
    deconv_ct_cleaned_norm <- do_normalisation(deconv_ct_cleaned, norm_type = 'deseq2_vst_onestep')
  } else if(deconv_norm_type %in%  c('q3_norm', 'log_norm', 'deseq2', 'libsize_log')){ 
    
    deconv_ct_cleaned_norm <- do_normalisation(deconv_ct_cleaned, norm_type = deconv_norm_type,
                                               meta_data = meta_dt, 
                                               aoi_segment_var, main_experimental_condition)
    
  } else{ 
    stop(print('choose either deseq2, deseq2_vst, q3_norm, log_norm or libsize_log as deconv_norm_type'))
  }
  
  # remove aois with only inf/nan normalised counts
  deconv_ct_cleaned_norm <- deconv_ct_cleaned_norm[, which(!is.na(colSums(deconv_ct_cleaned_norm)))]
  
  plot_expr_distribution(deconv_ct_cleaned_norm, paste(ct_name, deconv_norm_type),
                         file.path(output_dir, 'deconvolution', 'bayes_prism',
                                   scrna_anno, 'hist', paste0('expr_hist_', ct_name, '_', deconv_norm_type, '.png')), is_log = norm_is_log)
  
  return(deconv_ct_cleaned_norm)
})

names(deconv_ct_norm_list) <- ct_names

# clean list from ct for which norm was not computed
deconv_ct_norm_list[sapply(deconv_ct_norm_list, is.null)] <- NULL

saveRDS(deconv_ct_norm_list, file = file.path(output_dir,'deconvolution', 'bayes_prism',
                                     paste0('bp_res_', scrna_anno, '_expr_mtx_cleaned_', deconv_norm_type, '.RDS')))

# # do batch effect correction ----------------------------------------------

deconv_batch_rm_list <- lapply(names(deconv_ct_norm_list), function(ct_name){
  print(ct_name)

  if(!norm_is_log){
    # make log2 transformed normalised counts if norm_type not in log scale
    deconv_ct_norm <- log2(deconv_ct_norm_list[[ct_name]] + 1)
  } else{
    deconv_ct_norm <- deconv_ct_norm_list[[ct_name]]
  }

  # filter meta if some ROI does not contain given ct
  mata_dt_ct <- meta_dt[meta_dt[[aoi_id]] %in% colnames(deconv_ct_norm), ]

  if(batch_rm_type == 'harmony'){

    deconv_batch_rm <- t(HarmonyMatrix(deconv_ct_norm,
                                          meta_data = mata_dt_ct,
                                          vars_use = c(primary_batch_var, secondary_batch_var)))
  } else if(batch_rm_type == 'limma'){

    design <- model.matrix(exp_design, data = mata_dt_ct)
    prim_batch <-  mata_dt_ct[[primary_batch_var]]

    if(!is.null(secondary_batch_var)){
      second_batch <- mata_dt_ct[[secondary_batch_var]]
    }else{second_batch <- NULL}

    deconv_batch_rm <- limma::removeBatchEffect(deconv_ct_norm, batch = prim_batch,
                                                 batch2 = second_batch, design = design)

  } else{ stop(print('choose either harmony or limma as batch_rm_type'))}


  plot_expr_distribution(deconv_batch_rm, paste(ct_name, batch_rm_type, ' batch effect corr'),
                         file.path(output_dir, 'deconvolution', 'bayes_prism',
                                   scrna_anno, 'hist', paste0('expr_hist_', ct_name, '_',
                                                              deconv_norm_type, '_', batch_rm_type, '_corr.png')),
                         is_log = T)

  return(deconv_batch_rm)
})

names(deconv_batch_rm_list) <- names(deconv_ct_norm_list)

saveRDS(deconv_batch_rm_list, file = file.path(output_dir,'deconvolution', 'bayes_prism',
                                                    paste0('bp_res_', scrna_anno,'_expr_mtx_cleaned_', deconv_norm_type, '_', batch_rm_type, '_corr.RDS')))




# pull all deconv cells together and norm + batch corr --------------------


# TODO move to cellchat
# TODO adjust - now for deseq2 normal metadata
# deconv_norm_all <- create_norm_pseudosc_from_deconv(bp_res_path, bp_ct_frac_path, scrna_anno, cell_frac_cutoff, 
#                                                     bp_pseudosc_path, 
#                                                     norm_type = deconv_norm_type,
#                                                     meta_data = meta_dt, 
#                                                     aoi_segment_var, 
#                                                     main_experimental_condition)
# 
# # make metadata per dcc-ct to match deconv colnames
# dcc_ct <- data.frame('dcc_filename' = gsub('_.*', '', colnames(deconv_norm_all)), 
#                      'dcc_ct' = colnames(deconv_norm_all))
# meta_data_ct <- left_join(dcc_ct, meta_dt, by = 'dcc_filename')
# meta_data_ct$ct_label <- gsub('^[^_]*', '', meta_data_ct$dcc_ct)
# meta_data_ct$ct_label <- gsub('^_', '', meta_data_ct$ct_label)


# batch effect correction 

# if(!norm_is_log){
#   # make log2 transformed normalised counts if norm_type not in log scale
#   deconv_norm_all <- log2(deconv_norm_all + 1)
# }
# 
# # do batch effect correction
# if(batch_rm_type == 'harmony'){
#   
#   deconv_batch_rm_all <- t(HarmonyMatrix(deconv_norm_all,
#                                      meta_data = meta_data_ct,
#                                      vars_use = c(primary_batch_var, secondary_batch_var)))
# } else if(batch_rm_type == 'limma'){
#   
#   design <- model.matrix(exp_design, data = meta_data_ct)
#   prim_batch <-  meta_data_ct[[primary_batch_var]]
#   
#   if(!is.null(secondary_batch_var)){
#     second_batch <- meta_data_ct[[secondary_batch_var]]
#   }else{second_batch <- NULL}
#   
#   deconv_batch_rm_all <- limma::removeBatchEffect(deconv_norm_all, batch = prim_batch,
#                                               batch2 = second_batch, design = design)
#   
# } else{ stop(print('choose either harmony or limma as batch_rm_type'))}
# 
# plot_expr_distribution(deconv_batch_rm_all, paste(batch_rm_type, ' batch effect corr'),
#                        file.path(output_dir, 'deconvolution', 'bayes_prism',
#                                  paste0('expr_hist_bp_res_pseudosc_', 
#                                           scrna_anno, 'ct_frac_', cell_frac_cutoff, '_', 
#                                           deconv_norm_type, '_', batch_rm_type, '_corr.png')),
#                        is_log = T)
# 
# 
# fwrite(deconv_batch_rm_all, file = file.path(output_dir,'deconvolution', 'bayes_prism',
#                                                paste0('bp_res_pseudosc_', scrna_anno, 'ct_frac_',
#                                                       cell_frac_cutoff, '_', deconv_norm_type, '_', batch_rm_type, '_corr.csv')),
#        row.names = T)


# prepare data for SpatialDecon -------------------------------------------
# from
# https://bioconductor.org/packages/release/bioc/vignettes/SpatialDecon/inst/doc/SpatialDecon_vignette_NSCLC.html

# load cleaned scrna and geomx
geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)
scrna_ref_obj <- readRDS(scrna_ref_cleaned_path)
scrna_ref_raw_counts_mtx <- GetAssayData(object = scrna_ref_obj[["RNA"]], layer = raw_counts_layer)

geomx_filtered_pc <-   remove_low_complex_and_noncoding_genes(geomx_obj, scrna_ref_obj, raw_counts_layer = 'counts')

featureType(geomx_filtered_pc) <- "Target"
sampleNames(geomx_filtered_pc) <- sData(geomx_filtered_pc)[[aoi_id]]

dim(geomx_obj)
dim(geomx_filtered_pc)

# prepare cell profile matrix from reference scRNAseq

# format annotations
scrna_anno_dt <- scrna_ref_obj@meta.data[, c(cell_id, scrna_anno)]
rownames(scrna_anno_dt) <- NULL
colnames(scrna_anno_dt) <- c('cell_name', 'cell_type')

# TODO examine scalingFactor: 1 or 5 or what?
custom_oc_mtx <- create_profile_matrix(mtx = scrna_ref_raw_counts_mtx,            # cell x gene count matrix
                                       cellAnnots = scrna_anno_dt,  # cell annotations with cell type and cell name as columns
                                       cellTypeCol = "cell_type",  # column containing cell type
                                       cellNameCol = "cell_name",           # column containing cell ID/name
                                       matrixName = "oc_scrnaseq_ref_cell_type_filt_pc_norm", # name of final profile matrix
                                       outDir = output_dir,                    # path to desired output directory, set to NULL if matrix should not be written
                                       normalize = TRUE,                # Should data be normalized?
                                       minCellNum = ct_nr_thr,                   # minimum number of cells of one type needed to create profile, exclusive
                                       minGenes = 10,                    # minimum number of genes expressed in a cell, exclusive
                                       scalingFactor = 1,                # what should all values be multiplied by for final matrix
                                       discardCellTypes = TRUE)          # should cell types be filtered for types like mitotic, doublet, low quality, unknown, etc.

# run extended SpatialDecon with custom oc mtx ----------------------------
# TODO run with nuclei_counts when it will be counted reliably from cycif 

# estimate bcg for every segment based on neg probes (under the hood)
patient_nr <- length(unique(sData(geomx_filtered_pc)$Patient))

# sd_res_custom <- runspatialdecon(object = geomx_filtered_pc,
#                                     norm_elt = sd_norm_type,                # normalized data
#                                     raw_elt = "exprs",                    
#                                     X = custom_oc_mtx)
#                                     #cell_counts = geomx_obj$Nuclei,      # nuclei counts, used to estimate total cells
#                                     #is_pure_tumor = geomx_obj$istumor,   # identities of the Tumor segments/observations
#                                     #n_tumor_clusters = patient_nr)               # how many distinct tumor profiles to append to safeTME

# from inside of runspatialdecon()  to work also for ROI-merged
# some stats are not computed
# estimate background
bg <- derive_GeoMx_background(norm = as.matrix(Biobase::assayDataElement(geomx_filtered_pc , elt = sd_norm_type)),
                              # access the probe pool information from the feature metadata
                              probepool = Biobase::fData(geomx_filtered_pc)$Module,
                              # access the names of the negative control probes
                              negnames = Biobase::fData(geomx_filtered_pc)$TargetName[Biobase::fData(geomx_filtered_pc)$Negative])

sd_res_custom <- spatialdecon(norm = as.matrix(Biobase::assayDataElement(geomx_filtered_pc , elt = sd_norm_type)),
                              raw = as.matrix(Biobase::assayDataElement(geomx_filtered_pc , elt = "exprs")),
                              bg = bg,
                              X = custom_oc_mtx)

saveRDS(sd_res_custom, file = file.path(output_dir, 'deconvolution', 'spatial_decon', 
                                        paste0('sd_res_', scrna_anno, '_geomxfiltpc_rawresults.RDS')))


# extract and save ct fractions 
# ct_frac_st <- rownames_to_column(data.frame(pData(sd_res_custom)[, 'prop_of_all']), aoi_id)

ct_frac_st <- rownames_to_column(data.frame(t(sd_res_custom$prop_of_all)), aoi_id)
ct_frac_st <- left_join(ct_frac_st, sData(geomx_obj)[, meta_names],
                     by = aoi_id)


fwrite(ct_frac_st, file.path(output_dir,'deconvolution', 'spatial_decon', 
                             paste0('sd_res_', scrna_anno, 
                                    '_geomxfiltpc_ct_fraction.csv')))

# write logs --------------------------------------------------------------

# save logs
writeLines(c('deconvolution logs:',
             '; spatial decon input normalisation type : ', sd_norm_type,
             '; minimum cell type number : ', ct_nr_thr,
             '; cell type annotation  : ', scrna_anno,
             '; bp normalisation  : ', deconv_norm_type,
             '; limma primary batch effect variable : ', primary_batch_var,
             '; limma secondary batch effect variable : ', secondary_batch_var,
             '; limma experimental design : ', as.character(exp_design)[2]),
           deconv_logs_path)
