# README: script for performing deconvolution using SpatialDecon (output: cell_types fractions)
# and BayesPrism (output: cell types fractions and cell-type-specific transcriptional profiles
# which undergo vst normalisation and batch effect correction)

# recommended usage for scRNAseq reference dataset is raw counts
# although not log transofrmation of both sc and bulk is also ok

# TODO don't add NA samples for missing ones in deconv
# TODO improve adjust synonym genes function..
# TODO calculate q3 normalised dataframes and choose normalisation method

# define variables --------------------------------------------------------

ct_nr_thr <- 50 # min recommended is 20 - rmv cell states lower than thr in scrnaseq
tumor_ct_name <- 'Epithelial cells' # tumor ct label in scrna_anno
adjust_synonym_gene_names <- T # whether or not to adjust synonymical gene names between scRNAsea and GeoMX
# that help rescue typically around 300 genes with synonym names, but sometimes Ensembl not work
raw_counts_layer <- "counts" # raw counts slot (layer) name in reference scRNAseq dataset

deconv_norm_type <- 'q3_norm' # c('q3_norm', 'deseq2_vst') which norm should be used for bayesprism results
sd_norm_type <- 'q3_norm' # suggested normalisation type for SpatialDecon (CANNOT BE IN LOG FORM), BayesPrism uses raw counts
batch_rm_type <- 'harmony' # c('harmony', 'limma')

# variables to merge the final csv with
meta_names <- c(aoi_id, roi_id, aoi_segment_var, sample_name, main_experimental_condition, 
                main_roi_label, other_vars_bio)

# main experimental conditions for limma batch eff rmv
exp_design <- as.formula(paste('~', aoi_segment_var, '+', main_experimental_condition))

# make dirs and set additional vars ---------------------------------------

dir.create(file.path(output_dir, 'deconvolution'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'deconvolution', 'spatial_decon', scrna_anno), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'deconvolution', 'bayes_prism', scrna_anno), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'deconvolution', 'bayes_prism', scrna_anno, 'hist'), showWarnings = T, recursive = T)

scrna_ref_cleaned_path <- file.path(output_dir, 'deconvolution', gsub('.RDS', '_cleaned_for_deconv.RDS', basename(scrna_ref_path)))

# prepare scrnaseq reference dataset --------------------------------------

# TODO move to functions script
adjust_scrna_ref <- function(scrna_ref_path, scrna_ref_cleaned_path, raw_counts_layer = "counts", 
                             tumor_ct_name = 'Epithelial cells', ct_nr_thr = 50,
                             pt_colname = 'publication_patient_code_final', 
                             adjust_synonym_gene_names = F, geomx_norm_batch_eff_rm_path = NULL){
  
  scrna_ref_obj <- readRDS(scrna_ref_path)
  raw_counts_mtx <- GetAssayData(object = scrna_ref_obj[["RNA"]], layer = raw_counts_layer)
  
  # check genes outliers
  scrna_stat <- plot.scRNA.outlier(
    input=t(raw_counts_mtx), #make sure the colnames are gene symbol or ENSMEBL ID
    cell.type.labels=scrna_ref_obj@meta.data$cell_type,
    species="hs", 
    return.raw=TRUE, #return the data used for plotting.
    pdf.prefix= gsub('.RDS', '', scrna_ref_cleaned_path) # specify pdf.prefix if need to output to pdf
  )
  
  # filter out outlier genes
  # this mtx is t()
  scrna_filt <- cleanup.genes (input=t(raw_counts_mtx),
                               input.type="count.matrix",
                               species="hs", 
                               gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                               exp.cells=5)
  
  # geomx doesn't have to be filtered since later on they took only intersection of genes
  
  # subset to protein coding genes and t() back
  scrna_filt_pc <-  t(select.gene.type(scrna_filt, gene.type = "protein_coding"))
  
  dim(raw_counts_mtx)
  dim(t(scrna_filt))
  dim(scrna_filt_pc)
  
  # subset initial object to filtered genes (before names adjustment)
  scrna_ref_obj_filt <- subset(scrna_ref_obj, features = rownames(scrna_filt_pc))
  
  if(adjust_synonym_gene_names){
    
    geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)
    # repair synonymuous gene names
    length(rownames(geomx_obj@assayData$exprs))
    length(rownames(scrna_filt_pc))
    
    adjusted_genes_rna <- adjust_synonym_genes(rownames(geomx_obj@assayData$exprs), rownames(scrna_filt_pc))

    print(paste0('initial number of intersected genes between scRNAseq and geomx_object: ', 
                 length(intersect(rownames(geomx_obj@assayData$exprs), rownames(scrna_filt_pc)))))
    print(paste0('number of intersected genes after adjustment: ', 
                 length(intersect(rownames(geomx_obj@assayData$exprs), adjusted_genes_rna))))
    print(which(!(rownames(scrna_filt_pc) == adjusted_genes_rna)))
    
    # change rownames to adjusted gene names
    rownames(scrna_filt_pc) <- adjusted_genes_rna
  }
  
  which(!(rownames(scrna_filt_pc) == rownames(scrna_ref_obj_filt))) # double check
  
  # make new counts slot with filtered and adjusted mtx (Warning is expected since gene names are adjusted)
  scrna_ref_obj_filt <- SetAssayData(
    object = scrna_ref_obj_filt,
    layer = raw_counts_layer,
    new.data = scrna_filt_pc,
    assay = "RNA"
  )
  
  # adjust gene names in the main object
  rownames(scrna_ref_obj_filt) <- rownames(scrna_filt_pc)
  
  # change cell_type name to tumor
  scrna_ref_obj_filt@meta.data$cell_type <- ifelse(scrna_ref_obj_filt@meta.data$cell_type == tumor_ct_name, 
                                              'tumor', scrna_ref_obj_filt@meta.data$cell_type)
  
  # cell states - clustering tumor cells by patient
  scrna_ref_obj_filt@meta.data$cell_state <- ifelse(scrna_ref_obj_filt@meta.data$cell_type == 'tumor', 
                                               paste0('tumor_', scrna_ref_obj_filt@meta.data[[pt_colname]]), 
                                               scrna_ref_obj_filt@meta.data$cell_type)
  
  # remove cells from cell states with nr < thr 
  ct_freq <- as.data.frame(table(scrna_ref_obj_filt@meta.data$cell_state))
  cells_above_ct_thr <- scrna_ref_obj_filt@meta.data$cell_name[scrna_ref_obj_filt@meta.data$cell_state %in% 
                                                            as.character(ct_freq$Var1[ct_freq$Freq > ct_nr_thr])]
  
  scrna_ref_obj_filt <- subset(scrna_ref_obj_filt, cells = cells_above_ct_thr)
  
  # TODO find solution not connected with seurat version
  # set up correct dimnames in slot (bug in seurat v5 and its lost)
  scrna_ref_obj_filt@assays$RNA@layers[[raw_counts_layer]]@Dimnames <- dimnames(scrna_ref_obj_filt)
  
  print(paste0('dim of the original scRNAseq reference dataset: ', dim(scrna_ref_obj)))
  print(paste0('dim of the final cleaned scRNAseq reference dataset: ', dim(scrna_ref_obj_filt)))
  
  # save adjusted scRNAseq file
  saveRDS(scrna_ref_obj_filt, file = scrna_ref_cleaned_path)
  
  print(paste0("cleaned reference scRNAseq dataset made from ", scrna_ref_path,
               " have been saved under the path: ", scrna_ref_cleaned_path))
  
  ########################################
  # QC of cell states
  
  # plot.cor.phi (input=t(scrna_ref_obj@assays$RNA@data),
  #               input.labels=scrna_ref_obj@meta.data$cell_state,
  #               title="cell state correlation",
  #               #specify pdf.prefix if need to output to pdf
  #               #pdf.prefix="gbm.cor.cs",
  #               cexRow=0.6, cexCol=0.6,
  #               margins=c(6,6))
  # 
  # dev.off()
  # 
  # plot.cor.phi (input=t(scrna_ref_obj@assays$RNA@data),
  #               input.labels=scrna_ref_obj@meta.data$mid_lvl_ct,
  #               title="cell type correlation",
  #               #specify pdf.prefix if need to output to pdf
  #               #pdf.prefix="gbm.cor.ct",
  #               cexRow=0.5, cexCol=0.5,
  # )
  # 
  # dev.off()
}



if(!file.exists(scrna_ref_cleaned_path)){
  
  adjust_scrna_ref(scrna_ref_path, scrna_ref_cleaned_path, raw_counts_layer, 
                   tumor_ct_name, ct_nr_thr, pt_colname = 'publication_patient_code_final',
                   adjust_synonym_gene_names, geomx_norm_batch_eff_rm_path)
} else{
  print(paste0("cleaned reference scRNAseq dataset made from ", scrna_ref_path,
               " already exists under the path: ", scrna_ref_cleaned_path))
}

# deconvolution by bayesprism ---------------------------------------------
# https://github.com/Danko-Lab/BayesPrism/blob/main/tutorial_deconvolution.html

# load cleaned scrna and geomx
geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)
scrna_ref_obj <- readRDS(scrna_ref_cleaned_path)
scrna_ref_raw_counts_mtx <- GetAssayData(object = scrna_ref_obj[["RNA"]], layer = raw_counts_layer)

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
saveRDS(bprism_res, file = file.path(output_dir,'deconvolution', 'bayes_prism', 
                                     paste0('bp_res_', scrna_anno, '.RDS')))


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

ct_frac <- rownames_to_column(as.data.frame(ct_frac), 'dcc_filename')
ct_frac <- left_join(ct_frac, sData(geomx_obj)[, meta_names],
                     by = 'dcc_filename')

fwrite(ct_frac, file.path(output_dir,'deconvolution', 'bayes_prism', 
                          paste0('bp_res_', scrna_anno, '_ct_fraction.csv')))


# normalise deconvolution expr mtx  ---------------------------------------

deconv_ct_norm_list <- lapply(ct_names, function(ct_name){
  print(ct_name)
  cell_frac_cv <- as.data.frame(bprism_res@posterior.theta_f@theta.cv)
  # mask ct_frac results if cv > 0.2-0.5 (0.1 thr for bulk, 0.5 for Visium, GeoMx should be in the middle)
  # histogram from batch 1 suggests 0.2 as thr
  cell_to_rm <- rownames(cell_frac_cv)[cell_frac_cv[[ct_name]] > 0.2]
  
  deconv_ct <- BayesPrism::get.exp(bp=bprism_res,
                                   state.or.type="type",
                                   cell.name=ct_name)
  
  deconv_ct_cleaned <- deconv_ct[!(rownames(deconv_ct) %in% cell_to_rm), ] # here genes in cols, dcc in rows
  
  plot_expr_distribution(deconv_ct, paste0(ct_name, '_raw'), 
                         file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                   scrna_anno, 'hist', paste0('expr_hist_', ct_name, '_raw.png')), is_log = F)
  
  if(deconv_norm_type == 'deseq2_vst'){
    # do vst normalisation
    # before it was wrapped in trycatch, may be needed to comeback
    deconv_ct_cleaned_norm <- varianceStabilizingTransformation(round(t(deconv_ct_cleaned)))
    islog <- T
    
    # TODO may be needed on the downstream analysis steps
    # remove genes with 0 variance across whole dataset (artifact from deconv + vst)
    # per_gene_variance <- apply(deconv_ct_cleaned_vst, 1, var)
    # genes_var0 <- names(per_gene_variance)[which(per_gene_variance == 0)]
    # deconv_ct_cleaned_vst <- deconv_ct_cleaned_vst[!(rownames(deconv_ct_cleaned_vst) %in% genes_var0),]
  } else if(deconv_norm_type == 'q3_norm'){
    # upper quartile normalisation from sourcecode of GeoMxTools normalize() function
    qs <- apply(t(deconv_ct_cleaned), 2, function(x) stats::quantile(x, 0.75))
    deconv_ct_cleaned_norm <- sweep(t(deconv_ct_cleaned), 2L, qs / ngeoMean(qs), FUN = "/")
    islog <- F
    
    #qs_norm_factors <- qs / ngeoMean(qs)     # may be needed for additional validation
  } else{ stop(print('choose either deseq2_vst or q3_norm as deconv_norm_type'))}
  
  plot_expr_distribution(deconv_ct_cleaned_norm, paste(ct_name, deconv_norm_type), 
                         file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                   scrna_anno, 'hist', paste0('expr_hist_', ct_name, '_', deconv_norm_type, '.png')), is_log = islog)

  return(deconv_ct_cleaned_norm)
})

names(deconv_ct_norm_list) <- ct_names

# clean list from ct for which vst was not computed 
deconv_ct_norm_list[sapply(deconv_ct_norm_list, is.null)] <- NULL

saveRDS(deconv_ct_norm_list, file = file.path(output_dir,'deconvolution', 'bayes_prism', 
                                     paste0('bp_res_', scrna_anno, '_expr_mtx_cleaned_', deconv_norm_type, '.RDS')))



# do batch effect correction ----------------------------------------------

# make metadata for batch effect correction
meta_dt <- sData(geomx_obj)[, c(meta_names, primary_batch_var, secondary_batch_var)]

# do batch effect removal with harmony
deconv_batch_rm_list <- lapply(names(deconv_ct_norm_list), function(ct_name){
  print(ct_name)
  deconv_norm <- deconv_ct_norm_list[[ct_name]]
  
  # filter meta if some ROI does not contain given ct
  mata_dt_ct <- meta_dt[meta_dt[[aoi_id]] %in% colnames(deconv_norm), ]
  
  deconv_harmony_res <- t(HarmonyMatrix(deconv_norm, 
                                        meta_data = mata_dt_ct,
                                        vars_use = c(primary_batch_var, secondary_batch_var)))
  
  plot_expr_distribution(deconv_harmony_res, paste0(ct_name, '_harmony_corr'), 
                         file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                   scrna_anno, 'hist', paste0('expr_hist_', ct_name, '_', deconv_norm_type, '_harmony_corr.png')), is_log = islog)
  
  return(deconv_harmony_res)
})

names(deconv_batch_rm_harm_list) <- names(deconv_ct_list)


saveRDS(deconv_batch_rm_harm_list, file = file.path(output_dir,'deconvolution', 'bayes_prism', 
                                                    paste0('bp_res_', scrna_anno,'_expr_mtx_cleaned_', deconv_norm_type, '_harmony_batch_corr.RDS')))

# remove batch effect with limma

deconv_batch_rm_limma_list <- lapply(names(deconv_ct_list), function(ct_name){
  print(ct_name)
  deconv_vst <- deconv_ct_list[[ct_name]]
  
  # filter meta if some ROI does not contain given ct
  mata_dt_ct <- meta_dt[meta_dt[[aoi_id]] %in% colnames(deconv_vst), ]
  
  design <- model.matrix(exp_design, data = mata_dt_ct)
  
  prim_batch <-  mata_dt_ct[[primary_batch_var]]
  
  if(!is.null(secondary_batch_var)){
    second_batch <- mata_dt_ct[[secondary_batch_var]]
  }else{second_batch <- NULL} 
  
  if(!is.null(cov_design)){
    cov <- model.matrix(cov_design, data = mata_dt_ct)
  }else{cov <- NULL} 

  
  deconv_limma_res <- limma::removeBatchEffect(deconv_vst, batch = prim_batch, batch2 = second_batch,
                                        covariates = cov, design = design)
  

  
  plot_expr_distribution(deconv_limma_res, paste0(ct_name, '_limma_corr'), 
                         file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                   scrna_anno, 'hist',
                                   paste0('expr_hist_', ct_name, '_limma_batch_corr_', 
                                          primary_batch_var, secondary_batch_var,'.png')), is_log = T)


  return(deconv_limma_res)
})

names(deconv_batch_rm_limma_list) <- names(deconv_ct_list)

saveRDS(deconv_batch_rm_limma_list, file = file.path(output_dir,'deconvolution', 'bayes_prism', 
                                                     paste0('bp_res_', scrna_anno, '_expr_mtx_cleaned_vst_limma_batch_corr_', 
                                                            primary_batch_var, secondary_batch_var, '.RDS')))


# prepare data for SpatialDecon -------------------------------------------
# from
# https://bioconductor.org/packages/release/bioc/vignettes/SpatialDecon/inst/doc/SpatialDecon_vignette_NSCLC.html

# load cleaned scrna and geomx
geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)
scrna_ref_obj <- readRDS(scrna_ref_cleaned_path)
scrna_ref_raw_counts_mtx <- GetAssayData(object = scrna_ref_obj[["RNA"]], layer = raw_counts_layer)

scrna_mtx_name <- ifelse('RNA_common_genes' %in% colnames(scrna_ref_obj@meta.data), 'RNA_common_genes', 'RNA')

# filter geomx object from low complexity genes
geomx_stat <- plot.bulk.outlier(
  bulk.input=t(geomx_obj@assayData$exprs),#make sure the colnames are gene symbol or ENSMEBL ID
  sc.input=t(scrna_ref_raw_counts_mtx), #make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels=scrna_ref_obj@meta.data$cell_type,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE,
  pdf.prefix= file.path(output_dir, 'deconvolution', 'spatial_decon', 
                        paste0('sd_res_', scrna_anno)) #specify pdf.prefix if need to output to pdf
)


geomx_stat_to_rm <- geomx_stat[ rowSums(geomx_stat[, -c(1,2)]) >= 1, ]
geomx_filtered <- geomx_obj[!(rownames(geomx_obj) %in% rownames(geomx_stat_to_rm)),  ]

# subset to protein coding genes (neg probe have to be added for bg modelling)
geomx_pc <-  colnames(select.gene.type(t(geomx_filtered@assayData$exprs), gene.type = "protein_coding"))
geomx_filtered_pc <- geomx_filtered[rownames(geomx_filtered) %in% c(geomx_pc, "NegProbe-WTX"),  ]

featureType(geomx_filtered_pc) <- "Target"
sampleNames(geomx_filtered_pc) <- sData(geomx_filtered_pc)[['dcc_filename']]

dim(geomx_obj)
dim(geomx_filtered)
dim(geomx_filtered_pc)

# prepare cell profile matrix from reference scRNAseq

# format annotations
scrna_anno_dt <- scrna_ref_obj@meta.data[, c('cell_name', scrna_anno)]
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

# run spatial decon with bg estimated genes
# estimate bcg for every segment based on neg probes (under the hood)
patient_nr <- length(unique(sData(geomx_filtered_pc)$Patient))

sd_res_custom <- runspatialdecon(object = geomx_filtered_pc,
                                    norm_elt = sd_norm_type,                # normalized data
                                    raw_elt = "exprs",                    
                                    X = custom_oc_mtx,                            
                                    #cell_counts = geomx_obj$Nuclei,      # nuclei counts, used to estimate total cells
                                    #is_pure_tumor = geomx_obj$istumor,   # identities of the Tumor segments/observations
                                    n_tumor_clusters = patient_nr)               # how many distinct tumor profiles to append to safeTME


saveRDS(sd_res_custom, file = file.path(output_dir, 'deconvolution', 'spatial_decon', 
                                        paste0('sd_res_', scrna_anno, '_geomxfiltpc.RDS')))


# extract and save ct fractions 
ct_frac_st <- rownames_to_column(data.frame(pData(sd_res_custom)[, 'prop_of_all']), 'dcc_filename')
ct_frac_st <- left_join(ct_frac_st, sData(geomx_obj)[, meta_names],
                     by = 'dcc_filename')

fwrite(ct_frac_st, file.path(output_dir,'deconvolution', 'spatial_decon', 
                             paste0('sd_res_', scrna_anno, 
                                    '_geomxfiltpc_ct_fraction.csv')))

# write logs --------------------------------------------------------------

# save logs
writeLines(c('deconvolution logs:',
             '; spatial decon normalisation type : ', norm_type,
             '; minimum cell type number : ', ct_nr_thr,
             '; cell type annotation  : ', scrna_anno,
             '; limma primary batch effect variable : ', primary_batch_var,
             '; limma secondary batch effect variable : ', secondary_batch_var,
             '; limma experimental design : ', as.character(exp_design)[2],
             '; limma covariate : ', as.character(cov_design)[2]), deconv_logs_path)
