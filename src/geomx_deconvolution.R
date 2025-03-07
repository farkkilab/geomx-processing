
# recommended usage is raw counts, although not log transofrmation of both sc and bulk is also ok
#TODO normalise scRNAseq with deseq2norm and compare to raw


# define variables --------------------------------------------------------

norm_type <- 'q3_norm' # quantile is best for sd, bp works on raw counts, for sd norm cannot be in the log scale
ct_nr_thr <- 45 # best 45 for batch1 and 2 - to rmv cell states not abundant enough in scrnaseq

tumor_ct_name <- 'Epithelial cells' # tumor ct label in scrna_anno
adjust_synonym_gene_names <- F # whether or not to adjust synonymical gene names between scRNAsea and GeoMX
# that help rescue typically around 300 genes with synonym names, but sometimes Ensembl not work

meta_names <- c('dcc_filename', 'Patient', 'Segment', 'Sample', 'NACT_status', 'Annotation_cell', 'Site')

# main cause of the batch effect, from 1st PVCA plot
# should be the same as in batch effect rm script
main_batch_var <- 'batch_nr'
secondary_batch_var <- NULL

# main experimental conditions for limma batch eff rmv
exp_design <- formula(~ Segment + NACT_status)

# biological covariates which effect should be ignored by limma have to be in meta_names
# if NULL no cov are added to limma rmv batch eff
# TODO check if this is beneficial 
cov_design <- formula(~ Patient + Site) 
covname <- 'patient_site'

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'deconvolution'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'deconvolution', 'spatial_decon', scrna_anno), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'deconvolution', 'bayes_prism', scrna_anno), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'deconvolution', 'bayes_prism', scrna_anno, 'hist'), showWarnings = T, recursive = T)

scrna_ref_cleaned_path <- file.path(output_dir, 'deconvolution', gsub('.RDS', '_cleaned_for_deconv.RDS', basename(scrna_ref_path)))

# prepare scrnaseq reference dataset --------------------------------------

#TODO can be moved to some other script

if(!file.exists(scrna_ref_cleaned_path)){
  geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)
  scrna_ref_obj <- readRDS(scrna_ref_path)
  
  
  if(adjust_synonym_gene_names){
    # repair synonymuous gene names
    length(rownames(geomx_obj@assayData$exprs))
    length(rownames(scrna_ref_obj@assays$RNA@data))
    length(intersect(rownames(geomx_obj@assayData$exprs), rownames(scrna_ref_obj@assays$RNA@data)))
    
    adjusted_genes_rna <- adjust_synonym_genes(rownames(geomx_obj@assayData$exprs), rownames(scrna_ref_obj@assays$RNA@data))
    
    #  make a new assay with renamed genes
    RNA_common_genes <- scrna_ref_obj@assays$RNA
    RNA_common_genes@counts@Dimnames[[1]] <- adjusted_genes_rna
    RNA_common_genes@data@Dimnames[[1]] <- adjusted_genes_rna
    scrna_ref_obj@assays$RNA_common_genes <- RNA_common_genes
    
    length(intersect(rownames(geomx_obj@assayData$exprs), rownames(scrna_ref_obj@assays$RNA@data)))
    length(intersect(rownames(geomx_obj@assayData$exprs), rownames(scrna_ref_obj@assays$RNA_common_genes@data)))
    
    rna_mtx_touse <- 'RNA_common_genes'
  } else{
    rna_mtx_touse <- 'RNA'
  }
  
  # clean cell labels
  scrna_ref_obj@meta.data$cell_type <- ifelse(scrna_ref_obj@meta.data$cell_type == tumor_ct_name, 
                                              'tumor', scrna_ref_obj@meta.data$cell_type)
  
  # cell states - clustering tumor cells by patient
  scrna_ref_obj@meta.data$cell_state <- ifelse(scrna_ref_obj@meta.data$cell_type == 'tumor', 
                                               paste0('tumor_', scrna_ref_obj@meta.data$patient), 
                                               scrna_ref_obj@meta.data$cell_type)
  
  ########################################
  # QC of cell states
  # TODO think of changing labels for mast cells, Th17, tumor_H103
  
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
  #################################
  
  # check genes outliers
  scrna_stat <- plot.scRNA.outlier(
    input=t(scrna_ref_obj@assays[[rna_mtx_touse]]@data), #make sure the colnames are gene symbol or ENSMEBL ID
    cell.type.labels=scrna_ref_obj@meta.data$cell_type,
    species="hs", 
    return.raw=TRUE, #return the data used for plotting.
    pdf.prefix= gsub('.RDS', '', scrna_ref_cleaned_path) # specify pdf.prefix if need to output to pdf
  )
  
  
  # filter out outlier genes
  scrna_filt <- cleanup.genes (input=t(scrna_ref_obj@assays[[rna_mtx_touse]]@data),
                               input.type="count.matrix",
                               species="hs", 
                               gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                               exp.cells=5)
  
  dim(t(scrna_ref_obj@assays[[rna_mtx_touse]]@data))
  dim(scrna_filt)
  
  # geomx doesn't have to be filtered since later on they took only intersection of genes
  
  # subset to protein coding genes
  scrna_filt_pc <-  select.gene.type(scrna_filt, gene.type = "protein_coding")
  
  #  make a new assay with filtered genes
  RNA_filt_pc <- scrna_ref_obj@assays[[rna_mtx_touse]]
  RNA_filt_pc@counts <- RNA_filt_pc@counts[rownames(RNA_filt_pc@counts) %in% colnames(scrna_filt_pc),  ]
  RNA_filt_pc@data <- RNA_filt_pc@data[rownames(RNA_filt_pc@data) %in% colnames(scrna_filt_pc),  ]
  scrna_ref_obj@assays[[paste0(rna_mtx_touse, '_filt_pc')]] <- RNA_filt_pc
  
  # save adjusted scRNAseq file
  saveRDS(scrna_ref_obj, file = scrna_ref_cleaned_path)
  
  
  ###########################
  # TODO takes > 64G of memory, if needed have to be run on linux machine 
  # subset to signature genes (differentially expressed trough cell types)
  # diff_exp_stat <- get.exp.stat(sc.dat=scrna_raw[,colSums(scrna_raw>0)>3],# filter genes to reduce memory use
  #                               cell.type.labels=scrna_ref_obj@meta.data$cell_type,
  #                               cell.state.labels=scrna_ref_obj@meta.data$cell_state,
  #                               pseudo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
  #                               cell.count.cutoff=20, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
  #                               n.cores=8 #number of threads
  # )
  
  # scrna_filt_pc_sig <- select.marker (sc.dat=scrna_filt_pc,
  #                                          stat=diff_exp_stat,
  #                                          pval.max=0.01,
  #                                          lfc.min=0.1)
  
  # dim(scrna_filt_pc_sig)
  ##############################
  
} else{
  print(paste0("cleaned reference scRNAseq dataset made from ", scrna_ref_path,
  " already exists under the path: ", scrna_ref_cleaned_path))
}



# deconvolution by bayesprism ---------------------------------------------
# https://github.com/Danko-Lab/BayesPrism/blob/main/tutorial_deconvolution.html

# load cleaned scrna and geomx
geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)
scrna_ref_obj <- readRDS(scrna_ref_cleaned_path)

# remove cells from cell states with < thr cells in ref scrnaseq
ct_freq <- as.data.frame(table(scrna_ref_obj@meta.data$cell_state))
low_ct_cells <- scrna_ref_obj@meta.data$cell_name[scrna_ref_obj@meta.data$cell_state %in% 
                                                    as.character(ct_freq$Var1[ct_freq$Freq < ct_nr_thr])]

scrna_ref_obj <- scrna_ref_obj[, !colnames(scrna_ref_obj) %in% low_ct_cells]

# make a prism object
prism_obj <- new.prism(
  reference=t(scrna_ref_obj@assays$RNA_filt_pc@data), 
  mixture=t(geomx_obj@assayData$exprs),
  input.type="count.matrix", 
  cell.type.labels = scrna_ref_obj@meta.data[[scrna_anno]], 
  cell.state.labels = scrna_ref_obj@meta.data$cell_state,
  key="tumor",
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

# run bayesprism
bprism_res <- run.prism(prism = prism_obj, n.cores=18)

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

deconv_ct_list <- lapply(ct_names, function(ct_name){
  print(ct_name)
  cell_frac_cv <- as.data.frame(bprism_res@posterior.theta_f@theta.cv)
  # mask ct_frac results if cv > 0.2-0.5 (0.1 thr for bulk, 0.5 for Visium, GeoMx should be in the middle)
  # histogram from batch 1 suggests 0.2 as thr
  cell_to_rm <- rownames(cell_frac_cv)[cell_frac_cv[[ct_name]] > 0.2]
  
  deconv_ct <- BayesPrism::get.exp(bp=bprism_res,
                                   state.or.type="type",
                                   cell.name=ct_name)
  
  deconv_ct_cleaned <- deconv_ct[!(rownames(deconv_ct) %in% cell_to_rm), ]
  
  plot_expr_distribution(deconv_ct, paste0(ct_name, '_raw'), 
                         file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                   scrna_anno, 'hist', paste0('expr_hist_', ct_name, '_raw.png')), log = F)
  
  deconv_ct_cleaned_vst <- tryCatch({
    # do vst normalisation
    deconv_ct_cleaned_vst <- varianceStabilizingTransformation(round(t(deconv_ct_cleaned)))
    
    plot_expr_distribution(deconv_ct_cleaned_vst, paste0(ct_name, '_vst'), 
                           file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                     scrna_anno, 'hist', paste0('expr_hist_', ct_name, '_vst.png')), log = F)
    
    return(ct_name = deconv_ct_cleaned_vst)  # Return the result 
  }, error = function(e) {
    print('not enough AOIs with trustable predictions to perform vst. cell type is removed')
    return()
  })
})

names(deconv_ct_list) <- ct_names

# clean list from ct for which vst was not computed 
deconv_ct_list[sapply(deconv_ct_list, is.null)] <- NULL

saveRDS(deconv_ct_list, file = file.path(output_dir,'deconvolution', 'bayes_prism', 
                                     paste0('bp_res_', scrna_anno, '_expr_mtx_cleaned_vst.RDS')))


# do batch effect correction ----------------------------------------------

# make metadata for batch effect correction
meta_dt <- pData(geomx_obj)[, c(meta_names, main_batch_var, secondary_batch_var)]

# do batch effect removal with harmony
deconv_batch_rm_harm_list <- lapply(names(deconv_ct_list), function(ct_name){
  print(ct_name)
  deconv_vst <- deconv_ct_list[[ct_name]]
  
  # filter meta if some ROI does not contain given ct
  mata_dt_ct <- meta_dt[meta_dt$dcc_filename %in% colnames(deconv_vst), ]
  
  deconv_harmony_res <- t(HarmonyMatrix(deconv_vst, 
                                        meta_data = mata_dt_ct,
                                        vars_use = c(main_batch_var, secondary_batch_var)))
  
  plot_expr_distribution(deconv_harmony_res, paste0(ct_name, '_harmony_corr'), 
                         file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                   scrna_anno, 'hist', paste0('expr_hist_', ct_name, '_harmony_corr.png')), log = F)
  
  return(deconv_harmony_res)
})

names(deconv_batch_rm_harm_list) <- names(deconv_ct_list)

saveRDS(deconv_batch_rm_harm_list, file = deconv_bp_harm_path)

# remove batch effect with limma

deconv_batch_rm_limma_list <- lapply(names(deconv_ct_list), function(ct_name){
  print(ct_name)
  deconv_vst <- deconv_ct_list[[ct_name]]
  
  # filter meta if some ROI does not contain given ct
  mata_dt_ct <- meta_dt[meta_dt$dcc_filename %in% colnames(deconv_vst), ]
  
  design <- model.matrix(exp_design, data = mata_dt_ct)
  batch <-  mata_dt_ct[[main_batch_var]]
  
  batch2 <- switch((!is.null(secondary_batch_var)), 
                   mata_dt_ct[[secondary_batch_var]], NULL)
  
  cov <- switch((!is.null(cov_design)), 
                model.matrix(cov_design, data = mata_dt_ct), NULL)
  
  deconv_limma_res <- limma::removeBatchEffect(deconv_vst, batch = batch, batch2 = batch2,
                                        covariates = NULL, design = design)
  
  # TODO rmv after assessing what is better
  deconv_limma_res_cov <- limma::removeBatchEffect(deconv_vst, batch = batch, batch2 = batch2,
                                            covariates = cov, design = design)
  
  plot_expr_distribution(deconv_limma_res, paste0(ct_name, '_limma_corr'), 
                         file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                   scrna_anno, 'hist', paste0('expr_hist_', ct_name, '_limma_corr.png')), log = F)
  
  plot_expr_distribution(deconv_limma_res_cov, paste0(ct_name, '_limma_cov_corr'), 
                         file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                   scrna_anno, 'hist', paste0('expr_hist_', ct_name, '_limma_cov_corr.png')), log = F)

  return(list(limma = deconv_limma_res, limma_cov = deconv_limma_res_cov))
})

deconv_batch_rm_limma <- lapply(deconv_batch_rm_limma_list, `[[`, 1)
deconv_batch_rm_limma_cov <- lapply(deconv_batch_rm_limma_list, `[[`, 2)

names(deconv_batch_rm_limma) <- names(deconv_ct_list)
names(deconv_batch_rm_limma_cov) <- names(deconv_ct_list)

saveRDS(deconv_batch_rm_limma, file = gsub('harmony', 'limma', deconv_bp_harm_path))
saveRDS(deconv_batch_rm_limma_cov, file = gsub('harmony', 'limma_cov', deconv_bp_harm_path))

# prepare data for SpatialDecon -------------------------------------------
# from
# https://bioconductor.org/packages/release/bioc/vignettes/SpatialDecon/inst/doc/SpatialDecon_vignette_NSCLC.html

# load cleaned scrna and geomx
geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)
scrna_ref_obj <- readRDS(scrna_ref_cleaned_path)

scrna_mtx_name <- ifelse('RNA_common_genes' %in% colnames(scrna_ref_obj@meta.data), 'RNA_common_genes', 'RNA')

# filter geomx object from low complexity genes
geomx_stat <- plot.bulk.outlier(
  bulk.input=t(geomx_obj@assayData$exprs),#make sure the colnames are gene symbol or ENSMEBL ID
  sc.input=t(scrna_ref_obj@assays[[scrna_mtx_name]]@data), #make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels=scrna_ref_obj@meta.data$cell_type,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE,
  pdf.prefix= file.path(output_dir, 'deconvolution', 'spatial_decon', 
                        paste0('sd_res_', scrna_anno)) #specify pdf.prefix if need to output to pdf
)

geomx_stat_to_rm <- geomx_stat[ rowSums(geomx_stat[, -c(1,2)]) >= 1, ]
geomx_filtered <- geomx_obj[!(rownames(geomx_obj) %in% rownames(geomx_stat_to_rm)),  ]


featureType(geomx_obj) <- "Target"
sampleNames(geomx_obj) <- sData(geomx_obj)[['dcc_filename']]

featureType(geomx_filtered) <- "Target"
sampleNames(geomx_filtered) <- sData(geomx_filtered)[['dcc_filename']]

# prepare cell profile matrix from reference scRNAseq

# format annotations
scrna_anno_dt <- scrna_ref_obj@meta.data[, c('cell_name', scrna_anno)]
rownames(scrna_anno_dt) <- NULL
colnames(scrna_anno_dt) <- c('cell_name', 'cell_type')

# TODO examine scalingFactor: 1 or 5 or what?
custom_oc_mtx <- create_profile_matrix(mtx = scrna_ref_obj@assays$RNA_filt_pc@data,            # cell x gene count matrix
                                       cellAnnots = scrna_anno_dt,  # cell annotations with cell type and cell name as columns
                                       cellTypeCol = "cell_type",  # column containing cell type
                                       cellNameCol = "cell_name",           # column containing cell ID/name
                                       matrixName = "oc_scrnaseq_ref_cell_type_filt_pc", # name of final profile matrix
                                       outDir = output_dir,                    # path to desired output directory, set to NULL if matrix should not be written
                                       normalize = FALSE,                # Should data be normalized?
                                       minCellNum = ct_nr_thr,                   # minimum number of cells of one type needed to create profile, exclusive
                                       minGenes = 10,                    # minimum number of genes expressed in a cell, exclusive
                                       scalingFactor = 1,                # what should all values be multiplied by for final matrix
                                       discardCellTypes = TRUE)          # should cell types be filtered for types like mitotic, doublet, low quality, unknown, etc.

# run extended SpatialDecon with custom oc mtx ----------------------------

# TODO code repetition - rmv after checking if bg or no-bg is better 
# TODO run with nuclei_counts when it will be counted reliably from cycif 

sd_res_custom <- runspatialdecon(object = geomx_filtered,
                                    norm_elt = norm_type,                # normalized data
                                    raw_elt = "exprs",                    
                                    X = custom_oc_mtx,                            
                                    #cell_counts = geomx_obj$Nuclei,      # nuclei counts, used to estimate total cells
                                    #is_pure_tumor = geomx_obj$istumor,   # identities of the Tumor segments/observations
                                    n_tumor_clusters = 5)               # how many distinct tumor profiles to append to safeTME

# run spatial decon with bg estimated genes
# estimate bcg for every segment based on neg probes 
# TODO re-check if >1 module
negativeProbefData <- subset(fData(geomx_filtered), CodeClass == "Negative")
geomx_bg <- derive_GeoMx_background(norm = geomx_filtered@assayData[[norm_type]],
                                    probepool = fData(geomx_filtered)$Module,
                                    negnames = negativeProbefData$TargetName)

sd_res_custom_bg <- spatialdecon(norm = geomx_filtered@assayData[[norm_type]],                # normalized data
                                 bg = geomx_bg, # expected background counts for every data point in norm
                                 raw = geomx_filtered@assayData$exprs,                      
                                 X = custom_oc_mtx,                            
                                 #cell_counts = geomx_obj$Nuclei,      # nuclei counts, used to estimate total cells
                                 #is_pure_tumor = geomx_obj$istumor,   # identities of the Tumor segments/observations
                                 n_tumor_clusters = 5)               # how many distinct tumor profiles to append to safeTME



saveRDS(sd_res_custom, file = file.path(output_dir, 'deconvolution', 'spatial_decon', 
                                        paste0('sd_res_', scrna_anno, '_geomxfilt.RDS')))

saveRDS(sd_res_custom_bg, file = file.path(output_dir, 'deconvolution', 'spatial_decon', 
                                        paste0('sd_res_bg', scrna_anno, '_geomxfilt.RDS')))


# extract and save ct fractions 
ct_frac_st <- rownames_to_column(data.frame(pData(sd_res_custom)[, 'prop_of_all']), 'dcc_filename')
ct_frac_st <- left_join(ct_frac_st, sData(geomx_obj)[, meta_names],
                     by = 'dcc_filename')

fwrite(ct_frac_st, file.path(output_dir,'deconvolution', 'spatial_decon', 
                             paste0('sd_res_', scrna_anno, 
                                    '_geomxfilt_ct_fraction.RDS')))

ct_frac_st_bg <- rownames_to_column(data.frame(t(sd_res_custom_bg$prop_of_all)), 'dcc_filename')
ct_frac_st_bg <- left_join(ct_frac_st_bg, sData(geomx_obj)[, meta_names],
                        by = 'dcc_filename')

fwrite(ct_frac_st_bg, deconv_sd_path)
