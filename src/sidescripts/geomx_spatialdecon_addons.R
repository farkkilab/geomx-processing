######################################################################
######################################################################
# alternatives (probably less usefull) to run SpatialDecon

######################################################
######################################################
# all bells and whistles from the vignette
# prepare data for SpatialDecon -------------------------------------------
# from
# https://bioconductor.org/packages/release/bioc/vignettes/SpatialDecon/inst/doc/SpatialDecon_vignette_NSCLC.html

#TODO check geomx_obj and geomx_obj_filtered 

featureType(geomx_obj) <- "Target"

sampleNames(geomx_obj) <- sData(geomx_obj)[['dcc_filename']]

# get negative probes (aggregated to features already) names
negativeProbefData <- subset(fData(geomx_obj), CodeClass == "Negative")

# estimate bcg for every segment based on neg probes 
# TODO re-check if >1 module
geomx_bg <- derive_GeoMx_background(norm = geomx_obj@assayData[[norm_type]],
                                    probepool = fData(geomx_obj)$Module,
                                    negnames = negativeProbefData$TargetName)

# load pre-defined TME cell profile matrix
tme_mtx <- download_profile_matrix(species = "Human",
                                   age_group = "Adult", 
                                   matrixname = "ImmuneTumor_safeTME")

data("safeTME")
data("safeTME.matches")

# prepare cell profile matrix from reference scRNAseq ---------------------

scrna_ref_obj <- readRDS(scrna_ref_path)

# format annotations
scrna_anno_dt <- scrna_ref_obj@meta.data[, c('cell_name', scrna_anno)]
rownames(scrna_anno_dt) <- NULL
colnames(scrna_anno_dt) <- c('cell_name', 'cell_type')

# TODO examine scalingFactor: 1 or 5 or what?
custom_oc_mtx <- create_profile_matrix(mtx = scrna_ref_obj@assays$SCT@data,            # cell x gene count matrix
                                       cellAnnots = scrna_anno_dt,  # cell annotations with cell type and cell name as columns
                                       cellTypeCol = "cell_type",  # column containing cell type
                                       cellNameCol = "cell_name",           # column containing cell ID/name
                                       matrixName = "oc_scrnaseq_ref_cell_type_sct", # name of final profile matrix
                                       outDir = output_dir,                    # path to desired output directory, set to NULL if matrix should not be written
                                       normalize = FALSE,                # Should data be normalized?
                                       minCellNum = 50,                   # minimum number of cells of one type needed to create profile, exclusive
                                       minGenes = 10,                    # minimum number of genes expressed in a cell, exclusive
                                       scalingFactor = 1,                # what should all values be multiplied by for final matrix
                                       discardCellTypes = TRUE)          # should cell types be filtered for types like mitotic, doublet, low quality, unknown, etc.



# run basic SpatialDecon --------------------------------------------------
decon_res <-  runspatialdecon(object = geomx_obj,
                              norm_elt = norm_type,
                              raw_elt = "exprs",
                              X = tme_mtx,
                              align_genes = TRUE)

#heatmap(t(decon_res$beta), cexCol = 0.5, cexRow = 0.7, margins = c(10,7))

# run extended SpatialDecon -----------------------------------------------

#give info about tumor
geomx_obj$istumor = geomx_obj$Segment == "tumor"

# TODO nuclei counts from geomx are unreliable - match with info from cycif
# TODO examine if istumor should be used - it's not pure in our case
# TODO examine n_tumor_clusters param with different n

decon_res_ext <- runspatialdecon(object = geomx_obj,
                                 norm_elt = norm_type,                # normalized data
                                 raw_elt = "exprs",                      # expected background counts for every data point in norm
                                 X = safeTME,                            # safeTME matrix, used by default
                                 cellmerges = safeTME.matches,           # safeTME.matches object, used by default
                                 #cell_counts = geomx_obj$Nuclei,      # nuclei counts, used to estimate total cells
                                 #is_pure_tumor = geomx_obj$istumor,   # identities of the Tumor segments/observations
                                 n_tumor_clusters = 5)               # how many distinct tumor profiles to append to safeTME

heatmap(sweep(decon_res_ext@experimentData@other$SpatialDeconMatrix, 1, apply(decon_res_ext@experimentData@other$SpatialDeconMatrix, 1, max), "/"),
        labRow = NA, margins = c(10, 5))

# run extended SpatialDecon with custom oc mtx ----------------------------

decon_res_custom <- runspatialdecon(object = geomx_obj,
                                    norm_elt = norm_type,                # normalized data
                                    raw_elt = "exprs",                      # expected background counts for every data point in norm
                                    X = custom_oc_mtx,                            # safeTME matrix, used by default
                                    #cell_counts = geomx_obj$Nuclei,      # nuclei counts, used to estimate total cells
                                    #is_pure_tumor = geomx_obj$istumor,   # identities of the Tumor segments/observations
                                    n_tumor_clusters = 5)               # how many distinct tumor profiles to append to safeTME

# run extended SpatialDecon with custom oc mtx and estimated bg------------

#TODO to use bg geomx obj have to be converted to seurat
decon_res_custom_bg <- runspatialdecon(object = geomx_obj,
                                       bg = geomx_bg,                      # expected background counts for every data point in norm
                                       X = custom_oc_mtx,                            # safeTME matrix, used by default
                                       #cell_counts = geomx_obj$Nuclei,      # nuclei counts, used to estimate total cells
                                       #is_pure_tumor = geomx_obj$istumor,   # identities of the Tumor segments/observations
                                       n_tumor_clusters = 5)               # how many distinct tumor profiles to append to safeTME


# save spatialdecon results -----------------------------------------------

# colnames from output
res_cols <- c("beta", "p", "t", "se", "prop_of_all", "prop_of_nontumor")

res_ext <- pData(decon_res_ext)[, c(res_cols, "sigma")]
saveRDS(res_ext, file = file.path(output_dir, 'deconvolution', 'spatial_decon', scrna_anno, 
                                  'spat_dec_res_ext.rds'))

res_custom <- pData(decon_res_custom)[, c(res_cols, "sigmas")]
saveRDS(res_ext, file = file.path(output_dir, 'deconvolution', 'spatial_decon', scrna_anno, 
                                  'spat_dec_res_custom.rds'))


