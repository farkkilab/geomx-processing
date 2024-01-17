# TODO from the previous script
output_rds_path <- file.path(output_dir, 'geomx_qc.RDS')

geomx_qc <- readRDS(output_rds_path)

# Q3 normalisation --------------------------------------------------------

plot_q3_stats(geomx_obj, main_var, file.path(output_dir, 'qc/q3_stats.png'))


geomx_obj <- normalize(geomx_obj ,
                       norm_method = "quant", 
                       desiredQuantile = .75,
                       toElt = "q3_norm")

# quantile normalisation --------------------------------------------------

norm.quantile = normalize.quantiles(as.matrix(geomx_obj@assayData$exprs))
dimnames(norm.quantile) = dimnames(geomx_obj@assayData$exprs)

# hacking GeoMx class object 
# TODO this is experimental - newassay is not identical and it may cause problems
# if so, store this in another mtx and use when needed
newassay <- new.env(parent=geomx_obj@assayData)
newassay$exprs <- geomx_obj@assayData$exprs
newassay$q3_norm <- geomx_obj@assayData$q3_norm
newassay$quant_norm <- norm.quantile

geomx_obj@assayData <- newassay


# plot effects of normalisation -------------------------------------------

plot_norm_effect(exprs(geomx_obj)[,1:10], 'Raw Counts', file.path(output_dir, 'qc/norm_raw.png'))


plot_norm_effect(assayDataElement(geomx_obj[,1:10], elt = "q3_norm"),
                 'Q3 normalised', file.path(output_dir, 'qc/norm_q3.png'))

# TODO I don't like sth with this plot, why all outliers are the same in each segment?
plot_norm_effect(assayDataElement(geomx_obj[,1:10], elt = "quant_norm"),
                 'Quantile normalised', file.path(output_dir, 'qc/norm_quant.png'))


# make UMAP and t-SNE -----------------------------------------------------

#TODO change for any segment type
# divide for tumor and stroma and do dimentionality reduction for all
geomx_obj_tumor <- geomx_obj[, geomx_obj@phenoData@data$Segment == "tumor"]
geomx_obj_stroma <- geomx_obj[, geomx_obj@phenoData@data$Segment == "stroma"]

geomx_list <- list(all = geomx_obj, tumor = geomx_obj_tumor, stroma = geomx_obj_stroma)

geomx_list_dim_red <- lapply(1:length(geomx_list), function(n){
  
  geomx <- geomx_list[[n]]
  
  # run UMAP and tSNE on Q3 and quantile norm
  for(norm in c('q3_norm', 'quant_norm')){
    # update defaults for umap to contain a stable random_state (seed)
    custom_umap <- umap::umap.defaults
    custom_umap$random_state <- 42
    
    umap_out <-
      umap(t(log2(assayDataElement(geomx , elt = norm))),  
           config = custom_umap)
    
    # save UMAP1 and 2 results to pData
    pData(geomx)[, c(paste0("UMAP1_", norm), paste0("UMAP2_", norm))] <- umap_out$layout[, c(1,2)]
    
    # set the seed for tSNE as well
    set.seed(42) 
    tsne_out <-
      Rtsne(t(log2(assayDataElement(geomx , elt = norm))),
            perplexity = ncol(geomx)*.15)
    
    # save tSNE1 and 2 results to pData
    pData(geomx)[, c(paste0("tSNE1_", norm), paste0("tSNE2_", norm))] <- tsne_out$Y[, c(1,2)]
  }
  
  # generate umap and tsne plots and color by variables
  for(method in c('UMAP', 'tSNE')){
    for(norm in c('q3', 'quant')){
      for(color_var in c('Annotation_cell', 'Patient', 'NACT status', 'PFS', 'Site', 'Sample')){
        plot_umap_tsne(pData(geomx), method_type = method, 
                       norm_type = norm, color_var = color_var,
                       output_name = file.path(output_dir, 'umap_tsne2', names(geomx_list)[n], 
                                               paste0(method, '_', norm, '_', color_var, '.pdf')))
      }
    }
  }
  
  return(geomx)
})

# update objects
geomx_obj <- geomx_list_dim_red[[1]]
geomx_obj_tumor <- geomx_list_dim_red[[2]]
geomx_obj_stroma <- geomx_list_dim_red[[3]]

rm(geomx_list)
rm(geomx_list_dim_red)



