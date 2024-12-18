library(standR)
library(limma)


geomx_obj <- readRDS(geomx_norm_path)

neg_control_subset <- negativeControlSubset(geomx_obj)


geomx_batch_correct <- geomxBatchCorrection(geomx_obj@assayData$q3_norm, k = 0,
                                            factors = 'Segment', NCGs = NULL,
                                batch = sData(geomx_obj)$`Slide Name`, method = "Limma",
                                design = model.matrix(~Segment, data = sData(geomx_obj)))


  assay(spe, "logcounts") <- (function(.) .[rownames(spe), 
                                            colnames(spe)])(limma::removeBatchEffect(assay(spe, 
                                                                                           n_assay), batch = batch, batch2 = batch2, covariates = covariates, 
                                                                                     design = design)
                                                            
    
kk <- limma:removeBatchEffect(geomx_obj@assayData$q3_norm)