library(standR)
library(limma)
library(pvca)
library(harmony)
library(dplyr)

#########################
proj_dir <- '~/Documents/phd/st'
output_dir <- file.path(proj_dir, 'geomx-processing', 'results', 'batch2')
geomx_norm_path <<- file.path(output_dir, 'geomx_qc_norm.RDS')


geomx_obj <- readRDS(geomx_norm_path)

# make log2 transformed normalised counts
expr_norm_log <- log2(geomx_obj@assayData$deseq2_norm + 1)

# make expression sets for PVCA
phenoData <- new("AnnotatedDataFrame", data=geomx_obj@phenoData@data, 
                 varMetadata=geomx_obj@phenoData@varMetadata)
featureData <- new("AnnotatedDataFrame", data=geomx_obj@featureData@data, 
                 varMetadata=geomx_obj@featureData@varMetadata)

exprset_deseq2_norm_log <- ExpressionSet(assayData=expr_norm_log, 
                                         phenoData = phenoData,
                                         featureData = featureData)

########################
# check batch effect with PVCA

pData(geomx_obj) <- rename(pData(geomx_obj), Slide_Name = `Slide Name`)

# change vars into factors
batch_factors <- c('Slide_Name', 'Sample', 'Patient', 'Site', 'batch_nr')
for(colname in batch_factors){
  pData(geomx_obj)[[paste0(colname, '_factor')]] <- as.factor(pData(geomx_obj)[[colname]])
}
batch_factors_names <- paste0(batch_factors, '_factor')


pct_threshold <- 0.6

pvcaObj <- pvcaBatchAssess(exprset_deseq2_norm_log, batch_factors_names, pct_threshold) 


bp <- barplot(pvcaObj$dat,  xlab = "Effects",
              ylab = "Weighted average proportion variance", ylim= c(0,1.1),
              col = c("blue"), las=2, main="PVCA estimation bar chart")
axis(1, at = bp, labels = pvcaObj$label, xlab = "Effects", cex.axis = 0.5, las=2)
values = pvcaObj$dat
new_values = round(values , 3)
text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.8) 

##########################
# remove batch effect with limma
batch <-  sData(geomx_obj)$`Slide Name`
batch2 <-  sData(geomx_obj)$batch_nr
design <- model.matrix(~Segment, data = sData(geomx_obj))
cov <- NULL

limma_res1 <- limma::removeBatchEffect(expr_norm_log, batch = batch, covariates = cov,
                               design = design)

# better
limma_res2 <- limma::removeBatchEffect(expr_norm_log, batch = batch2, covariates = cov,
                                                design = design)

limma_res3 <- limma::removeBatchEffect(expr_norm_log, batch = batch2, batch2 = batch, covariates = cov,
                                                design = design)

##########################
# remove batch effect with harmony
meta_dt <- pData(geomx_obj)[, c(batch_factors, 'dcc_filename', 'Segment', 'NACT status', 'Annotation_cell')]

# harmony res have to be flipped
#ok
harmony_res1 <- t(HarmonyMatrix(expr_norm_log, 
                                 meta_data = meta_dt,
                                 vars_use = 'Slide_Name'))


harmony_res2 <- t(HarmonyMatrix(expr_norm_log, 
                                 meta_data = meta_dt,
                                 vars_use = 'batch_nr'))
# best
harmony_res3 <- t(HarmonyMatrix(expr_norm_log, 
                                meta_data = meta_dt,
                                vars_use = c('Slide_Name', 'batch_nr')))


##########################
# check PVCA after batch effect removal

expr_after_batch_corr <- harmony_res3

exprset_after_batch_corr <- ExpressionSet(assayData=expr_after_batch_corr, 
                                         phenoData = phenoData,
                                         featureData = featureData)

pvcaObj2 <- pvcaBatchAssess(exprset_after_batch_corr, batch_factors_names, pct_threshold) 


bp2 <- barplot(pvcaObj2$dat,  xlab = "Effects",
              ylab = "Weighted average proportion variance", ylim= c(0,1.1),
              col = c("blue"), las=2, main="PVCA estimation bar chart")
axis(1, at = bp2, labels = pvcaObj2$label, xlab = "Effects", cex.axis = 0.5, las=2)
values = pvcaObj2$dat
new_values = round(values , 3)
text(bp2,pvcaObj2$dat,labels = new_values, pos=3, cex = 0.7) 


####################################
####################################
#############################
# check scatter plots per PCA
#TODO rewrite from source
plotPairPCA(spe_lrb, assay = 2, color = Type, title = "Limma removeBatch")

p1 <- do_scatter(harmony_results, meta_dt, 'Slide_Name') + 
  labs(title = 'Colored by Slide_Name')

print(p1)

p2 <- do_scatter(harmony_results2, meta_dt, 'batch_nr') + 
  labs(title = 'Colored by batch_nr')

print(p2)


#######################
#######################
#######################
#####################
# do_scatter function from harmony

do_scatter <- function(xy, meta_data, label_name, base_size = 12) {    
  palette_use <- c(`jurkat` = '#810F7C', `t293` = '#D09E2D',`half` = '#006D2C')
  xy <- xy[, 1:2]
  colnames(xy) <- c('X1', 'X2')
  plt_df <- xy %>% data.frame() %>% cbind(meta_data)
  plt <- ggplot(plt_df, aes(X1, X2, col = !!rlang::sym(label_name), fill = !!rlang::sym(label_name))) + 
    theme_test(base_size = base_size) +
    guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1,
                                                    shape = 16, size = 4))) +
    scale_color_manual(values = palette_use) +
    scale_fill_manual(values = palette_use) +
    theme(plot.title = element_text(hjust = .5)) +
    labs(x = "PC 1", y = "PC 2") +
    theme(legend.position = "none") +
    geom_point(shape = '.')
  
  ## Add labels
  data_labels <- plt_df %>%
    dplyr::group_by(!!rlang::sym(label_name)) %>%
    dplyr::summarise(X1 = mean(X1), X2 = mean(X2)) %>%
    dplyr::ungroup()
  plt + geom_label(data = data_labels, aes(label = !!rlang::sym(label_name)), 
                   color = "white", size = 4)
}
