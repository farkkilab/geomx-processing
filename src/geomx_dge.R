
# WARNING: DGE with mixed model will take around 30G RAM
# best to run in >10 cores

# problems with Matrix package - i has to be lower that 1.7 to work with lmer
#devtools::install_version("Matrix","1.6.4")

library(data.table)
#library(clusterProfiler)
# library(msigdbr)
# library(progeny)
# library(biomaRt)
# library(GSVA)
# library(ggpubr)
# library(topGO)
# library(fgsea)
# library(fpc)
# library(dbscan)

#TODO add visualsation to DGE script
#TODO clean calling variables
#TODO adjust for deconvoluted data as well


# get variables -----------------------------------------------------------

# comparison_type <- 'within' 
# 'within' when you compare different ROI types within sample
# 'between' - comparisons between slides

# cofounder_name <- 'Sample' # don't change it
# then 'Sample' is added as a cofounder (random intercept in LLM model)

# main_var_name <- "Annotation_cell" 
# main_var - main variable to make comparison 

# main_var_is_bin <- TRUE 
# if main_var_is_bin is True, main_var_main_val will be compared 
# with all other categories in main_var
# if False - each category in main_var will be compared with every other one
# main_var_main_val is set to NULL

# main_var_main_val <- 'CD4_CD8_CD11_Iba1' 
# main_var_main_val value among main_var which needs to be compare against all other vals
# or part of the value eg 'CD8' within values for comparison


# dge_categories <- c('Segment', 'NACT status')
# dge_categories - all conditions for which we want to make DGE separately

################
# comparison_type <- 'within'
# cofounder_name <- 'Sample'
# main_var_name <- "Annotation_cell" 
# main_var_is_bin <- TRUE 
# main_var_main_val <- 'CD4_CD8_CD11_Iba1' 
# dge_categories <- c('Segment', 'NACT status')
################
# comparison_type <- 'between'
# cofounder_name <- 'Sample'
# main_var_name <- "NACT status"
# main_var_is_bin <- FALSE
# dge_categories <- c('Segment', 'Annotation_cell')
###############

norm_type <- 'limma_batch_corr' # best on batch-effect corrected data: 'limma_batch_corr' or 'harmony_batch_corr'
norm_is_log <- TRUE # if normalised expr matrix is in the log scale, both limma and harmony batch corr are

multicore = TRUE # if Linux or macOS, for Windows multicore = FALSE

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'dge'), showWarnings = T, recursive = T)

#source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')

# load geomx obj from rds -------------------------------------------------

geomx_obj <- readRDS(geomx_now_path)

if(!norm_is_log){
  # convert normalized counts to log scale
  assayDataElement(object = geomx_obj, elt = paste0("log_", norm_type)) <-
    assayDataApply(geomx_obj, 2, FUN = log, base = 2, elt = norm_type)
  norm_type <- paste0("log_", norm_type)
}


# DGE with main variable comparison ---------------------------------------

if(main_var_is_bin){
  # make binary vector - either main variable has the desired value or not
  pData(geomx_obj)$main_var <- ifelse(grepl(main_var_main_val, pData(geomx_obj)[, main_var_name]),
                                          main_var_main_val, 'other_roi_type')
} else{
  pData(geomx_obj)$main_var <- pData(geomx_obj)[, main_var_name]
}


# convert test variables to factors
for(col in c(dge_categories, 'main_var')){
  pData(geomx_obj)[[paste0(col, "_factor")]] <- factor(pData(geomx_obj)[[col]])
}

pData(geomx_obj)$cofounder_factor <- factor(pData(geomx_obj)[[cofounder_name]])

# make variable with all dge categories
pData(geomx_obj)$dge_group <- apply(pData(geomx_obj), 1, function(row){
  group <- sapply(dge_categories, function(var){
    paste(row[var])
  })
  group <- paste(group, collapse = '_')
  return(group)
})

# create formula for the LLM model:
# Sample is used as a mixed effect (cofounder)
if(comparison_type == 'within'){
  # within slide analysis - with random slope in LLM
  model_formula <- ~ main_var_factor + (1 + main_var_factor | cofounder_factor) # random slope + random intercept
} else if(comparison_type == 'between'){
  model_formula <- ~ main_var_factor + (1 | cofounder_factor) # random intercept
} else{stop('comparison type can be either "within" or "between"')}


# run LMM:
# formula follows conventions defined by the lme4 package

dge_results <- c()
for(data_group in unique(pData(geomx_obj)[, 'dge_group'])){
  
  print(data_group)
  ind <- geomx_obj@phenoData@data$dge_group == data_group
  
  mixed_result <- tryCatch({
    mixedOutmc <- mixedModelDE(
      geomx_obj[, ind],
      elt = norm_type,
      modelFormula = model_formula, 
      groupVar = 'main_var_factor',
      nCores = (parallel::detectCores() - 2),
      multiCore = multicore
    )
    mixedOutmc  # Return the result of mixedModelDE
  }, error = function(e) {
    # Return an empty dataframe if an error occurs eg to little ROIs
    data.frame()
    
  })
  
  gc()
  
  if(nrow(mixed_result) > 1){
    # format results as data.frame
    r_test <- do.call(rbind, mixed_result["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests
    
    # use lapply in case you have multiple levels of your test factor to
    # correctly associate gene name with it's row in the results table
    r_test$Gene <-
      unlist(lapply(colnames(mixed_result),
                    rep, nrow(mixed_result["lsmeans", ][[1]])))
    r_test$data_group <- data_group
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
    r_test <- r_test[, c("Gene", "data_group",  "Contrast", "Estimate",
                         "Pr(>|t|)", "FDR")]
    dge_results <- rbind(dge_results, r_test)
  } else{
    print('error while computing dge. probably too little ROI for comparison')
    dge_results <- dge_results
  }
  
}


# write results table -----------------------------------------------------

fwrite(dge_results, file.path(output_dir, 'dge', 
                          paste('dge_', comparison_type, '_slide_', main_var_name, 
                                 '_bin_', main_var_is_bin, '_', paste0(dge_categories, collapse = '_'), 
                                '.csv')))


# write logs with parameters ----------------------------------------------

writeLines(c('DGE logs:',
             'Comparison type: ', comparison_type, 
             '; ', main_var_name, ' bin ', main_var_is_bin, '; categories to compare: ', dge_categories), dge_logs_path)
close(dge_logs_path)