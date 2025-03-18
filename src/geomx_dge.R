# WARNING: DGE with mixed model will take around 30G RAM
# best to run in >10 cores

# get variables -----------------------------------------------------------

################
# comparison_type <- 'within'
# cofounder_name <- 'Sample'
# main_var_name <- "Annotation_cell" 
# main_var_is_bin <- TRUE 
# main_var_main_val <- 'CD4_CD8_CD11_Iba1' # 'CD8_.*Iba1'
# dge_categories <- c('Segment', 'NACT status')
################
# comparison_type <- 'between'
# cofounder_name <- 'Sample'
# main_var_name <- "NACT_status"
# main_var_is_bin <- FALSE
# dge_categories <- c('Segment', 'Annotation_cell')
###############

# best on batch-effect corrected data: 'limma_batch_corr' or 'harmony_batch_corr' (both log)
# q3 also ok but its not batch corrected
norm_type <- 'harmony_batch_corr' 

cofounder_name <- 'Sample' # better don't change - is added as a cofounder (random intercept in LLM model)

# path to cleaned scrna which should be calculated in deconvolution step
scrna_ref_cleaned_path <- file.path(output_dir, 'deconvolution', gsub('.RDS', '_cleaned_for_deconv.RDS', basename(scrna_ref_path)))

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'dge'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'dge', dge_name), showWarnings = T, recursive = T)

norm_is_log <- ifelse(norm_type %in% c('exprs', 'q3_norm', 'deseq2_norm'), FALSE, TRUE)
deconv_bp_path <- ifelse(grepl('harmony', norm_type), deconv_bp_harm_path, deconv_bp_limma_path)

# load geomx obj from rds -------------------------------------------------

geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)

low_complex_rmv <- ifelse(file.exists(scrna_ref_cleaned_path), TRUE, FALSE)
norm_name <- ifelse(norm_is_log, norm_type, paste0("log_", norm_type))

expr_list <- list()

if('all' %in% dge_inp_data_type){
  # remove low complexity genes and change to log if needed
  expr_mtx <- prepare_expr_mtx(geomx_norm_batch_eff_rm_path, norm_type, norm_is_log, 
                               scrna_ref_cleaned_path)
  
  expr_list[[length(expr_list) + 1]] <- expr_mtx
  names(expr_list) <- 'dge_all'
}


# load deconvoluted signal ------------------------------------------------

if('bp' %in% dge_inp_data_type){
  
  deconv_ct_list <- readRDS(deconv_bp_path)
  
  # artificially add missing AOIs to prevent issues with geomx object
  deconv_ct_list_padded <- lapply(deconv_ct_list, function(expr){

    if(!identical(colnames(expr), colnames(geomx_obj@assayData[[norm_type]]))){
      expr <- as.data.frame(expr)
      # add empty columns
      expr[, setdiff(colnames(geomx_obj@assayData[[norm_type]]), colnames(expr))] <- NA

      # merge with expr and ensure order
      expr <- as.matrix(expr[, colnames(geomx_obj@assayData[[norm_type]])])

      stopifnot(identical(colnames(expr), colnames(geomx_obj@assayData[[norm_type]])))
    }

    return(expr)
  })

  names(deconv_ct_list_padded) <- paste0('dge_deconv_', names(deconv_ct_list))
  
  expr_list <- c(expr_list, deconv_ct_list_padded)
}

# DGE with main variable comparison ---------------------------------------

# create formula for the LLM model:
# Sample is used as a mixed effect (cofounder)
if(comparison_type == 'within'){
  # within slide analysis - with random slope in LLM
  model_formula <- ~ main_var_factor + (1 + main_var_factor | cofounder_factor) # random slope + random intercept
  reduced_model_formula <- ~ (1 + main_var_factor | cofounder_factor) # for testing if model add any information
} else if(comparison_type == 'between'){
  model_formula <- ~ main_var_factor + (1 | cofounder_factor) # random intercept
  reduced_model_formula <- ~ (1 | cofounder_factor)
} else{stop('comparison type can be either "within" or "between"')}


# iterate through all + deconv matrices
lapply(names(expr_list), function(expr_name){
  print(expr_name)
  
  # hacking GeoMx class object 
  newassay <- new.env(parent=geomx_obj@assayData)
  newassay[[expr_name]] <- expr_list[[expr_name]]
  
  geomx_obj_dge <- geomx_obj
  geomx_obj_dge@assayData <- newassay
  
  pData(geomx_obj_dge) <- prepare_dge_metadata(pData(geomx_obj), main_var_name, main_var_is_bin, main_var_main_val,
                                               dge_categories, cofounder_name) 
  
  dge_results <- c()
  
  # iterate through data groups
  for(data_group in unique(pData(geomx_obj_dge)[, 'dge_group'])){
    
    print(data_group)
    
    # filter to group
    ind <- pData(geomx_obj_dge)$dge_group == data_group
    geomx_obj_dge_group <- geomx_obj_dge[, ind]
    
    # remove samples with <2 nr of each AOI comparison group (not enough to compare, only adds noise)
    #TODO thr 2 or 1?
    min_aoi_nr <- 2
    geomx_obj_dge_group_cleaned <- rm_too_small_groups(geomx_obj_dge_group, min_aoi_nr, main_var_is_bin, comparison_type)
    
    ###########################
    # TODO likelihood ratio test - anova(full model, reduced model)
    # check if main_var significantly improved the effect
    # but it have to be done per gene - maybe worth to check afterwards for interesting genes
    
    # genename <- 'ACAT1'
    # dat <- data.frame(expr = geomx_obj_dge_group_cleaned@assayData[[expr_name]][genename, ], 
    #                   pData(geomx_obj_dge_group_cleaned)[, c('main_var_factor', 'cofounder_factor')])
    # 
    # full_model <- lmerTest::lmer(formula(expr ~ main_var_factor + (1 + main_var_factor | cofounder_factor)), dat)
    # reduced_model <- lmerTest::lmer(formula(expr ~ (1 + main_var_factor | cofounder_factor)), dat)
    # anova_results <- anova(reduced_model, full_model)
    #########################
    

    # run LMM:
    # formula follows conventions defined by the lme4 package
    mixed_result <- tryCatch({
      mixedOutmc <- mixedModelDE(
        geomx_obj_dge_group_cleaned,
        elt = expr_name,
        modelFormula = model_formula, 
        groupVar = 'main_var_factor',
        nCores = (parallel::detectCores() - 2),
        multiCore = unname(ifelse(Sys.info()['sysname'] == 'Windows', FALSE, TRUE))
      )
      mixedOutmc  # Return the result of mixedModelDE
    }, error = function(e) {
      # Return an empty dataframe if an error occurs eg to little ROIs
      print('wtf')
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
  
  # write results table 
  out_path <- file.path(output_dir, 'dge', dge_name, paste0(expr_name, '_', dge_name, '.csv'))
  if(nrow(dge_results) > 1){fwrite(dge_results, out_path)}
  
  print(paste0('results saved for ', expr_name))
  
  # make volcano plots for visualisation ------------------------------------
  
  dir.create(file.path(output_dir, 'dge', dge_name, expr_name))
  
  for(dt_group in unique(dge_results$data_group)){
    print(dt_group)
    dge_results_group <- dge_results[dge_results$data_group == dt_group, ]
    
    for(cont in unique(dge_results_group$Contrast)){
      dge_results_group_cont <- dge_results_group[dge_results_group$Contrast == cont, ]
      groups <- strsplit(cont, split = ' - ', fixed = T)
      
      plot_volcano_deg(dge_results_group_cont, dt_group, 20, groups[[1]][1], groups[[1]][2],
                       file.path(output_dir, 'dge', dge_name, expr_name))
    }
  }

})

# write logs with parameters ----------------------------------------------

writeLines(c('DGE logs:',
             'Comparison type: ', comparison_type, 
             '; ', main_var_name, ' bin ', main_var_is_bin, 
             '; main var value: ', main_var_main_val,
             '; categories to compare: ', dge_categories,
             '; normalisation type: ', norm_name,
             '; low complex gene removed : ', low_complex_rmv), dge_logs_path)
#close(dge_logs_path)
