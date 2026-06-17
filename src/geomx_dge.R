#README: script for performing DGE analysis between selected groups (both binary and multi-group options)
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
norm_type <- 'harmony_batch_corr_deseq2_vst' 

cofounder_name <- sample_name # better don't change - is added as a cofounder (random intercept in LLM model)
# remove samples with <2 nr of each AOI comparison group (not enough to compare, only adds noise)
min_aoi_nr <- 1

# for deconvolution DGE: remove AOIs with lower ct fraction - too unstable
min_ct_fraction <- 0.01

# if main_var_is_bin - wheter it should be used as regex or as it is
main_val_use_regex <- F

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'dge'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'dge', dge_name), showWarnings = T, recursive = T)

#TODO fix this - causing issues for dge_all
norm_is_log <- ifelse(norm_type %in% c('exprs', 'q3_norm', 'deseq2_norm'), FALSE, TRUE)
norm_name <- ifelse(norm_is_log, norm_type, paste0("log_", norm_type)) # TODO needed?

# path to cleaned scrna which should be calculated in deconvolution step
scrna_ref_cleaned_path <- file.path(output_dir, 'deconvolution', gsub('.RDS', '_cleaned_for_deconv.RDS', basename(scrna_ref_path)))

# path to deconvolution mtx
#TODO adjust to new naming and parse with deconv norm and batch corr
deconv_bp_path <- ifelse(grepl('harmony', norm_type), 
                         file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                   paste0('bp_res_', scrna_anno, '_expr_mtx_cleaned_deseq2_vst_harmony_corr.RDS')), 
                         file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                   paste0('bp_res_', scrna_anno, '_expr_mtx_cleaned_deseq2_vst_limma_corr_', 
                                          primary_batch_var, secondary_batch_var,
                                          '_cov_', covname, '.RDS'))) 

# path to deconvolution ct fractions 
deconv_bp_cellfrac_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                     paste0('bp_res_', scrna_anno, '_ct_fraction.csv'))
# deconv_bp_pulled_path <- file.path(output_dir, 'deconvolution', 'bayes_prism',
#                             paste0('bp_res_pseudosc_mid_lvl_ct_updated_ct_frac_0.005_deseq2_vst_harmony.csv'))

# load geomx obj from rds -------------------------------------------------

geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)

# merge geomx metadata with custom metadata
if(!is.null(custom_metadt_path)){
  custom_metadt <- fread(custom_metadt_path)
  if(aoi_id %in% colnames(custom_metadt)){
    pData(geomx_obj) <- left_join(pData(geomx_obj), custom_metadt, by = aoi_id, suffix = c("_orig", ""))
    colnames(geomx_obj) <- pData(geomx_obj)[[aoi_id]] # returning lost colnames
    pData(geomx_obj)[['Roi_geomx_factor']] <- as.factor(pData(geomx_obj)[[roi_id]])
    pData(geomx_obj)[['sample_roi_factor']] <- as.factor(pData(geomx_obj)[['sample_roi']])
  } else{
    stop('custom_metadt have to contain aoi_id column to be merged with metadata')
  }
}

expr_list <- list()

if('all' %in% dge_inp_data_type){
  
  # remove low complexity genes and change to log if needed
  if(file.exists(scrna_ref_cleaned_path)){
    scrna_ref_obj <- readRDS(scrna_ref_cleaned_path)
    geomx_filt <- remove_low_complex_and_noncoding_genes(geomx_obj, scrna_ref_obj, raw_counts_layer = 'counts')
    expr_mtx <- geomx_filt@assayData[[norm_type]]
  } else{
    expr_mtx <- geomx_obj@assayData[[norm_type]]
  }
  
  if(!norm_is_log){
    expr_mtx <- log2(expr_mtx + 1)
  }
  
  expr_list[[length(expr_list) + 1]] <- expr_mtx
  names(expr_list) <- 'dge_all'
}


# load deconvoluted signal ------------------------------------------------

if('bp' %in% dge_inp_data_type){
  
  # load deconvoluted profiles and ct fractions
  deconv_ct_list <- readRDS(deconv_bp_path)
  deconv_ct_frac <- fread(deconv_bp_cellfrac_path, select = c(aoi_id, ct_of_interest))
  deconv_ct_frac[is.na(deconv_ct_frac)] <- 0
  colnames(deconv_ct_frac) <- c(aoi_id, paste0(ct_of_interest, '_aoi_ct_frac'))
  
  # merge with metadata
  pData(geomx_obj) <- left_join(pData(geomx_obj), deconv_ct_frac, by = aoi_id, suffix = c("_orig", ""))
  colnames(geomx_obj) <- pData(geomx_obj)[[aoi_id]] # returning lost colnames
  
  # filter to cell types of interest
  # filter to AOIs with > min fraction of given cell
  deconv_ct_list_filt <- lapply(ct_of_interest, function(ct_name){
    deconv_ct <- deconv_ct_list[[ct_name]]
    deconv_ct <- deconv_ct[, deconv_ct_frac$dcc_filename[deconv_ct_frac[[paste0(ct_name, '_aoi_ct_frac')]] >= min_ct_fraction]]
    return(deconv_ct)
  })
  
  # TODO add error if name not in names from deconv list
  # TODO parse if norm is not log
  # TODO may throw an error if 0 AOI remain after filtering
  
  names(deconv_ct_list_filt) <- paste0('dge_deconv_', ct_of_interest)
  expr_list <- c(expr_list, deconv_ct_list_filt)
} 

# DGE with main variable comparison ---------------------------------------

# create formula for the LLM model:
# Sample is used as a mixed effect (cofounder)
if(comparison_type == 'within'){
  # within slide analysis - with random slope in LMM
  model_formula <- ~ main_var_factor + (1 + main_var_factor | cofounder_factor) # random slope + random intercept
  # correction for ct fraction in deconv
  #model_formula_ctfrac_corr <- ~ main_var_factor + ct_fraction + (1 + main_var_factor | cofounder_factor) 
  #reduced_model_formula <- ~ (1 + main_var_factor | cofounder_factor) # for testing if model add any information
} else if(comparison_type == 'between'){
  model_formula <- ~ main_var_factor + (1 | cofounder_factor) # random intercept
  # correction for ct fraction in deconv
  #model_formula_ctfrac_corr <- ~ main_var_factor + ct_fraction + (1 | cofounder_factor) 
  #reduced_model_formula <- ~ (1 | cofounder_factor)
} else{stop('comparison type can be either "within" or "between"')}

# stores information about eg comparion groups and samples
runlogs <- c()

# iterate through all + deconv matrices
lapply(names(expr_list), function(expr_name){
  print(paste0('########## ', expr_name, ' ###########'))

  # hacking GeoMx class object 
  newassay <- new.env(parent=geomx_obj@assayData)
  newassay[[expr_name]] <- expr_list[[expr_name]]
  
  geomx_obj_dge <- geomx_obj
  geomx_obj_dge@assayData <- newassay
  
  # filter to cells and genes which remained after deconvolution
  geomx_obj_dge <- geomx_obj_dge[rownames(expr_list[[expr_name]]),  colnames(expr_list[[expr_name]])]
  
  # add back pData which was removed during filtering.. geomxtools bug
  pData(geomx_obj_dge) <- pData(geomx_obj)[pData(geomx_obj)[[aoi_id]] %in% colnames(expr_list[[expr_name]]), ]
  
  pData(geomx_obj_dge) <- prepare_dge_metadata(pData(geomx_obj_dge), main_var_name, main_var_is_bin, main_var_main_val,
                                               dge_categories, cofounder_name, main_val_use_regex) 
  
  
  dge_results <- data.frame()
  
  # iterate through data groups
  for(data_group in unique(pData(geomx_obj_dge)[['dge_group']])){
    
    print(data_group)
    runlogs <- c(runlogs, paste('dge run info for', expr_name, data_group, ':'))
    
    # filter to group
    ind <- pData(geomx_obj_dge)$dge_group == data_group
    geomx_obj_dge_group <- geomx_obj_dge[, ind]
    
    # TODO optimise this logic
    # for 'within' comparison remove samples which desn't contain enough nr of AOI per group
    cleaned_dt <- rm_too_small_groups(geomx_obj_dge_group, min_aoi_nr, main_var_is_bin, comparison_type)
    geomx_obj_dge_group_cleaned <- cleaned_dt$geomx_obj
    runlogs <- c(runlogs, cleaned_dt$logs)
    
    # if >1 main_variable value present in cleaned_dt, make dge
    if(length(unique(pData(geomx_obj_dge_group_cleaned)$main_var_factor)) > 1){
      
      # run LMM:
      # formula follows conventions defined by the lme4 package
      mixed_result <- tryCatch({
        mixedOutmc <- mixedModelDE2(
          geomx_obj_dge_group_cleaned,
          elt = expr_name,
          modelFormula = model_formula, 
          groupVar = 'main_var_factor',
          nCores = (parallel::detectCores() - 2),
          multiCore = unname(ifelse(Sys.info()['sysname'] == 'Windows', FALSE, TRUE))
        )
        mixedOutmc  # Return the result of mixedModelDE
      }, error = function(e) {
        # Return an empty dataframe if an error occurs eg to little AOIs
        print('not enough AOI for comparison!')
        data.frame()
        
      })
    } else{
      mixed_result <- data.frame()
    }
    
    gc()
    
    if(nrow(mixed_result) > 1){
      # format results as data.frame
      r_test <- do.call(rbind, mixed_result["lsmeans", ])
      tests <- rownames(r_test)
      r_test <- as.data.frame(r_test)
      r_test$Contrast <- tests
      
      # use lapply in case you have multiple levels of your test factor to
      # correctly associate gene name with it's row in the results table
      r_test$Gene <- unlist(lapply(colnames(mixed_result), function(x){
        lsmeans_res <- mixed_result["lsmeans", ][[x]]
        
        if(nrow(lsmeans_res) == 1 & all(is.na(lsmeans_res))){
          gene_rep <- NA
        } else{
          gene_rep <- rep(x, nrow(lsmeans_res))
        }
      }))
      
      r_test <- r_test[!is.na(r_test$Gene), ]
      
      # don't merge if NA results for all genes
      if(nrow(r_test) >0){
        r_test$data_group <- data_group
        r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
        r_test <- r_test[, c("Gene", "data_group",  "Contrast", "Estimate",
                             "Pr(>|t|)", "FDR")]
        dge_results <- rbind(dge_results, r_test)
      } else{
        err <- paste('no computed dge for', expr_name, data_group, 'probable reason: too little signal from given cell type')
        print(err)
        runlogs <- c(runlogs, err)
        dge_results <- dge_results
      }

    } else{
      err <- paste('error while computing dge for', expr_name, data_group, 'probably too little AOI for comparison. Check the comparison groups!!')
      print(err)
      runlogs <- c(runlogs, err)
      dge_results <- dge_results
    }
    
  }
  
  # write results table 
  out_path <- file.path(output_dir, 'dge', dge_name, paste0(expr_name, '_', dge_name, '.csv'))
  if(nrow(dge_results) > 1){
    fwrite(dge_results, out_path)
    print(out_path)
    print(paste0('results saved for ', expr_name))}
  
  
  
  # make volcano plots for visualisation ------------------------------------

  if(nrow(dge_results) > 1){  
    dir.create(file.path(output_dir, 'dge', dge_name, expr_name))
    
    for(dt_group in unique(dge_results$data_group)){
      print(dt_group)
      dge_results_group <- dge_results[dge_results$data_group == dt_group, ]
      
      for(cont in unique(dge_results_group$Contrast)){
        
        if(cont != ''){ # artifact from returning NA when lmm fails
          dge_results_group_cont <- dge_results_group[dge_results_group$Contrast == cont, ]
          groups <- strsplit(cont, split = ' - ', fixed = T)
          
          plot_volcano_deg(dge_results_group_cont, dt_group, 20, groups[[1]][1], groups[[1]][2],
                           file.path(output_dir, 'dge', dge_name, expr_name))
        }
      }
    }
  }
  
  })

# write logs with parameters ----------------------------------------------

writeLines(c('DGE logs:',
             'Comparison type: ', comparison_type, 
             '; ', main_var_name, ' bin ', main_var_is_bin, 
             '; main var value: ', main_var_main_val,
             '; categories to compare: ', dge_categories,
             '; min number of categories to compare within slide: ', min_aoi_nr,
             '; normalisation type: ', norm_name,
             '; low complex gene removed : ', ifelse(file.exists(scrna_ref_cleaned_path), 'TRUE', 'FALSE'),
             '; deconv mtx used : ', deconv_bp_path,
             'runlogs: ', runlogs), dge_logs_path)
