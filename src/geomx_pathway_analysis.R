# README: script to perform ssgsea/gsva 

# TODO low complex removal should be done soewhere else and saved to the main geomx object
# TODO better output_names - deconv and all may have diff norm types!
# get variables -----------------------------------------------------------

# variables to merge the final csv with
meta_names <- c(aoi_id, roi_id, aoi_segment_var, sample_name, main_experimental_condition, 
               main_roi_label, other_vars_bio)

# best to use batch effect corrected or at least vst data (all in log form) 
norm_type <- 'harmony_batch_corr_q3_norm' # from geomx assays
deconv_norm_type <- 'deseq2_vst' # c('q3_norm', 'deseq2_vst') which norm should be used for bayesprism results
deconv_batch_rm_type <- 'harmony' # c('harmony', 'limma')

# whethr or not rmv low complexity and non-coding genes from full signal geomx obj  (as for bp deconvolution)
low_complex_rmv <- TRUE 

adj_synonym <- T # whether or not adjust synonyms genes
# around 300 genes can be rescued this way but ensembl does not always work
# if there are issues, turn it off
min_sign_gene_nr <- 5 # signatures with less nr of genes will be removed, 5 is min in msigdb

compute_hallmark <- T
# should GSEA for msigdb hallmark be computed

msigdb_subcat <- c('CP:BIOCARTA', 'CP:KEGG','CP:KEGG_MEDICUS', 'GO:BP')
# subcategories ('gs_subcat') of msigdb database for GSEA calculation

# make dirs and set additional vars ---------------------------------------

dir.create(file.path(output_dir, 'pathway_analysis'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'pathway_analysis', 'gsea'), showWarnings = T, recursive = T)

norm_is_log <- ifelse(norm_type %in% c('exprs', 'q3_norm', 'deseq2_norm'), FALSE, TRUE)
norm_name <- ifelse(norm_is_log, norm_type, paste0("log_", norm_type)) #TODO is it needed?


# path to deconvolution mtx
deconv_bp_path <- file.path(output_dir,'deconvolution', 'bayes_prism', 
                            paste0('bp_res_', scrna_anno,'_expr_mtx_cleaned_', 
                                   deconv_norm_type, '_', deconv_batch_rm_type, '_corr.RDS'))

#TODO move somewhere else - before normalisation? 
# path to cleaned scrna which should be calculated in deconvolution step - for low complex gene rmv 
scrna_ref_cleaned_path <- file.path(output_dir, 'deconvolution', gsub('.RDS', '_cleaned_for_deconv.RDS', basename(scrna_ref_path)))
raw_counts_layer <- 'counts'

# load geomx obj from rds -------------------------------------------------

geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)

if(low_complex_rmv){
  # removing low complexity genes as it was done before bp deconvolution
  scrna_ref_obj <- readRDS(scrna_ref_cleaned_path)
  geomx_obj <-   remove_low_complex_and_noncoding_genes(geomx_obj, scrna_ref_obj, raw_counts_layer = 'counts')
} 

expr_list <- list()

if('all' %in% pathway_inp_data_type){
  # make log expression mtx if needed
  if(!norm_is_log){
    # make log2 transformed normalised counts if norm_type not in log scale
    expr_mtx <- log2(geomx_obj@assayData[[norm_type]] + 1)
  } else{
    expr_mtx <- geomx_obj@assayData[[norm_type]]
  }
  
  expr_list[[length(expr_list) + 1]] <- expr_mtx
  names(expr_list) <- 'all'
}

# load deconvoluted signal ------------------------------------------------

if('bp' %in% pathway_inp_data_type){
  deconv_ct_list <- readRDS(deconv_bp_path)
  
  # filter to cell types of interest
  if(!is.null(ct_of_interest)){
    deconv_ct_list <- deconv_ct_list[ct_of_interest]
  }
  
  names(deconv_ct_list) <- paste0('deconv_', names(deconv_ct_list))
  
  expr_list <- c(expr_list, deconv_ct_list)
}

# prepare signatures list -------------------------------------------------

if(signature_type == 'msigdb'){
  # signatures from all Hallmark + selected CP from msigDB 
  sign_list <- prepare_msigdb_sign_list(adjust_synonym = adj_synonym, geomx_obj = geomx_obj, hal = compute_hallmark, 
                                               db_subcat_list = msigdb_subcat)
  out_name <- 'msigdb'
} else if(signature_type == 'custom'){
  # signatures from custom file
  sign_list <- prepare_custom_sign_list(fread(custom_sign_path), adjust_synonym = adj_synonym,
                                               geomx_obj = geomx_obj)
  out_name <- paste0('custom_', gsub('//.csv', '', basename(custom_sign_path)))
} else{
  stop("signature_type parameter can only be 'msigb' or 'custom'")
}

sign_list <- sign_list[sapply(sign_list, length) >= min_sign_gene_nr]

# calculate gsea ----------------------------------------------------------

gsva_list_long <- lapply(1:length(expr_list), function(x){

  if(gsea_type == 'gsva'){
    # do gsva
    gsea <- gsva(gsvaParam(expr_list[[x]], sign_list, kcdf="Gaussian", minSize = 5))
  } else if(gsea_type == 'ssgsea'){
    # do ssgsea
    gsea <- gsva(ssgseaParam(expr_list[[x]], sign_list, minSize = 5, normalize = T))
  } else{
    stop("gsea_type parameter can only be 'gsva' or 'ssgsea'")
  }
  
  # adjust df and save
  gsea_long <- melt(gsea)
  colnames(gsea_long) <- c('pathway','dcc_filename', paste0(gsea_type, '_score'))
  gsea_long$expr_signal <- names(expr_list)[x]
  gsea_long <- left_join(gsea_long, sData(geomx_obj)[meta_names])
  
  #TODO better names - deconv and all may have diff norm types!
  fwrite(gsea_long, file.path(output_dir,'pathway_analysis', 'gsea', 
                              paste0(gsea_type, '_norm_', norm_name, '_',
                                     names(expr_list)[x], '_', out_name,  '.csv')))
  
  return(gsea_long)
})

writeLines(c('GSEA logs:',
             'GSEA type: ', gsea_type, 
             '; normalisation type : ', norm_name,
             '; signature type : ', signature_type,
             '; low complex gene removed : ', low_complex_rmv,
             '; synonym genes adjusted : ', adj_synonym,
             '; deconv mtx used : ', deconv_bp_path), gsea_logs_path)


# calculate limma rotation gene set test ----------------------------------

#TODO
#input matrix in log
# limma::fry(v, index = x, design = design, contrast = contr.matrix[,1], robust = TRUE)


# PROGENy scores ----------------------------------------------------------

#TODO progeny needs an updated Matrix package, while >1.7 does not work for DGE
#TODO use pathway significance info from prog_perm[[2]] (nulldist)

progeny_type <- NULL 
# whether 'perm' or 'nonperm' - some quirks in progeny algorithm, results similar but perm is preferred
# if NULL <- no progeny calculation

if(!is.null(progeny_type)){
  
  dir.create(file.path(output_dir, 'progeny'), showWarnings = T, recursive = T)
  
  if(progeny_type == 'perm'){
    prog_res <- progeny(
      expr_list$all,
      organism = "Human",
      top = 100,
      perm = 10,
      z_scores = FALSE,
      get_nulldist = TRUE
    )
    
  } else if(progeny_type == 'noperm'){
    prog_res <- progeny(
      geomx_obj@assayData[[norm_type]],
      scale = TRUE,
      organism = "Human",
      top = 100,
    )
  } else{
    stop("progeny_type can be either 'perm' or 'noperm'")
  }
  
  # adjust the table and save
  prog_res[[1]] <- t(prog_res[[1]])
  rownames(prog_res[[1]]) <- gsub('\\.', '\\-', rownames(prog_res[[1]]))
  rownames(prog_res[[1]]) <- gsub('\\-dcc', '\\.dcc', rownames(prog_res[[1]]))
  
  prog_res <- melt(prog_res)
  colnames(prog_res) <- c('dcc_filename', 'progeny_path', 'progeny_score')
  prog_res <- left_join(prog_res, sData(geomx_obj)[meta_names])
  
  fwrite(prog_res, file.path(output_dir, 'progeny', paste0('progeny_', progeny_type, '.csv')))
}

