
# get variables -----------------------------------------------------------

# TODO put it somewhere in the main script and reuse through scripts
imp_vars <- c("Segment", "Annotation_cell", "NACT_status", "PFS", "PFS_months", "Sample") 
gsva_vars <- c(imp_vars, 'dcc_filename', 'Patient') 

# best to use batch effect corrected or at least vst data in log form 
norm_type <- 'limma_batch_corr' # or harmony_batch_corr or deseq2_vst_scaled
norm_is_log <- TRUE # limma and harmony batch eff corr are in log scale, vst is similar to log

adj_synonym <- T # whether or not adjust synonyms genes
# around 300 genes can be rescued this way but ensembl does not always work
# if there are issues, turn it off
min_sign_gene_nr <- 5 # signatures with less nr of genes will be removed, 5 is min in msigdb

progeny_type <- NULL 
# whether 'perm' or 'nonperm' - some quirks in progeny algorithm, results similar but perm is preferred
# if NULL <- no progeny calculation

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'pathway_analysis'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'pathway_analysis', 'gsea'), showWarnings = T, recursive = T)

# load geomx obj from rds -------------------------------------------------

geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)

# read expression mtx
expr_mtx <- assayDataElement(geomx_obj, elt = norm_type)

# make log expression mtx if needed
if(!norm_is_log){
  # convert normalized counts to log scale
  assayDataElement(object = geomx_obj, elt = paste0("log_", norm_type)) <-
    assayDataApply(geomx_obj, 2, FUN = log, base = 2, elt = norm_type)
  norm_type <- paste0("log_", norm_type)
}

# load deconvoluted signal ------------------------------------------------

# TODO load from 
# saveRDS(deconv_ct_list, file = file.path(output_dir,'deconvolution', 'bayes_prism', 
# paste0('bp_res_', scrna_anno, '_', ct_nr_thr, '_expr_mtx_cleaned_norm.RDS')))


# prepare signatures list -------------------------------------------------

if(signature_type == 'msigdb'){
  # signatures from all Hallmark + selected CP from msigDB 
  sign_list <- prepare_msigdb_sign_list(adjust_synonym = adj_synonym, geomx_obj = geomx_obj, hal = T, 
                                               db_subcat_list = c('CP:BIOCARTA', 'CP:KEGG','GO:BP'))
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

#TODO for deconvolution
# expr_list <- deconv_ct_list
# expr_list[[length(expr_list) + 1]] <- expr_mtx
# names(expr_list) <- c(paste0('deconv_', ct_names, '_', deconv_type), 'all')

expr_list <- list(geomx_obj@assayData[[norm_type]])
names(expr_list) <- 'all'


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
  gsea_long <- left_join(gsea_long, pData(geomx_obj)[gsva_vars])
  
  fwrite(gsea_long, file.path(output_dir,'pathway_analysis', 'gsea', 
                              paste0(gsea_type, '_norm_', norm_type, '_',
                                     names(expr_list)[x], '_', out_name,  '.csv')))
  
  return(gsea_long)
})


# calculate limma rotation gene set test ----------------------------------

#TODO
#input matrix in log
# limma::fry(v, index = x, design = design, contrast = contr.matrix[,1], robust = TRUE)


# PROGENy scores ----------------------------------------------------------

#TODO progeny needs an updated Matrix package, while >1.7 does not work for DGE
#TODO use pathway significance info from prog_perm[[2]] (nulldist)

if(!is.null(progeny_type)){
  
  dir.create(file.path(output_dir, 'progeny'), showWarnings = T, recursive = T)
  
  if(progeny_type == 'perm'){
    prog_res <- progeny(
      geomx_obj@assayData[[norm_type]],
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
  prog_res <- left_join(prog_res, pData(geomx_obj)[gsva_vars])
  
  fwrite(prog_res, file.path(output_dir, 'progeny', paste0('progeny_', progeny_type, '.csv')))
}

