# TODO check if all packages are needed
# library(NanoStringNCTools)
# library(GeomxTools)
# library(GeoMxWorkflows)
# library(GSVA)
# library(plyr)
# library(dplyr)
# library(data.table)
# library(biomaRt)
# library(DESeq2)
# library(msigdbr)
# library(tibble)
#library(GeoDiff)
# library(ggplot2)
# library(ggforce)
# 
# library(cowplot)
# library(preprocessCore)
# library(Biobase)
# library(reshape2)
# 
# library(clusterProfiler)
# library(progeny)
# library(reshape2)
# library(ggpubr)


# get variables -----------------------------------------------------------
# data_dir <- '/home/iganiemi/Documents/phd/st/data/geomx/geomx_batch1_nact/'
# output_dir <- '/home/iganiemi/Documents/phd/st/geomx-processing/results/nact'
# 
# input_rds_path <- file.path(output_dir, 'geomx_qc_norm.RDS')
# 
# input_bp_deconv_path <- file.path(output_dir, 'deconvolution', 'bp', 'bp_res_mid_lvl_ct_45.RDS')
# deconv_type <- 'mid_lvl_ct' # either mid_lvl_ct or cell_type
# 
# input_sd_deconv_path <- file.path(output_dir, 'deconvolution', 'sd', 'sd_res_mid_lvl_ct_nofilt.rds')
# 
# sig_additional_path <- '/home/iganiemi/Documents/phd/st/geomx-processing/data/signatures/stromal_cell_subtype_signatures_symbols_ensembl_ids_revised.csv'
#sig_additional_path <- '/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/additional_signatures_macro_tcells_msigdb_filt.csv'
# sig_path_macro <- '/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/additional_signatures_macro.csv'
# sig_path_tcell <- '/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/additional_signatures_tcells.csv'
#sig_name <- 'additional_macro'

# if not doing all hall_cp
# selected_sig_path <- file.path('/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/immune_signatures_selected_names.csv')

######
# TODO put it somewhere in the main script and reuse through scripts
imp_vars <- c("Segment", "Annotation_cell", "NACT status", "PFS", "PFS_months", "Sample") 
gsva_vars <- c(imp_vars, 'dcc_filename', 'Patient') 

#TODO check if this is the best norm type or can be used with batch correction
norm_type <- 'q3_norm' # deseq2_norm / vst / limma_batch_corr

adj_synonym <- T # whether or not adjust synonyms genes
# around 300 genes can be rescued this way but ensembl does not always work
# if there are issues, turn it off

progeny_type <- 'perm' 
# whether 'perm' or 'nonperm' - some quirks in progeny algorithm, results similar but perm is preferred
# if NULL <- no progeny calculation

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'pathway_analysis'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'gsea'), showWarnings = T, recursive = T)

# load geomx obj from rds -------------------------------------------------

geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)

# read expression mtx
expr_mtx <- assayDataElement(geomx_obj, elt = norm_type)


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

