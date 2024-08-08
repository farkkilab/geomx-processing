# TODO check if all packages are needed
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(GSVA)
library(plyr)
library(dplyr)
library(data.table)
library(biomaRt)
library(DESeq2)
library(msigdbr)
library(tibble)
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
data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'
output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2'

input_rds_path <- file.path(output_dir, 'geomx_qc_norm.RDS')

input_bp_deconv_path <- file.path(output_dir, 'deconvolution', 'bp', 'bp_res_mid_lvl_ct_45.RDS')
deconv_type <- 'mid_lvl_ct' # either mid_lvl_ct or cell_type

input_sd_deconv_path <- file.path(output_dir, 'deconvolution', 'sd', 'sd_res_mid_lvl_ct_nofilt.rds')

sig_additional_path <- '/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/additional_signatures_macro_tcells_msigdb_filt.csv'

# sig_path_macro <- '/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/additional_signatures_macro.csv'
# sig_path_tcell <- '/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/additional_signatures_tcells.csv'
#sig_name <- 'additional_macro'

# if not doing all hall_cp
selected_sig_path <- file.path('/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/immune_signatures_selected_names.csv')

######
imp_vars <- c("Segment", "Annotation_cell", "NACT status", "PFS") # vals used for sankey, detection rate plots, 
gsva_vars <- c(imp_vars, 'dcc_filename', 'Patient') #TODO add 'Sample

norm_type <- 'q3_norm' # either 'q3_norm' or 'quant_norm'

do_gsva_hal_cp_all <- FALSE # whether do gsva on all hallmark and cp paths from msigdb

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'gsva'), showWarnings = T, recursive = T)
#dir.create(file.path(output_dir, 'progeny'), showWarnings = T, recursive = T)

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')

# load geomx obj from rds -------------------------------------------------

geomx_obj <- readRDS(input_rds_path)

# read expression mtx
expr_mtx <- assayDataElement(geomx_obj, elt = norm_type)

# load deconvoluted signal for macrophages and tcells ---------------------
# TODO move theta cv and normalisation to deconvolution script

if(deconv_type == 'mid_lvl_ct'){
  ct_names <- ct_names <- c("Tcells", "Bcells", "Fibroblasts", "NKcells", "Macrophages", "DCs", "tumor", "Endothelial cells")
} else if(deconv_type == 'cell_type'){

  ct_names <- c("Fibroblasts","Macrophages", "tumor", "Endothelial cells", "Classical monocytes",
                "Tem/Trm cytotoxic Tcells", "CD16+ NK cells", "pDC", "Plasma cells", "DC1", "Regulatory Tcells",
                "Naive B cells", "Migratory DCs", "NK cells", "Type 17 helper T cells", "Memory B cells",
                "Tcm/Naive helper Tcells", "CD16- NK cells")
}

deconv_res <- readRDS(input_bp_deconv_path)

# extract coeff of variation per cell type
# mask ct_frac results if cv > 0.2-0.5 (0.1 thr for bulk, 0.5 for Visium, GeoMx should be in the middle)
# histogram suggests 0.2 as thr
# TODO move it to deconvolution script
deconv_ct_list <- lapply(ct_names, function(ct_name){
  cell_frac_cv <- as.data.frame(deconv_res@posterior.theta_f@theta.cv)
  cell_to_rm <- rownames(cell_frac_cv)[cell_frac_cv[[ct_name]] > 0.2]
  
  deconv_ct <- BayesPrism::get.exp(bp=deconv_res,
                                      state.or.type="type",
                                      cell.name=ct_name)
  
  deconv_ct <- deconv_ct[!(rownames(deconv_ct) %in% cell_to_rm), ]
  deconv_ct <- varianceStabilizingTransformation(round(t(deconv_ct))) # normalisation
  
  return(deconv_ct)
})

# GSVA on all Hallmark + CP + additional -------------------------------------
# do GSVA on all Hallmark + CP from msigDB or on selected Hal + CP + additional pathways

# prepare msigdb signatures list
msigdb_df <- msigdbr(species = "Homo sapiens")
msigdb_df <- filter(msigdb_df, gs_cat == 'H' | 
                      gs_subcat %in% c('CP:BIOCARTA', 'CP:KEGG', 'CP:REACTOME', 'CP:PID', 'CP:WIKIPATHWAYS', 'GO:BP'))

msigdb_df$gene_symbol_adj <- adjust_synonym_genes(rownames(geomx_obj), msigdb_df$gene_symbol)

hal_cp_list <- lapply(unique(msigdb_df$gs_name), function(x){
  gs <- filter(msigdb_df, gs_name == x)
  gs_genes <- unique(gs$gene_symbol_adj)
})

names(hal_cp_list) <- unique(msigdb_df$gs_name)

# filter list to selected pathways
if(!do_gsva_hal_cp_all){
  selected_sig <- fread(selected_sig_path)
  hal_cp_list <- hal_cp_list[names(hal_cp_list) %in% selected_sig$pathway]
  out_name <- 'selected_and_additional'
} else{
  out_name <- 'hal_cp_full_and_additional'
}

# prepare additional signatures list
sig_list_additional <- as.list(fread(sig_additional_path))
sig_list_additional <- lapply(sig_list_additional, function(l){l[l !=""]})
sig_list_additional <- lapply(sig_list_additional, function(x){
  adjust_synonym_genes(rownames(geomx_obj), x)})

sig_list_all <- c(hal_cp_list, sig_list_additional)



# do ssgsea
# ssgsea_hal_cp <- gsva(expr_mtx, hal_cp_list, method = 'ssgsea', kcdf="Gaussian", min.sz = 5)

expr_list <- deconv_ct_list
expr_list[[length(expr_list) + 1]] <- expr_mtx
names(expr_list) <- c(paste0('deconv_', ct_names, '_', deconv_type), 'all')

gsva_list_long <- lapply(1:length(expr_list), function(x){
  # do gsva
  gsva <- gsva(gsvaParam(expr_list[[x]], sig_list_all, kcdf="Gaussian", minSize = 5))
  
  # adjust df and save
  gsva_long <- melt(gsva)
  colnames(gsva_long) <- c('pathway','dcc_filename', 'gsva_score')
  gsva_long$expr_signal <- names(expr_list)[x]
  gsva_long <- left_join(gsva_long, pData(geomx_obj)[gsva_vars])
  
  fwrite(gsva_long, file.path(output_dir, 'gsva', paste0('gsva_', names(expr_list)[x], '_', out_name,  '.csv')))
  
  return(gsva_long)
})

###############################################################################
###############################################################################
# adjusting gsva scores from full signal for spatialdecon cell freq -------
dir.create(file.path(output_dir, 'gsva', 'sd_lm_additional_msigdb_filt'), showWarnings = T, recursive = T)

sd_deconv <- readRDS(input_sd_deconv_path)
sd_deconv <- data.frame(pData(sd_deconv)[, 'prop_of_all'])
sd_deconv <- dplyr::select(sd_deconv, -Mast.cells, -other)
colnames(sd_deconv) <- paste0('deconv_', colnames(sd_deconv))
deconv_names <- colnames(sd_deconv)
sd_deconv <- tibble::rownames_to_column(sd_deconv, 'dcc_filename')

#gsva_all_long <- fread(file.path(output_dir, 'gsva', 'gsva_all_selected.csv'))
gsva_all_long <- gsva_list_long[[1]]
gsva_all_long <- left_join(gsva_all_long, sd_deconv)

selected_sig <- fread(selected_sig_path)
gsva_df <- filter(gsva_all_long, pathway %in% selected_sig$pathway)

##############
# calculate lm for each pathway vs deconvolution cell type
gsva_lm <- lapply(unique(as.vector(gsva_df$pathway)), function(path_name){
  gsva_path <- gsva_df[gsva_df$pathway == path_name, ]
  
  lapply(deconv_names, function(ct){
    
    # fit lm with ct fraction as explanatory var
    lm_res <- lm(gsva_score~get(ct),data=gsva_path)
    lm_coef <- summary(lm_res)$coefficients[2]
    lm_rsq <- summary(lm_res)$adj.r.squared
    
    gsva_lm_res <- list('pathway' = path_name, 'deconv_ct' = ct,
                              lm_coef = lm_coef, lm_rsq = lm_rsq)
    
    if(!do_gsva_hal_cp_all){
      png(file = file.path(output_dir, 'gsva', 'sd_lm_additional_msigdb_filt', paste0('scatter_', path_name, '_', ct, '.png')))
      plot(gsva_path[[ct]], gsva_path$gsva_score, xlab = path_name, ylab = ct)
      abline(lm(gsva_score~get(ct),data=gsva_path),col='red')
      dev.off()
    }
    return(gsva_lm_res)
  })
})

gsva_lm <- unlist(gsva_lm, recursive = F)
gsva_lm_df <- rbindlist(gsva_lm, fill=TRUE)

fwrite(gsva_lm_df, file.path(output_dir, 'gsva', paste0('sd_lm_gsva_', out_name, '.csv')))


##########################################
# correct gsva scores based on lm parameters
# adjusted for lin reg
# y = a + xb
# y = a // -xb

rsq_thr <- 0.4

gsva_lm_adj_tcell <- lapply(unique(gsva_df$pathway), function(path_name){
  gsva_df_path <- gsva_df[gsva_df$pathway == path_name, ]

  lm_tcell <- gsva_lm_df[gsva_lm_df$pathway == path_name & grepl(cd8_ct, gsva_lm_df$deconv_ct), ]

  if(lm_tcell$lm_rsq > rsq_thr){
    gsva_df_path$gsva_score <- gsva_df_path$gsva_score - (gsva_df_path[[paste0('deconv_', cd8_ct)]] * lm_tcell$lm_coef)
  }
  
  return(gsva_df_path)
})

gsva_lm_adj_macro <- lapply(unique(gsva_df$pathway), function(path_name){
  gsva_df_path <- gsva_df[gsva_df$pathway == path_name, ]
  
  lm_macro <- gsva_lm_df[gsva_lm_df$pathway == path_name & grepl(macro_ct, gsva_lm_df$deconv_ct), ]
  
  if(lm_macro$lm_rsq > rsq_thr){
    gsva_df_path$gsva_score <- gsva_df_path$gsva_score - (gsva_df_path[[paste0('deconv_', macro_ct)]] * lm_macro$lm_coef)
  }
  return(gsva_df_path)
})

gsva_lm_adj_tcell <- do.call(rbind, gsva_lm_adj_tcell)
gsva_lm_adj_macro <- do.call(rbind, gsva_lm_adj_macro)

fwrite(gsva_lm_adj_tcell, file.path(output_dir, 'gsva', 
                                 paste0('gsva_all_', out_name, '_sd_lm_adjusted_tcell_additional_msigdb_filt_', as.character(rsq_thr), '.csv')))

fwrite(gsva_lm_adj_macro, file.path(output_dir, 'gsva', 
                                    paste0('gsva_all_', out_name, '_sd_lm_adjusted_macro_additional_msigdb_filt_', as.character(rsq_thr), '.csv')))



# VISUALISATION IN  GSVA_VISUALISATION.R  

# PROGENy scores ----------------------------------------------------------

# prog_noperm <- progeny(
#   geomx_obj@assayData[[norm_type]],
#   scale = TRUE,
#   organism = "Human",
#   top = 100,
# )
# 
# prog_perm <- progeny(
#   geomx_obj@assayData[[norm_type]],
#   organism = "Human",
#   top = 100,
#   perm = 10,
#   z_scores = FALSE,
#   get_nulldist = TRUE
# )
# 
# prog_perm[[1]] <- t(prog_perm[[1]])
# rownames(prog_perm[[1]]) <- gsub('\\.', '\\-', rownames(prog_perm[[1]]))
# rownames(prog_perm[[1]]) <- gsub('\\-dcc', '\\.dcc', rownames(prog_perm[[1]]))
# 
# #TODO use pathway significance infor from prog_perm[[2]] (nulldist)
# 
# # adjust and save dfs
# prog_list <- list(noperm = prog_noperm, perm = prog_perm[[1]])
# 
# sapply(1:length(prog_list), function(x){
#   prog_df <- prog_list[[x]]
#   prog_name <- names(prog_list)[x]
#   
#   # adjust df
#   prog_long <- melt(prog_df)
#   colnames(prog_long) <- c('dcc_filename', 'progeny_path', 'progeny_score')
#   prog_long <- left_join(prog_long, pData(geomx_obj)[gsva_vars])
#   
#   fwrite(prog_long, file.path(output_dir, 'progeny', paste0('progeny_', prog_name, '.csv')))
# })
# 

 # UCell scores for selected pathways --------------------------------------
# 
# 
# # changing signatures to HGNC symbols
# 
# # sig_list <- lapply(sig_all, function(x){
# #   sig <- fread(x)
# #   sig_list <- as.list(sig)
# #   return(sig_list)
# # })
# # 
# # sig_list <- unlist(sig_list, recursive = F)
# # 
# # sig_list_names <- lapply(sig_list, function(x){
# #   print('XXXXX')
# #   print(x)
# #   x <- x[!is.na(x)]
# #   print(x)
# #   x <- unlist(gene_2names(x, conv='entrez'))
# #   print(x)
# #   })
# # 
# # library(qpcR)
# # sig_df <- do.call(qpcR:::cbind.na, sig_list_names)
# # fwrite(sig_df, '/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/texh_macro_mhc.csv')
# 
# 
