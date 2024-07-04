library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(DESeq2)
library(plyr)
library(dplyr)
library(data.table)
library(biomaRt)
library(tibble)
library(ggplot2)
library(ggpubr)

# get variables -----------------------------------------------------------
data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'
output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2'

input_rds_path <- file.path(output_dir, 'geomx_qc_norm.RDS')
input_deconv_path <- file.path(output_dir, 'deconvolution', 'bp', 'bp_res_mid_lvl_ct_45.RDS')

macro_ct <- 'Macrophages'
cd8_ct <- 'Tcells'

######
imp_vars <- c("Segment", "Annotation_cell", "NACT status", "PFS") # vals used for sankey, detection rate plots, 
gsva_vars <- c(imp_vars, 'dcc_filename', 'Patient') #TODO add 'Sample

norm_type <- 'q3_norm' # either 'q3_norm' or 'quant_norm'

# goi <- c('NECTIN2', 'TIGIT', 'CD96', 'LAG3', 'CD226', 'CXCL12', 'CXCR4', 'HGF',
#          'CD44', 'PDCD1', 'HAVCR2', 'CXCL9', 'CXCL10', 'CXCR3')
# HAVCR2 = TIM3

goi_macro <- c('NECTIN2', 'CD163')
goi_tcell <- c('TIGIT', 'LAG3', 'HAVCR2', 'PDCD1', 'CTLA4')

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'ind-genes'), showWarnings = T, recursive = T)
#dir.create(file.path(output_dir, 'progeny'), showWarnings = T, recursive = T)

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')

# load geomx obj from rds -------------------------------------------------

geomx_obj <- readRDS(input_rds_path)

# read expression mtx
expr_mtx_all <- assayDataElement(geomx_obj, elt = norm_type)


goi_macro <- adjust_synonym_genes(rownames(expr_mtx_all), goi_macro)
goi_tcell <- adjust_synonym_genes(rownames(expr_mtx_all), goi_tcell)

expr_mtx_all <- expr_mtx_all[intersect(rownames(expr_mtx_all), c(goi_macro, goi_tcell)), ]

# read deconv results
deconv_res <- readRDS(input_deconv_path)
deconv_tcell <- BayesPrism::get.exp(bp=deconv_res,
                                    state.or.type="type",
                                    cell.name=cd8_ct)

deconv_tcell <- varianceStabilizingTransformation(round(t(deconv_tcell))) # normalisation
deconv_tcell <- deconv_tcell[intersect(rownames(deconv_tcell), c(goi_macro, goi_tcell)), ]


deconv_macro <- BayesPrism::get.exp(bp=deconv_res,
                                    state.or.type="type",
                                    cell.name=macro_ct)

deconv_macro <- varianceStabilizingTransformation(round(t(deconv_macro))) # normalisation
deconv_macro <- deconv_macro[intersect(rownames(deconv_macro), c(goi_macro, goi_tcell)), ]

# individual genes expression ---------------------------------------------

exp_mtx_list <- list('all' = expr_mtx_all, 'deconv_macro' = deconv_macro, 'deconv_tcell' = deconv_tcell)

sapply(1:length(exp_mtx_list), function(x){
  exp_mtx <- exp_mtx_list[[x]]
  exp_name <- names(exp_mtx_list)[x]
  
  goi_long <- melt(exp_mtx)
  colnames(goi_long) <- c('gene','dcc_filename', 'expr')
  goi_long <- left_join(goi_long, pData(geomx_obj)[c('dcc_filename', 'Segment', 'Annotation_cell', 'NACT status', 'PFS', 'Patient')])
  fwrite(goi_long, file.path(output_dir, 'ind-genes', paste0('ind_genes_expr_', exp_name, '.csv')))
  
  # boxplot visualisation ---------------------------------------------------
  
  pathway_boxplot(goi_long, 'gene', 'expr', 'Annotation_cell', c('Segment'), 'individual gene expression ',
                  file.path(output_dir, 'ind-genes', paste0('box_ind_', exp_name, '_anno.pdf')))
  
  pathway_boxplot(goi_long, 'gene', 'expr', 'NACT status', c('Segment'), 'individual gene expression ',
                  file.path(output_dir, 'ind-genes', paste0('box_ind_', exp_name, '_nact_all.pdf')))
  
  pathway_boxplot(goi_long, 'gene', 'expr', 'NACT status', c('Annotation_cell', 'Segment'), 'individual gene expression ',
                  file.path(output_dir, 'ind-genes', paste0('box_ind_', exp_name, '_nact_peranno.pdf')))
  
  pathway_boxplot(goi_long, 'gene', 'expr', 'PFS', c('Segment'), 'individual gene expression ',
                  file.path(output_dir, 'ind-genes', paste0('box_ind_', exp_name, '_pfs_all.pdf')))
  
  pathway_boxplot(goi_long, 'gene', 'expr', 'PFS', c('Annotation_cell', 'Segment'), 'individual gene expression ',
                  file.path(output_dir, 'ind-genes', paste0('box_ind_', exp_name, '_pfs_peranno.pdf')))
  
  gc()
})








