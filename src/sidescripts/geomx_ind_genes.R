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
library(BayesPrism)

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
goi_other <- c('CD226', 'CD96')

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

expr_mtx_all <- expr_mtx_all[intersect(rownames(expr_mtx_all), c(goi_macro, goi_tcell, goi_other)), ]

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



# ind gene expr vs gsva and cell fraq -------------------------------------
# NECTIN2 vs TIGIT, CD96, CD226
# NECTIN2 vs CD163/ITGAX DEGS signatures
# NECTIN2 vs tumor and macro fractions
dir.create(file.path(output_dir, 'ind-genes', 'NECTIN2'), showWarnings = T, recursive = T)

gsva_paths <- c("DEG_CD163_VS_ITGAX_NEG_05FC", "DEG_CD163_VS_ITGAX_NEG_ALL", 
                "DEG_CD163_VS_ITGAX_POS_05FC", "DEG_CD163_VS_ITGAX_POS_ALL")
deconv_cells <- c('Tcells', 'tumor', 'Macrophages')

exp_mtx_dt <- as.data.frame(t(expr_mtx_all)) %>%
  rownames_to_column('dcc_filename')

gsva_df <- fread(file.path(output_dir, 'gsva', 'gsva_all_selected.csv'))
gsva_df <- dcast(gsva_df, dcc_filename + Segment + `NACT status` + PFS ~ pathway, value.var = 'gsva_score')
gsva_df <- gsva_df[, c('dcc_filename', 'Segment', 'NACT status', 'PFS', gsva_paths)]

deconv_frac <- BayesPrism::get.fraction(bp=deconv_res,
                         which.theta="final",
                         state.or.type="type")

deconv_frac <- rownames_to_column(as.data.frame(deconv_frac), 'dcc_filename')
deconv_frac <- deconv_frac[, c('dcc_filename', deconv_cells)]

dt_all <- left_join(gsva_df, exp_mtx_dt) %>%
  left_join(deconv_frac)

#####################
coln1 <- 'NECTIN2'


for(coln2 in c('TIGIT', 'CD96', 'CD226', gsva_paths, deconv_cells)){
  lm_rsq <- summary(lm(get(coln1)~get(coln2), data=dt_all))$adj.r.squared
  cor_pe <- cor.test(dt_all[[coln1]], dt_all[[coln2]], method = 'pearson')$estimate
  
  png(file = file.path(output_dir, 'ind-genes', 'NECTIN2', paste0('scatter_', coln2, '.png')))
  plot(dt_all[[coln1]], dt_all[[coln2]], xlab = coln1, ylab = coln2, sub = paste0('R2 = ', round(lm_rsq, 2), 
                                                                                  ' cor = ', round(cor_pe, 2)))
  abline(lm(get(coln1)~get(coln2), data=dt_all),col='red')
  dev.off()
}




