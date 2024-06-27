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
library(ggplot2)
library(ggforce)

library(cowplot)
library(preprocessCore)
library(Biobase)
library(reshape2)

library(umap)
library(Rtsne)

library(clusterProfiler)
library(progeny)
library(reshape2)
library(ggpubr)


# get variables -----------------------------------------------------------
data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'
output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2'

input_rds_path <- file.path(output_dir, 'geomx_qc_norm.RDS')

input_bp_deconv_path <- file.path(output_dir, 'deconvolution', 'bp', 'bp_res_mid_lvl_ct_45.RDS')
deconv_type <- 'mid_lvl_ct' # either mid_lvl_ct or cell_type

input_sd_deconv_path <- file.path(output_dir, 'deconvolution', 'sd', 'sd_res_mid_lvl_ct_nofilt.rds')

sig_path_macro <- '/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/additional_signatures_macro.csv'
sig_path_tcell <- '/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/additional_signatures_tcells.csv'
#sig_name <- 'additional_macro'

hal_cp_selected_path <- file.path('/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/hal_cp_immune_pathways_selected.csv')

######
imp_vars <- c("Segment", "Annotation_cell", "NACT status", "PFS") # vals used for sankey, detection rate plots, 
gsva_vars <- c(imp_vars, 'dcc_filename', 'Patient') #TODO add 'Sample

norm_type <- 'q3_norm' # either 'q3_norm' or 'quant_norm'

do_gsva_hal_cp <- FALSE # whether do gsva on all hallmark and cp paths from msigdb

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'gsva'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'progeny'), showWarnings = T, recursive = T)

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')

# load geomx obj from rds -------------------------------------------------

geomx_obj <- readRDS(input_rds_path)


# load deconvoluted signal for macrophages and tcells ---------------------
# TODO remove redundancy tcell macro

if(deconv_type == 'mid_lvl_ct'){
  macro_ct <- 'Macrophages'
  cd8_ct <- 'Tcells'
} else if(deconv_type == 'cell_type'){
  macro_ct <- 'Macrophages' #TODO and 'Classical monocytes' ??
  cd8_ct <- 'Tem/Trm cytotoxic T cells'
}

deconv_res <- readRDS(input_bp_deconv_path)

# extract coeff of variation per cell type
# mask ct_frac results if cv > 0.2-0.5 (0.1 thr for bulk, 0.5 for Visium, GeoMx should be in the middle)
# histogram suggests 0.2 as thr
cell_frac_cv <- as.data.frame(deconv_res@posterior.theta_f@theta.cv)

hist(cell_frac_cv[[cd8_ct]], breaks = 1000)
hist(cell_frac_cv[[macro_ct]], breaks = 1000)

tcell_to_rm <- rownames(cell_frac_cv)[cell_frac_cv[[cd8_ct]] > 0.2]
macro_to_rm <- rownames(cell_frac_cv)[cell_frac_cv[[macro_ct]] > 0.2]

deconv_tcell <- get.exp (bp=deconv_res,
                       state.or.type="type",
                       cell.name=cd8_ct)

deconv_tcell <- deconv_tcell[!(rownames(deconv_tcell) %in% tcell_to_rm), ]
deconv_tcell <- varianceStabilizingTransformation(round(t(deconv_tcell))) # normalisation


deconv_macro <- get.exp (bp=deconv_res,
                       state.or.type="type",
                       cell.name=macro_ct)

deconv_macro <- deconv_macro[!(rownames(deconv_macro) %in% macro_to_rm), ]
deconv_macro <- varianceStabilizingTransformation(round(t(deconv_macro))) # normalisation

# GSVA and ssGSEA on macro + tcells pathways ------------------------------

# TODO Correlating Z (after normalization using vst or from bp.res@reference.update@psi_mal) 
# with theta to understand how gene expression of each gene (in malignant cells) 
#correlates with the cell type fraction of non-malignant cells in tumor micro-environment, 
# followed by gene set enrichment analysis (as done in BayesPrism paper).

# read selected pathways
sig_list_macro <- as.list(fread(sig_path_macro))
sig_list_macro <- lapply(sig_list_macro, function(l){l[l !=""]})
sig_list_macro <- lapply(sig_list_macro, function(x){
  adjust_synonym_genes(rownames(geomx_obj), x)})

sig_list_tcell <- as.list(fread(sig_path_tcell))
sig_list_tcell <- lapply(sig_list_tcell, function(l){l[l !=""]})
sig_list_tcell <- lapply(sig_list_tcell, function(x){
  adjust_synonym_genes(rownames(geomx_obj), x)})

# do gsva on pathway lists
gsva_sig_macro <- gsva(gsvaParam(deconv_macro, sig_list_macro, kcdf="Gaussian", minSize = 5))
gsva_sig_tcell <- gsva(gsvaParam(deconv_tcell, sig_list_tcell, kcdf="Gaussian", minSize = 5))

# do ssgsea on pathway lists
# ssgsea_sel_sig <- gsva(expr_mtx, sig_list, method = 'ssgsea', kcdf="Poisson", min.sz = 5)

# adjust df
gsva_sig_macro_long <- melt(gsva_sig_macro)
colnames(gsva_sig_macro_long) <- c('pathway','dcc_filename', 'gsva_score')
gsva_sig_macro_long <- left_join(gsva_sig_macro_long, pData(geomx_obj)[gsva_vars])
fwrite(gsva_sig_macro_long, file.path(output_dir, 'gsva', paste0('gsva_deconv_macro_additional_', deconv_type, '.csv')))

gsva_sig_tcell_long <- melt(gsva_sig_tcell)
colnames(gsva_sig_tcell_long) <- c('pathway','dcc_filename', 'gsva_score')
gsva_sig_tcell_long <- left_join(gsva_sig_tcell_long, pData(geomx_obj)[gsva_vars])
fwrite(gsva_sig_tcell_long, file.path(output_dir, 'gsva', paste0('gsva_deconv_tcell_additional_', deconv_type, '.csv')))

# GSVA on all Hallmark + CP + Go:BP ---------------------------------------
# do GSVA on all Hallmark + CP from msigDB
# read expression mtx
expr_mtx <- assayDataElement(geomx_obj, elt = norm_type)

if(do_gsva_hal_cp){
  # all Hallmark + CP from msigDB
  msigdb_df <- msigdbr(species = "Homo sapiens")
  msigdb_df <- filter(msigdb_df, gs_cat == 'H' | 
                        gs_subcat %in% c('CP:BIOCARTA', 'CP:KEGG', 'CP:REACTOME', 'CP:PID', 'CP:WIKIPATHWAYS', 'GO:BP'))
  
  msigdb_df$gene_symbol_adj <- adjust_synonym_genes(rownames(geomx_obj), msigdb_df$gene_symbol)
  
  hal_cp_list <- lapply(unique(msigdb_df$gs_name), function(x){
    gs <- filter(msigdb_df, gs_name == x)
    gs_genes <- unique(gs$gene_symbol_adj)
  })
  
  names(hal_cp_list) <- unique(msigdb_df$gs_name)
  
  #saveRDS(hal_cp_list, file.path(output_dir, 'gsva', 'hal_cp_adj_names.rds'))
  
  # do gsva
  gsva_hal_cp_all <- gsva(gsvaParam(expr_mtx, hal_cp_list, kcdf="Gaussian", minSize = 5))
  gsva_hal_cp_macro <- gsva(gsvaParam(deconv_macro, hal_cp_list, kcdf="Gaussian", minSize = 5))
  gsva_hal_cp_tcell <- gsva(gsvaParam(deconv_tcell, hal_cp_list, kcdf="Gaussian", minSize = 5))
  
  # do ssgsea
  # ssgsea_hal_cp <- gsva(expr_mtx, hal_cp_list, method = 'ssgsea', kcdf="Poisson", min.sz = 5)
  gsva_list <- list('all' = gsva_hal_cp_all, paste0('deconv_macro_', deconv_type) = gsva_hal_cp_macro, 
                    paste0('deconv_tcell_', deconv_type) = gsva_hal_cp_tcell)
  
  sapply(1:length(gsva_list), function(x){
    # adjust df
    gsva_hal_cp_long <- melt(gsva_list[x])
    colnames(gsva_hal_cp_long) <- c('pathway','dcc_filename', 'gsva_score', 'expr_signal')
    gsva_hal_cp_long <- left_join(gsva_hal_cp_long, pData(geomx_obj)[gsva_vars])
    
    fwrite(gsva_hal_cp_long, file.path(output_dir, 'gsva', paste0('gsva_hal_cp_', names(gsva_list)[x], '.csv')))
  })

}

###################
# do GSVA on macro + tcell list on the whole signal
gsva_macro_tcell_all <- gsva(gsvaParam(expr_mtx, c(sig_list_macro, sig_list_tcell), kcdf="Gaussian", minSize = 5))

gsva_macro_tcell_all_long <- melt(gsva_macro_tcell_all)
colnames(gsva_macro_tcell_all_long) <- c('pathway','dcc_filename', 'gsva_score')
gsva_macro_tcell_all_long <- left_join(gsva_macro_tcell_all_long, pData(geomx_obj)[gsva_vars])
fwrite(gsva_macro_tcell_all_long, file.path(output_dir, 'gsva', 'gsva_macro_tcell_additional_all.csv'))

###############################################################################
###############################################################################
# adjusting gsva scores from full signal for spatialdecon cell freq -------

sd_deconv <- readRDS(input_sd_deconv_path)
sd_deconv <- data.frame(pData(sd_deconv)[, 'prop_of_all'])
sd_deconv <- dplyr::select(sd_deconv, -Mast.cells, -other)
colnames(sd_deconv) <- paste0('deconv_', colnames(sd_deconv))
deconv_names <- colnames(sd_deconv)
sd_deconv <- tibble::rownames_to_column(sd_deconv, 'dcc_filename')

gsva_all_long <- fread(file.path(output_dir, 'gsva', 'gsva_hal_cp_all.csv'))
gsva_all_long <- left_join(gsva_all_long, sd_deconv)

gsva_macro_tcell_all_long <- left_join(gsva_macro_tcell_all_long, sd_deconv)

hal_cp_selected <- fread(hal_cp_selected_path)

###############
gsva_df <- gsva_all_long

gsva_df <- filter(gsva_df, pathway %in% hal_cp_selected$pathway)

##############
gsva_lm <- lapply(unique(as.vector(gsva_df$pathway)), function(path_name){
  gsva_path <- gsva_df[gsva_df$pathway == path_name, ]
  
  lapply(deconv_names, function(ct){
    
    # fit lm with ct fraction as explanatory var
    lm_res <- lm(gsva_score~get(ct),data=gsva_path)
    lm_coef <- summary(lm_res)$coefficients[2]
    lm_rsq <- summary(lm_res)$adj.r.squared
    
    gsva_lm_res <- list('pathway' = path_name, 'deconv_ct' = ct,
                              lm_coef = lm_coef, lm_rsq = lm_rsq)
    
    png(file = file.path(output_dir, 'gsva', 'sd_lm_hal_cp_selected', paste0('scatter_', path_name, '_', ct, '.png')))
    plot(gsva_path[[ct]], gsva_path$gsva_score, xlab = path_name, ylab = ct)
    abline(lm(gsva_score~get(ct),data=gsva_path),col='red')
    dev.off()
    
    return(gsva_lm_res)
  })
})

gsva_lm <- unlist(gsva_lm, recursive = F)
gsva_lm_df <- rbindlist(gsva_lm, fill=TRUE)

fwrite(gsva_lm_df, file.path(output_dir, 'gsva', 'gsva_hal_cp_sd_lm.csv'))


###############
# testing
path_name <- unique(as.vector(gsva_df$pathway))[4]

gsva_path <- gsva_df[gsva_df$pathway == path_name, ]

plot(gsva_path$deconv_Macrophages, gsva_path$gsva_score,xlab = path_name)
abline(lm(gsva_score~deconv_Macrophages,data=gsva_path),col='red') 
lm_res <- lm(gsva_score~deconv_Macrophages,data=gsva_path)
summary(lm_res)

# adjusted for lin reg
# y = a + xb
# y = a // -xb
lm_coef <- summary(lm_res)$coefficients[2]
lm_rsq <- summary(lm_res)$adj.r.squared
gsva_path$gsva_adj <- (gsva_path$gsva_score) - (gsva_path$deconv_Macrophages * 2.39702)
gsva_path$macro_adj <- gsva_path$deconv_Macrophages

plot(gsva_path$macro_adj, gsva_path$gsva_adj,xlab = path_name)
abline(lm(gsva_adj~macro_adj,data=gsva_path),col='red') 
lm_res <- lm(gsva_adj~macro_adj,data=gsva_path)
summary(lm_res)

# VISUALISATION IN  GSVA_VISUALISATION.R  

# PROGENy scores ----------------------------------------------------------

prog_noperm <- progeny(
  geomx_obj@assayData[[norm_type]],
  scale = TRUE,
  organism = "Human",
  top = 100,
)

prog_perm <- progeny(
  geomx_obj@assayData[[norm_type]],
  organism = "Human",
  top = 100,
  perm = 10,
  z_scores = FALSE,
  get_nulldist = TRUE
)

prog_perm[[1]] <- t(prog_perm[[1]])
rownames(prog_perm[[1]]) <- gsub('\\.', '\\-', rownames(prog_perm[[1]]))
rownames(prog_perm[[1]]) <- gsub('\\-dcc', '\\.dcc', rownames(prog_perm[[1]]))

#TODO use pathway significance infor from prog_perm[[2]] (nulldist)

# adjust and save dfs
prog_list <- list(noperm = prog_noperm, perm = prog_perm[[1]])

sapply(1:length(prog_list), function(x){
  prog_df <- prog_list[[x]]
  prog_name <- names(prog_list)[x]
  
  # adjust df
  prog_long <- melt(prog_df)
  colnames(prog_long) <- c('dcc_filename', 'progeny_path', 'progeny_score')
  prog_long <- left_join(prog_long, pData(geomx_obj)[gsva_vars])
  
  fwrite(prog_long, file.path(output_dir, 'progeny', paste0('progeny_', prog_name, '.csv')))
})



##############
# make boxplots

# per Anno cell type
# prog_boxpl <- ggplot(data = prog_long, aes(x = progeny_path, y = progeny_score, color = Annotation_cell)) +
#   geom_boxplot() +
#   facet_wrap(~Segment, scales = "fixed", dir="v") +
#   geom_pwc(method = "t_test", label = "p.signif", hide.ns = TRUE) +
#   theme(axis.text.x = element_text(angle=45, hjust=1)) +
#   ggtitle(paste('progeny', prog_name, 'scores'))
# 
# plot(prog_boxpl)
# ggsave(file.path(output_dir, 'progeny', paste0('box_progeny_anno_', prog_name, '.png')), 
#        height = 2000, width = 3000, unit = 'px')
# 
# ################################
# # per NACT status
# 
# prog_boxpl_nact_all <- ggplot(data = prog_long, aes(x = progeny_path, y = progeny_score, color = `NACT status`)) +
#   geom_boxplot() +
#   facet_wrap(~Segment, scales = "fixed", dir="v") +
#   geom_pwc(method = "t_test", label = "p.signif", hide.ns = TRUE) +
#   theme(axis.text.x = element_text(angle=45, hjust=1))+
#   ggtitle(paste('progeny', prog_name, 'scores'))
# 
# plot(prog_boxpl_nact_all)
# ggsave(file.path(output_dir, 'progeny', paste0('box_progeny_prepost_all_', prog_name, '.png')), 
#        height = 2000, width = 3000, unit = 'px')
# 
# prog_boxpl_nact_peranno <- ggplot(data = prog_long, aes(x = progeny_path, y = progeny_score, color = `NACT status`)) +
#   geom_boxplot() +
#   facet_wrap(Segment~Annotation_cell, scales = "fixed", ncol=4, nrow=2) +
#   geom_pwc(method = "t_test", label = "p.signif", hide.ns = TRUE) +
#   theme(axis.text.x = element_text(angle=45, hjust=1)) +
#   ggtitle(paste('progeny', prog_name, 'scores'))
# 
# plot(prog_boxpl_nact_peranno)
# ggsave(file.path(output_dir, 'progeny', paste0('box_progeny_prepost_anno_', prog_name, '.png')), 
#        height = 2000, width = 4000, unit = 'px')
# 
# 
# ################################
# # per PFS in post samples
# prog_long_post <- filter(prog_long, `NACT status` == 'post')
# 
# prog_boxpl_pfs_all <- ggplot(data = prog_long_post, aes(x = progeny_path, y = progeny_score, color = PFS)) +
#   geom_boxplot() +
#   facet_wrap(~Segment, scales = "fixed", dir="v") +
#   geom_pwc(method = "t_test", label = "p.signif", hide.ns = TRUE) +
#   theme(axis.text.x = element_text(angle=45, hjust=1))+
#   ggtitle(paste('progeny', prog_name, 'scores'))
# 
# plot(prog_boxpl_pfs_all)
# ggsave(file.path(output_dir, 'progeny', paste0('box_progeny_pfs_all_', prog_name, '.png')), height = 2000, width = 3000, unit = 'px')
# 
# prog_boxpl_pfs_peranno <- ggplot(data = prog_long_post, aes(x = progeny_path, y = progeny_score, color = PFS)) +
#   geom_boxplot() +
#   facet_wrap(Segment~Annotation_cell, scales = "fixed", ncol=4, nrow=2) +
#   geom_pwc(method = "t_test", label = "p.signif", hide.ns = TRUE) +
#   theme(axis.text.x = element_text(angle=45, hjust=1))+
#   ggtitle(paste('progeny', prog_name, 'scores'))
# 
# plot(prog_boxpl_pfs_peranno)
# ggsave(file.path(output_dir, 'progeny', paste0('box_progeny_pfs_anno_', prog_name, '.png')), 
#        height = 2000, width = 4000, unit = 'px')




########################
#######################
#######################
# OLD MESSY PART OF A SCRIPT

# -------------------------------------------------------------------------



# # signature files
# sig_t_exhaustion <- '/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/t_cell_exhaustion.csv'
# sig_macro <- '/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/macrophages.csv'
# sig_mhc <- '/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/MHC.csv'
# 
# sig_all <- list(sig_t_exhaustion, sig_macro, sig_mhc)
# 
# ##########################################################################
# ##########################################################################
# ##########################################################################
# # TODO PREPROCESSING ENDS HERE UPPER PART OF A SCRIPT SEPARATED FROM THE ANALYSIS BELOW
# # TODO SAVE GEOMX_OBJ AS EXTERNAJ OBJ AND LOAD FOR THE NEXT SCRIPTS
# # make DGE between selected ROI groups ------------------------------------
# 
# # within slide analysis - with random slope in LLM
# # comparison between ++ (posCD8_posIBA1) and other groups
# 
# # convert test variables to factors
# pData(geomx_obj)[["Annotation_cell_factor"]] <- factor(pData(geomx_obj)[["Annotation_cell"]])
# pData(geomx_obj)[["Sample_factor"]] <- factor(pData(geomx_obj)[["Sample"]])
# pData(geomx_obj)[["NACT_status_factor"]] <- factor(pData(geomx_obj)[["NACT status"]])
# pData(geomx_obj)[["PFS_factor"]] <- factor(pData(geomx_obj)[["PFS"]])
# 
# # convert normalized counts to log scale
# assayDataElement(object = geomx_obj, elt = "log_q3_norm") <-
#   assayDataApply(geomx_obj, 2, FUN = log, base = 2, elt = "q3_norm")
# 
# assayDataElement(object = geomx_obj, elt = "log_quant_norm") <-
#   assayDataApply(geomx_obj, 2, FUN = log, base = 2, elt = "quant_norm")
# 
# # TODO same for quant norm
# # run LMM:
# # formula follows conventions defined by the lme4 package
# results <- c()
# for(segment in c("tumor", "stroma")){
#   # careful! this is from all table with umap made for all ROIs - should be change to avoid confusion
#   geomx_segment <- geomx_obj[, geomx_obj@phenoData@data$Segment == segment]
#   for(status in c("pre", "post")) {
#     ind <- geomx_segment@phenoData@data$`NACT status` == status
#     
#     mixedOutmc <-
#       mixedModelDE(geomx_segment[, ind],
#                    elt = "log_q3_norm",
#                    modelFormula = ~ Annotation_cell_factor + (1 + Annotation_cell_factor | Sample_factor), # random slope
#                    groupVar = "Annotation_cell_factor",
#                    nCores = (parallel::detectCores() - 1),
#                    multiCore = FALSE)
#     
#     
#     # format results as data.frame
#     r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
#     tests <- rownames(r_test)
#     r_test <- as.data.frame(r_test)
#     r_test$Contrast <- tests
#     
#     # use lapply in case you have multiple levels of your test factor to
#     # correctly associate gene name with it's row in the results table
#     r_test$Gene <- 
#       unlist(lapply(colnames(mixedOutmc),
#                     rep, nrow(mixedOutmc["lsmeans", ][[1]])))
#     r_test$Subset <- status
#     r_test$Segment <- segment
#     r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
#     r_test <- r_test[, c("Gene", "Subset", "Segment",  "Contrast", "Estimate", 
#                          "Pr(>|t|)", "FDR")]
#     results <- rbind(results, r_test)
#   }
# }
# 
# fwrite(results, file.path(output_dir, 'dge/dge_annotation_cell_pre_post_separately.csv'))
# 
# # without differentiation to pre and post
# 
# results2 <- c()
# for(segment in c("tumor", "stroma")){
#   # careful! this is from all table with umap made for all ROIs - should be change to avoid confusion
#   geomx_segment <- geomx_obj[, geomx_obj@phenoData@data$Segment == segment]
#   
#   mixedOutmc <-
#     mixedModelDE(geomx_segment,
#                  elt = "log_q3_norm",
#                  modelFormula = ~ Annotation_cell_factor + (1 + Annotation_cell_factor | Sample_factor), # random slope
#                  groupVar = "Annotation_cell_factor",
#                  nCores = (parallel::detectCores() - 1),
#                  multiCore = FALSE)
#   
#   
#   # format results as data.frame
#   r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
#   tests <- rownames(r_test)
#   r_test <- as.data.frame(r_test)
#   r_test$Contrast <- tests
#   
#   # use lapply in case you have multiple levels of your test factor to
#   # correctly associate gene name with it's row in the results table
#   r_test$Gene <- 
#     unlist(lapply(colnames(mixedOutmc),
#                   rep, nrow(mixedOutmc["lsmeans", ][[1]])))
#   r_test$Segment <- segment
#   r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
#   r_test <- r_test[, c("Gene", "Segment",  "Contrast", "Estimate", 
#                        "Pr(>|t|)", "FDR")]
#   results2 <- rbind(results2, r_test)
# }
# 
# fwrite(results2, file.path(output_dir, 'dge/dge_annotation_cell_all.csv'))
# 
# 
# #####################################
# results_signif <- results[results$FDR <= 0.05, ]
# results2_signif <- results2[results2$FDR <= 0.05, ] 
# 
# fwrite(results_signif, file.path(output_dir, 'dge/dge_annotation_cell_pre_post_separately_signif.csv'))
# fwrite(results2_signif, file.path(output_dir, 'dge/dge_annotation_cell_all_signif.csv'))
# 
# # TODO redo for 1group vs 3groups all together
# 
# ######################################
# # BETWEEN SLIDES COMPARISON
# 
# # run LMM without random slope:
# # formula follows conventions defined by the lme4 package
# # results_pre_post_per_cell <- c()
# # for(segment in c("tumor", "stroma")){
# #   # careful! this is from all table with umap made for all ROIs - should be change to avoid confusion
# #   geomx_segment <- geomx_obj[, geomx_obj@phenoData@data$Segment == segment]
# #   for(anno_cell in unique(geomx_obj@phenoData@data$Annotation_cell)) {
# #     ind <- geomx_segment@phenoData@data$Annotation_cell == anno_cell
# #     
# #     # TODO trycatch if too litle nr of ROIs - return empty frame
# #     mixedOutmc <-
# #       mixedModelDE(geomx_segment[, ind],
# #                    elt = "log_q3_norm",
# #                    modelFormula = ~ NACT_status_factor + (1 | Sample_factor), # random slope
# #                    groupVar = "NACT_status_factor",
# #                    nCores = (parallel::detectCores() - 1),
# #                    multiCore = FALSE)
# #     
# #     
# #     # format results as data.frame
# #     r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
# #     tests <- rownames(r_test)
# #     r_test <- as.data.frame(r_test)
# #     r_test$Contrast <- tests
# #     
# #     # use lapply in case you have multiple levels of your test factor to
# #     # correctly associate gene name with it's row in the results table
# #     r_test$Gene <- 
# #       unlist(lapply(colnames(mixedOutmc),
# #                     rep, nrow(mixedOutmc["lsmeans", ][[1]])))
# #     r_test$Annotation_cell <- anno_cell
# #     r_test$Segment <- segment
# #     r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
# #     r_test <- r_test[, c("Gene", "Annotation_cell", "Segment",  "Contrast", "Estimate", 
# #                          "Pr(>|t|)", "FDR")]
# #     results_pre_post_per_cell <- rbind(results_pre_post_per_cell, r_test)
# #   }
# # }
# # 
# # fwrite(results_pre_post_per_cell, file.path(output_dir, 'dge/dge_pre_post_per_cell_separately.csv'))
# 
# #########################
# ##########################
# # all cell anno mixed together doesnt give any meaningful results!!!
# 
# results_pre_post_doublepos <- c()
# for(segment in c("tumor", "stroma")){
#   # careful! this is from all table with umap made for all ROIs - should be change to avoid confusion
#   geomx_segment <- geomx_obj[, geomx_obj@phenoData@data$Segment == segment]
#   geomx_segment_doublepos <- geomx_segment[, geomx_segment@phenoData@data$Annotation_cell == "posCD8_posIBA1"]
#   
#   mixedOutmc <-
#     mixedModelDE(geomx_segment_doublepos,
#                  elt = "log_q3_norm",
#                  modelFormula = ~ NACT_status_factor + (1 | Sample_factor), # random slope
#                  groupVar = "NACT_status_factor",
#                  nCores = (parallel::detectCores() - 1),
#                  multiCore = FALSE)
#   
#   
#   # format results as data.frame
#   r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
#   tests <- rownames(r_test)
#   r_test <- as.data.frame(r_test)
#   r_test$Contrast <- tests
#   
#   # use lapply in case you have multiple levels of your test factor to
#   # correctly associate gene name with it's row in the results table
#   r_test$Gene <- 
#     unlist(lapply(colnames(mixedOutmc),
#                   rep, nrow(mixedOutmc["lsmeans", ][[1]])))
#   r_test$Segment <- segment
#   r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
#   r_test <- r_test[, c("Gene", "Segment",  "Contrast", "Estimate", 
#                        "Pr(>|t|)", "FDR")]
#   results_pre_post_doublepos <- rbind(results_pre_post_doublepos, r_test)
# }
# 
# fwrite(results_pre_post_doublepos, file.path(output_dir, 'dge/dge_pre_post_doublepos.csv'))
# results_pre_post_doublepos_signif <- results_pre_post_doublepos[results_pre_post_doublepos$FDR <= 0.05, ]
# fwrite(results_pre_post_doublepos_signif, file.path(output_dir, 'dge/dge_pre_post_doublepos_signif.csv'))
# 
# ################################
# # for alltogether post samples long vs short pfs
# geomx_post <- geomx_obj[, geomx_obj@phenoData@data$`NACT status` == 'post']
# 
# results_pfs_doublepos <- c()
# for(segment in c("tumor", "stroma")){
#   # careful! this is from all table with umap made for all ROIs - should be change to avoid confusion
#   geomx_segment <- geomx_post[, geomx_post@phenoData@data$Segment == segment]
#   geomx_segment_doublepos <- geomx_segment[, geomx_segment@phenoData@data$Annotation_cell == "posCD8_posIBA1"]
#   
#   mixedOutmc <-
#     mixedModelDE(geomx_segment_doublepos,
#                  elt = "log_q3_norm",
#                  modelFormula = ~ PFS_factor + (1 | Sample_factor), # random slope
#                  groupVar = "PFS_factor",
#                  nCores = (parallel::detectCores() - 1),
#                  multiCore = FALSE)
#   
#   
#   # format results as data.frame
#   r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
#   tests <- rownames(r_test)
#   r_test <- as.data.frame(r_test)
#   r_test$Contrast <- tests
#   
#   # use lapply in case you have multiple levels of your test factor to
#   # correctly associate gene name with it's row in the results table
#   r_test$Gene <- 
#     unlist(lapply(colnames(mixedOutmc),
#                   rep, nrow(mixedOutmc["lsmeans", ][[1]])))
#   r_test$Segment <- segment
#   r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
#   r_test <- r_test[, c("Gene", "Segment",  "Contrast", "Estimate", 
#                        "Pr(>|t|)", "FDR")]
#   results_pfs_doublepos <- rbind(results_pre_post_doublepos, r_test)
# }
# 
# fwrite(results_pfs_doublepos, file.path(output_dir, 'dge/dge_pfs_doublepos.csv'))
# results_pfs_doublepos_signif <- results_pfs_doublepos[results_pfs_doublepos$FDR <= 0.05, ]
# fwrite(results_pfs_doublepos_signif, file.path(output_dir, 'dge/dge_pfs_doublepos_signif.csv'))
# 
# 
# #############################################################
# #############################################################
# # separately per each patient pre and post, all 
# paired_patients <- c("S015", "S027", "S032", "S084", "S139")
# 
# results_patient_pairs <- c()
# for(segment in c("tumor", "stroma")){
#   # careful! this is from all table with umap made for all ROIs - should be change to avoid confusion
#   geomx_segment <- geomx_obj[, geomx_obj@phenoData@data$Segment == segment]
#   
#   for(patient in paired_patients){
#     geomx_patient <- geomx_segment[, geomx_segment@phenoData@data$Patient == patient]
#     
#     mixedOutmc <-
#       mixedModelDE(geomx_patient,
#                    elt = "log_q3_norm",
#                    modelFormula = ~ NACT_status_factor + (1 | Sample_factor), # random slope
#                    groupVar = "NACT_status_factor",
#                    nCores = (parallel::detectCores() - 1),
#                    multiCore = FALSE)
#     
#     
#     # format results as data.frame
#     r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
#     tests <- rownames(r_test)
#     r_test <- as.data.frame(r_test)
#     r_test$Contrast <- tests
#     
#     # use lapply in case you have multiple levels of your test factor to
#     # correctly associate gene name with it's row in the results table
#     r_test$Gene <- 
#       unlist(lapply(colnames(mixedOutmc),
#                     rep, nrow(mixedOutmc["lsmeans", ][[1]])))
#     r_test$Segment <- segment
#     r_test$Patient <- patient
#     r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
#     r_test <- r_test[, c("Gene", "Segment", 'Patient',  "Contrast", "Estimate", 
#                          "Pr(>|t|)", "FDR")]
#     results_patient_pairs <- rbind(results_patient_pairs, r_test)
#   }
# }
# 
# fwrite(results_patient_pairs, file.path(output_dir, 'dge/dge_patient_pairs.csv'))
# results_patient_pairs_signif <- results_patient_pairs[results_patient_pairs$FDR <= 0.05, ]
# fwrite(results_patient_pairs_signif, file.path(output_dir, 'dge/dge_patient_pairs_signif.csv'))
# 
# ############
# # separately for patients, only doublepos
# results_patient_pairs_doublepos <- c()
# for(segment in c("tumor", "stroma")){
#   # careful! this is from all table with umap made for all ROIs - should be change to avoid confusion
#   geomx_segment <- geomx_obj[, geomx_obj@phenoData@data$Segment == segment]
#   geomx_segment_doublepos <- geomx_segment[, geomx_segment@phenoData@data$Annotation_cell == "posCD8_posIBA1"]
#   
#   for(patient in paired_patients){
#     geomx_patient <- geomx_segment_doublepos[, geomx_segment_doublepos@phenoData@data$Patient == patient]
#     
#     mixedOutmc <-
#       mixedModelDE(geomx_patient,
#                    elt = "log_q3_norm",
#                    modelFormula = ~ NACT_status_factor + (1 | Sample_factor), # random slope
#                    groupVar = "NACT_status_factor",
#                    nCores = (parallel::detectCores() - 1),
#                    multiCore = FALSE)
#     
#     
#     # format results as data.frame
#     r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
#     tests <- rownames(r_test)
#     r_test <- as.data.frame(r_test)
#     r_test$Contrast <- tests
#     
#     # use lapply in case you have multiple levels of your test factor to
#     # correctly associate gene name with it's row in the results table
#     r_test$Gene <- 
#       unlist(lapply(colnames(mixedOutmc),
#                     rep, nrow(mixedOutmc["lsmeans", ][[1]])))
#     r_test$Segment <- segment
#     r_test$Patient <- patient
#     r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
#     r_test <- r_test[, c("Gene", "Segment", 'Patient',  "Contrast", "Estimate", 
#                          "Pr(>|t|)", "FDR")]
#     results_patient_pairs_doublepos <- rbind(results_patient_pairs_doublepos, r_test)
#   }
# }
# 
# fwrite(results_patient_pairs_doublepos, file.path(output_dir, 'dge/dge_patient_pairs_doublepos.csv'))
# results_patient_pairs_doublepos_signif <- results_patient_pairs_doublepos[results_patient_pairs_doublepos$FDR <= 0.05, ]
# fwrite(results_patient_pairs_doublepos_signif, file.path(output_dir, 'dge/dge_patient_pairs_doublepos_signif.csv'))
# 
# 
# ###########################################
# # check DEG results
# list.files(file.path(output_dir, 'dge/old'))
# 
# dge_anno <- fread(file.path(output_dir, 'dge/old/dge_annotation_cell_all_signif.csv'))
# dge_anno_prepost <- fread(file.path(output_dir, 'dge/old/dge_annotation_cell_pre_post_separately_signif.csv'))
# 
# # length = 0
# #dge_pfs <- fread(file.path(output_dir, 'dge/old/dge_pfs_all_signif.csv'))
# #dge_pfs_doublepos <- fread(file.path(output_dir, 'dge/old/dge_pfs_doublepos_signif.csv'))
# 
# #dge_prepost <- fread(file.path(output_dir, 'dge/old/dge_pre_post_all_signif.csv'))
# #dge_prepost_doublepos <- fread(file.path(output_dir, 'dge/old/dge_pre_post_doublepos_signif.csv'))
# 
# dge_patient_pairs <- fread(file.path(output_dir, 'dge/old/dge_patient_pairs_signif.csv'))
# dge_patient_pairs_doublepos <- fread(file.path(output_dir, 'dge/old/dge_patient_pairs_doublepos_signif.csv'))
# 
# dge_list <- list(dge_anno, dge_anno_prepost, dge_patient_pairs, dge_patient_pairs_doublepos)
# 
# ##################################################
# ###################################################
# # prepare volcano plots
# 
# library(ggrepel) 
# 
# ####
# # anno
# res_anno <- fread(file.path(output_dir, 'dge/old/dge_annotation_cell_all.csv'))
# 
# for(contrast in unique(res_anno$Contrast)){
#   
#   cont_pos <- unlist(strsplit(contrast, ' - '))[1]
#   cont_neg <- unlist(strsplit(contrast, ' - '))[2]
#   
#   plot_volcano_deg(filter(res_anno, Contrast == contrast),
#                    'anno', 10, cont_pos, cont_neg)
# }
# 
# #####
# # anno pre-post
# res_anno_prepost <- fread(file.path(output_dir, 'dge/old/dge_annotation_cell_pre_post_separately.csv'))
# 
# for(subset in unique(res_anno_prepost$Subset)){
#   for(contrast in unique(res_anno_prepost$Contrast)){
#     
#     cont_pos <- unlist(strsplit(contrast, ' - '))[1]
#     cont_neg <- unlist(strsplit(contrast, ' - '))[2]
#     
#     plot_volcano_deg(filter(res_anno_prepost, Contrast == contrast & Subset == subset),
#                      paste0('anno_', subset), 10, cont_pos, cont_neg)
#   }
# }
# 
# ######
# # patient_pairs all
# res_patient_pairs <- fread(file.path(output_dir, 'dge/old/dge_patient_pairs.csv'))
# 
# for(patient in unique(res_patient_pairs$Patient)){
#   for(contrast in unique(res_patient_pairs$Contrast)){
#     
#     cont_pos <- unlist(strsplit(contrast, ' - '))[1]
#     cont_neg <- unlist(strsplit(contrast, ' - '))[2]
#     
#     plot_volcano_deg(filter(res_patient_pairs, Contrast == contrast & Patient == patient),
#                      paste0('patient_pairs_all_', patient), 10, cont_pos, cont_neg)
#   }
# }
# 
# ######
# # patient_pairs doublepos
# res_patient_pairs_doublepos <- fread(file.path(output_dir, 'dge/old/dge_patient_pairs_doublepos.csv'))
# 
# for(patient in unique(res_patient_pairs_doublepos$Patient)){
#   for(contrast in unique(res_patient_pairs_doublepos$Contrast)){
#     
#     cont_pos <- unlist(strsplit(contrast, ' - '))[1]
#     cont_neg <- unlist(strsplit(contrast, ' - '))[2]
#     
#     plot_volcano_deg(filter(res_patient_pairs_doublepos, Contrast == contrast & Patient == patient),
#                      paste0('patient_pairs_doublepos_', patient), 10, cont_pos, cont_neg)
#   }
# }
# 
# #######################################
# #######################################
# 
# # ORA on all Hallmarks + CP -----------------------------------------------
# 
# 
# 
# # get bcg genes - all genes in dataset
# bcg_genes <- rownames(geomx_obj)
# 
# # prepare mdigdb
# msigdb_df <- msigdbr(species = "Homo sapiens")
# msigdb_df <- filter(msigdb_df, gs_cat %in% c("H", "C2") & gs_subcat != "CGP")
# 
# ################################
# # ora for dge_anno
# ora_anno <- data.frame()
# for(segment in unique(dge_anno$Segment)){
#   for(contrast in unique(dge_anno$Contrast)){
#     dge_pos <- filter(dge_anno, Segment == segment & Contrast == contrast & Estimate > 0)
#     dge_neg <- filter(dge_anno, Segment == segment & Contrast == contrast & Estimate < 0)
#     
#     ora_pos <- calculate_ora(dge_pos$Gene, bcg_genes, msigdb_df, padj = 0.05)
#     ora_neg <- calculate_ora(dge_neg$Gene, bcg_genes, msigdb_df, padj = 0.05)
#     
#     if(nrow(ora_pos > 0)){
#       ora_pos$direction <- 'up'
#     }
#     
#     if(nrow(ora_neg > 0)){
#       ora_neg$direction <- 'down'
#     }
#     
#     ora <- rbind(ora_pos, ora_neg)
#     
#     if(nrow(ora) > 0){
#       ora$contrast <- contrast
#       ora$segment <- segment
#       ora$GeneRatio_perc <- as.numeric(gsub("\\/[0-9]*", "", ora$GeneRatio))/
#         as.numeric(gsub("[0-9]*\\/", "", ora$GeneRatio))
#       
#       ora_anno <- rbind(ora_anno, ora)
#     }
#     
#   }
# }
# 
# fwrite(ora_anno, file.path(output_dir, 'dge/ora/ora_anno.csv'))
# #######################################################
# # ora for dge_anno_prepost
# 
# ora_anno_prepost <- data.frame()
# for(segment in unique(dge_anno_prepost$Segment)){
#   for(subset in unique(dge_anno_prepost$Subset)){
#     for(contrast in unique(dge_anno_prepost$Contrast)){
#       dge_pos <- filter(dge_anno_prepost, Segment == segment & Subset == subset & 
#                           Contrast == contrast & Estimate > 0)
#       dge_neg <- filter(dge_anno_prepost, Segment == segment & Subset == subset &
#                           Contrast == contrast & Estimate < 0)
#       
#       ora_pos <- calculate_ora(dge_pos$Gene, bcg_genes, msigdb_df, padj = 0.05)
#       ora_neg <- calculate_ora(dge_neg$Gene, bcg_genes, msigdb_df, padj = 0.05)
#       
#       if(nrow(ora_pos > 0)){
#         ora_pos$direction <- 'up'
#       }
#       
#       if(nrow(ora_neg > 0)){
#         ora_neg$direction <- 'down'
#       }
#       
#       ora <- rbind(ora_pos, ora_neg)
#       
#       if(nrow(ora) > 0){
#         ora$contrast <- contrast
#         ora$segment <- segment
#         ora$subset <- subset
#         ora$GeneRatio_perc <- as.numeric(gsub("\\/[0-9]*", "", ora$GeneRatio))/
#           as.numeric(gsub("[0-9]*\\/", "", ora$GeneRatio))
#         
#         ora_anno_prepost <- rbind(ora_anno_prepost, ora)
#       }
#     }
#   }
# }
# fwrite(ora_anno_prepost, file.path(output_dir, 'dge/ora/ora_anno_prepost.csv'))
# #######################################################
# # ora for dge_patient_pairs
# 
# ora_patient_pairs <- data.frame()
# for(segment in unique(dge_patient_pairs$Segment)){
#   for(patient in unique(dge_patient_pairs$Patient)){
#     for(contrast in unique(dge_patient_pairs$Contrast)){
#       dge_pos <- filter(dge_patient_pairs, Segment == segment & Patient == patient & 
#                           Contrast == contrast & Estimate > 0)
#       dge_neg <- filter(dge_patient_pairs, Segment == segment & Patient == patient &
#                           Contrast == contrast & Estimate < 0)
#       
#       ora_pos <- calculate_ora(dge_pos$Gene, bcg_genes, msigdb_df, padj = 0.05)
#       ora_neg <- calculate_ora(dge_neg$Gene, bcg_genes, msigdb_df, padj = 0.05)
#       
#       if(nrow(ora_pos > 0)){
#         ora_pos$direction <- 'up'
#       }
#       
#       if(nrow(ora_neg > 0)){
#         ora_neg$direction <- 'down'
#       }
#       
#       ora <- rbind(ora_pos, ora_neg)
#       
#       if(nrow(ora) > 0){
#         ora$contrast <- contrast
#         ora$segment <- segment
#         ora$patient <- patient
#         ora$GeneRatio_perc <- as.numeric(gsub("\\/[0-9]*", "", ora$GeneRatio))/
#           as.numeric(gsub("[0-9]*\\/", "", ora$GeneRatio))
#         
#         ora_patient_pairs <- rbind(ora_patient_pairs, ora)
#       }
#     }
#   }
# }
# 
# fwrite(ora_patient_pairs, file.path(output_dir, 'dge/ora/ora_patient_pairs.csv'))
# #######################################################
# # ora for dge_patient_pairs_doublepos
# 
# ora_patient_pairs_doublepos <- data.frame()
# for(segment in unique(dge_patient_pairs_doublepos$Segment)){
#   for(patient in unique(dge_patient_pairs_doublepos$Patient)){
#     for(contrast in unique(dge_patient_pairs_doublepos$Contrast)){
#       dge_pos <- filter(dge_patient_pairs_doublepos, Segment == segment & Patient == patient & 
#                           Contrast == contrast & Estimate > 0)
#       dge_neg <- filter(dge_patient_pairs_doublepos, Segment == segment & Patient == patient &
#                           Contrast == contrast & Estimate < 0)
#       
#       ora_pos <- calculate_ora(dge_pos$Gene, bcg_genes, msigdb_df, padj = 0.05)
#       ora_neg <- calculate_ora(dge_neg$Gene, bcg_genes, msigdb_df, padj = 0.05)
#       
#       if(nrow(ora_pos > 0)){
#         ora_pos$direction <- 'up'
#       }
#       
#       if(nrow(ora_neg > 0)){
#         ora_neg$direction <- 'down'
#       }
#       
#       ora <- rbind(ora_pos, ora_neg)
#       
#       if(nrow(ora) > 0){
#         ora$contrast <- contrast
#         ora$segment <- segment
#         ora$patient <- patient
#         ora$GeneRatio_perc <- as.numeric(gsub("\\/[0-9]*", "", ora$GeneRatio))/
#           as.numeric(gsub("[0-9]*\\/", "", ora$GeneRatio))
#         
#         ora_patient_pairs_doublepos <- rbind(ora_patient_pairs_doublepos, ora)
#       }
#     }
#   }
# }
# 
# fwrite(ora_patient_pairs_doublepos, file.path(output_dir, 'dge/ora/ora_patient_pairs_doublepos.csv'))
# 
# 
# # PROGENy scores ----------------------------------------------------------
# 
# 
# dim(geomx_obj)
# colnames(geomx_obj)[1:10]
# rownames(geomx_obj)[1:10]
# 
# prog_noperm <- progeny(
#   geomx_obj@assayData$q3_norm,
#   scale = TRUE,
#   organism = "Human",
#   top = 100,
#   perm = 1
# )
# 
# 
# prog_perm <- progeny(
#   geomx_obj@assayData$q3_norm,
#   organism = "Human",
#   top = 100,
#   perm = 10,
#   z_scores = FALSE,
#   get_nulldist = FALSE
# )
# 
# rownames(prog_perm) <- gsub('\\.', '\\-', rownames(prog_perm))
# rownames(prog_perm) <- gsub('\\-dcc', '\\.dcc', rownames(prog_perm))
# 
# 
# #####
# #TODO do it in loop
# prog_df <- prog_perm
# prog_name <- 'perm'
# 
# # adjust df
# prog_long <- melt(prog_df)
# colnames(prog_long) <- c('dcc_filename', 'progeny_path', 'progeny_score')
# prog_long <- left_join(prog_long, pData(geomx_obj)[c('dcc_filename', 'Segment', 'Annotation_cell', 'NACT status', 'PFS', 'Patient')])
# 
# fwrite(prog_long, file.path(output_dir, 'progeny', paste0('progeny_', prog_name, '.csv')))
# 
# ##############
# # make boxplots
# 
# # per Anno cell type
# prog_boxpl <- ggplot(data = prog_long, aes(x = progeny_path, y = progeny_score, color = Annotation_cell)) +
#   geom_boxplot() +
#   facet_wrap(~Segment, scales = "fixed", dir="v") +
#   geom_pwc(method = "t_test", label = "p.signif", hide.ns = TRUE) +
#   theme(axis.text.x = element_text(angle=45, hjust=1)) +
#   ggtitle(paste('progeny', prog_name, 'scores'))
# 
# plot(prog_boxpl)
# ggsave(file.path(output_dir, 'progeny', paste0('box_progeny_anno_', prog_name, '.png')), 
#        height = 2000, width = 3000, unit = 'px')
# 
# ################################
# # per NACT status
# 
# prog_boxpl_nact_all <- ggplot(data = prog_long, aes(x = progeny_path, y = progeny_score, color = `NACT status`)) +
#   geom_boxplot() +
#   facet_wrap(~Segment, scales = "fixed", dir="v") +
#   geom_pwc(method = "t_test", label = "p.signif", hide.ns = TRUE) +
#   theme(axis.text.x = element_text(angle=45, hjust=1))+
#   ggtitle(paste('progeny', prog_name, 'scores'))
# 
# plot(prog_boxpl_nact_all)
# ggsave(file.path(output_dir, 'progeny', paste0('box_progeny_prepost_all_', prog_name, '.png')), 
#        height = 2000, width = 3000, unit = 'px')
# 
# prog_boxpl_nact_peranno <- ggplot(data = prog_long, aes(x = progeny_path, y = progeny_score, color = `NACT status`)) +
#   geom_boxplot() +
#   facet_wrap(Segment~Annotation_cell, scales = "fixed", ncol=4, nrow=2) +
#   geom_pwc(method = "t_test", label = "p.signif", hide.ns = TRUE) +
#   theme(axis.text.x = element_text(angle=45, hjust=1)) +
#   ggtitle(paste('progeny', prog_name, 'scores'))
# 
# plot(prog_boxpl_nact_peranno)
# ggsave(file.path(output_dir, 'progeny', paste0('box_progeny_prepost_anno_', prog_name, '.png')), 
#        height = 2000, width = 4000, unit = 'px')
# 
# 
# ################################
# # per PFS in post samples
# prog_long_post <- filter(prog_long, `NACT status` == 'post')
# 
# prog_boxpl_pfs_all <- ggplot(data = prog_long_post, aes(x = progeny_path, y = progeny_score, color = PFS)) +
#   geom_boxplot() +
#   facet_wrap(~Segment, scales = "fixed", dir="v") +
#   geom_pwc(method = "t_test", label = "p.signif", hide.ns = TRUE) +
#   theme(axis.text.x = element_text(angle=45, hjust=1))+
#   ggtitle(paste('progeny', prog_name, 'scores'))
# 
# plot(prog_boxpl_pfs_all)
# ggsave(file.path(output_dir, 'progeny', paste0('box_progeny_pfs_all_', prog_name, '.png')), height = 2000, width = 3000, unit = 'px')
# 
# prog_boxpl_pfs_peranno <- ggplot(data = prog_long_post, aes(x = progeny_path, y = progeny_score, color = PFS)) +
#   geom_boxplot() +
#   facet_wrap(Segment~Annotation_cell, scales = "fixed", ncol=4, nrow=2) +
#   geom_pwc(method = "t_test", label = "p.signif", hide.ns = TRUE) +
#   theme(axis.text.x = element_text(angle=45, hjust=1))+
#   ggtitle(paste('progeny', prog_name, 'scores'))
# 
# plot(prog_boxpl_pfs_peranno)
# ggsave(file.path(output_dir, 'progeny', paste0('box_progeny_pfs_anno_', prog_name, '.png')), 
#        height = 2000, width = 4000, unit = 'px')
# 
# 
# 
# 
# # UCell scores for selected pathways --------------------------------------
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
# # GSVA and ssGSEA ---------------------------------------------------------
# library(GSVA)
# 
# # selected pathways
# #sig_texh_macro_mhc_list <- as.list(fread('/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/texh_macro_mhc.csv'))
# sig_texh_macro_mhc_list <- as.list(fread('/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/texh_macro_mhc_ifng_forpaper.csv'))
# sig_caf_revised_list <- as.list(fread('/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/stromal_cell_subtype_signatures_symbols_ensembl_ids_revised.csv'))
# 
# # pathways from Glasgow
# sig_pycr1 <- fread('~/Downloads/Gs_markers_02_04_10_23.csv')
# sig_pycr1_05 <- filter(sig_pycr1, avg_log2FC >= 0.5)
# pycr_list <- list('pycr1' = sig_pycr1$gene, 'pycr1_05' = sig_pycr1_05$gene)
# 
# pycr_sig <- data.frame('pycr1' = sig_pycr1$gene, 'pycr1_05' = sig_pycr1_05$gene)
# 
# # all Hallmark + CP from msigDB
# msigdb_df <- msigdbr(species = "Homo sapiens")
# msigdb_df <- filter(msigdb_df, gs_cat == 'H' | gs_subcat %in% c('CP:BIOCARTA', 'CP:KEGG', 'CP:REACTOME'))
# 
# hal_cp_list <- lapply(unique(msigdb_df$gs_name), function(x){
#   gs <- filter(msigdb_df, gs_name == x)
#   gs_genes <- unique(gs$human_gene_symbol)
# })
# 
# names(hal_cp_list) <- unique(msigdb_df$gs_name)
# 
# #####
# 
# expr_mtx <- assayDataElement(geomx_obj, elt = "q3_norm")
# View(expr_mtx[1:10, 1:10])
# 
# # do gsva on pathway lists
# gsva_hal_cp <- gsva(expr_mtx, hal_cp_list, method = 'gsva', kcdf="Poisson", min.sz = 5)
# gsva_texh_macro_mhc <- gsva(expr_mtx, sig_texh_macro_mhc_list, method = 'gsva', kcdf="Poisson", min.sz = 5)
# gsva_caf <- gsva(expr_mtx, sig_caf_revised_list, method = 'gsva', kcdf="Poisson", min.sz = 5)
# gsva_pycr1 <- gsva(expr_mtx, pycr_list, method = 'gsva', kcdf="Poisson", min.sz = 5)
# 
# # do ssgsea on selected lists
# # ssgsea_hal_cp <- gsva(expr_mtx, hal_cp_list, method = 'ssgsea', kcdf="Poisson", min.sz = 5)
# # ssgsea_texh_macro_mhc <- gsva(expr_mtx, sig_texh_macro_mhc_list, method = 'ssgsea', kcdf="Poisson", min.sz = 5)
# # ssgsea_caf <- gsva(expr_mtx, sig_caf_revised_list, method = 'ssgsea', kcdf="Poisson", min.sz = 5)
# 
# ######
# 
# # adjust df
# 
# # TODO do the same for hallmark+cp
# 
# gsva_texh_macro_mhc_long <- melt(gsva_texh_macro_mhc)
# colnames(gsva_texh_macro_mhc_long) <- c('pathway','dcc_filename', 'gsva_score')
# gsva_texh_macro_mhc_long <- left_join(gsva_texh_macro_mhc_long, pData(geomx_obj)[c('dcc_filename', 'Segment', 'Annotation_cell', 'NACT status', 'PFS', 'Patient')])
# 
# fwrite(gsva_texh_macro_mhc_long, file.path(output_dir, 'gsva', paste0('gsva_texh_macro_mhc_forpaper.csv')))
# 
# gsva_caf_long <- melt(gsva_caf)
# colnames(gsva_caf_long) <- c('pathway','dcc_filename', 'gsva_score')
# gsva_caf_long <- left_join(gsva_caf_long, pData(geomx_obj)[c('dcc_filename', 'Segment', 'Annotation_cell', 'NACT status', 'PFS')])
# 
# fwrite(gsva_caf_long, file.path(output_dir, 'gsva', paste0('gsva_caf.csv')))
# 
# gsva_list <- list(gsva_texh_macro_mhc_long, gsva_caf_long)
# names(gsva_list) <- c('texh_macro_mhc', 'caf')
# 
# gsva_pycr1_long <- melt(gsva_pycr1)
# colnames(gsva_pycr1_long) <- c('pathway','dcc_filename', 'gsva_score')
# gsva_pycr1_long <- left_join(gsva_pycr1_long, pData(geomx_obj)[c('dcc_filename', 'Segment', 'Annotation_cell', 'NACT status', 'PFS')])
# 
# fwrite(gsva_pycr1_long, file.path(output_dir, 'gsva', paste0('gsva_pycr1.csv')))
# 
# # VISUALISATION IN  GSVA_VISUALISATION.R  
# 
# 
# # individual genes expression ---------------------------------------------
# NECTIN2, TIGIT, CD96, LAG3, CD226, CXCL12, CXCR4, HGF, CD44
# PD1, TIM3
# CXCL9, CXCL10, CXCR3

# dir.create(file.path(output_dir, 'ind_genes'), showWarnings = T, recursive = T)
# 
# goi <- c('NECTIN2', 'TIGIT', 'CD96', 'LAG3', 'CD226', 'CXCL12', 'CXCR4', 'HGF', 
#          'CD44', 'PDCD1', 'HAVCR2', 'CXCL9', 'CXCL10', 'CXCR3')
# 
# expr_mtx <- assayDataElement(geomx_obj, elt = "q3_norm")
# expr_mtx_goi <- expr_mtx[goi, ]
# 
# expr_mtx_goi_long <- melt(expr_mtx_goi)
# colnames(expr_mtx_goi_long) <- c('gene','dcc_filename', 'expr')
# expr_mtx_goi_long <- left_join(expr_mtx_goi_long, pData(geomx_obj)[c('dcc_filename', 'Segment', 'Annotation_cell', 'NACT status', 'PFS', 'Patient')])
# fwrite(expr_mtx_goi_long, file.path(output_dir, 'ind_genes', paste0('ind_genes_exh_lr_fin.csv')))
