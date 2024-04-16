library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(GeoDiff)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggforce)
library(data.table)
library(cowplot)
library(preprocessCore)
library(Biobase)
library(reshape2)

library(umap)
library(Rtsne)

library(clusterProfiler)
library(msigdbr)
library(progeny)
library(reshape2)
library(biomaRt)
library(GSVA)
library(ggpubr)
library(topGO)
library(fgsea)


# get variables -----------------------------------------------------------
data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'
output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2'

input_rds_path <- file.path(output_dir, 'geomx_qc_norm.RDS')

imp_vars <- c("Segment", "Annotation_cell", "NACT status", "PFS") 
gsva_vars <- c(imp_vars, 'dcc_filename', 'Patient')

norm_type <- 'q3_norm' # either 'q3_norm' or 'quant_norm'

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'dge'), showWarnings = T, recursive = T)

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')

# load geomx obj from rds -------------------------------------------------

geomx_obj <- readRDS(input_rds_path)

# # make DGE between selected ROI groups ------------------------------------
# 
# # within slide analysis - with random slope in LLM
# # comparison between ++ (posCD8_posIBA1) and other groups
# 
# # convert test variables to factors
# for(col in c(imp_vars, 'Sample')){
#   pData(geomx_obj)[[paste0(col, "_factor")]] <- factor(pData(geomx_obj)[[col]])
# }
# 
# # convert normalized counts to log scale
# assayDataElement(object = geomx_obj, elt = paste0("log_", norm_type)) <-
#   assayDataApply(geomx_obj, 2, FUN = log, base = 2, elt = norm_type)
# 
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



# make DGE between selected ROI groups ------------------------------------

# within slide analysis - with random slope in LLM
# comparison between ++ (posCD8_posIBA1) and other groups

# convert test variables to factors
for(col in c(imp_vars, 'Sample')){
  pData(geomx_obj)[[paste0(col, "_factor")]] <- factor(pData(geomx_obj)[[col]])
}

# convert normalized counts to log scale
assayDataElement(object = geomx_obj, elt = paste0("log_", norm_type)) <-
  assayDataApply(geomx_obj, 2, FUN = log, base = 2, elt = norm_type)

# run LMM:
# formula follows conventions defined by the lme4 package
results <- c()
for(segment in c("tumor", "stroma")){
  # careful! this is from all table with umap made for all ROIs - should be change to avoid confusion
  geomx_segment <- geomx_obj[, geomx_obj@phenoData@data$Segment == segment]
  for(status in c("pre", "post")) {
    ind <- geomx_segment@phenoData@data$`NACT status` == status

    mixedOutmc <-
      mixedModelDE(geomx_segment[, ind],
                   elt = "log_q3_norm",
                   modelFormula = ~ Annotation_cell_factor + (1 + Annotation_cell_factor | Sample_factor), # random slope
                   groupVar = "Annotation_cell_factor",
                   nCores = (parallel::detectCores() - 1),
                   multiCore = FALSE)


    # format results as data.frame
    r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests

    # use lapply in case you have multiple levels of your test factor to
    # correctly associate gene name with it's row in the results table
    r_test$Gene <-
      unlist(lapply(colnames(mixedOutmc),
                    rep, nrow(mixedOutmc["lsmeans", ][[1]])))
    r_test$Subset <- status
    r_test$Segment <- segment
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
    r_test <- r_test[, c("Gene", "Subset", "Segment",  "Contrast", "Estimate",
                         "Pr(>|t|)", "FDR")]
    results <- rbind(results, r_test)
  }
}

fwrite(results, file.path(output_dir, 'dge/dge_annotation_cell_pre_post_separately.csv'))


# enrichment on DGE -------------------------------------------------------
# TODO examine FDR calculation on the lvl of DGE
dge_data_dir <- file.path("/media/iganiemi/T7-iga/st/geomx-processing/results/nact/dge")

#list.files(dge_data_dir)

#dge_df <- fread(file.path(dge_data_dir, "dge_annotation_cell_all.csv"))
dge_df <- fread(file.path(dge_data_dir, "old",  "dge_annotation_cell_pre_post_separately.csv"))
#dge_df <- fread(file.path(dge_data_dir, "dge_patient_pairs_doublepos.csv"))

# select thr
fc_thr <- 1
pval_thr <- 0.05
# filter to post and significant

dge_df_sig <- filter(dge_df, Subset == 'post' & FDR <= pval_thr & (Estimate >= fc_thr | Estimate <= -fc_thr))

table(dge_df$Contrast)

# subset to anno group
dge_sub <- dge_df_sig[dge_df_sig$Contrast == 'negCD8_negIBA1 - posCD8_posIBA1', ]

# convert to Entrez ID
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

# bcg genes - all from geomx dataset
gene_entrez_universe <- getBM(attributes=c('external_gene_name', 'entrezgene_id'),
                    filters = 'external_gene_name',
                    values = rownames(geomx_obj),
                    mart = ensembl)

#rmv duplicates
gene_entrez_universe <- gene_entrez_universe[!duplicated(gene_entrez_universe$external_gene_name),]

dge_sub <- left_join(dge_sub, gene_entrez_universe, by = c('Gene' = 'external_gene_name'))


# filter msigdb -----------------------------------------------------------

# Hallmark +CP
msigdb_df <- msigdbr(species = "Homo sapiens")
msigdb_df <- filter(msigdb_df, gs_cat %in% c("H", "C2", "C5") & !(gs_subcat %in% c("CGP", "GO:CC", "GO:MF", "HPO", "CP")))

msigdb_list <- lapply(unique(msigdb_df$gs_name), function(x){
  ls <- msigdb_df$gene_symbol[msigdb_df$gs_name == x]
})

names(msigdb_list) <- unique(msigdb_df$gs_name)


# ORA ---------------------------------------------------------------------

# ORA on GO
go_res <- as.data.frame(enrichGO(gene = as.character(unlist(dge_sub$entrezgene_id)),
                   ont = "BP",
                   OrgDb ="org.Hs.eg.db",
                   universe = as.character(unlist(gene_entrez_universe$entrezgene_id)),
                   readable=TRUE,
                   pvalueCutoff = pval_thr))


# ORA on H+CP
msigdb_res <- enricher(
  gene = as.character(unlist(dge_sub$entrezgene_id)),
  pvalueCutoff = pval_thr, # Can choose a FDR cutoff
  pAdjustMethod = "BH", 
  universe = as.character(unlist(gene_entrez_universe$entrezgene_id)), 
  TERM2GENE = dplyr::select(msigdb_df, gs_name, entrez_gene)
)

msigdb_res <- data.frame(msigdb_res@result) %>%
  filter(p.adjust <= pval_thr)


# GO enrichment with topGO ------------------------------------------------


# sampleGOdata <- new("topGOdata",
#                     + description = "Simple session", ontology = "BP",
#                     + allGenes = geneList, geneSel = topDiffGenes,
#                     + nodeSize = 10,
#                     + annot = annFUN.db, affyLib = affyLib)


# GSEA with fgsea ---------------------------------------------------------

rank_type <- 'rank_fdr_fcval' # c('rank_p', 'rank_p_fcval', 'rank_fdr_fcval')

gsea_all <- lapply(unique(dge_df$Contrast), function(contrast){
  lapply(unique(dge_df$Segment), function(segment){
    
    dge_all_sub <- filter(dge_df, Subset == 'post' &
                            Contrast == contrast &
                            Segment == segment)
    
    # TODO why is it different than the one saved from DGE calculation??
    dge_all_sub$FDR2 <- p.adjust(dge_all_sub$`Pr(>|t|)`, method = "fdr") 
    
     dge_all_sub$rank_p <- sign(dge_all_sub$Estimate)*(-log10(dge_all_sub$`Pr(>|t|)`))
     dge_all_sub$rank_p_fcval <- dge_all_sub$Estimate*(-log10(dge_all_sub$`Pr(>|t|)`))
     dge_all_sub$rank_fdr_fcval <- dge_all_sub$Estimate*(-log10(dge_all_sub$FDR2))
    
     # fdr rank alnone makes no sense - too many ties
     # plot(dge_all_sub$rank_p, dge_all_sub$rank_fdr) #linear + plateau - may be used interchangeably
     # plot(dge_all_sub$rank_p, dge_all_sub$rank_p_fcval) #linear + plateau - may be used interchangeably
     # plot(dge_all_sub$rank_fdr, dge_all_sub$rank_fdr_fcval) # linear
    
     dge_sub_rank <- dge_all_sub[[rank_type]]
     names(dge_sub_rank) <- dge_all_sub$Gene
     dge_sub_rank <- sort(dge_sub_rank, decreasing = T)
    
     plot(dge_sub_rank)
    
     # fix infinite ranks if needed
     # Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
     max_ranking <- max(dge_sub_rank[is.finite(dge_sub_rank)])
     min_ranking <- min(dge_sub_rank[is.finite(dge_sub_rank)])
     dge_sub_rank <- replace(dge_sub_rank, dge_sub_rank > max_ranking, max_ranking * 10)
     dge_sub_rank <- replace(dge_sub_rank, dge_sub_rank < min_ranking, min_ranking * 10)
     dge_sub_rank <- sort(dge_sub_rank, decreasing = TRUE) # sort genes by ranking
    
    
     gsea_res <- fgsea(pathways = msigdb_list, # List of gene sets to check
                       stats = dge_sub_rank,
                       scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                       minSize = 10,
                      maxSize = 500,
                       nproc = 18) # for parallelisation
    
    
     gsea_res <- arrange(gsea_res, padj) %>%
       filter(padj <= 0.01)
    
     # Select only independent pathways, removing redundancies/similar pathways
     collapsedPathways <- collapsePathways(gsea_res, msigdb_list, dge_sub_rank)
     mainPathways <- gsea_res[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
     gsea_res_main <- filter(gsea_res, pathway %in% mainPathways)
    
     gsea_res_main$Segment <- segment
     gsea_res_main$Contrast <- contrast

     # plotEnrichment(msigdb_list[[head(gsea_res[order(padj), ], 1)$pathway]],
     #                dge_sub_rank) +
     #   labs(title = head(gsea_res[order(padj), ], 1)$pathway)

    print(paste(contrast, segment, 'pass'))
    return(gsea_res_main)
  })
})


gsea_all <- do.call(rbind, unlist(gsea_all, recursive=FALSE))

fwrite(gsea_all, file.path(dge_data_dir, 'gsea', paste0('gsea_post_', rank_type, '.csv')))


# inspect GSEA results ----------------------------------------------------

gsea_res_path <- "/media/iganiemi/T7-iga/st/geomx-processing/results/nact/dge/gsea/gsea_post_rank_p_fcval.csv"

gsea_res <- fread(gsea_res_path)

# filter only for doublepositive comparisons
gsea_res <- gsea_res[grepl('posCD8_posIBA1', gsea_res$Contrast), ]
gsea_up <- gsea_res[gsea_res$NES > 0, ]
gsea_down <- gsea_res[gsea_res$NES < 0, ]

length(unique(gsea_up$pathway))
sort(table(gsea_up$pathway), decreasing = T)

gsea_up_stroma <- gsea_up[gsea_up$Segment == 'stroma', ]
gsea_up_tumor <- gsea_up[gsea_up$Segment == 'tumor', ]
