library(data.table)
library(dplyr)
library(msigdbr)


output_signs_path <- '~/Documents/phd/st/geomx-processing/results/batch123-2808/lr_interactions/multi_niche_netr/myelonets_signatures_list.RDS'
######
ct_name1 <- 'Tcells_CD8'
ct_name1 <- 'Macrophages_Monocytes'
ct_names_rec <- c('Macrophages_Monocytes', 'Tcells_CD8', 'Fibroblasts_Mesothelial')

group_name <- 'stroma_post_Macro_domin'
group_name_dge <- 'stroma_post'
######

nn_res_dir1 <- '~/Documents/phd/st/geomx-processing/results/batch123-2808/lr_interactions/multi_niche_netr/Macro_domin_vs_others_stroma_post_Macro_CD8_Fibro_onlyposdge_priorupdown'
nn_res_dir2 <- '~/Documents/phd/st/geomx-processing/results/batch123-2808/lr_interactions/multi_niche_netr/Macro_domin_vs_others_stroma_post_Macro_CD8_Fibro_priorupdown'


dge_path_macro <- "/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2808/dge/dge_within_slide_roi_cluster_label_gmm_bin_TRUE_Macro_domin_Segment_NACT_status/dge_deconv_Macrophages_Monocytes_dge_within_slide_roi_cluster_label_gmm_bin_TRUE_Macro_domin_Segment_NACT_status.csv"
cc_res_path <- "/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2808/lr_interactions/cell_chat/stroma_post_roi_cluster_gmm_deseq2log_harmony/CellChat_df_Segment_NACT_status_roi_cluster_label_gmm_tumor_Bcells_Tcells_CD4_Tcells_other_Tcells_CD8_Fibroblasts_Mesothelial_Macrophages_Monocytes_DCs_lr.csv"

nn_res_act_path1 <- file.path(nn_res_dir1,'plots_and_csv_files', paste0('stroma_post_Macro_domin-stroma_post_other_roi_type_', ct_name2, '_', ct_name1, '_LR_pairs_dfplot_ligand_activity.csv'))
nn_res_expr_path1 <- file.path(nn_res_dir1,'plots_and_csv_files', paste0('stroma_post_Macro_domin-stroma_post_other_roi_type_', ct_name2, '_', ct_name1, '_LR_pairs_dfplot_median_bulk_expr.csv'))
nn_res_act_path2 <- file.path(nn_res_dir2,'plots_and_csv_files', paste0('stroma_post_Macro_domin-stroma_post_other_roi_type_', ct_name2, '_', ct_name1, '_LR_pairs_dfplot_ligand_activity.csv'))
nn_res_expr_path2 <- file.path(nn_res_dir2,'plots_and_csv_files', paste0('stroma_post_Macro_domin-stroma_post_other_roi_type_', ct_name2, '_', ct_name1, '_LR_pairs_dfplot_median_bulk_expr.csv'))

cc_thr <- 0.1

# load files --------------------------------------------------------------

dge_res <- fread(dge_path_macro)
dge_res_group <- dge_res[dge_res$data_group == group_name_dge]
dge_res_group <- dge_res_group[dge_res_group$FDR <= 0.05, ]

cc_res <- fread(cc_res_path)
cc_res_group <- cc_res[cc_res$source == ct_name1 & cc_res$target == ct_name2 & cc_res$group == group_name, ]
cc_res_group <- cc_res_group[cc_res_group$prob >= cc_thr, ]

nn_act1 <- fread(nn_res_act_path1)
nn_expr1 <- fread(nn_res_expr_path1) %>% dplyr::rename(diff_median_contrast = 'group')
nn_res1 <- left_join(nn_act1, nn_expr1[, c('lr_interaction', 'diff_median', 'diff_median_contrast')])

nn_act2 <- fread(nn_res_act_path2)
nn_expr2 <- fread(nn_res_expr_path2) %>% dplyr::rename(diff_median_contrast = 'group')
nn_res2 <- left_join(nn_act2, nn_expr2[, c('lr_interaction', 'diff_median', 'diff_median_contrast')])


# filter nn and cc results
# remove ligands where both up and down activity scaled is negative
nn_res_group1 <- nn_res1 %>%
  filter(group == !!group_name) %>%
  group_by(lr_interaction) %>%
  filter(max(activity_scaled) > 0) %>%
  ungroup() %>%
  filter(prioritization_score >= 0.75)

nn_res_group2 <- nn_res2 %>%
  filter(group == !!group_name) %>%
  group_by(lr_interaction) %>%
  filter(max(activity_scaled) > 0) %>%
  ungroup() %>%
  filter(prioritization_score >= 0.75)


lr <- sort(unique(c(nn_res_group1$ligand, nn_res_group1$receptor)))
lr2 <- sort(unique(c(nn_res_group2$ligand, nn_res_group2$receptor)))
lr_cc <- sort(unique(c(cc_res_group$ligand, cc_res_group$receptor)))

length(intersect(lr, lr2))
setdiff(lr, lr2)
setdiff(lr2, lr)
length(intersect(lr_cc, lr))
length(intersect(lr_cc, lr2))

# save for copying

fwrite(nn_res_group1[nn_res_group1$direction_regulation == 'up', ], '~/Downloads/nn_res1.csv')
fwrite(cc_res_group, '~/Downloads/cc_res.csv')


# myelonets signature -----------------------------------------------------

# top 50 from most expressed genes + top 10 macro-macro, macro-cd8, macro-fibro and macro-dc ?
dge_res_group_overexpr <- dplyr::arrange(dge_res_group, -Estimate) %>%
  filter(Estimate >= 0.5)

# top LR for all receivers
ct_of_interest <- 'Macrophages_Monocytes'

ct_names_other <- c('Macrophages_Monocytes', 'Tcells_CD8', 'Fibroblasts_Mesothelial')

nn_res_all_rec <- lapply(ct_names_rec, function(ct_name2){
  ct_name1 <- ct_of_interest
  
  nn_res_act_path1 <- file.path(nn_res_dir1,'plots_and_csv_files', paste0('stroma_post_Macro_domin-stroma_post_other_roi_type_', ct_name2, '_', ct_name1, '_LR_pairs_dfplot_ligand_activity.csv'))
  nn_res_expr_path1 <- file.path(nn_res_dir1,'plots_and_csv_files', paste0('stroma_post_Macro_domin-stroma_post_other_roi_type_', ct_name2, '_', ct_name1, '_LR_pairs_dfplot_median_bulk_expr.csv'))
  
  nn_act1 <- fread(nn_res_act_path1)
  nn_expr1 <- fread(nn_res_expr_path1) %>% dplyr::rename(diff_median_contrast = 'group')
  nn_res1 <- left_join(nn_act1, nn_expr1[, c('lr_interaction', 'diff_median', 'diff_median_contrast')])
  
  return(nn_res1)
})


nn_res_all_send <- lapply(ct_names_rec, function(ct_name1){
  ct_name2 <- ct_of_interest
  
  nn_res_act_path1 <- file.path(nn_res_dir1,'plots_and_csv_files', paste0('stroma_post_Macro_domin-stroma_post_other_roi_type_', ct_name2, '_', ct_name1, '_LR_pairs_dfplot_ligand_activity.csv'))
  nn_res_expr_path1 <- file.path(nn_res_dir1,'plots_and_csv_files', paste0('stroma_post_Macro_domin-stroma_post_other_roi_type_', ct_name2, '_', ct_name1, '_LR_pairs_dfplot_median_bulk_expr.csv'))
  
  nn_act1 <- fread(nn_res_act_path1)
  nn_expr1 <- fread(nn_res_expr_path1) %>% dplyr::rename(diff_median_contrast = 'group')
  nn_res1 <- left_join(nn_act1, nn_expr1[, c('lr_interaction', 'diff_median', 'diff_median_contrast')])
  
  return(nn_res1)
})

nn_res_all_rec <- do.call(rbind, nn_res_all_rec)
nn_res_all_send <- do.call(rbind, nn_res_all_send)

nn_res_all <- rbind(nn_res_all_rec, nn_res_all_send) %>%
  distinct()

# loop through scores and top_n
# activity_thr <- 0.5
# diff_median_thr <- 0.5
# top_dge <- 10
# top_lr <- 10
# 

signs <- c()
sign_names <- c()

for(activity_thr in c(0.5, 1)){
  for(diff_median_thr in c(0.5, 1)){
    for(top_dge in c(10, 20, 50, 100)){
      for(top_lr in c(10)){
        
        # filter results:
        nn_res_all_filt <- nn_res_all %>%
          filter(group == !!group_name) %>%
          group_by(sender_receiver, lr_interaction) %>%
          filter(max(abs(activity_scaled)) >= activity_thr) %>%
          ungroup() %>%
          filter(diff_median >= diff_median_thr) %>%
          filter(prioritization_score >= 0.75) %>%
          group_by(sender_receiver, lr_interaction) %>%
          filter(activity_scaled == max(abs(activity_scaled))) %>% # get one LR pair per direction
          ungroup()
        
        # filter best activity and top 10 per send-rec pair
        nn_top_lig <- nn_res_all_filt %>%
          filter(sender == ct_of_interest) %>%
          group_by(sender_receiver, ligand) %>%
          filter(prioritization_score == max(prioritization_score)) %>% # choose best ligand per pair
          ungroup() %>%
          group_by(sender_receiver) %>%
          slice_max(order_by = prioritization_score, n = 10)
        
        nn_top_rec <- nn_res_all_filt %>%
          filter(receiver == ct_of_interest) %>%
          group_by(sender_receiver, receptor) %>%
          filter(prioritization_score == max(prioritization_score)) %>% # choose best receptor per pair
          ungroup() %>%
          group_by(sender_receiver) %>%
          slice_max(order_by = prioritization_score, n = 10)
        
        fin_signature <- unique(c(nn_top_lig$ligand, nn_top_rec$receptor, dge_res_group_overexpr$Gene[1:top_dge]))
        sign_name <- paste0('dgetop', as.character(top_dge), '_lrtop', as.character(top_lr), '_lractiv', as.character(activity_thr), '_lrexprdiff', as.character(diff_median_thr))
        
        # append 
        signs <- c(signs, list(fin_signature))
        sign_names <- c(sign_names, sign_name)
      }
    }
  }
}

names(signs) <- sign_names

saveRDS(signs, output_signs_path)

#######################################################################################
# signatures based on transcriptomics programmes
output_myeprog_path <- '~/Documents/phd/st/geomx-processing/results/batch123-2808/downstream/myelonets_signatures_validation/myelonets_pathways_selected.RDS'


additional_signs_csv <- fread('/home/iganiemi/Documents/phd/st/geomx-processing/data/signatures/eyemt_immune/additional_signatures_immune.csv')
msigdb_df <- msigdbr(species = "Homo sapiens")
msigdb_df <- filter(msigdb_df, gs_subcollection %in% c('CP:BIOCARTA', 'CP:REACTOME', 'GO:BP'))

msigdb_list <- lapply(unique(msigdb_df$gs_name), function(x){
  gs <- filter(msigdb_df, gs_name == x)
  gs_genes <- unique(gs$gene_symbol)
})
names(msigdb_list) <- unique(msigdb_df$gs_name)

additional_sign_list <- as.list(additional_signs_csv)
additional_sign_list <- lapply(additional_sign_list, function(l){l[l !=""]})

all_sign_list <- c(msigdb_list, additional_sign_list)

###################################
# create myeloid signatures
myelonets_lipid <- c('REACTOME_INTERLEUKIN_10_SIGNALING',
                          'REACTOME_ROS_AND_RNS_PRODUCTION_IN_PHAGOCYTES',
                          'GOBP_ACUTE_INFLAMMATORY_RESPONSE',
                          'ADDITIONAL_WANG_M2',
                          'GOBP_CHOLESTEROL_EFFLUX',
                          'GOBP_CHOLESTEROL_STORAGE',
                          'GOBP_FOAM_CELL_DIFFERENTIATION'                     
                        )

myelonets_il1 <- c('REACTOME_INTERLEUKIN_1_SIGNALING',
                          'REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION',
                          'REACTOME_SCAVENGING_BY_CLASS_A_RECEPTORS',
                          'REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY'
                          )

myelonets_tnfa <- c('REACTOME_TNF_SIGNALING',
                    'BIOCARTA_TGFB_PATHWAY',
                    'REACTOME_TNFR1_INDUCED_PROAPOPTOTIC_SIGNALING',
                    'REACTOME_TNFR1_INDUCED_NF_KAPPA_B_SIGNALING_PATHWAY'
                    
)

myelonets_csfil2vegf <- c('REACTOME_INTERLEUKIN_2_SIGNALING',
                    'REACTOME_SIGNALING_BY_CSF1_M_CSF_IN_MYELOID_CELLS',
                    'BIOCARTA_VEGF_PATHWAY'
)


myelonets_il1tnfcsf <- c('REACTOME_INTERLEUKIN_1_SIGNALING',
                         'REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION',
                         'REACTOME_SCAVENGING_BY_CLASS_A_RECEPTORS',
                         'REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY',
                         'REACTOME_TNF_SIGNALING',
                         'BIOCARTA_TGFB_PATHWAY',
                         'REACTOME_TNFR1_INDUCED_PROAPOPTOTIC_SIGNALING',
                         'REACTOME_TNFR1_INDUCED_NF_KAPPA_B_SIGNALING_PATHWAY',
                         'REACTOME_SIGNALING_BY_CSF1_M_CSF_IN_MYELOID_CELLS'
)

mye_prog_list <- list(myelonets_lipid, myelonets_il1, myelonets_tnfa, myelonets_csfil2vegf, myelonets_il1tnfcsf)
names(mye_prog_list) <- c('myelonets_lipid','myelonets_il1', 'myelonets_tnfa', 'myelonets_csfil2vegf', 'myelonets_il1tnfcsf')
########################################
# loading and selecting dge
dge_path_macro <- "/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2808/dge/dge_within_slide_roi_cluster_label_gmm_bin_TRUE_Macro_domin_Segment_NACT_status/dge_deconv_Macrophages_Monocytes_dge_within_slide_roi_cluster_label_gmm_bin_TRUE_Macro_domin_Segment_NACT_status.csv"
dge_res <- fread(dge_path_macro)
dge_res_group <- dge_res[dge_res$data_group == 'stroma_post']
dge_res_group <- dge_res_group[dge_res_group$FDR <= 0.05 & dge_res_group$Estimate >= 0.5, ]

##########################
# collecting genes for each programme
mye_prog_genes_list <- sapply(mye_prog_list, function(prog){
  prog_genes <- all_sign_list[names(all_sign_list) %in% prog]
  prog_genes <- unique(c(unlist(prog_genes)))
  return(prog_genes)
})
names(mye_prog_genes_list) <- names(mye_prog_list)

# intersection with dge
mye_prog_genes_list_dge <- sapply(mye_prog_genes_list, function(prog_genes){
  prog_genes_dge <- intersect(dge_res_group$Gene, prog_genes)
  return(prog_genes_dge)
})

names(mye_prog_genes_list_dge) <- paste0(names(mye_prog_list), '_dge')

mye_prog_genes_all <- c(mye_prog_genes_list, mye_prog_genes_list_dge)

saveRDS(mye_prog_genes_all, output_myeprog_path)

###########################################3
# just selected signatures one by one
mye_sel_sigs_names <- unique(c(myelonets_lipid, myelonets_il1, myelonets_tnfa, myelonets_csfil2vegf, myelonets_il1tnfcsf))

mye_sel_sigs_genes <- sapply(mye_sel_sigs_names, function(sig){
  sig_genes <- all_sign_list[names(all_sign_list) == sig]
  return(sig_genes)
})

saveRDS(mye_sel_sigs_genes, '~/Documents/phd/st/geomx-processing/results/batch123-2808/downstream/myelonets_signatures_validation/myelonets_signs_selected.RDS')

# intersection with dge
mye_sel_sigs_genes_dge <- sapply(mye_sel_sigs_genes, function(sig_genes){
  sig_genes_dge <- intersect(dge_res_group$Gene, sig_genes)
  if(length(sig_genes_dge) > 0){
    return(sig_genes_dge)
  }
})

mye_sel_sigs_genes_dge[sapply(mye_sel_sigs_genes_dge, is.null)] <- NULL

saveRDS(mye_sel_sigs_genes_dge, '~/Documents/phd/st/geomx-processing/results/batch123-2808/downstream/myelonets_signatures_validation/myelonets_signs_selected_dge.RDS')
