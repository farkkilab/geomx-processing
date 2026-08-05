library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ComplexHeatmap)
library(viridis)
library(scales)
library(tidyr)
library(circlize)
library(paletteer)
library(car)
library(ggpubr)
library(GeomxTools)
library(forcats)
library(ggpmisc)
library(umap)
library(PCAtools)


# load files --------------------------------------------------------------

proj_dir <<- '~/Documents/phd/st'
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808')
out_dir <- file.path(output_dir, 'downstream', 'figures_for_manuscript') 

#input files
geomx_path <- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS')
# clusters in metadata from geomx_relabel_roi_2nd_approach.R
metadt_path <- file.path(proj_dir, 'data/geomx/metadata_full_SENSITIVE.csv')

# deconv ct fractions per roi and aoi
# from geomx_roi_hubs_integration.R
# TODO confirm from github bcs it changed the output now per each deconv-cycif check
deconv_ct_count_path_4_imm_roi <- file.path(output_dir, 'deconvolution', 'ct_frac_deconv_roi_mid_lvl_ct_updated_4mainimmune.csv')
deconv_ct_count_path_4_imm_aoi <- file.path(output_dir, 'deconvolution', 'ct_frac_deconv_mid_lvl_ct_updated_4mainimmune.csv')

deconv_ct_count_path_roi <- file.path(output_dir, 'deconvolution', 'ct_frac_deconv_roi_mid_lvl_ct_updated.csv')
deconv_ct_count_path_aoi <- file.path(output_dir, 'deconvolution', 'ct_frac_deconv_mid_lvl_ct_updated.csv')

# ct fractions from geomx_roi_hubs_integration.R
ct_frac_all_path <- file.path(output_dir, 'cycif_integration/ct_frac_comparison_batch123tls', 'ct_frac_all_roi_batch123tls.csv')

#ssgsea scores from eyemt_downstream_analysis_gsea.R
gsea_res_path <- file.path(output_dir, 'pathway_analysis', 'gsea', 'gsea_all_cells_immune_reactome.csv')

# ssgsea_scores dir for msigdb and additional
# gsea_out_dir <- file.path(output_dir, 'pathway_analysis', 'gsea')
# 
# # names of pathways of interest
# sigs_path <- file.path(proj_dir, 'geomx-processing', 'data', 'signatures', 'eyemt_immune', 'immune_pathways_names.csv')
# sigs_path <- file.path(proj_dir, 'geomx-processing', 'data', 'signatures', 'eyemt_immune', 'immune_pathways_names_reactome.csv')
# # RDS object with list of pathways and their genes
# signs_genes_list_path <- file.path(proj_dir, 'geomx-processing', 'data', 'signatures', 'eyemt_immune', 'sign_list_immune.RDS')
# 
# # source and create output dir
# source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))
# dir.create(out_dir, recursive = T, showWarnings = F)

dir.create(out_dir)
# input params ------------------------------------------------------------

ct_all <- c("tumor", "Macrophages_Monocytes", "Tcells_CD8", "Tcells_CD4", "DCs", "Bcells", 
            "Tcells_other", "NKcells", "Endothelial_cells", "Fibroblasts_Mesothelial")
ct_of_interest <- c("tumor", "Macrophages_Monocytes", "Tcells_CD8", "Tcells_CD4", "DCs", "Bcells")
ct_names_immune <- c("Tcells_CD4", "Tcells_CD8", "DCs", "Macrophages_Monocytes", "Bcells", "NKcells")
ct_names_myeloids <- c("DCs", "Macrophages_Monocytes")
ct_names_lymphoids <- c("Tcells_CD4", "Tcells_CD8", "Bcells")

min_frac <- 0.01 # gsea scores computed for dcc with smaller fraction, will be removed

#clust_type <- paste0(clust_types[1], '_label')
clust_type <- 'roi_cluster_label_gmm'
clust_type_name <- 'gmm' # for plotting

#cols_segment_nact <- c('#FED18C','#9c5f01', '#A7AEF2', '#15208e') # pre-postnact
cols_segment_nact <- c('#508791', '#305157', '#D33F49','#861f26')
cols_nact <- c('#FED18C', '#A7AEF2')
cols_cell_types_imm <- c('#7d28d7', '#00BA38', '#0095B6', '#F8766D', '#ffff41', '#fa7828') # Bcells, DC, Macro, other_imm, CD4, CD8
cols_segm <- c("#508791", "#ffff41",  "#D33F49")
# Colormap for annotations
# color_map = {'pre':'#FED18C', 'post':'#A7AEF2',                                      # NACT status   'pre':'gold', 'post':'royalblue',
#   'HRP':'#508791','HRD':'#D33F49',                                  # HRD status    'HRP':'darkturquoise','HRD':'tomato',
#   'short':'#D33F49','mid':'#C5DEF7','long':'#59A14F',         # survival
#   'BRCAwt':'lightblue','gBRCA1':'lightpink','sBRCA1':'darkred',          # BRCA
#   'yes':'#b7c78f',#'#98ad60',#'#8eaa44',                               # WGD status, processed info
#   'missense_variant':'cornflowerblue','stop_gained':'firebrick','frameshift_variant':'palegoldenrod','inframe_deletion':'pink', 'splice_acceptor_variant':'orange',                                    # TP53,
#   np.nan:'white',
#   'unknown':'white',
#   'no':'white',
#   'present':'darksalmon',
#   'PD':'#D33F49','SD':'#FED18C','PR':'#A7AEF2','CR':"#508791",           # primary treatment response
#   'inside_tumor': '#D33F49', 'outside_tumor':'#508791',
#   'Tumor':'#D33F49', 'Stroma':'#508791'
# }

# Define the order and color
# celltype_colors = {
#   'Stromal': '#fed18c',
#   
#   'Cancer': '#dc19dc',
#   'CD8.T.cells': '#fa7828',
#   'CD4.T.cells': '#ffff41',
#   
#   'Myeloid': '#0095B6',
#   'B.cells': '#7d28d7',
#   'Other.immune': '#508791',
#   
#   'Other': '#D9DEE4'
# }

# load metadata and filter to dcc after qc --------------------------------

geomx_dcc <- colnames(readRDS(geomx_path))
metadt <- as.data.frame(fread(metadt_path))
metadt <- metadt[metadt$dcc_filename %in% geomx_dcc, ]
rownames(metadt) <- NULL

deconv_ct_count_roi <- fread(deconv_ct_count_path_roi)
deconv_ct_count_aoi <- fread(deconv_ct_count_path_aoi)

ct_frac_all <- fread(ct_frac_all_path)

gsea_res <- fread(gsea_res_path)

# deconv ct fractions -----------------------------------------------------

# each ct freq in tumor/stroma segment across nact status + freq cluster
aoi_ctfrac_long <- deconv_ct_count_aoi %>%
  left_join(metadt[, c('dcc_filename', 'Segment', 'NACT_status', clust_type)]) %>%
  mutate(segment_nact = paste0(Segment, '_', NACT_status)) %>%
  mutate(segment_nact_order = factor(segment_nact, 
                                     levels=c('stroma_pre', 'stroma_post', 'tumor_pre', 'tumor_post'), 
                                     ordered = T)) %>%
  filter(cell_type %in% !!ct_names_immune)


# boxpl all cells at once, color by segment_nact
ggplot(aoi_ctfrac_long, aes(x = factor(cell_type), y = ct_frac_sd, fill = factor(segment_nact_order))) +
  geom_boxplot() +
  labs(title = 'immune cell type frequencies aross segment and NACT status', 
       x = "cell type", y = "cell type fraction in AOI", fill = "segment and NACT status") +
  scale_fill_manual(values = cols_segment_nact) +
  geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = TRUE, size = 0.2, label.size = 2.8) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
  

ggsave(file.path(out_dir, paste0('ct_frac_aoi_segment_nact.pdf')), height = 8, width = 10, unit = 'in')
ggsave(file.path(out_dir, paste0('ct_frac_aoi_segment_nact.svg')), height = 8, width = 10, unit = 'in')


# roi celltype clustering -------------------------------------------------
# TODO find + check calcuations of immunefractions (from geomx_deconvolution_count_unify)
# get immunefraction per roi, transform to long

immunefrac_cluster_roi <- metadt %>%
  select(sample_roi, !!clust_type, starts_with('ct_immunefrac_sd_roi')) %>%
  pivot_longer(cols = starts_with('ct_immunefrac_sd_roi'), names_to = 'cell_type') %>%
  mutate(cell_type = gsub('ct_immunefrac_sd_roi_', '', cell_type)) %>%
  distinct()

# stacked bar plot with all samples per cluster
ggplot(immunefrac_cluster_roi, aes(x = sample_roi, y = value, fill = cell_type)) +
  geom_bar(stat = "identity") +
  labs(title = paste0("immune cell type fractions across clustered ROIs"), x = "ROIs", y = "immune cell type fraction") +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  facet_wrap(~ get(clust_type), scales = "free", ncol = 2) +
  scale_fill_manual(values = cols_cell_types_imm)

ggsave(file.path(out_dir, paste0('barplot_clust_immunefrac_', clust_type,  '.pdf')),height = 6, width = 8, unit = 'in')
ggsave(file.path(out_dir, paste0('barplot_clust_immunefrac_', clust_type,  '.svg')),height = 6, width = 8, unit = 'in')


ct_frac_clust_mean <- immunefrac_cluster_roi %>%
  dplyr::group_by_at(c(clust_type, 'cell_type')) %>%
  dplyr::summarise(mean_ct_frac_of_immune = mean(value))

# stacked barplot for mean ct fraction per cluster
ggplot(ct_frac_clust_mean, aes(x = get(clust_type), y = mean_ct_frac_of_immune, fill = cell_type)) +
  geom_bar(stat = "identity") +
  labs(title = paste0("mean cell type fraction of immune per ROI cluster"), x = "ROI cluster", y = "mean cell type fraction") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6)) +
  scale_fill_manual(values = cols_cell_types_imm)

ggsave(file.path(out_dir, paste0('barplot_mean_clust_immunefrac_', clust_type,  '.pdf')), height = 6, width = 8, unit = 'in')
ggsave(file.path(out_dir, paste0('barplot_mean_clust_immunefrac_', clust_type,  '.svg')), height = 6, width = 8, unit = 'in')


# ROI clusters distributions ----------------------------------------------

# nr of ROI clusters in different sample groups
roi_clustnr <- metadt %>%
  select(sample_roi, !!clust_type, NACT_status, HRP_status) %>%
  distinct() %>%
  mutate(nact_order = factor(NACT_status, levels=c('pre', 'post'), ordered = T))
  
# barplots with nr of clusters
ggplot(roi_clustnr, aes(x = get(clust_type), fill = nact_order)) +
  geom_bar(position = position_dodge()) +
  labs(title = paste0("numbers of ROI clusters across ", 'NACT status'), x = "ROI cluster", fill= "NACT status") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 8)) +
  scale_fill_manual(values = cols_nact)

ggsave(file.path(out_dir, paste0('roi_nr_cluster_', clust_type_name, '_boxpl_color_', 'NACT_status', '.pdf')))
ggsave(file.path(out_dir, paste0('roi_nr_cluster_', clust_type_name, '_boxpl_color_', 'NACT_status', '.svg')))

# TODO old plots taken from eyemt_downstream_analysis_gsea.R old version from github

# ct frac clusters distribution across samples (pre, post, HRD, PFS, OS)
vars_labels <- c('NACT_status', 'HRP_status', 'BRCA_status') #, 'HRP_status', 'primary_treatment_response'
vars_cont <- c('PFS_days', 'OS_days') #c('TMB', 'ovaHRDscar_score', 'PFS_days', 'OS_days')

# count ROI label frequency per sample (per ROI not AOIs!)
cluster_freqs_per_sample <- metadt[, c('Sample', 'sample_roi', clust_type)] %>%
  distinct() %>%
  group_by(Sample, get(clust_type)) %>%
  mutate(clust_nr_per_sample = n()) %>%
  ungroup() %>%
  group_by(Sample) %>%
  mutate(clust_freq_per_sample = clust_nr_per_sample/n()) %>%
  select(Sample, !!clust_type, clust_nr_per_sample, clust_freq_per_sample) %>%
  distinct() %>%
  full_join(clust_allcombs, by = c('Sample', clust_type)) %>% # join with all combs to get 0
  replace(is.na(.), 0) %>%
  left_join(distinct(metadt[, c('Sample', vars_labels, vars_cont)])) %>%
  mutate(nact_order = factor(NACT_status, levels=c('pre', 'post'), ordered = T))

# frequencies of ROIs clusters - boxplots
ggplot(cluster_freqs_per_sample, aes(x = get(clust_type), y = clust_freq_per_sample, fill = nact_order)) +
  geom_boxplot() +
  geom_point(position= position_jitterdodge(dodge.width = 1, jitter.width= .3, jitter.height = 0),
             size= 0.5, alpha = 0.6) +
  geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = TRUE, size = 0.2, label.size = 2.8) +
  labs(title = paste0("frequencies of ROI clusters per sample across ", 'NACT_status'), x = "ROI cluster", y = "ROI cluster frequency", fill= paste0(label_var)) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 8)) +
  ylim(0, 1.1) +
  scale_fill_manual(values = cols_nact)

ggsave(file.path(out_dir, paste0('roi_freq_cluster_', clust_type_name, '_boxpl_color_', 'NACT_status', '.pdf')))
ggsave(file.path(out_dir, paste0('roi_freq_cluster_', clust_type_name, '_boxpl_color_', 'NACT_status', '.svg')))

# correlation between Macrophages_domin ROI fraction and PFS
# cont_var <- 'PFS_days'
# 
# cluster2 <- cluster_freqs_per_sample %>%
#   filter(NACT_status == 'post' & roi_cluster_label_gmm == 'Macro_domin' & HRP_status == 'HRP' & BRCA_status == 'BRCAwt') #
# cor(cluster2$clust_freq_per_sample, cluster2$OS_days)
# 
# # scatterplots for continuous vars
# for(cont_var in vars_cont){
# 
#   ggplot(cluster2, aes(x = clust_freq_per_sample, y = get(cont_var), color = HRP_status, shape = BRCA_status)) +
#     geom_point(size = 3) +
#     xlab("ROI cluster frequency in sample") +
#     ylab(cont_var) +
#     scale_color_discrete(name = 'HRD status') +
#     ggtitle(paste0('all corr ', as.character(cor(cluster2$clust_freq_per_sample, cluster2$PFS_days)))) +
#     theme_bw()
# 
#   ggsave(file.path(out_dir, paste0('roi_freq_cluster_', clust_type_name, '_scatter_color_', cont_var, '.png')))
# }

# cell type fractions and roi labels comparison ---------------------------

# bp vs sd per aoi
all_scatter_aoi <- ggplot(data = deconv_ct_count_aoi, aes(x = ct_frac_bp, y = ct_frac_sd)) +
  geom_point(aes(color = Segment_geomx)) +
  facet_wrap(~ cell_type) +
  geom_smooth(method='lm', formula= y~x) +
  stat_correlation(method = 'pearson', output.type = 'text') +
  scale_color_manual(values = cols_segm) +
  labs(title = paste('cell type fraction bp vs sd'), x = paste0('cell type fraction bp'), y = paste0('cell type fraction sd'))

ggsave(file.path(out_dir, paste0('scatter_all_bp_vs_sd_per_aoi.pdf')), width = 8, height = 8, unit = 'in')
ggsave(file.path(out_dir, paste0('scatter_all_bp_vs_sd_per_aoi.svg')), width = 8, height = 8, unit = 'in')

# from geomx_roi_hubs_integration
ct_frac_all[ct_frac_all==""]<- "nolabel"
ct_frac_all[is.na(ct_frac_all)]<- "nolabel"

# scatterplot with geomx vs cycif ct fractions
comp_type <- 'ct_frac'

for(value_comb in c('bp_sd', 'bp_cycif', 'sd_cycif')){
  vals <- unlist(strsplit(value_comb, split = '_'))
  print(value_comb)

  # for all faceted by ct
  all_scatter <- ggplot(data = ct_frac_all, aes(x = get(paste0(comp_type, '_', vals[1])), y = get(paste0(comp_type, '_', vals[2])))) +
    geom_point(aes(color = Segment_geomx)) +
    facet_wrap(~ cell_type) +
    geom_smooth(method='lm', formula= y~x) +
    stat_correlation(method = 'pearson', output.type = 'text') +
    scale_color_manual(values = cols_segm) +
    labs(title = paste('cell type fraction', vals[1], 'vs', vals[2]), x = paste0('cell type fraction', '_', vals[1]), y = paste0('cell type fraction', '_', vals[2]))
  
  ggsave(file.path(out_dir, paste0('scatter_all_', vals[1], '_', vals[2], '_', comp_type, '_per_roi.pdf')),
         width = 8, height = 8, unit = 'in')
  ggsave(file.path(out_dir, paste0('scatter_all_', vals[1], '_', vals[2], '_', comp_type, '_per_roi.svg')),
         width = 8, height = 8, unit = 'in')
}

# comparison between geomx and cycif roi labels
# TODO rerun basic code with new clustering
lab_name1 <- 'community_cluster_label_manualnames_freq0.05'
lab_name2 <- 'roi_cluster_label_gmm'

labs_all <- ct_frac_all %>%
  select(sample_roi, lab_name1, lab_name2) %>%
  distinct()

labs_all[labs_all==""]<- "nolabel"
labs_all[is.na(labs_all)]<- "nolabel"

# compute cross-frequencies of different labels
labs_cross <- table(labs_all[[lab_name1]], labs_all[[lab_name2]])
labs_cross <- matrix(labs_cross, ncol=ncol(labs_cross), dimnames=dimnames(labs_cross))
labs_cross <- labs_cross[, order(colnames(labs_cross))]
labs_cross <- labs_cross[order(rownames(labs_cross)), ]

# do heatmap
col_fun = colorRamp2(c(0, max(labs_cross)), c("white", "red"))

ht <- Heatmap(labs_cross, col = col_fun, show_heatmap_legend = TRUE,
              cluster_rows = FALSE, cluster_columns = FALSE, row_title = lab_name1, column_title = lab_name2)

pdf(file = file.path(out_dir, paste0('heatmap_labels_', lab_name1, '_', lab_name2, '.pdf')), width=10, height=8)
draw(ht, heatmap_legend_side="bottom")
dev.off()


# ssGSEA correlation with ct ----------------------------------------------

gsea_res_macro <- gsea_res[gsea_res$expr_signal == 'deconv_Macrophages_Monocytes', ]
#unique(gsea_res_macro$pathway)

myelonets_paths <- c('REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION',
                     'REACTOME_COMPLEMENT_CASCADE',
                     'REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL',
                     'REACTOME_INTERLEUKIN_10_SIGNALING',
                     'REACTOME_INTERLEUKIN_1_SIGNALING',
                     'REACTOME_INTERLEUKIN_2_FAMILY_SIGNALING',
                     'REACTOME_INTERLEUKIN_2_SIGNALING',
                     'REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING',
                     'REACTOME_INTERLEUKIN_6_FAMILY_SIGNALING',
                     'REACTOME_INTERLEUKIN_6_SIGNALING',
                     'REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION',
                     'REACTOME_ROS_AND_RNS_PRODUCTION_IN_PHAGOCYTES',
                     'REACTOME_SIGNALING_BY_CSF1_M_CSF_IN_MYELOID_CELLS',
                     'REACTOME_TNF_SIGNALING',
                     'REACTOME_TOLL_LIKE_RECEPTOR_CASCADES',
                     'BIOCARTA_IL2_PATHWAY',
                     'BIOCARTA_TGFB_PATHWAY',
                     'BIOCARTA_TNFR1_PATHWAY',
                     'BIOCARTA_TNFR2_PATHWAY',
                     'BIOCARTA_VEGF_PATHWAY',
                     'GOBP_ACUTE_INFLAMMATORY_RESPONSE',
                     'GOBP_INTERLEUKIN_2_MEDIATED_SIGNALING_PATHWAY',
                     'GOBP_INTERLEUKIN_2_PRODUCTION',
                     'GOBP_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY',
                     'ADDITIONAL_WANG_M1',
                     'ADDITIONAL_WANG_M2',
                     'REACTOME_NR1H3_NR1H2_REGULATE_GENE_EXPRESSION_LINKED_TO_CHOLESTEROL_TRANSPORT_AND_EFFLUX',
                     'REACTOME_LYSOSPHINGOLIPID_AND_LPA_RECEPTORS',
                     'GOBP_CHOLESTEROL_EFFLUX',
                     'GOBP_CHOLESTEROL_IMPORT',
                     'GOBP_CHOLESTEROL_STORAGE',
                     'GOBP_FOAM_CELL_DIFFERENTIATION',
                     'REACTOME_SIGNALING_BY_PTK6',
                     'REACTOME_RUNX3_REGULATES_IMMUNE_RESPONSE_AND_CELL_MIGRATION',
                     'REACTOME_SCAVENGING_BY_CLASS_A_RECEPTORS',
                     'REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY',
                     'REACTOME_TNFR1_INDUCED_PROAPOPTOTIC_SIGNALING',
                     'REACTOME_TNFR1_INDUCED_NF_KAPPA_B_SIGNALING_PATHWAY'
                     )

myelonets_paths_less <- c('REACTOME_INTERLEUKIN_10_SIGNALING',
                     'REACTOME_INTERLEUKIN_1_SIGNALING',
                     'REACTOME_INTERLEUKIN_2_SIGNALING',
                     'REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION',
                     'REACTOME_ROS_AND_RNS_PRODUCTION_IN_PHAGOCYTES',
                     'REACTOME_SIGNALING_BY_CSF1_M_CSF_IN_MYELOID_CELLS',
                     'REACTOME_TNF_SIGNALING',
                     'BIOCARTA_TGFB_PATHWAY',
                     'BIOCARTA_VEGF_PATHWAY',
                     'GOBP_ACUTE_INFLAMMATORY_RESPONSE',
                     'ADDITIONAL_WANG_M2',
                     'GOBP_CHOLESTEROL_EFFLUX',
                     'GOBP_CHOLESTEROL_IMPORT',
                     'GOBP_CHOLESTEROL_STORAGE',
                     'GOBP_FOAM_CELL_DIFFERENTIATION',                     
                     'REACTOME_RUNX3_REGULATES_IMMUNE_RESPONSE_AND_CELL_MIGRATION',
                     'REACTOME_SCAVENGING_BY_CLASS_A_RECEPTORS',
                     'REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY',
                     'REACTOME_TNFR1_INDUCED_PROAPOPTOTIC_SIGNALING',
                     'REACTOME_TNFR1_INDUCED_NF_KAPPA_B_SIGNALING_PATHWAY'
                     
)

outname <- 'stroma_post_myelonets' # depend on filtering

metadt_sel <- metadt %>%
  filter(roi_cluster_label_gmm == 'Macro_domin', NACT_status == 'post' & Segment == 'stroma') %>% # roi_cluster_label_gmm == 'Macro_domin' & 
  #filter(NACT_status == 'post' & Segment == 'stroma') %>%
  #filter(roi_cluster_label_gmm != 'Macro_domin', Segment == 'stroma', NACT_status == 'post') %>%
  mutate(PFS_group = ifelse(PFS_days <= 350, 'short', ifelse(PFS_days >= 602, 'long', 'mid')))

dcc_annots <- c('Segment_geomx', 'HRP_status', 'BRCA_status', 'primary_treatment_response', 'PFS_group', 'roi_cluster_label_gmm')

# row annotations based on metadata
dcc_annot <- metadt_sel %>%
  select('dcc_filename', !!dcc_annots) %>%
  column_to_rownames(var="dcc_filename") %>%
  mutate_all(as.factor)

# select dcc and pathways
gsea_sel <- gsea_res_macro %>%
  select(dcc_filename, pathway, ssgsea_score) %>%
  filter(dcc_filename %in% metadt_sel$dcc_filename & pathway %in% myelonets_paths_less) %>% #
  distinct()

# transform to wide
gsea_sel_wide <- gsea_sel %>%
  pivot_wider(id_cols = 'dcc_filename', names_from = 'pathway', values_from = 'ssgsea_score') %>%
  column_to_rownames(var="dcc_filename")

# make heatmap with annotations
# filter annots and ensure ordering
dcc_annot_toplot <- dcc_annot[match(rownames(gsea_sel_wide), rownames(dcc_annot)),] 


identical(rownames(gsea_sel_wide), rownames(dcc_annot_toplot))

# for(path in colnames(gsea_sel_wide)){
#   corval <- cor(gsea_sel_wide[[path]], as.numeric(dcc_annot_toplot$PFS_days))
#   if(corval >= 0.3 | corval <= -0.3){
#     print(path)
#     print(corval)
#   }
# }


# set up annotations
dcc_ha = HeatmapAnnotation(df = dcc_annot_toplot, which = 'row', na_col = "white")

# define colors
scalemin <- ifelse(min(gsea_sel_wide) < 0, min(gsea_sel_wide), -0.01) # to keep the color scale b-r
col_fun <- colorRamp2(c(scalemin, 0, max(gsea_sel_wide)), c("blue", "white", "red"))

# do the hmap
pdf(file=file.path(out_dir, paste0('hmap_scores_', outname, '.pdf')), width=8, height=8)

gsea_heat <- Heatmap(gsea_sel_wide, name = "ssGSEA score", 
                     left_annotation = dcc_ha, 
                     col = col_fun,
                     heatmap_legend_param = list(
                       legend_direction = "vertical", 
                       legend_width = unit(1, "in")),
                     show_row_names = F,
                     column_names_gp = gpar(fontsize = 6),
                     column_names_max_height = unit(12, "in")
)

draw(gsea_heat, 
     heatmap_legend_side="right", 
     annotation_legend_side="right",
     merge_legend = TRUE)
dev.off()

#######################################
#######################################
# correlation between pathways
gsea_sel_wide_corr <- as.matrix(cor(gsea_sel_wide, method = 'pearson'))
corr_thr <- 0.5

if(!is.null(corr_thr)){
  gsea_sel_wide_corr[gsea_sel_wide_corr >= -corr_thr & gsea_sel_wide_corr <= corr_thr] <- 0
  cor_abovethr <- (colSums(gsea_sel_wide_corr, na.rm=T) != 1) # cor btw same pathway will be 1
  gsea_sel_wide_corr <- gsea_sel_wide_corr[cor_abovethr, cor_abovethr]
}


col_fun_corr <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# do the heatmap
pdf(file=file.path(out_dir, paste0('hmap_scores_corr_',outname, '.pdf')), width=11, height=7)

gsea_heat <- Heatmap(gsea_sel_wide_corr, name = "ssGSEA scores corr", 
                     col = col_fun_corr,
                     row_names_gp = gpar(fontsize = 6),
                     column_names_gp = gpar(fontsize = 6),
                     heatmap_legend_param = list(
                       legend_direction = "horizontal", 
                       legend_width = unit(1, "in")),
                     column_names_max_height = unit(8, "in"),
                     row_names_max_width = unit(8, "in")
)

draw(gsea_heat, heatmap_legend_side="top", )
dev.off()

########################
########################
# UMAP based on myelonets pathways
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42

# make umap
umap_out <- umap::umap(gsea_sel_wide, config = custom_umap)
umap_out_res <- umap_out$layout[, c(1,2)]
colnames(umap_out_res) <- c('UMAP_1', 'UMAP_2')
# merge with metadata
identical(rownames(umap_out_res), rownames(dcc_annot_toplot))

# make PCA
pca_obj <- pca(t(gsea_sel_wide), scale = F)
pca_res <- t(-1*pca_obj$rotated) # reverse the signs of eigen vectors
pca_loads <- -1*pca_obj$loadings
#pca_vars <- pca_obj$variance
pca_resdf <- t(pca_res)[, c(1,2)]


umap_annot <- cbind(dcc_annot_toplot, umap_out_res, pca_resdf)

# make plot
ggplot(umap_annot,
       aes(x = PC1, y = PC2, 
           #color = as.numeric(PFS_days), 
           #color = PFS_group,
           #color = HRP_status, 
           color = roi_cluster_label_gmm
           )) +
  geom_point(size = 3) +
  xlab(paste0('PC1 (',as.character(round(pca_obj$variance[1], 1)), '% of variance)')) +
  ylab(paste0('PC2 (',as.character(round(pca_obj$variance[2], 1)), '% of variance)')) +
  scale_color_discrete(name = 'ROI cluster') + 
  #scale_shape_discrete(name = 'ROI cluster type') + 
  theme_bw()

ggsave(file.path(out_dir, paste0('pca_gsea_myelonets_paths_',outname, '.pdf')), width = 5, height = 5, device=pdf)

#######################################################3
########################################################
# correlation between pathways and immune fractions
# wide df with immunefractions
# roi_immunefrac <- metadt %>%
#   filter(dcc_filename %in% rownames(gsea_sel_wide)) %>%
#   dplyr::select('dcc_filename', starts_with('ct_immunefrac_sd_roi')) %>%
#   column_to_rownames('dcc_filename')
corr_thr_ctfrac <- 0.3

metadt_sel <- metadt# %>%
  #filter(roi_cluster_label_gmm == 'Macro_domin' & NACT_status == 'post' & Segment == 'stroma') # 

# transform to wide
gsea_sel_wide_all <- gsea_res_macro %>%
  filter(dcc_filename %in% metadt_sel$dcc_filename) %>%
  dplyr::select(dcc_filename, pathway, ssgsea_score) %>%
  distinct() %>%
  pivot_wider(id_cols = 'dcc_filename', names_from = 'pathway', values_from = 'ssgsea_score') %>%
  column_to_rownames(var="dcc_filename")

aoi_ctfrac <- metadt %>%
  filter(dcc_filename %in% rownames(gsea_sel_wide_all)) %>%
  dplyr::select('dcc_filename', starts_with('ct_frac_sd_aoi')) %>%
  column_to_rownames('dcc_filename') %>%
  select(ends_with(ct_all))

colnames(aoi_ctfrac) <- gsub('ct_frac_sd_aoi_', '', colnames(aoi_ctfrac))

# correlation between pathways and immune ct fractions
gsea_seg_corr_aoi_ctfrac <- cor(gsea_sel_wide_all, aoi_ctfrac, method = 'spearman')
gsea_seg_corr_aoi_ctfrac[is.na(gsea_seg_corr_aoi_ctfrac)] <- 0

gsea_seg_corr_aoi_ctfrac[gsea_seg_corr_aoi_ctfrac >= -corr_thr_ctfrac & gsea_seg_corr_aoi_ctfrac <= corr_thr_ctfrac] <- 0
cor_abovethr_aoi_ctfrac_row <- (rowSums(gsea_seg_corr_aoi_ctfrac, na.rm=T) != 0) 
cor_abovethr_aoi_ctfrac_col <- (colSums(gsea_seg_corr_aoi_ctfrac, na.rm=T) != 0) 
gsea_seg_corr_aoi_ctfrac <- as.matrix(gsea_seg_corr_aoi_ctfrac[cor_abovethr_aoi_ctfrac_row, cor_abovethr_aoi_ctfrac_col])

pdf(file=file.path(out_dir, paste0('hmap_scores_sdaoictfrac_corr_all.pdf')), width=11, height=7)

gsea_heat <- Heatmap(gsea_seg_corr_aoi_ctfrac, name = "correlation between ssGSEA scores", 
                     col = col_fun_corr,
                     row_names_gp = gpar(fontsize = 8),
                     column_names_gp = gpar(fontsize = 8)
)

draw(gsea_heat, heatmap_legend_side="bottom")
dev.off()


#############