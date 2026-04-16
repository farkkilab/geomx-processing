library(plyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(GeomxTools)
library(biomaRt)
library(ggpubr)
library(viridis)
library(ggpmisc)
library(reshape2)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(PCAtools)
library(umap)
library(Rtsne)
library(tidyr)

# PROBLEM: BP predicts a lot more tumor in stromal AOIs, 
# while sd predicts a lot more Tcells 
# maybe binning all Tcells together created highly variable cluster which confuses the algorithm

# careful for genes with variance 0 removed from deconvolution

# TODO split into 2 scripts - for sparsity, gdr etc and deconvolution
# TODO gsea on canonical_markers (for scRNAseq cell typing) for deconv assesment

batch <- 'batch123'

###########################

proj_dir <<- '~/Documents/phd/st'

if(batch == 'batch1'){
  output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch1-1903') # batch1
} else if(batch == 'batch2'){
  output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch2-1903') # batch2
} else if(batch == 'batch3'){
  output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch3-2606') # batch3
} else if(batch == 'batch12'){
  # output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch12-1004') # batch12
  output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch12-1205-no-counts-shift2') # batch12
} else if(batch == 'batch23'){
  output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch23-2706') # batch23
} else if(batch == 'batch123'){
  output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808') # batch123
}else{
  stop('wrong batch nr')
}


#clinical_dt_path <- file.path(proj_dir, 'geomx-processing/data/b12_dcc_clinical_data.csv')
#TODO adjust script to this - clinical should be merged with dcc
clinical_dt_path <<- file.path(proj_dir, 'data/geomx/clinical_data/9_eyemt_patient_clinical_data_SENSITIVE_upd_0426.csv')

geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS')
#geomx_norm_batch_eff_rm_path2 <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm_deseq2.RDS') # without vst
metadata_orig_path <- file.path(proj_dir, 'data/geomx/batch123/metadata/dcc_metadata_batch123_no_tls_cleaned.csv')


#################
# files and params for deconvolution

#used for deconv
scrna_anno <- 'mid_lvl_ct_updated' # 'mid_lvl_ct' / 'mid_lvl_ct_updated' / 'low_lvl_ct'


if(scrna_anno == 'mid_lvl_ct'){
  ct_names <- c('Bcells', 'DCs', 'Endothelial cells', 'Fibroblasts', 'Macrophages', 'NKcells', 'Tcells', 'tumor')
  cells_immune <- c('Bcells', 'DCs', 'Macrophages', 'NKcells', 'Tcells')
} else if(scrna_anno == 'mid_lvl_ct_updated'){
  # ct_names <- c("Tcells_reg","Tcells_CD8","Tcells_CD4", "Tcells_other", "Bcells", "NKcells",
  #               "Macrophages_Monocytes", "DCs", "Mast_cells", "Fibroblasts", "Endothelial_cells", "tumor")
  # cells_immune <- c("Tcells_reg","Tcells_CD8","Tcells_CD4", "Tcells_other", "Bcells", 
  #                   "Macrophages", "Monocytes", "DCs", "Mast_cells")
  ct_names <- c("Tcells_other","Tcells_CD8","Tcells_CD4", "Bcells", 'NKcells', 'Mast_cells',
                "Macrophages_Monocytes", "DCs", "Fibroblasts_Mesothelial", "Endothelial_cells", "tumor")
  cells_immune <- c("Tcells_other","Tcells_CD8","Tcells_CD4", "Bcells", 'NKcells', 'Mast_cells',
                "Macrophages_Monocytes", "DCs")
  cells_stroma <- c("Fibroblasts_Mesothelial", "Endothelial_cells")
  cells_with_ct_sign <- c("tumor", "Fibroblasts_Mesothelial", "Endothelial_cells", 
                          "Tcells_CD4","Tcells_CD8", "Tcells_other", "NKcells",
                          "Bcells", "Macrophages_Monocytes")
  
  # ct_names <- c("tumor", "Fibroblasts", "Mesothelial", "Endothelial", 
  #               "Tcells_CD4","Tcells_CD8", "Tcells_Treg","Tcells_other", "NK", "ILC",
  #               "Plasma_cells", "B_cells", "pDC",
  #               "DC",  "Macrophages",  "Mast_cells")
  # cells_immune <- c("Tcells_CD4","Tcells_CD8", "Tcells_Treg","Tcells_other", "NK", "ILC",
  #                   "Plasma_cells", "B_cells", "DC",  "Macrophages",  "Mast_cells")
  # cells_stroma <- c("Fibroblasts", "Mesothelial", "Endothelial")
  # 
  # cells_with_ct_sign <- c("tumor", "Fibroblasts", "Endothelial", 
  #                         "Tcells_CD4","Tcells_CD8", "Tcells_Treg","Tcells_other", "NK",
  #                         "B_cells", "Macrophages") #TODO add for DC
  
} else if(scrna_anno == 'low_lvl_ct'){
  ct_names <- c("Tcells_NK", "Bcells", "Myeloids","Mast_cells", "Fibroblasts_Endothelial", "tumor")
  cells_immune <- c("Tcells_NK", "Bcells", "Myeloids","Mast_cells")
  ct_to_pathway <- list(Tcells_NK = c('Tcells', 'NKcells', 'Tcells_CD8'), Bcells = 'Bcells',
                        Myeloids = 'Macrophages', Mast_cells = 'Mast_cells', 
                        Fibroblasts_Endothelial = c('Fibroblasts', 'Endothelial_cells'),
                        tumor = c('tumor', 'tumor_old'), stroma = 'stroma')
} else{
  stop('wrong annotation')
}



ct_markers_path <- file.path(proj_dir, 'geomx-processing', 'data', 'signatures', 'ct_markers.csv')

ct_gsea_all_path <- file.path(output_dir, 'pathway_analysis', 'gsea', 'ssgsea_norm_harmony_q3_norm_all_custom_ct_markers.csv.csv')


bp_cellcounts_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', paste0('bp_res_', scrna_anno, '_ct_fraction.csv'))
sd_cellcounts_path <- file.path(output_dir, 'deconvolution', 'spatial_decon', paste0('sd_res_', scrna_anno, '_geomxfiltpc_ct_fraction.csv'))


deconv_bp_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', 'bp_res_mid_lvl_ct_updated_expr_mtx_cleaned_deseq2_vst_harmony_corr.RDS')
# TODO not used atm
# deconv_raw_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', paste0('bp_res_', scrna_anno, '.RDS'))
# deconv_harmony_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', paste0('bp_res_', scrna_anno, '_expr_mtx_cleaned_vst_harmony_batch_corr.RDS'))
# deconv_limma_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', paste0('bp_res_', scrna_anno, '_expr_mtx_cleaned_vst_limma_batch_corr_main_batch_nrbatch_nr_cov_no.RDS'))


# ct fractions per aoi from geomx_roi_hubs_integration
ct_frac_deconv_aoi_path <- file.path(output_dir, 'cycif_integration', 'b123_ct_frac_deconv_bcells.csv')

# ct fractions of immune per roi (used for clustering) from geomx_relabel_roi_2nd_approach
ct_frac_deconv_roi_path <- file.path(output_dir, 'deconvolution', 'relabel-roi-deconv-dimred', 'sd_mye_lymph_b_ct_fractions_of_immune.csv')

# ROI clusters based on deconvolution ct fractions from geomx_relabel_roi_2nd_approach
ct_frac_clust_path <- file.path(output_dir, 'deconvolution', 'relabel-roi-deconv-dimred', 'sd_mye_lymph_b_all_clustering_results.csv')
clust_types <- c('clusters_gmm_clustnr_5', 'clusters_hclust_cut2')
# descriptive labels for clusters - IN THIS CASE BOTH METHODS HAS THE SAME CLUSTERS DESCRIPTION
clust_labels <- list(CD8_Macro_domin = 1, mixed_w_CD4 = 2, mixed_w_others = 3, Macro_domin = 4, Bcell_domin = 5)

# set up metadata variables names -----------------------------------------

#aoi_id <<- 'dcc_filename'
aoi_id <<- 'dcc_filename'
roi_id <<- 'Roi'

main_batch_var <- 'main_batch_nr'
batch_var <<- 'batch_nr'


aoi_segment_var <<- "Segment"
main_roi_label <<- "Annotation_cell" 
main_experimental_condition <<- 'NACT_status'
sample_name <<- 'Sample'

other_vars_bio <<- c("Segment_geomx", "Patient", "Site") # 'PFS_months', 'PFS'
other_vars_tech <<- c('Slide_Name')

primary_batch_var <<- ifelse(batch %in% c('batch1', 'batch2', 'batch3'), batch_var, main_batch_var)
if(batch %in% c('batch1', 'batch2', 'batch3')){secondary_batch_var <<- NULL} else{secondary_batch_var <<- batch_var}

important_metadt <- unique(c(primary_batch_var, secondary_batch_var,  
                aoi_segment_var, sample_name, main_experimental_condition, 
                other_vars_bio, 'Segment_geomx'))

important_clindt_labs <- c('stage', 'primary_surgery_residual', 'HRP_status', 'BRCA_status')
important_clindt_cont <- c('PFS_days', 'OS_days', 'TMB', 'age_at_diagnosis')

# load util functions and create dirs -------------------------------------

source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))

dir.create(file.path(output_dir, 'sanity_check'), recursive = T, showWarnings = F)
dir.create(file.path(output_dir, 'sanity_check', 'deconv'), recursive = T, showWarnings = F)

#################################################################
################################################################
# basic sanity check, gene coverage etc

geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)

sort(colnames(sData(geomx_obj)))
dim(geomx_obj@assayData$exprs)

expr_list <- list(raw = geomx_obj@assayData$exprs, deseq_norm = geomx_obj@assayData$deseq2_norm,
                  q3_norm = geomx_obj@assayData$q3_norm, vst = geomx_obj@assayData$deseq2_vst_norm, 
                  limma_deseq2_vst_batch_corr = geomx_obj@assayData$limma_deseq2_vst_norm,
                  harmony_deseq2_vst_batch_corr = geomx_obj@assayData$harmony_deseq2_vst_norm,
                  limma_q3_norm_batch_corr = geomx_obj@assayData$limma_q3_norm,
                  harmony_q3_norm_batch_corr = geomx_obj@assayData$harmony_q3_norm)

# proportion of 0 reads ---------------------------------------------------

sapply(names(expr_list), function(x){
  expr <- expr_list[[x]]
  
  # limma and harmony batch effect corrected data are in a log scale so 0 --> 1
  #null_nr <- ifelse(x %in% c('limma_batch_corr', 'harmony_batch_corr'), 1, 0)
  
  prop_0 <- apply(expr, 2, function(x){
    length(which(x == 0))/length(x)
  })
  
  prop_0 <- cbind(sData(geomx_obj)[, c('dcc_filename', 'Segment', 'main_batch_nr')], prop_0)
  
  ggplot(data = prop_0) +
    geom_histogram(aes(x= prop_0, fill = Segment), alpha = 0.5, bins = 200) +
    #facet_wrap(~Segment, scales = "fixed", dir="v") +
    ylab('AOI counts') +
    ggtitle(paste('proportions of 0 reads in ', x, ' expr mtx'), 
            subtitle = paste0('prop of 0 in whole expr mtx: ', round(length(which(expr == 0))/length(expr), 5)))
  
  ggsave(file.path(output_dir, 'sanity_check', paste0('prop_of_0_', x, '.png')), 
         width = 2000, height = 1000, unit = 'px')
})

# nr of reads per nuclei --------------------------------------------------

reads_sum <- lapply(expr_list, colSums)
reads_sum <- as.data.frame(do.call(cbind, reads_sum))
reads_sum <- cbind(sData(geomx_obj)[, c('dcc_filename', 'nuclei', 'Segment', 'main_batch_nr')], reads_sum)
reads_sum$raw <- round(reads_sum$raw) # reverts weird behaviour of colsums to return floats for sum of integers

for(coln in names(expr_list)){
  ggplot(data = reads_sum) +
    geom_point(aes(x= nuclei, y = get(coln), shape = Segment, color = as.factor(main_batch_nr))) +
    ylab(coln) +
    ggtitle(paste(coln, 'reads nr per nuclei nr'))
  
  ggsave(file.path(output_dir, 'sanity_check', paste0('reads_nr_per_nuclei_', coln, '.png')), width = 2000, height = 1000, unit = 'px')
  
}


# gene coverage -----------------------------------------------------------

# what is the % of samples covered by each gene?

# for all data
fdt <- fData(geomx_obj)

# plot gene detection rate per gene
ggplot(data = fdt) +
  geom_histogram(aes(x = fdt$DetectionRate), bins = 100) +
  xlim(0, 1) +
  ggtitle('gene detection rate per gene')

ggsave(file.path(output_dir, 'sanity_check', paste0('gene_detection_rate_per_gene.png')))

# what is the gene coverage in each AOI?
# plot gene detection rate per AOI
ggplot(data = pData(geomx_obj)) +
  geom_histogram(aes(x = GeneDetectionRate, fill = Segment), bins = 100) +
  xlim(0, 1) +
  ggtitle('gene detection rate per AOI')

ggsave(file.path(output_dir, 'sanity_check', paste0('gene_detection_rate_per_aoi.png')))


# stromal/tumor markers ---------------------------------------------------
# check tum/str markers activity in tum/str AOIs

ct_markers <- fread(ct_markers_path)
all_markers <-   as.list(ct_markers)
all_markers <- lapply(all_markers, function(l){l[l !=""]})

# not in exprs:
# CCL4 IGLC2 IGLC3 KLRC1 TPSB2 XCL2
# SMA ? ACTA2

ct_gsea_all <- fread(ct_gsea_all_path)
ct_gsea_all$pathway <- gsub(' ', '_', ct_gsea_all$pathway)
ct_gsea_all$pathway <- ifelse(ct_gsea_all$pathway == 'Nkcells', 'NKcells', ct_gsea_all$pathway)
ct_gsea_all$pathway <- ifelse(ct_gsea_all$pathway == 'CD8_Tcells', 'Tcells_CD8', ct_gsea_all$pathway)
ct_gsea_all$pathway <- ifelse(ct_gsea_all$pathway == 'Mast cells', 'Mast_cells', ct_gsea_all$pathway)
ct_gsea_all$pathway <- ifelse(ct_gsea_all$pathway == 'Dcs', 'DC', ct_gsea_all$pathway)


# tum/stromal markers in tum/stromal AOIs
ct_boxpl <- pathway_boxplot(ct_gsea_all, 'pathway', 'ssgsea_score', aoi_segment_var, NULL, 
                            'cell type markers ssgsea score per segment',
                            file.path(output_dir, 'sanity_check', 'ct_markers_per_segment.pdf'))

# macro + Tcell paths activity should be higher in ++ (or single+)
ct_gsea_imm <- filter(ct_gsea_all, pathway %in% c('Tcells', 'Macrophages', 'CD8_Tcells'))

ct_boxpl_anno <- pathway_boxplot(ct_gsea_imm,'pathway', 'ssgsea_score', 'Annotation_cell', facet_var = aoi_segment_var, 
                            'cell type markers ssgsea score per annotation',
                            file.path(output_dir, 'sanity_check', 'ct_markers_per_annotation.pdf'),
                            manual_colours = viridis(17))

# DCs were removed bcs only 3 genes - find sth else to check
# T-cells work as expected (especially batch1, both new labels as well)
# macrophages - not big changes in the markers pathway activity


# deconvolution sanity check ----------------------------------------------
# cell markers: http://117.50.127.228/CellMarker/CellMarkerBrowse.jsp
# human, ovary


###########################
# ct markers activity should be higher in given ct
ct_gsea_deconv <- list.files(file.path(output_dir, 'pathway_analysis', 'gsea', scrna_anno), pattern = 'ct_markers.*csv', full.names = T)
#ct_gsea_deconv <- ct_gsea_deconv[-1]

ct_gsea_deconv <- lapply(ct_gsea_deconv, fread)
ct_gsea_deconv <- do.call(rbind, ct_gsea_deconv)
ct_gsea_deconv$expr_signal <- gsub('deconv_', '', ct_gsea_deconv$expr_signal)

ct_boxpl_deconv <- pathway_boxplot(ct_gsea_deconv,'pathway', 'ssgsea_score', 'expr_signal', facet_var = aoi_segment_var, 
                                 'cell type markers ssgsea score in each deconvoluted data',
                                 file.path(output_dir, 'sanity_check','deconv', 'ct_markers_in_deconv.pdf'),
                                 manual_colours = viridis(17))

########################
# ct markers activity in full signal correlated with nr of cells sd-bp
cell_fraq <- list(bp = fread(bp_cellcounts_path), sd = fread(sd_cellcounts_path))
#colnames(cell_fraq$sd) <- gsub('.', ' ', colnames(cell_fraq$sd), fixed = T)

sapply(1:length(cell_fraq), function(x){
  print(x)
  cell_fraq_res <- as.data.frame(cell_fraq[[x]])
  cell_fraq_res$stroma <- rowSums(cell_fraq_res[, cells_stroma])
  
  ct_gsea_all_fraq <- left_join(ct_gsea_all, cell_fraq_res[, c(aoi_id, c(ct_names, 'stroma'))])
  
  for(ct_name in c(cells_with_ct_sign)){ #TODO find markers for all ct 
    print(ct_name)
    
    #ct_gsea_all_fraq_ct <- ct_gsea_all_fraq[ct_gsea_all_fraq$pathway %in% ct_to_pathway[[ct_name]], ]
    ct_gsea_all_fraq_ct <- ct_gsea_all_fraq[grepl(unlist(strsplit(ct_name, '_'))[1], ct_gsea_all_fraq$pathway), ] # find matching name in ct markers 
    
    ct_scatter <- ggplot(data = ct_gsea_all_fraq_ct, aes(x = ssgsea_score, y = get(ct_name), color = get(aoi_segment_var))) +
      geom_point(aes(shape = get(aoi_segment_var))) + 
      ggtitle(paste0(ct_name, ' ssgsea marker activity vs cell count')) + 
      xlab(paste0(ct_name, ' markers ssgsea score')) +
      ylab(paste0(ct_name, ' cell count')) +
      geom_smooth(method='lm', formula= y~exp(x)) +
      stat_poly_eq(use_label(c("R2", "p"))) +
      facet_wrap(~pathway, scales = "fixed", dir="v") 
    
    ggsave(file.path(output_dir,'sanity_check','deconv', paste0(names(cell_fraq)[x], '_', ct_name, '_ssgsea_vs_cell_count.png')),
           width = 2000, height = 2000, unit = 'px')
  }
})

####################################
# pairwise bp-sd comparison

cell_fraq_bp <- as.data.frame(cell_fraq$bp)
#cell_fraq_bp$other <- cell_fraq_bp$`Mast cells` + cell_fraq_bp$other
cell_fraq_bp <- cell_fraq_bp[, c(aoi_id, aoi_segment_var, 'Annotation_cell', ct_names)]
colnames(cell_fraq_bp) <- c(aoi_id, aoi_segment_var, 'Annotation_cell', paste0(ct_names, '_bp'))       

cell_fraq_sd <- as.data.frame(cell_fraq$sd)
#cell_fraq_sd$other <- cell_fraq_sd$`Mast cells` + cell_fraq_bp$other
cell_fraq_sd <- cell_fraq_sd[, c(aoi_id,ct_names)]
colnames(cell_fraq_sd) <- c(aoi_id, paste0(ct_names, '_sd'))   

cell_fraq_both <- left_join(cell_fraq_bp, cell_fraq_sd)

# cell_fraq_both$stroma_bp <- cell_fraq_both$Fibroblasts_bp + cell_fraq_both$Mesothelial_bp + cell_fraq_both$Endothelial_bp
# cell_fraq_both$stroma_sd <- cell_fraq_both$Fibroblasts_sd + cell_fraq_both$Mesothelial_sd + cell_fraq_both$Endothelial_sd

cell_fraq_both$stroma_bp <- cell_fraq_both$Fibroblasts_Mesothelial_bp + cell_fraq_both$Endothelial_cells_bp
cell_fraq_both$stroma_sd <- cell_fraq_both$Fibroblasts_Mesothelial_sd + cell_fraq_both$Endothelial_cells_sd

cell_fraq_both_long <- melt(cell_fraq_both, id.vars = c(aoi_id, aoi_segment_var, 'Annotation_cell'),
                            variable.name = 'cell_type', value.name = 'fraction')

cell_fraq_both_long$deconv_type <- ifelse(grepl('bp', cell_fraq_both_long$cell_type), 'bp', 'sd')
cell_fraq_both_long$cell_type <- gsub('_bp|_sd', '', cell_fraq_both_long$cell_type)


# make a paired plot
for(ct in c(ct_names)){
  
  cell_fraq_both_long_ct <- cell_fraq_both_long[cell_fraq_both_long$cell_type == ct, ]
  
  ggplot(cell_fraq_both_long_ct, aes(x = deconv_type, y = fraction)) + 
    geom_boxplot(aes(fill = deconv_type), alpha = .2) +
    geom_line(aes(group = get(aoi_id)), size = 0.2, alpha = 0.8) + 
    geom_point(size = 0.2) + 
    ggtitle(paste0(ct, ' fractions bp vs sd')) +
    facet_wrap(~ get(aoi_segment_var))
  
  ggsave(file.path(output_dir, 'sanity_check','deconv', paste0('deconv_comparison_', ct, '_bp_sd.png')),
         width = 1500, height = 1000, unit = 'px')
}

# for all
ggplot(cell_fraq_both_long, aes(x = deconv_type, y = fraction)) +
  geom_boxplot(aes(fill = deconv_type), alpha = .2) +
  geom_line(aes(group = get(aoi_id)), size = 0.2, alpha = 0.8) +
  geom_point(size = 0.2) +
  ggtitle('ct fractions bp vs sd') +
  facet_wrap(~ cell_type)

ggsave(file.path(output_dir, 'sanity_check','deconv', paste0('deconv_comparison_all_bp_sd.png')),
       width = 1500, height = 1000, unit = 'px')

#########################
# make scatters bp vs sd
for(ct_name in c(ct_names, 'stroma')){
  
  print(ct_name)
  bp_sd_scatter <- ggplot(data = cell_fraq_both, aes(x = get(paste0(ct_name, '_bp')),
                                                            y = get(paste0(ct_name, '_sd')))) +
    geom_point(aes(color = Annotation_cell)) +
    ggtitle(paste0(ct_name, ' fractions bp vs sd'),
            subtitle = paste('overall pearson cor: ', round(stats::cor(cell_fraq_both[[paste0(ct_name, '_bp')]],
                                                        cell_fraq_both[[paste0(ct_name, '_sd')]], 
                                                        use = "complete.obs"), 2))) +
    xlab('bp fraq') +
    ylab('sd fraq') +
    facet_wrap(~ Segment) +
    geom_smooth(method='lm', formula= y~x) +
    # stat_poly_eq(use_label(c("R2", "p"))) +
    stat_correlation(method = 'pearson')
  
  ggsave(file.path(output_dir,'sanity_check', 'deconv', paste0('deconv_comparison_scatter_', ct_name, '.png')),
         width = 2000, height = 2000, unit = 'px')
}

################################################
# do UMAP on deconvoluted data
# TODO move to deconv script
top_var <- NULL # it doesn't matter if we take top PCA the differnc eis non=visible
top_pca <- NULL

ct_of_interest <- c('Tcells_CD8', 'Tcells_CD4', 'Bcells', 'Macrophages_Monocytes', 'DCs')
cluster_type <- 'clusters_gmm_clustnr_5'

deconv_bp <- readRDS(deconv_bp_path)

#TODO move it somewhere - it's important!
# create metadata
aois <- sData(readRDS(geomx_norm_batch_eff_rm_path))[[aoi_id]] #QCed AOIs from dataset
clindt <- fread(clinical_dt_path, select = c('Patient', important_clindt_cont, important_clindt_labs)) %>%
  distinct()

ct_freq <- fread(ct_frac_deconv_aoi_path, select = c('dcc_filename', 'cell_type', 'ct_frac_sd')) %>%
  spread(key = 'cell_type', value = 'ct_frac_sd') %>%
  dplyr::select(dcc_filename, !!ct_of_interest)

roi_clust <- fread(ct_frac_clust_path, select = c('sample_roi', cluster_type))

metadt <- fread(metadata_orig_path, select = c('dcc_filename', 'Roi_geomx', important_metadt)) %>%
  filter(dcc_filename %in% aois) %>%
  mutate(sample_roi = paste0(Sample, '_', Roi_geomx)) %>%
  left_join(clindt) %>%
  left_join(ct_freq) %>%
  left_join(roi_clust) %>%
  as.data.frame()

metadt$roi_cluster_label <- mapvalues(metadt[[cluster_type]], 
                                      from=c(unname(unlist(clust_labels))),
                                      to=c(names(clust_labels)))

##########3
ct_name <- ct_names[1]

# iterate through all cell types
sapply(ct_names[1], function(ct_name){
  
  print(ct_name)
  # get deconv df and filter metadata
  deconv_ct <- deconv_bp[[ct_name]]
  metadt_ct <- metadt[metadt$dcc_filename %in% colnames(deconv_ct),]
  
  # iterate through all + different segments
  seg_types <- c('all', unique(metadt[, aoi_segment_var]))
  
  sapply(seg_types, function(seg){
    
    print(seg)
    dir.create(file.path(output_dir, 'sanity_check', paste0('deconv_umap_tsne_', seg)), showWarnings = T, recursive = T)
    
    if(seg != 'all'){
      deconv_seg <- deconv_ct[, metadt_ct$dcc_filename[metadt_ct[[aoi_segment_var]] == seg]]
      metadt_seg <- metadt_ct[metadt_ct$dcc_filename %in% colnames(deconv_seg),]
    } else{
      deconv_seg <- deconv_ct
      metadt_seg <- metadt_ct
    }
    
    print(dim(deconv_seg))
    
    # run UMAP and tSNE 
    ###########################
    # get top N variable genes
    if(!is.null(top_var)){
      per_gene_variance <- apply(deconv_seg, 1, stats::var)
      top_var_genes <- names(sort(per_gene_variance, decreasing = T)[1:top_var])
      
      deconv_seg <- deconv_seg[rownames(deconv_seg) %in% top_var_genes, ]
    }
    
    # do PCA
    pca_obj <- pca(deconv_seg, scale = T)
    pca_res <- t(-1*pca_obj$rotated) # reverse the signs of eigen vectors
    pca_loads <- -1*pca_obj$loadings
    #pca_vars <- pca_obj$variance
    metadt_seg[, c("PCA1_","PCA2_")] <- t(pca_res)[, c(1,2)]
    
    if(!is.null(top_pca)){
      deconv_seg <- pca_res[1:top_pca, ]
    }

    # make umap
    custom_umap <- umap::umap.defaults
    custom_umap$random_state <- 42
    umap_out <- umap(t(deconv_seg), config = custom_umap)
    metadt_seg[, c("UMAP1_","UMAP2_")] <- umap_out$layout[, c(1,2)]
    
    # make tsne
    set.seed(42) 
    tsne_out <- Rtsne(t(deconv_seg), perplexity = ncol(deconv_seg)*.15)
    metadt_seg[, c("tSNE1_","tSNE2_")] <- tsne_out$Y[, c(1,2)]
    
    for(method in c('UMAP', 'tSNE', 'PCA')){
      # for discrete labels
      for(color_var in c(important_metadt, important_clindt_labs, important_clindt_cont, ct_of_interest, 'roi_cluster_label')){
        print(color_var)
        
        sub <- ifelse(method == 'PCA', paste0('% of variance explained: PC1= ', as.character(round(pca_obj$variance[1], 2)),
                                              ' PC2= ', as.character(round(pca_obj$variance[2], 2))), '')
        
        plot_umap_tsne(metadt_seg, method_type = method, 
                       assay_name = "", color_var = color_var,
                       subtitle = sub, 
                       output_name = file.path(output_dir, 'sanity_check', paste0('deconv_umap_tsne_', seg), 
                                               paste0(ct_name, '_', method, 
                                                      '_topvargenes_', ifelse(is.null(top_var), 'NULL', as.character(top_var)),
                                                      '_toppca_', ifelse(is.null(top_pca), 'NULL', as.character(top_pca)), 
                                                      '_', color_var, '.png')),
                       output_type = 'png')
      }
    }
    # save PCA res df
    fwrite(pca_loads[, 1:10], file.path(output_dir, 'sanity_check', paste0('deconv_umap_tsne_', seg), 
                                        paste0('pca_loads_', ct_name, '.csv')), row.names = F)
  })
})


################################333
###################################
# up here works for now

#####################
# compare macro vs Tcells - scatter per annotation cell (same per sample?)
for(deconv_type in c('bp', 'sd')){
  macro_tcells_scatter <- ggplot(data = cell_fraq_both, aes(x = get(paste0('Macrophages_', deconv_type)),
                                                            y = get(paste0('Tcells_', deconv_type)), 
                                                            color = Annotation_cell)) +
    geom_point() + 
    ggtitle(paste0(deconv_type, 'Macrophages vs Tcells ct fractions')) + 
    xlab('Macrophages fraq') +
    ylab('Tcells fraq') +
    facet_wrap(~ Segment)
    # geom_smooth(method='lm', formula= y~x) +
    # stat_poly_eq(use_label(c("R2", "p")))
  
  ggsave(file.path(output_dir,'sanity_check', paste0('macro_vs_tcell_', deconv_type, '.png')),
         width = 2000, height = 2000, unit = 'px')
}


# heatmaps and stacked barplots -------------------------------------------

# TODO not really visible - after making better lables - rerun for each Annotation_cell in facet
sapply(c('bp', 'sd'), function(deconv_type){
  sapply(c('stroma', 'tumor'), function(seg){
    
    cell_fraq_both_long_deconv_seg <- cell_fraq_both_long[cell_fraq_both_long$Segment == seg &
                                                            cell_fraq_both_long$deconv_type == deconv_type, ]

    cell_fraq_both_long_deconv_seg <- arrange(cell_fraq_both_long_deconv_seg, Annotation_cell) %>%
      rowid_to_column()
    cell_fraq_both_long_deconv_seg <- cell_fraq_both_long_deconv_seg[!cell_fraq_both_long_deconv_seg$cell_type == 'stroma', ]
    
    anno_color_mapping <- setNames(viridis(length(unique(cell_fraq_both_long_deconv_seg$Annotation_cell))),
                                   unique(cell_fraq_both_long_deconv_seg$Annotation_cell))
    
    ggplot(data = cell_fraq_both_long_deconv_seg, aes(x = reorder(dcc_filename, rowid), 
                                                      y = fraction, fill = cell_type)) +
      geom_bar(position="fill", stat="identity") +
      ggtitle(seg) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 3,
                                       color = anno_color_mapping[cell_fraq_both_long_deconv_seg$Annotation_cell])) +
      scale_x_discrete(labels = cell_fraq_both_long_deconv_seg$Annotation_cell) 
    
    ggsave(file.path(output_dir, 'sanity_check', paste0('barplot_',deconv_type, '_', seg, '.png')),
           width = 1500, height = 1000, unit = 'px')

    cell_fraq_both_long_deconv_seg_imm <- filter(cell_fraq_both_long_deconv_seg, cell_type %in% cells_immune)

    ggplot(data = cell_fraq_both_long_deconv_seg_imm, aes(x = reorder(dcc_filename, rowid), y = fraction, fill = cell_type)) +
      geom_bar(position="fill", stat="identity") +
      ggtitle(seg) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 3,
                                       color = anno_color_mapping[cell_fraq_both_long_deconv_seg$Annotation_cell])) +
      scale_x_discrete(labels = cell_fraq_both_long_deconv_seg$Annotation_cell)

    ggsave(file.path(output_dir, 'sanity_check', paste0('barplot_',deconv_type, '_', seg, '_immune.png')),
           width = 1500, height = 1000, unit = 'px')
  })
})

#########################33
# heatmaps

col_fun = viridis(100)


sapply(c('bp', 'sd'), function(deconv_type){
  
  png(filename = file.path(output_dir, 'sanity_check', paste0('hmap_',deconv_type, '_all.png')), width=1000, height=750)
  Heatmap(t(as.matrix(cell_fraq_both[, paste0(c(ct_names, 'other'), '_', deconv_type)])), col = col_fun) %v%
    HeatmapAnnotation(segment = cell_fraq_both$Segment, 
                      col = list(segment = c("stroma" = "green", "tumor" = "blue")))
  dev.off()
  
  sapply(c('stroma', 'tumor'), function(seg){
    ct_frac_seg <- cell_fraq_both[cell_fraq_both$Segment == seg, ]

    png(filename = file.path(output_dir, 'sanity_check', paste0('hmap_',deconv_type, '_', seg, '.png')), width=1000, height=750)
    print(Heatmap(t(as.matrix(ct_frac_seg[, paste0(c(ct_names, 'other'), '_', deconv_type)])), col = col_fun) %v%
            HeatmapAnnotation(roi_type = ct_frac_seg$Annotation_cell))
    dev.off()

    png(filename = file.path(output_dir, 'sanity_check', paste0('hmap_',deconv_type, '_', seg, '_immune.png')), width=1000, height=750)
    print(Heatmap(t(as.matrix(ct_frac_seg[, paste0(cells_immune, '_', deconv_type)])), col = col_fun) %v%
            HeatmapAnnotation(roi_type = ct_frac_seg$Annotation_cell))
    dev.off()
  })
})




################################################
# this is not important
# for sisana
demo <- fread('/home/iganiemi/Documents/phd/st/gene-regulatory-networks/sisana/sisana/example_input/BRCA_TCGA_20_LumA_LumB_samps_5000_genes_exp.tsv')

geomx_norm_batch_eff_rm_path <- '/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2706/geomx_qc_norm_batch_eff_rm.RDS'
geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)

geomx_harm <- data.frame(geomx_obj@assayData$harmony_batch_corr)
geomx_harm <- rownames_to_column(geomx_harm, var = 'Target')
write_tsv(geomx_harm, '/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2706/sisana/geomx_input/geomx_batch123_harmony_corr_expr_mtx.tsv',
          col_names = T)

meta_segment <- sData(geomx_obj)[, c('dcc_filename', 'Segment')]
meta_segment$dcc_filename <- gsub('\\-', '\\.', meta_segment$dcc_filename)
rownames(meta_segment) <- NULL
colnames(meta_segment) <- NULL

write_csv(meta_segment, '/home/iganiemi/Documents/phd/st/gene-regulatory-networks/sisana/sisana/geomx_input/geomx_meta_segment.csv', 
          col_names = F)


meta <- pData(geomx_obj)
meta <- meta[, c('dcc_filename', 'Sample')]

clin <- fread('/home/iganiemi/Documents/phd/st/data/geomx/clinical_data/9_eyemt_patient_clinical_data.csv')
clin <- clin[, c(1, 10:21)]

meta <- left_join(meta, clin)
fwrite(meta, '/home/iganiemi/Documents/phd/st/geomx-processing/data/b12_dcc_clinical_data.csv')

########################
# deconv b1 sanity check

deconv_list <- readRDS('/home/iganiemi/Documents/phd/st/geomx-processing/results/batch1-1903/deconvolution/bayes_prism/bp_res_mid_lvl_ct.RDS')

deconv_tcell <-   deconv_ct <- BayesPrism::get.exp(bp=deconv_list,
                                                   state.or.type="type",
                                                   cell.name='Tcells')

deconv_macro <-   deconv_ct <- BayesPrism::get.exp(bp=deconv_list,
                                                   state.or.type="type",
                                                   cell.name='Macrophages')

tcell_genes <- c('CD8A', 'CD4', 'CTLA4', 'CXCR6', 'IL16')
macro_genes <- c('CD80', 'CD86', 'CXCL16','HLA-DQA2', 'HLA-DPB1')

deconv_tcell <- deconv_tcell[, which(colnames(deconv_tcell) %in% c(tcell_genes, macro_genes))]
colnames(deconv_tcell) <- paste0(colnames(deconv_tcell), "_tcell")

deconv_macro <- deconv_macro[, which(colnames(deconv_macro) %in% c(tcell_genes, macro_genes))]
colnames(deconv_macro) <- paste0(colnames(deconv_macro), "_macro")

identical(rownames(deconv_tcell), rownames(deconv_macro))
deconv_both <- as.data.frame(cbind(deconv_tcell, deconv_macro))
deconv_both <- rownames_to_column(deconv_both, 'dcc')

library(tibble)
library(tidyr)
deconv_both2 <- gather(deconv_both, gene, value, -dcc)

boxpl <- ggplot(data = deconv_both2) +
  geom_boxplot(aes(x = gene, y = value)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylim(0, 150)

print(boxpl)

data_combined <- left_join(data_combined, features_info, by = c('ligand'))
data_filtered <- data_combined_df[data_combined_df$source != data_combined_df$cluster, ]


mean(deconv_tcell$CD4)
