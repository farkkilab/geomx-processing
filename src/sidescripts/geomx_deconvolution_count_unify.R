# count cell nrs and frequencies from sd and bp deconv for AOI and ROI
# TODO finish taking this part out from geomx_roi_hubs_integration
# or is it needed?
library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(tibble)
library(ComplexHeatmap)
library(GeomxTools)
library(readxl)
library(ggpmisc)
library(circlize)
library(tidyr)
library(gtools)
library(ggpubr)

# define vars -------------------------------------------------------------

# cell fraction from deconv below that lvl will be changed to 0 
# max nr of cells = 300 so 0.005 cell fraction is 1 cell/200 cells 1,5 cell/300 cells
min_frac <- 0.01 
#scrna_anno <<- 'mid_lvl_ct_updated' # deconvolution lvl
scrna_anno <- 'mid_lvl_ct_updated'

# input files
proj_dir <<- '~/Documents/phd/st'
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808') # batch123

metadt_path <- file.path(output_dir, 'metadata_full_SENSITIVE.csv')
geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS') 

bp_cellcounts_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', paste0('bp_res_', scrna_anno, '_ct_fraction.csv'))
sd_cellcounts_path <- file.path(output_dir, 'deconvolution', 'spatial_decon', paste0('sd_res_', scrna_anno, '_geomxfiltpc_ct_fraction.csv'))

deconv_list <- list(bp = bp_cellcounts_path, sd = sd_cellcounts_path)
#deconv_list <- list(sd = sd_cellcounts_path)

# output files
source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))

output_ct_frac_deconv_path <- file.path(output_dir, 'deconvolution', paste0('ct_frac_deconv_', scrna_anno, '_4mainimmune.csv'))
output_ct_frac_deconv_roi_path <- file.path(output_dir, 'deconvolution', paste0('ct_frac_deconv_roi_', scrna_anno, '_4mainimmune.csv'))
output_ct_frac_deconv_both_wide_path <- file.path(output_dir, 'deconvolution', paste0('ct_frac_deconv_aoiroi_wide_', scrna_anno, '_4mainimmune.csv'))

##############
# setting up ct names for counting
meta_names <- c('dcc_filename', 'Sample', 'Segment_geomx',  'Segment', 'NACT_status') 


if(scrna_anno == 'mid_lvl_ct_updated'){
  # all names of deconvoluted cells
  ct_names_all <- c("tumor", "Bcells", "Tcells_CD4", "Tcells_other", "Tcells_CD8", 
                    "Fibroblasts_Mesothelial", "Macrophages_Monocytes", "Mast_cells",
                    "NKcells", "Endothelial_cells", "DCs")
  
  # ct from deconv counted as stroma in cycif
  ct_names_stroma <- c("Fibroblasts_Mesothelial", "Endothelial_cells")
  
  # main immune cells from deconv - also counted in cycif phenotyping
  ct_names_immune <- c("Tcells_CD4", "Tcells_CD8", "DCs", "Macrophages_Monocytes") # , , no Bcells in basic phenotyping eg b1b2 
  ct_names_myeloids <- c("Macrophages_Monocytes", "DCs")
  ct_names_lymphoids <- c("Tcells_CD4", "Tcells_CD8")
  
  #TODO organise it better to work both here and in comparison with deconv
  # additional cells from deconv not counted in phenotyping and should be treated as 'other'
  ct_names_other <- c("Tcells_other", "Mast_cells", "Bcells", "NKcells") # "Bcells" goes here when basic phenotyping b1b2 and 'consensus_label_clean'
  
} else if(scrna_anno == 'cell_type'){
  
  # vaharautio
  # ct_names_all <- c("tumor", "Fibroblasts", 'Endothelial.cells', "Mesothelial.cells",
  #                   'Plasma.cells', 'Naive.B.cells', "Tcm.Naive.helper.T.cells",
  #                   "Regulatory.T.cells", "Tem.Trm.cytotoxic.T.cells", "Macrophages", "Classical.monocytes",
  #                   "Migratory.DCs", "pDC", "CD16..NK.cells.1")
  # 
  # # stroma
  # ct_names_stroma <- c("Fibroblasts", 'Endothelial.cells', "Mesothelial.cells")
  # 
  # # main immune cells
  # ct_names_immune <- c('Plasma.cells', 'Naive.B.cells', "Tcm.Naive.helper.T.cells",
  #                      "Regulatory.T.cells", "Tem.Trm.cytotoxic.T.cells", "Macrophages", "Classical.monocytes",
  #                      "Migratory.DCs", "pDC", "CD16..NK.cells.1") 
  # ct_names_myeloids <- c("Macrophages", "Classical.monocytes", "Migratory.DCs", "pDC")
  # ct_names_lymphoids <- c("Tcm.Naive.helper.T.cells", "Regulatory.T.cells", "Tem.Trm.cytotoxic.T.cells")
  
  # hautaniemi
  ct_names_all <- c("tumor", "CAF_3", "Mesothelial", "CAF_1", "CAF_2", "Endothelial",             
                    "Tcells_CD8_NaiveLike", "Tcells_Treg", "Plasma_cells", "Tcells_CD4_Tfh",
                    "DC_2", "B_cells", "Tcells_CD8_EffectorMemory", "Tcells_CD8_EarlyActiv",
                    "Macrophages", "pDC", "ILC")

  # stroma
  ct_names_stroma <- c("CAF_3", "Mesothelial", "CAF_1", "CAF_2", "Endothelial")

  # main immune cells
  ct_names_immune <- c("Tcells_CD8_NaiveLike", "Tcells_Treg", "Plasma_cells", "Tcells_CD4_Tfh",
                       "DC_2", "B_cells", "Tcells_CD8_EffectorMemory", "Tcells_CD8_EarlyActiv",
                       "Macrophages", "pDC", "ILC")
  ct_names_myeloids <- c("Macrophages", "pDC", "DC_2")
  ct_names_lymphoids <- c("Tcells_CD8_NaiveLike", "Tcells_Treg", "Tcells_CD4_Tfh", 
                          "Tcells_CD8_EffectorMemory", "Tcells_CD8_EarlyActiv")

}

# load geomx, merge with cleaned metadata ---------------------------------
# TODO run once again in 1811 with already cleaned metadata and just load meta from geomx
geomx_dcc <- colnames(readRDS(geomx_norm_batch_eff_rm_path))
metadt <- as.data.frame(fread(metadt_path))
metadt <- metadt[metadt$dcc_filename %in% geomx_dcc, ]
rownames(metadt) <- NULL

# load deconv and transform to long format --------------------------------

# loop through sd, bp results, convert to long
ct_frac_deconv_long <- lapply(1:length(deconv_list), function(n){
  ct_frac_deconv <- as.data.frame(fread(deconv_list[[n]], select = c('dcc_filename', ct_names_all)))
  
  # clean too small cell fractions
  ct_frac_deconv <- mutate_at(ct_frac_deconv, all_of(ct_names_all), funs(ifelse((is.na(.) | . < min_frac), 0, .)))
  
  ct_frac_deconv$stroma <- rowSums(ct_frac_deconv[, ct_names_stroma])
  ct_frac_deconv$immune <- rowSums(ct_frac_deconv[, ct_names_immune])
  ct_frac_deconv$other <- rowSums(ct_frac_deconv[, ct_names_other])
  ct_frac_deconv$immune_other <- rowSums(ct_frac_deconv[, c(ct_names_immune, ct_names_other)])
  ct_frac_deconv$myeloids <- rowSums(ct_frac_deconv[, ct_names_myeloids])
  ct_frac_deconv$lymphoids <- rowSums(ct_frac_deconv[, ct_names_lymphoids])
  
  ct_frac_long <- melt(setDT(ct_frac_deconv), id.vars = 'dcc_filename', variable.name = "cell_type")
  colnames(ct_frac_long)[which(colnames(ct_frac_long) == 'value')] <- paste0('ct_frac_', names(deconv_list[n]))
  
  return(ct_frac_long)
})

if(length(deconv_list) > 1){
  ct_frac_deconv_long <- do.call(left_join, ct_frac_deconv_long)
} else {
  ct_frac_deconv_long <- as.data.frame(ct_frac_deconv_long[[1]])
}


# add metadata
ct_frac_deconv_long <- left_join(ct_frac_deconv_long, metadt[, c(meta_names,'sample_roi', 'Nuclei')], by = 'dcc_filename') %>%
  dplyr::rename(total_cell_nr_geomx = Nuclei)

# calculate cell number 
if("bp" %in% names(deconv_list)){
  ct_frac_deconv_long$ct_nr_bp <- round(ct_frac_deconv_long$ct_frac_bp * ct_frac_deconv_long$total_cell_nr_geomx)
}

ct_frac_deconv_long$ct_nr_sd <- round(ct_frac_deconv_long$ct_frac_sd * ct_frac_deconv_long$total_cell_nr_geomx)

# before merging per ROI, remove AOI from tsi regions which lost their pair during eg QC
roi_incomplete <- ct_frac_deconv_long %>%
  select(dcc_filename, sample_roi, Segment_geomx) %>%
  distinct() %>%
  group_by(sample_roi) %>%
  mutate(nr_aoi = n()) %>%
  filter(Segment_geomx == 'tsi' & nr_aoi == 1)

file.path(output_dir, 'cycif_integration', 'b123_ct_frac_rois_incomplete.csv')

# merge per ROI for incomplete ROIs: 
# assume 50/50 ct number in both AOIs - double numbers
# cell fractions stay the same
# TODO remember about it
#TODO summarise depends if bp was there..
ct_frac_deconv_long_roi_incomplete <- ct_frac_deconv_long %>%
  filter(dcc_filename %in% roi_incomplete$dcc_filename) %>% # get 20 incomplete ROIs
  group_by(sample_roi, cell_type) %>%
  summarise(total_cell_nr_geomx = total_cell_nr_geomx*2, ct_nr_bp = ct_nr_bp*2, ct_nr_sd = ct_nr_sd*2,
            ct_frac_bp = ct_frac_bp, ct_frac_sd = ct_frac_sd) %>%
  #summarise(total_cell_nr_geomx = total_cell_nr_geomx*2, ct_nr_sd = ct_nr_sd*2, ct_frac_sd = ct_frac_sd) %>%
  ungroup()

# merge per ROI
#TODO summarise depends if bp was there..
ct_frac_deconv_long_roi <- ct_frac_deconv_long %>%
  filter(!(dcc_filename %in% roi_incomplete$dcc_filename)) %>% # rmv 20 AOIs from incomplete ROIs
  group_by(sample_roi, cell_type) %>%
  summarise(total_cell_nr_geomx = sum(total_cell_nr_geomx), ct_nr_bp = sum(ct_nr_bp), ct_nr_sd = sum(ct_nr_sd),
            ct_frac_bp = mean(ct_frac_bp), ct_frac_sd = mean(ct_frac_sd)) %>%
  #summarise(total_cell_nr_geomx = sum(total_cell_nr_geomx), ct_nr_sd = sum(ct_nr_sd), ct_frac_sd = mean(ct_frac_sd)) %>%
  ungroup()

ct_frac_deconv_long_roi <- rbind(ct_frac_deconv_long_roi, ct_frac_deconv_long_roi_incomplete)


fwrite(ct_frac_deconv_long, output_ct_frac_deconv_path)
fwrite(ct_frac_deconv_long_roi, output_ct_frac_deconv_roi_path)

#######################################
# additional transformations
# transform to long and save both
ct_frac_aoi <- ct_frac_deconv_long[, c('dcc_filename','sample_roi', 'cell_type', 'ct_frac_sd')] %>%
  spread(key = 'cell_type', value = 'ct_frac_sd')
colnames(ct_frac_aoi)[3:ncol(ct_frac_aoi)] <- paste0('ct_frac_sd_aoi_', colnames(ct_frac_aoi)[3:ncol(ct_frac_aoi)])

ct_frac_roi <- ct_frac_deconv_long_roi[, c('sample_roi', 'cell_type', 'ct_frac_sd')] %>%
  spread(key = 'cell_type', value = 'ct_frac_sd')
colnames(ct_frac_roi)[2:ncol(ct_frac_aoi)] <- paste0('ct_frac_sd_roi_', colnames(ct_frac_roi)[2:ncol(ct_frac_aoi)])

# calculate roi immunefraction
ct_immunefrac_roi <- ct_frac_deconv_long_roi[, c('sample_roi', 'cell_type', 'ct_frac_sd')] %>%
  filter(cell_type %in% c(!!ct_names_immune, 'immune')) %>%
  spread(key = 'cell_type', value = 'ct_frac_sd') %>%
  mutate_at(ct_names_immune, ~./immune) %>%
  select(-immune)
colnames(ct_immunefrac_roi)[2:ncol(ct_immunefrac_roi)] <- paste0('ct_immunefrac_sd_roi_', colnames(ct_immunefrac_roi)[2:ncol(ct_immunefrac_roi)])


ct_frac_both <- left_join(ct_frac_aoi, ct_frac_roi) %>%
  left_join(ct_immunefrac_roi)

fwrite(ct_frac_both, output_ct_frac_deconv_both_wide_path)

#####################################################################
# clusters distributions vs ct frac and clinical vars ---------------------

# each ct freq in tumor/stroma segment across nact status + freq cluster
#clust_type <- paste0(clust_types[1], '_label')
clust_type <- 'roi_cluster_label_gmm'
clust_type_name <- 'gmm' # for plotting

out_dir_clust <- file.path('~/Documents/phd/st/geomx-processing/results/batch123-2808/downstream/roi_clusters_freq_finegrained_deconv_hautaniemi', clust_type_name)
dir.create(out_dir_clust, recursive = T)

#######
ct_frac_deconv_long <- ct_frac_deconv_long %>%
  left_join(metadt[, c('dcc_filename', clust_type)]) %>%
  mutate(segment_nact = paste0(Segment, '_', NACT_status, '_NACT')) %>%
  filter(cell_type %in% !!ct_names_immune)

ct_frac_deconv_long_roi <- ct_frac_deconv_long_roi %>%
  left_join(distinct(metadt[, c('sample_roi','Segment_geomx', 'NACT_status', clust_type)])) %>%
  filter(cell_type %in% !!ct_names_immune)

# boxpl all cells at once, color by segment_nact
ggplot(ct_frac_deconv_long, aes(x = factor(cell_type), y = ct_frac_sd, fill = segment_nact)) +
  geom_boxplot() +
  geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = FALSE, size = 0.2, label.size = 2.8) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 8)) +
  labs(title = 'ct frequencies aross segment and NACT status', x = "cell type", y = "AOI ct fraction", fill = 'segment and NACT status')

ggsave(file.path(out_dir_clust, paste0('ct_frac_aoi_segment_nact_w_wilcox.png')),
       width = 12, height = 8, units = c("in"))

# boxpl faceted by cell, color by segment_nact all clusters at once
ggplot(ct_frac_deconv_long, aes(x = factor(get(clust_type)), y = ct_frac_sd, fill = factor(segment_nact))) +
  geom_boxplot() +
  labs(title = 'ct freq aross segment and nact status', x = "cell type", y = "ct frac sd in AOI") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  facet_wrap(~cell_type, dir="v")

ggsave(file.path(out_dir_clust, paste0('ct_frac_aoi_per_cluster_', clust_type_name, '_segment_nact.png')),
       width = 10, height = 8, units = c("in"))

# boxpl faceted by cell, color by segment_nact all clusters at once
ggplot(ct_frac_deconv_long_roi, aes(x = factor(get(clust_type)), y = ct_frac_sd, fill = cell_type)) +
  geom_boxplot() +
  labs(title = 'ct freq aross segment and nact status', x = "cell type", y = "ct frac sd in ROI") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

ggsave(file.path(out_dir_clust, paste0('ct_frac_roi_per_cluster_', clust_type_name, '_nact.png')),
       width = 10, height = 8, units = c("in"))

############################################################
# stacked oxplots with ct distribution across clusters

# immunefrac to long
ct_immunefrac_roi_long <- pivot_longer(ct_immunefrac_roi, cols = starts_with("ct_immunefrac_sd_roi"),
                                       names_to = 'cell_type', values_to = 'ct_immunefrac_sd_roi') %>%
  mutate(cell_type = gsub('ct_immunefrac_sd_roi_', '', cell_type)) %>%
  left_join(distinct(metadt[, c('sample_roi', clust_type)]))

# stacked bar plot with all samples per cluster
ggplot(ct_immunefrac_roi_long, aes(x = sample_roi, y = ct_immunefrac_sd_roi, fill = cell_type)) +
  geom_bar(stat = "identity") +
  labs(title = paste0("ROI ", clust_type, " ct frac of immune"), x = "ROIs", y = "ct_fraction") +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  facet_wrap(~ get(clust_type), scales = "free", ncol = 2)

ggsave(file.path(out_dir_clust, paste0('barplot_clust_immunefrac_', clust_type,  '.png')),
       width = 12, height = 12, units = c("in"))


ct_frac_clust_mean <- ct_immunefrac_roi_long %>%
  dplyr::group_by_at(c(clust_type, 'cell_type')) %>%
  dplyr::summarise(mean_ct_frac_of_immune = mean(ct_immunefrac_sd_roi))

# stacked barplot for mean ct fraction per cluster
ggplot(ct_frac_clust_mean, aes(x = get(clust_type), y = mean_ct_frac_of_immune, fill = cell_type)) +
  geom_bar(stat = "identity") +
  labs(title = paste0("mean ct frac of immune per clust in ", clust_type), x = "ROI clusters", y = "mean ct fraction") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 6))

ggsave(file.path(out_dir_clust, paste0('barplot_mean_clust_immunefrac_', clust_type,  '.png')))


###########################################
# ct frac clusters distribution across samples (pre, post, HRD, PFS, OS)
# TODO move to downstream analysis of clusters

vars_labels <- c('NACT_status', 'HRP_status', 'primary_treatment_response', 'PFS_group')
vars_cont <- c('TMB', 'ovaHRDscar_score', 'PFS_days', 'OS_days')

pfs_unique <- distinct(metadt, Patient, PFS_days)
pfs_lo <- unname(quantile(pfs_unique$PFS_days, 0.25))
pfs_hi <- unname(quantile(pfs_unique$PFS_days, 0.75))

metadt <- mutate(metadt, PFS_group = ifelse(PFS_days <= pfs_lo, 'low', ifelse(PFS_days >= pfs_hi, 'high', 'medium')))

# count total nr of clusters in dataset (for ordering)
cluster_labels_count <- metadt %>%
  select(dcc_filename, Sample, !!vars_labels, !!vars_cont, !!clust_type) %>%
  group_by(get(clust_type)) %>%   
  mutate(clust_name_occur_total = n()) %>%
  ungroup() %>%
  as.data.frame()


# all combinations
clust_allcombs <- tidyr::expand(metadt, Sample, get(clust_type))
colnames(clust_allcombs) <- c('Sample', clust_type) #fixing stupid names

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
  left_join(distinct(metadt[, c('Sample', vars_labels, vars_cont)]))


# stacked barplot for nr of clusters across samples faceted by discrete vars
for(label_var in vars_labels){
  
  # nrs of AOIs from given ROI cluster 
  ggplot(cluster_labels_count, aes(x = reorder(Sample, clust_name_occur_total), fill = get(clust_type))) +
    geom_bar(stat = "count") +
    labs(title = paste0("nr of AOIs per sample across ", label_var), x = "Sample", y = "AOI nr", fill='ROI cluster type') +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 8)) +
    facet_wrap(~get(label_var), dir="v", scales = "free")
  
  ggsave(file.path(out_dir_clust, paste0('aoi_nr_cluster_', clust_type_name, '_color_', label_var, '.png')))
  
  # frequencies of ROIs clusters - stacked barplots
  ggplot(cluster_freqs_per_sample, aes(x = Sample, y = clust_freq_per_sample, fill = get(clust_type))) +
    geom_bar(stat = "identity") +
    labs(title = paste0("frequencies of ROI clusters per sample across ", label_var), x = "Sample", y = "ROI cluster frequency", fill='ROI cluster type') +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 8)) +
    facet_wrap(~get(label_var), dir="v", scales = "free")
  
  ggsave(file.path(out_dir_clust, paste0('roi_freq_cluster_', clust_type_name, '_color_', label_var, '.png')))
  
  # frequencies of ROIs clusters - boxplots
  ggplot(cluster_freqs_per_sample, aes(x = get(clust_type), y = clust_freq_per_sample, fill = get(label_var))) +
    geom_boxplot() +
    geom_point(position= position_jitterdodge(dodge.width = 1, jitter.width= .3, jitter.height = 0),
               size= 0.5, alpha = 0.6) +
    geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = TRUE, size = 0.2, label.size = 2.8) +
    labs(title = paste0("frequencies of ROI clusters per sample across ", label_var), x = "ROI cluster", y = "ROI cluster frequency", fill= paste0(label_var)) +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 8)) 
  #+facet_wrap(~NACT_status, dir="v", scales = "free")
  
  ggsave(file.path(out_dir_clust, paste0('roi_freq_cluster_', clust_type_name, '_boxpl_color_', label_var, '.png')))
}

# scatterplots for continuous vars
for(cont_var in vars_cont){
  
  ggplot(cluster_freqs_per_sample, aes(x = clust_freq_per_sample, y = get(cont_var), color = get(clust_type))) +
    geom_point(size = 3) +
    xlab("ROI cluster frequency in sample") +
    ylab(cont_var) +
    scale_color_discrete(name = clust_type) +
    theme_bw() +
    facet_wrap(~get(clust_type), scales = "fixed", dir="v")
  
  ggsave(file.path(out_dir_clust, paste0('roi_freq_cluster_', clust_type_name, '_scatter_color_', cont_var, '.png')))
}


