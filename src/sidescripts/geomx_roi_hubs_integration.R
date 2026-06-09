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

# TODO ct_frac_cycif_roi 129 ROI, but 126 in the final ct_frac_all - probably without hub labs - to check
# TODO ensure # "S130_iOme_5" "S130_iOme_6" "S197_iOme_1" - everywhere in metadata - probably removed during QC
# TODO compare with deconv made per-ROI (not per-AOI)
# TODO copare with annotation cell from semi-automated method

# define vars -------------------------------------------------------------

# cell fraction from deconv below that lvl will be changed to 0 
# max nr of cells = 300 so 0.005 cell fraction is 1 cell/200 cells 1,5 cell/300 cells
min_frac <- 0.01 
min_label_frac <- 0.05 # min fraction of immune cells with given component/community label to give a label to whole ROI
# TODO 'consensus_label_clean_all_ct' when taking bcells!!, or final_label in new Elias files
cycif_main_ct_label <- 'final_label' 
hubs_labels_list <- c('component_label', 'community_cluster_label') # before with 'interaction_hub_type'

# main immune cells from deconv - also counted in cycif phenotyping
#TODO decide if B_cells should be counted as immune or as other
ct_names_immune <- c("Tcells_CD4", "Tcells_CD8", "DCs", "Macrophages_Monocytes", "Bcells") # , , no Bcells in basic phenotyping eg b1b2 
ct_names_myeloids <- c("Macrophages_Monocytes", "DCs")
ct_names_lymphoids <- c("Tcells_CD4", "Tcells_CD8")
# additional cells from deconv not counted in phenotyping and should be treated as 'other'
ct_names_other <- c("Tcells_other", "Mast_cells", "NKcells") # "Bcells" goes here when basic phenotyping b1b2 and 'consensus_label_clean'

###########
proj_dir <<- '~/Documents/phd/st'
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808') # batch123

metadt_path <- file.path(proj_dir,'data', 'geomx', 'metadata_full_SENSITIVE.csv')
geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS') 

bp_cellcounts_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', 'bp_res_mid_lvl_ct_updated_ct_fraction.csv')
sd_cellcounts_path <- file.path(output_dir, 'deconvolution', 'spatial_decon', 'sd_res_mid_lvl_ct_updated_geomxfiltpc_ct_fraction.csv')

######################
# all cells within ROIs with components/communities annotations computed with geomx_cycif_integration.R

# # combined myeloids, no bcells,  min 20 cells components, no mixing
hubs_inroi_path <- file.path(output_dir, "cycif_integration", "batch3tls_hubs_cells_inroi_combined_myeloids_min20cells_nomix_dt171715_ct10_dt300_new_polygons.csv")
hubs_outname <- "batch3tls_combined_myeloids_min20cells_nomix_new_polygons"
ct_frac_deconv_outname <- "b123_ct_frac_deconv_bcells"
clust_names <- paste0('cluster_', seq(0, 16))
clust_manualnames <- c('CD4_Macro',
                       'CD8_Macro', 'loCD8_Macro', 'CD8', 'CD4_CD8', 'Macro',
                       'loCD8_Macro', 'CD4', 'loCD4_CD8', 'loCD4_CD8_loMacro', 'CD8_loMacro',
                       'CD4_loMacro', 'loCD4_loCD8_Macro', 'loCD4_loCD8_loMacro', 'loCD4_Macro', 'loCD8_Macro',
                       'loCD8_Macro')

#####################
# # combined myeloids, no bcells,  min 2 components, no mixing
# hubs_inroi_path <- file.path(output_dir, "cycif_integration", "batch3tls_hubs_cells_inroi_combined_myeloids_min2comp_nomix_dt171715_ct10_dt300.csv")
# hubs_outname <- "batch3tls_combined_myeloids_min2comp_nomix"
# ct_frac_deconv_outname <- "b123_ct_frac_deconv_bcells"
# clust_names <- paste0('cluster_', seq(0, 16))
# clust_manualnames <- c('loCD8_Macro',
#                        'loCD4_CD8_loMacro', 'CD4', 'Macro', 'CD4_Macro', 'CD8_Macro',
#                        'CD8_Macro', 'CD8', 'loCD4_loCD8_Macro', 'CD8_Macro', 'CD4_CD8',
#                        'CD4_CD8_Macro', 'CD4_CD8_Macro', 'CD4_CD8_Macro', 'CD4_Macro', 'CD8_Macro',
#                        'CD4_CD8')

#####################
# combined myeloids, no bcells,  min 2 components, forced mixing
# hubs_inroi_path <- file.path(output_dir, "cycif_integration", "batch3tls_hubs_cells_inroi_combined_myeloids_min2comp_mixed_dt171715_ct10_dt300.csv")
# hubs_outname <- "batch3tls_combined_myeloids_min2comp_mixed"
# ct_frac_deconv_outname <- "b123_ct_frac_deconv_bcells"
# clust_names <- paste0('cluster_', seq(0, 16))
# clust_manualnames <- c('loCD8_Macro', 'CD8_loMacro', 'CD4_Macro', 'CD4_CD8', 'CD8_Macro', 'loCD4_Macro',
#                        'loCD4_loCD8_loMacro', 'loCD4_loCD8_Macro', 'loCD8_Macro', 'CD4_loMacro',
#                        'loCD4_loCD8_Macro', 'loCD4_CD8', 'loCD4_CD8_loMacro', 'loCD4_Macro',
#                        'loCD4_loCD8_loMacro', 'loCD8_Macro', 'CD4_loCD8')

#####################
# bcells + combined myeloids, min 3 components, old distances from centroids
# hubs_inroi_path <- file.path(output_dir, "cycif_integration", "batch3tls_hubs_cells_inroi_bcells_combined_myeloids_15151517_ct10_dt300.csv")
# hubs_outname <- "batch3tls_bcells_combined_myeloids_15151517_ct10"
# ct_frac_deconv_outname <- "b123_ct_frac_deconv_bcells" # bcells are counted separately not as 'other_immune'
# clust_names <- paste0('cluster_', seq(0, 19))
# clust_manualnames <- c('Macro', 'Macro_loCD4_loCD8', 'Macro_loCD8_loBcells', 'CD4', 'loMacro_Bcells',
#                        'Macro_CD8', 'Macro_loMixed', 'Macro_loCD8', 'loMacro_CD4', 'CD8', 'Macro_loBcells',
#                        'loCD4_CD8', 'Macro_loCD8', 'Bcells_loCD8', 'Macro_CD8', 'loMacro_loCD4_loBcells',
#                        'loMacro_CD8_loBcells', 'Macro_loCD4', 'Bcells', 'loMacro_Bcells')


##################################################
###################################################

# output files
source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))

outp_plot_dir <- file.path(output_dir, 'cycif_integration', paste0('ct_frac_comparison_', hubs_outname))
dir.create(outp_plot_dir, recursive = T)
dir.create(file.path(outp_plot_dir,'hmaps'), recursive = T)
dir.create(file.path(outp_plot_dir,'ct_distribution'), recursive = T)

output_ct_frac_deconv_path <- file.path(outp_plot_dir, paste0(ct_frac_deconv_outname, '.csv'))
output_ct_frac_deconv_roi_path <- file.path(outp_plot_dir, paste0(ct_frac_deconv_outname, '_roi.csv'))

output_ct_frac_cycif_roi_path <- file.path(outp_plot_dir,  paste0('ct_frac_cycif_roi_', hubs_outname, '.csv'))
output_ct_frac_all_roi_path <- file.path(outp_plot_dir, paste0('ct_frac_all_roi_', hubs_outname, '.csv'))
output_roi_labels_path <- file.path(outp_plot_dir, paste0('roi_labels_', hubs_outname, '.csv'))

##############
meta_names <- c('dcc_filename', 'Sample', 'Annotation_cell', 'Roi_geomx', "roi_cluster_label_gmm", "roi_cluster_label_hclust",
                'Segment_geomx',  'Segment', 'tCycIF_preselection_initial_label') 

meta_names_per_roi <- c('Sample', 'Segment_geomx', 'Annotation_cell', 'tCycIF_preselection_initial_label', 
                        "roi_cluster_label_gmm", "roi_cluster_label_hclust")

ct_names_all <- c("tumor", "Bcells", "Tcells_CD4", "Tcells_other", "Tcells_CD8", 
                  "Fibroblasts_Mesothelial", "Macrophages_Monocytes", "Mast_cells",
                  "NKcells", "Endothelial_cells", "DCs")

# ct from deconv counted as stroma in cycif
ct_names_stroma <- c("Fibroblasts_Mesothelial", "Endothelial_cells")

# load geomx, merge with cleaned metadata ---------------------------------
# TODO run once again in 1811 with already cleaned metadata and just load meta from geomx
geomx_dcc <- colnames(readRDS(geomx_norm_batch_eff_rm_path))
metadt <- as.data.frame(fread(metadt_path))
metadt <- metadt[metadt$dcc_filename %in% geomx_dcc, ]
rownames(metadt) <- NULL

# load deconv and transform to long format --------------------------------

deconv_list <- list(bp = bp_cellcounts_path, sd = sd_cellcounts_path)

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

ct_frac_deconv_long <- do.call(left_join, ct_frac_deconv_long)

# add metadata
ct_frac_deconv_long <- left_join(ct_frac_deconv_long, metadt[, c(meta_names, 'Nuclei', 'sample_roi')], by = 'dcc_filename') %>%
  dplyr::rename(total_cell_nr_geomx = Nuclei)

# calculate cell number 
ct_frac_deconv_long$ct_nr_bp <- round(ct_frac_deconv_long$ct_frac_bp * ct_frac_deconv_long$total_cell_nr_geomx)
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
ct_frac_deconv_long_roi_incomplete <- ct_frac_deconv_long %>%
  filter(dcc_filename %in% roi_incomplete$dcc_filename) %>% # get 20 incomplete ROIs
  group_by(sample_roi, cell_type) %>%
  summarise(total_cell_nr_geomx = total_cell_nr_geomx*2, ct_nr_bp = ct_nr_bp*2, ct_nr_sd = ct_nr_sd*2,
            ct_frac_bp = ct_frac_bp, ct_frac_sd = ct_frac_sd) %>%
  ungroup()

# merge per ROI
ct_frac_deconv_long_roi <- ct_frac_deconv_long %>%
  filter(!(dcc_filename %in% roi_incomplete$dcc_filename)) %>% # rmv 20 AOIs from incomplete ROIs
  group_by(sample_roi, cell_type) %>%
  summarise(total_cell_nr_geomx = sum(total_cell_nr_geomx), ct_nr_bp = sum(ct_nr_bp), ct_nr_sd = sum(ct_nr_sd),
            ct_frac_bp = mean(ct_frac_bp), ct_frac_sd = mean(ct_frac_sd)) %>%
  ungroup()

ct_frac_deconv_long_roi <- rbind(ct_frac_deconv_long_roi, ct_frac_deconv_long_roi_incomplete)


fwrite(ct_frac_deconv_long, output_ct_frac_deconv_path)
fwrite(ct_frac_deconv_long_roi, output_ct_frac_deconv_roi_path)

# calculate cell nr/fractions from phenotyped cells in cycif --------------

hubs_cells_inroi <- fread(hubs_inroi_path)

# rename cells to match deconvolution
hubs_cells_inroi$cell_type <-  mapvalues(hubs_cells_inroi[[cycif_main_ct_label]], 
                                         from = c("Macrophages", "Myeloids", "CD4_Tcells", "CD8_Tcells",
                                                  "Tumor", "CD11c", "Undefined", "NK", 
                                                  "Stroma", "CD31+_Endothelial", "HEV+_Endothelial", "Bcells"),
                                         to=c("Macrophages_Monocytes", "Macrophages_Monocytes", "Tcells_CD4", "Tcells_CD8",
                                              "tumor", "DCs", "other", "NKcells", "stroma", "stroma", "stroma", "Bcells"))

# count nr of cells per ROI and AOI
ct_frac_cycif_roi <- as.data.frame(dcast(hubs_cells_inroi, sample_roi ~ cell_type))
ct_frac_cycif_roi$total_cell_nr_cycif <- rowSums(ct_frac_cycif_roi[, -1])
ct_frac_cycif_roi$immune <- rowSums(ct_frac_cycif_roi[, intersect(ct_names_immune, colnames(ct_frac_cycif_roi))])
ct_frac_cycif_roi$immune_other <- rowSums(ct_frac_cycif_roi[, c(intersect(ct_names_immune, colnames(ct_frac_cycif_roi)), "other")])
ct_frac_cycif_roi$lymphoids <- rowSums(ct_frac_cycif_roi[, intersect(ct_names_lymphoids, colnames(ct_frac_cycif_roi))])
if(length(intersect(ct_names_myeloids, colnames(ct_frac_cycif_roi))) > 1){
  ct_frac_cycif_roi$myeloids <- rowSums(ct_frac_cycif_roi[, intersect(ct_names_myeloids, colnames(ct_frac_cycif_roi))])
} else{
  ct_frac_cycif_roi$myeloids <- ct_frac_cycif_roi[, intersect(ct_names_myeloids, colnames(ct_frac_cycif_roi))]
}

# transform to long
ct_frac_cycif_long_roi <- melt(setDT(ct_frac_cycif_roi), id.vars = 'sample_roi', variable.name = "cell_type", value.name = "ct_nr_cycif")

# add total cell nr as column, remove total and calculate ct frac
ct_frac_cycif_long_roi <- left_join(ct_frac_cycif_long_roi, ct_frac_cycif_roi[, c('sample_roi', 'total_cell_nr_cycif')]) %>%
  filter(cell_type != 'total_cell_nr_cycif') %>%
  mutate(ct_frac_cycif = ct_nr_cycif / total_cell_nr_cycif) %>%
  filter(sample_roi %in% metadt$sample_roi) # rmv roi not in metadata eg removed during qc

fwrite(ct_frac_cycif_long_roi, output_ct_frac_cycif_roi_path)

##########################################################
##########################################################
# experimental labels coverage checking
# count fractions of cells with given hub label ---------------------------

hubs_inroi_labs <- dplyr::select(hubs_cells_inroi, sample_roi, !!hubs_labels_list) %>%
  dplyr::mutate(across(c('component_label'), clean_labs)) %>%
  dplyr::mutate(across(hubs_labels_list, ~replace(., . ==  '' , 'notinhub')))

# count nr of each labels 
hubs_inroi_labs <- melt(setDT(hubs_inroi_labs), id.vars = 'sample_roi', variable.name = "label_type", value.name = "label")
hubs_inroi_labs<- group_by(hubs_inroi_labs, sample_roi, label_type, label) %>%
  count() %>%
  ungroup() %>%
  tidyr::spread(label_type, n)

# merge with nr of immune cells and compute fraction of immune cells labelled with given label
hubs_inroi_labs_immunefrac <- left_join(hubs_inroi_labs, ct_frac_cycif_roi[, c('sample_roi', 'immune')]) %>%
  filter(label != 'notinhub') %>%
  mutate(across(hubs_labels_list, ~./immune)) %>%
  select(-immune)

hubs_inroi_labs_immunefrac_dupl <- hubs_inroi_labs_immunefrac %>% 
  group_by(sample_roi) %>% 
  filter(n()>1) %>%
  ungroup()

# remove labels with insufficient fraction of immune cells labelled with a given label
hubs_inroi_labs_frequent <- hubs_inroi_labs_immunefrac %>%
  mutate(across(hubs_labels_list, ~ifelse(. >= min_label_frac, label, NA))) %>%
  filter(if_any(hubs_labels_list, ~!is.na(.))) %>%
  select(-label)

# merge network labels in roi
hubs_inroi_labs_frequent <- hubs_inroi_labs_frequent %>%
  group_by(sample_roi) %>%
  mutate(component_label = paste0(unique(na.omit(component_label)), collapse = "_")) %>%
  distinct() %>%
  mutate(across(hubs_labels_list[!hubs_labels_list == 'component_label'], ~paste0(unique(na.omit(.)), collapse = "|"))) %>%
  distinct() 

# relabel after merging to ensure ordering
hubs_inroi_labs_frequent$component_label <- clean_labs(hubs_inroi_labs_frequent$component_label)

colnames(hubs_inroi_labs_frequent) <- c('sample_roi', paste0(hubs_labels_list, '_freq', as.character(min_label_frac)))

# weird to check:
# b2 S069_pOme_roi-002 - all network_hub_type within CD4/CD8 while community label is CD11_Iba1

##########################################################
##########################################################

# clean hubs labs ---------------------------------------------------------
# TODO a bit of repetition with the previous section

# remove cells without labels
hubs_inroi <- dplyr::select(hubs_cells_inroi, sample_roi, !!hubs_labels_list) %>%
  dplyr::mutate(across(c('component_label'), clean_labs)) %>%
  dplyr::filter(if_any(hubs_labels_list, ~!. == '')) %>%
  distinct()

# group by rois and merge labels
hubs_inroi <- hubs_inroi %>%
  group_by(sample_roi) %>%
  mutate(component_label = paste0(unique(na.omit(component_label)), collapse = "_")) %>%
  distinct() %>%
  mutate(across(hubs_labels_list[!hubs_labels_list == 'component_label'], ~paste0(unique(na.omit(.)), collapse = "|"))) %>%
  distinct() 

hubs_inroi$component_label <- clean_labs(hubs_inroi$component_label) # relabel after merging to ensure ordering

# merge all information together ------------------------------------------

# immune_other = all immune + 'other'
label_vars <- c("Annotation_cell", hubs_labels_list, paste0(hubs_labels_list, '_freq', as.character(min_label_frac)),
                "tCycIF_preselection_initial_label_cleaned", "roi_cluster_label_gmm", "roi_cluster_label_hclust")

# join cycif and deconv ct fraction tables, add metadata and hubs info
ct_frac_all <- left_join(ct_frac_cycif_long_roi, ct_frac_deconv_long_roi, by = c('sample_roi', 'cell_type')) %>%
  left_join(distinct(metadt[, c(meta_names_per_roi, 'sample_roi')]), by = 'sample_roi') %>% # before metadt_filt
  left_join(hubs_inroi, by = 'sample_roi') %>%
  left_join(hubs_inroi_labs_frequent, by = 'sample_roi')

ct_frac_all[ct_frac_all == ""] <- NA

# clean labels
ct_frac_all$tCycIF_preselection_initial_label_cleaned <- clean_labs(ct_frac_all$tCycIF_preselection_initial_label)
ct_frac_all$tCycIF_preselection_initial_label_cleaned <- ifelse(ct_frac_all$tCycIF_preselection_initial_label_cleaned == '', 'otherlabel',
                                                                ct_frac_all$tCycIF_preselection_initial_label_cleaned)

ct_frac_all <- ct_frac_all %>%
  dplyr::mutate(across(label_vars, ~replace(., is.na(.) , 'nolabel')))

ct_frac_all$community_cluster_label_manualnames <- ct_frac_all$community_cluster_label
ct_frac_all[[paste0('community_cluster_label_manualnames_freq', as.character(min_label_frac))]] <- ct_frac_all[[paste0('community_cluster_label_freq', as.character(min_label_frac))]]


for(i in 1:length(clust_names)){
  ct_frac_all$community_cluster_label_manualnames <- gsub(paste0(clust_names[i], '$'), clust_manualnames[i], ct_frac_all$community_cluster_label_manualnames)
  ct_frac_all$community_cluster_label_manualnames <- gsub(paste0(clust_names[i], '\\|'), paste0(clust_manualnames[i], '|'), ct_frac_all$community_cluster_label_manualnames)
  
  ct_frac_all[[paste0('community_cluster_label_manualnames_freq', as.character(min_label_frac))]] <- gsub(paste0(clust_names[i], '$'), clust_manualnames[i], ct_frac_all[[paste0('community_cluster_label_manualnames_freq', as.character(min_label_frac))]])
  ct_frac_all[[paste0('community_cluster_label_manualnames_freq', as.character(min_label_frac))]] <- gsub(paste0(clust_names[i], '\\|'), paste0(clust_manualnames[i], '|'), ct_frac_all[[paste0('community_cluster_label_manualnames_freq', as.character(min_label_frac))]])
}

fwrite(ct_frac_all, output_ct_frac_all_roi_path)

labs_all <- ct_frac_all %>%
  select(sample_roi, !!label_vars, !!c('community_cluster_label_manualnames', 
                                       paste0('community_cluster_label_manualnames_freq', as.character(min_label_frac)))) %>%
  distinct()

labs_all[labs_all==""]<- "nolabel"
labs_all[is.na(labs_all)]<- "nolabel"

fwrite(labs_all, output_roi_labels_path)

# compare cell fractions --------------------------------------------------

label_vars_forplots <- c("Annotation_cell", "component_label", paste0('component_label_freq', as.character(min_label_frac)),
                         "roi_cluster_label_gmm", "roi_cluster_label_hclust", 'community_cluster_label_manualnames', 
                         paste0('community_cluster_label_manualnames_freq', as.character(min_label_frac)))

# scatterplot with geomx vs cycif total cell count
cell_count_roi <- select(ct_frac_all, sample_roi, Segment_geomx, total_cell_nr_cycif, total_cell_nr_geomx) %>%
  distinct()

cell_count_scatter <- ggplot(data = cell_count_roi, aes(x = total_cell_nr_cycif, y = total_cell_nr_geomx)) +
  geom_point(aes(color = Segment_geomx)) +
  ggtitle('total cell count per ROI cycif vs geomx') +
  geom_smooth(method='lm', formula= y~x) +
  stat_correlation(method = 'pearson')

ggsave(file.path(outp_plot_dir,'ct_distribution', paste0('total_cellnr_cycif_vs_geomx.png')),
       width = 2000, height = 2000, unit = 'px')


##########################################################
# compare ct fractions between deconv sd/bp and cycif phenotyping
comp_type <- 'ct_frac' # or ct_nr

# scatterplots with value comparisons between methods
for(value_comb in c('bp_sd', 'bp_cycif', 'sd_cycif')){
  vals <- unlist(strsplit(value_comb, split = '_'))
  
  print(value_comb)
  # per cell type
  for(color_var in label_vars_forplots){
    print(color_var)
    for(ct_name in unique(ct_frac_all$cell_type)){
      print(ct_name)
      ct_frac_ct <- ct_frac_all[ct_frac_all$cell_type == ct_name, ]

      ct_scatter <- ggplot(data = ct_frac_ct, aes(x = get(paste0(comp_type, '_', vals[1])), y = get(paste0(comp_type, '_', vals[2])))) +
        geom_point(aes(color = get(color_var), shape = Segment_geomx)) +
        geom_smooth(method='lm', formula= y~x) +
        stat_correlation(method = 'pearson', output.type = 'text') +
        labs(title = paste(ct_name, comp_type, vals[1], 'vs', vals[2]),
             x = paste0(comp_type, '_', vals[1]), y = paste0(comp_type, '_', vals[2]), color = color_var)

      ggsave(file.path(outp_plot_dir, 'ct_distribution', paste0('scatter_', ct_name, '_', vals[1], '_', vals[2], '_', comp_type, '_', color_var, '.png')),
             width = 2000, height = 2000, unit = 'px')
    }
  }
  
  # for all faceted by ct
  all_scatter <- ggplot(data = ct_frac_all, aes(x = get(paste0(comp_type, '_', vals[1])), y = get(paste0(comp_type, '_', vals[2])))) +
    geom_point(aes(color = Segment_geomx)) +
    facet_wrap(~ cell_type) +
    geom_smooth(method='lm', formula= y~x) +
    stat_correlation(method = 'pearson', output.type = 'text') +
    labs(title = paste(comp_type, vals[1], 'vs', vals[2]), x = paste0(comp_type, '_', vals[1]), y = paste0(comp_type, '_', vals[2]))
  
  ggsave(file.path(outp_plot_dir, 'ct_distribution', paste0('scatter_all_', vals[1], '_', vals[2], '_', comp_type, '.png')),
         width = 2000, height = 2000, unit = 'px')
}

#######################################################################
# boxplots with cell nr/fractions per different labels

# make long dataframe for sd/bp/cycif methods
ct_frac_all_method <- select(ct_frac_all, sample_roi, Segment_geomx, cell_type, !!label_vars_forplots, starts_with('ct_frac'), starts_with('ct_nr'))
ct_frac_all_method_long <- melt(setDT(ct_frac_all_method), id.vars = c('sample_roi', 'Segment_geomx', 'cell_type', label_vars_forplots),
                                variable.name = "method_type")

# for each ct faceted by method
for(comp_type in c('ct_frac', 'ct_nr')){
  ct_frac_all_method_comp <- ct_frac_all_method_long[grepl(comp_type, ct_frac_all_method_long$method_type), ]
  
  for(label_var in label_vars_forplots){
    for(ct_name in unique(ct_frac_all$cell_type)){
      
      ct_frac_all_method_comp_ct <- ct_frac_all_method_comp[ct_frac_all_method_comp$cell_type == ct_name, ]
      
      ggplot(ct_frac_all_method_comp_ct, aes(x = get(label_var), y = value)) + 
        geom_boxplot(alpha = .2) +
        geom_point(size = 0.2) + 
        facet_wrap(~ method_type, nrow = length(unique(ct_frac_all_method_comp_ct$method_type))) +
        labs(title = ct_name, x = label_var, y = comp_type) +
        scale_x_discrete(guide = guide_axis(angle = 45))
      
      
      ggsave(file.path(outp_plot_dir,'ct_distribution', paste0('boxpl_', ct_name, '_', label_var, '_', comp_type, '.png')),
             width = 1500, height = 2000, unit = 'px')
    }
  }
}

# for each method, faceted by label
for(method_name in unique(ct_frac_all_method_long$method_type)){
  ct_frac_all_method_sel <- ct_frac_all_method_long[ct_frac_all_method_long$method_type == method_name, ]
  
  for(label_var in label_vars_forplots){
    ggplot(ct_frac_all_method_sel, aes(x = cell_type, y = value)) + 
      geom_boxplot(alpha = .2) +
      geom_point(size = 0.2) + 
      facet_wrap(~ get(label_var), ncol = 4) +
      labs(title = paste(method_name, 'per', label_var), x = 'cell_type', y = method_name) +
      scale_x_discrete(guide = guide_axis(angle = 90))
    
    
    ggsave(file.path(outp_plot_dir,'ct_distribution', paste0('label_boxpl_allct_', label_var, '_', method_name, '.png')),
           width = 1500, height = 3000, unit = 'px')
    
    # filtered to immune
    ggplot(ct_frac_all_method_sel[ct_frac_all_method_sel$cell_type %in% ct_names_immune, ], aes(x = cell_type, y = value)) + 
      geom_boxplot(alpha = .2) +
      geom_point(size = 0.2) + 
      facet_wrap(~ get(label_var), ncol = 4) +
      labs(title = paste(method_name, 'per', label_var), x = 'cell_type', y = method_name) +
      scale_x_discrete(guide = guide_axis(angle = 90))
    
    
    ggsave(file.path(outp_plot_dir,'ct_distribution', paste0('label_boxpl_immune_', label_var, '_', method_name, '.png')),
           width = 1500, height = 3000, unit = 'px')
  }
}


# compare labels ----------------------------------------------------------

labs_comb <- combinations(length(label_vars_forplots), 2, label_vars_forplots)

apply(labs_comb, 1, function(x){
  lab_name1 <- x[1]
  lab_name2 <- x[2]
  
  # compute cross-frequencies of different labels
  labs_cross <- table(labs_all[[lab_name1]], labs_all[[lab_name2]])
  labs_cross <- matrix(labs_cross, ncol=ncol(labs_cross), dimnames=dimnames(labs_cross))
  labs_cross <- labs_cross[, order(colnames(labs_cross))]
  labs_cross <- labs_cross[order(rownames(labs_cross)), ]
  
  
  # do heatmap
  col_fun = colorRamp2(c(0, max(labs_cross)), c("white", "red"))
  
  ht <- Heatmap(labs_cross, col = col_fun, show_heatmap_legend = FALSE,
                cluster_rows = FALSE, cluster_columns = FALSE, row_title = lab_name1, column_title = lab_name2)
  
  at = seq(0, max(labs_cross), by = 1)
  lgd = Legend(at = at, title = "nr_of_matched_labels", legend_gp = gpar(fill = col_fun(at)))
  
  png(filename = file.path(outp_plot_dir, 'hmaps', paste0('heatmap_labels_', lab_name1, '_', lab_name2, '.png')), width=1000, height=750)
  draw(ht, heatmap_legend_list = lgd)
  dev.off()
})


################################################
# metadt <- fread(file.path(output_dir, 'metadata_full_SENSITIVE.csv'))
# # mixed
# # ct_frac_all <- fread('/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2808/cycif_integration/ct_frac_comparison_batch3tls_combined_myeloids_min2comp_mixed/ct_frac_all_roi_batch3tls_combined_myeloids_min2comp_mixed.csv')
# # roi_labels <- fread('/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2808/cycif_integration/ct_frac_comparison_batch3tls_combined_myeloids_min2comp_mixed/roi_labels_batch3tls_combined_myeloids_min2comp_mixed.csv')
# # roi_nolabel_ct_frac_outpath <- '/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2808/cycif_integration/ct_frac_comparison_batch3tls_combined_myeloids_min2comp_mixed/roi_nolabel_ct_nr_cycif.csv'
# # 
# 
# # nomix
# ct_frac_all <- fread('/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2808/cycif_integration/ct_frac_comparison_batch3tls_combined_myeloids_min20cells_nomix/ct_frac_all_roi_batch3tls_combined_myeloids_min20cells_nomix.csv')
# roi_labels <- fread('/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2808/cycif_integration/ct_frac_comparison_batch3tls_combined_myeloids_min20cells_nomix/roi_labels_batch3tls_combined_myeloids_min20cells_nomix.csv')
# roi_nolabel_ct_frac_outpath <- '/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2808/cycif_integration/ct_frac_comparison_batch3tls_combined_myeloids_min20cells_nomix/roi_nolabel_ct_nr_cycif.csv'
# 
# table(roi_labels$community_cluster_label_manualnames_freq0.05, roi_labels$roi_cluster_label_gmm)
# 
# ############################
# # ROIs without any labels
# # 47 ROIs wo labels - 6 has some component (but less than 0.5 freq), 1 - cluster undefined
# 
# roi_nolabel <- roi_labels[roi_labels$community_cluster_label_manualnames_freq0.05 == 'nolabel', ]
# 
# ct_frac_nolabel <- ct_frac_all[ct_frac_all$sample_roi %in% roi_nolabel$sample_roi[roi_nolabel$component_label == 'nolabel']]
# 
# ct_frac_nolabel_long <- spread(ct_frac_nolabel[, c('sample_roi', 'cell_type', 'ct_nr_cycif')], key = cell_type, value = ct_nr_cycif) %>%
#   mutate(total_nr = immune + stroma + tumor + other)
# 
# fwrite(ct_frac_nolabel_long, roi_nolabel_ct_frac_outpath)
# 
# # distribution of cells through rois without labels
# hist(distinct(ct_frac_nolabel[, c('sample_roi', 'total_cell_nr_cycif')])$total_cell_nr_cycif, breaks = 50)
# 
# ggplot(ct_frac_nolabel, aes(x=ct_nr_cycif, color = cell_type)) + 
#   geom_density() +
#   ylim(0, 0.15)
# 
# ############################
# # mixed_w_CD4 ROIs without community with CD4
# # 49 ROIs mixed_w_CD4, 10 wo label, 22 off-labelled, 17 with ok label
# table(roi_labels$roi_cluster_label_gmm)
# 
# CD4_off <- roi_labels[roi_labels$roi_cluster_label_gmm == 'mixed_w_CD4' & !grepl('CD4', roi_labels$community_cluster_label_manualnames_freq0.05)]
# CD4_off_comp <- roi_labels[roi_labels$roi_cluster_label_gmm == 'mixed_w_CD4' & !grepl('CD4', roi_labels$component_label_freq0.05)]
# 
# 
# CD4_off_wlab <- CD4_off[CD4_off$component_label != 'nolabel']
# 
# ct_frac_CD4off <- ct_frac_all[ct_frac_all$sample_roi %in% CD4_off_wlab$sample_roi]
# ct_frac_CD4off_long <- spread(ct_frac_CD4off[, c('sample_roi', 'cell_type', 'ct_nr_cycif')], key = cell_type, value = ct_nr_cycif)
