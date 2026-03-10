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

# TODO ct_frac_cycif_roi 129 ROI, but 126 in the final ct_frac_all - probably without hub labs - to check
# TODO ensure # "S130_iOme_5" "S130_iOme_6" "S197_iOme_1" - everywhere in metadata - probably removed during QC
# TODO compare with deconv made per-ROI (not per-AOI)
# TODO copare with annotation cell from semi-automated method

# define vars -------------------------------------------------------------

# cell fraction from deconv below that lvl will be changed to 0 
# max nr of cells = 300 so 0.005 cell fraction is 1 cell/200 cells 1,5 cell/300 cells
min_frac <- 0.005
min_label_frac <- 0.05
cycif_main_ct_label <- 'consensus_label_clean_all_ct'
hubs_labels_list <- c('network_hub_type', 'community_cluster_label') # before with 'interaction_hub_type'

batch_name <- 'batch3tls'
proj_dir <<- '~/Documents/phd/st'
data_dir <<- '~/Documents/phd/st/data/geomx/batch123/' # batch1 2 and 3
anno_path <<- file.path(data_dir, 'metadata', 'dcc_metadata_batch123_no_tls_cleaned.xlsx') #batch1 and 2 and 3
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808') # batch123
geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS') 

bp_cellcounts_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', 'bp_res_mid_lvl_ct_updated_ct_fraction.csv')
sd_cellcounts_path <- file.path(output_dir, 'deconvolution', 'spatial_decon', 'sd_res_mid_lvl_ct_updated_geomxfiltpc_ct_fraction.csv')

# all cells within ROIs with hubs annotations computed with geomx_cycif_integration.R
#hubs_inroi_path <- file.path(output_dir, "cycif_integration", "batch2_hubs_cells_inroi_dt15171517_ct15_dt300_res0015_.csv")
#hubs_inroi_path <- file.path(output_dir, "cycif_integration", paste0(batch_name, "_hubs_cells_inroi_bcells_dt1517151715_ct15_dt300.csv"))
hubs_inroi_path <- file.path(output_dir, "cycif_integration", paste0(batch_name, "_hubs_cells_inroi_bcells_dt1515151717_ct10_dt300.csv"))

##############
# output files
source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))

outp_plot_dir <- file.path(output_dir, 'cycif_integration', 'ct_frac_comparison_b3tls_bcells_ct10')
dir.create(outp_plot_dir, recursive = T)
dir.create(file.path(outp_plot_dir,'hmaps'), recursive = T)

output_ct_frac_deconv_path <- file.path(output_dir, 'cycif_integration', 'b123_ct_frac_deconv_bcells.csv')
output_ct_frac_deconv_roi_path <- file.path(output_dir, 'cycif_integration', 'b123_ct_frac_deconv_roi_bcells.csv')

output_ct_frac_cycif_roi_path <- file.path(output_dir, 'cycif_integration',  paste0(batch_name, '_ct_frac_cycif_roi_bcells.csv'))
output_ct_frac_all_roi_path <- file.path(output_dir, 'cycif_integration', paste0(batch_name, '_ct_frac_all_roi_bcells.csv'))

##############
meta_names <- c('dcc_filename', 'Sample', 'Annotation_cell', 'Roi_geomx', 
                'Segment_geomx',  'Segment', 'tCycIF_preselection_initial_label') 

meta_names_per_roi <- c('Sample', 'Segment_geomx', 'Annotation_cell', 'tCycIF_preselection_initial_label')

ct_names_all <- c("tumor", "Bcells", "Tcells_CD4", "Tcells_other", "Tcells_CD8", 
                  "Fibroblasts_Mesothelial", "Macrophages_Monocytes", "Mast_cells",
                  "NKcells", "Endothelial_cells", "DCs")

# ct from deconv counted as stroma in cycif
ct_names_stroma <- c("Fibroblasts_Mesothelial", "Endothelial_cells")

# main immune cells from deconv - also counted in cycif phenotyping
ct_names_immune <- c("Tcells_CD4", "Tcells_CD8", "Macrophages_Monocytes", "NKcells", "DCs", "Bcells") # no Bcells in basic phenotyping eg b1b2 
ct_names_myeloids <- c("Macrophages_Monocytes", "DCs")
ct_names_lymphoids <- c("Tcells_CD4", "Tcells_CD8")

# additional cells from deconv not counted in phenotyping and should be treated as 'other'
ct_names_other <- c("Tcells_other", "Mast_cells") # "Bcells" goes here when basic phenotyping b1b2

# load geomx, merge with cleaned metadata ---------------------------------
# TODO run once again in 1811 with already cleaned metadata and just load meta from geomx
meta_geomx <- pData(readRDS(geomx_norm_batch_eff_rm_path))
meta_cleaned <- read_excel(anno_path)

metadt <- left_join(meta_geomx[, c('dcc_filename', 'Slide_Name')], meta_cleaned) # join ensuring order
metadt$sample_roi <- paste0(metadt$Sample, '_', metadt$Roi_geomx)

rm(meta_geomx)
rm(meta_cleaned)

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
hubs_cells_inroi$sample_roi <- gsub(' ', '', hubs_cells_inroi$sample_roi)

# rename cells to match deconvolution
hubs_cells_inroi$cell_type <-  mapvalues(hubs_cells_inroi[[cycif_main_ct_label]], 
                                         from = c("Macrophages", "CD4Tcells", "CD8Tcells",
                                                  "Tumor", "CD11c", "Undefined", "NK", 
                                                  "Stroma", "CD31+_Endothelial", "HEV+_Endothelial", "Bcells"),
                                         to=c("Macrophages_Monocytes", "Tcells_CD4", "Tcells_CD8",
                                              "tumor", "DCs", "other", "NKcells", "stroma", "stroma", "stroma", "Bcells"))

# count nr of cells per ROI and AOI
ct_frac_cycif_roi <- as.data.frame(dcast(hubs_cells_inroi, sample_roi ~ cell_type))
ct_frac_cycif_roi$total_cell_nr_cycif <- rowSums(ct_frac_cycif_roi[, -1])
ct_frac_cycif_roi$immune <- rowSums(ct_frac_cycif_roi[, intersect(ct_names_immune, colnames(ct_frac_cycif_roi))])
ct_frac_cycif_roi$immune_other <- rowSums(ct_frac_cycif_roi[, c(intersect(ct_names_immune, colnames(ct_frac_cycif_roi)), "other")])
ct_frac_cycif_roi$myeloids <- rowSums(ct_frac_cycif_roi[, ct_names_myeloids])
ct_frac_cycif_roi$lymphoids <- rowSums(ct_frac_cycif_roi[, ct_names_lymphoids])

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
  dplyr::mutate(across(hubs_labels_list, clean_labs)) %>%
  #dplyr::mutate(across(c(network_hub_type), clean_labs)) %>% # if 'cluster_N' as names - should't be claaned
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
  mutate(network_hub_type = paste0(unique(na.omit(network_hub_type)), collapse = "_")) %>%
  distinct() %>%
  mutate(across(hubs_labels_list[!hubs_labels_list == 'network_hub_type'], ~paste0(unique(na.omit(.)), collapse = "|"))) %>%
  distinct() 

# relabel after merging to ensure ordering
hubs_inroi_labs_frequent$network_hub_type <- clean_labs(hubs_inroi_labs_frequent$network_hub_type)

colnames(hubs_inroi_labs_frequent) <- c('sample_roi', paste0(hubs_labels_list, '_freq', as.character(min_label_frac)))

# weird to check:
# b2 S069_pOme_roi-002 - all network_hub_type within CD4/CD8 while community label is CD11_Iba1

##########################################################
##########################################################

# clean hubs labs ---------------------------------------------------------
# TODO a bit of repetition with the previous section

# remove cells without labels
hubs_inroi <- dplyr::select(hubs_cells_inroi, sample_roi, !!hubs_labels_list) %>%
  dplyr::mutate(across(hubs_labels_list[!hubs_labels_list == 'network_hub_type'], clean_labs)) %>%
  dplyr::filter(if_any(hubs_labels_list, ~!. == '')) %>%
  distinct()

# group by rois and merge labels
hubs_inroi <- hubs_inroi %>%
  group_by(sample_roi) %>%
  mutate(network_hub_type = paste0(unique(na.omit(network_hub_type)), collapse = "_")) %>%
  distinct() %>%
  mutate(across(hubs_labels_list[!hubs_labels_list == 'network_hub_type'], ~paste0(unique(na.omit(.)), collapse = "|"))) %>%
  distinct() 

hubs_inroi$network_hub_type <- clean_labs(hubs_inroi$network_hub_type) # relabel after merging to ensure ordering

# merge all information together ------------------------------------------

# immune_other = all immune + 'other'

# join cycif and deconv ct fraction tables, add metadata and hubs info
ct_frac_all <- left_join(ct_frac_cycif_long_roi, ct_frac_deconv_long_roi, by = c('sample_roi', 'cell_type')) %>%
  left_join(distinct(metadt[, c(meta_names_per_roi, 'sample_roi')]), by = 'sample_roi') %>% # before metadt_filt
  left_join(hubs_inroi, by = 'sample_roi') %>%
  left_join(hubs_inroi_labs_frequent, by = 'sample_roi')

ct_frac_all[ct_frac_all == ""] <- NA
fwrite(ct_frac_all, output_ct_frac_all_roi_path)

# compare cell fractions --------------------------------------------------

# TODO just to compare 2 versions of community labels
# ct_frac_all300 <- fread(file.path(output_dir, 'cycif_integration', paste0(batch_name, '_ct_frac_all_roi_dt300.csv')))
# ct_frac_all500 <- fread(file.path(output_dir, 'cycif_integration', paste0(batch_name, '_ct_frac_all_roi_dt500.csv')))
# 
# ct_frac_all300 <- dplyr::rename(ct_frac_all300, community_cluster_300 = community_cluster_label, 
#                          community_cluster_freq0.05_300 = community_cluster_label_freq0.05)
# ct_frac_all500 <- dplyr::rename(ct_frac_all500, community_cluster_500 = community_cluster_label, 
#                                 community_cluster_freq0.05_500 = community_cluster_label_freq0.05)
# 
# ct_frac_all <- left_join(ct_frac_all300, 
#                          ct_frac_all500[, c('sample_roi','cell_type', 'community_cluster_500', 'community_cluster_freq0.05_500')],
#                          by = c('sample_roi', 'cell_type'))
# 
# ct_frac_all <- left_join(ct_frac_all, unique(metadt[, c('sample_roi', 'tCycIF_preselection_initial_label')])) %>%
#   mutate(tCycIF_preselection_initial_label = clean_labs(tCycIF_preselection_initial_label))
# 
# #label_names <- c('network_hub_type', 'community_cluster')
# 
# label_vars <- c("Annotation_cell", "interaction_hub_type", "network_hub_type",
#                 "community_cluster_300", paste0("interaction_hub_type_freq", as.character(min_label_frac)),
#                 paste0("network_hub_type_freq", as.character(min_label_frac)),
#                 paste0("community_cluster_freq", as.character(min_label_frac), "_300"), "community_cluster_500",
#                 paste0("community_cluster_freq", as.character(min_label_frac), "_500"),
#                 'tCycIF_preselection_initial_label')

ct_frac_all$tCycIF_preselection_initial_label_cleaned <- clean_labs(ct_frac_all$tCycIF_preselection_initial_label)
ct_frac_all$tCycIF_preselection_initial_label_cleaned <- ifelse(ct_frac_all$tCycIF_preselection_initial_label_cleaned == '', 'otherlabel',
                                                                ct_frac_all$tCycIF_preselection_initial_label_cleaned)

label_vars <- c("Annotation_cell", hubs_labels_list, paste0(hubs_labels_list, '_freq', as.character(min_label_frac)),
                "tCycIF_preselection_initial_label_cleaned")

##########################################################
# scatterplot with geomx vs cycif total cell count
cell_count_roi <- select(ct_frac_all, sample_roi, Segment_geomx, total_cell_nr_cycif, total_cell_nr_geomx) %>%
  distinct()

cell_count_scatter <- ggplot(data = cell_count_roi, aes(x = total_cell_nr_cycif, y = total_cell_nr_geomx)) +
  geom_point(aes(color = Segment_geomx)) +
  ggtitle('total cell count per ROI cycif vs geomx') +
  geom_smooth(method='lm', formula= y~x) +
  stat_correlation(method = 'pearson')

ggsave(file.path(outp_plot_dir, paste0('total_cellnr_cycif_vs_geomx.png')),
       width = 2000, height = 2000, unit = 'px')


##########################################################
# compare ct fractions between deconv sd/bp and cycif phenotyping
comp_type <- 'ct_frac' # or ct_nr

# scatterplots with value comparisons between methods
for(value_comb in c('bp_sd', 'bp_cycif', 'sd_cycif')){
  vals <- unlist(strsplit(value_comb, split = '_'))
  
  print(value_comb)
  # per cell type
  for(color_var in c('Annotation_cell')){
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

      ggsave(file.path(outp_plot_dir, paste0('scatter_', ct_name, '_', vals[1], '_', vals[2], '_', comp_type, '_', color_var, '.png')),
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
  
  ggsave(file.path(outp_plot_dir, paste0('scatter_all_', vals[1], '_', vals[2], '_', comp_type, '.png')),
         width = 2000, height = 2000, unit = 'px')
}

#######################################################################
# boxplots with cell nr/fractions per different labels

# make long dataframe for sd/bp/cycif methods
ct_frac_all_method <- select(ct_frac_all, sample_roi, Segment_geomx, cell_type, !!label_vars, starts_with('ct_frac'), starts_with('ct_nr'))
ct_frac_all_method_long <- melt(setDT(ct_frac_all_method), id.vars = c('sample_roi', 'Segment_geomx', 'cell_type', label_vars),
                                variable.name = "method_type")

# for each ct faceted by method
for(comp_type in c('ct_frac', 'ct_nr')){
  ct_frac_all_method_comp <- ct_frac_all_method_long[grepl(comp_type, ct_frac_all_method_long$method_type), ]
  
  for(label_var in label_vars){
    for(ct_name in unique(ct_frac_all$cell_type)){
      
      ct_frac_all_method_comp_ct <- ct_frac_all_method_comp[ct_frac_all_method_comp$cell_type == ct_name, ]
      
      ggplot(ct_frac_all_method_comp_ct, aes(x = get(label_var), y = value)) + 
        geom_boxplot(alpha = .2) +
        geom_point(size = 0.2) + 
        facet_wrap(~ method_type, nrow = length(unique(ct_frac_all_method_comp_ct$method_type))) +
        labs(title = ct_name, x = label_var, y = comp_type) +
        scale_x_discrete(guide = guide_axis(angle = 45))
      
      
      ggsave(file.path(outp_plot_dir, paste0('boxpl_', ct_name, '_', label_var, '_', comp_type, '.png')),
             width = 1500, height = 2000, unit = 'px')
    }
  }
}

# for each method, faceted by label
for(method_name in unique(ct_frac_all_method_long$method_type)){
  ct_frac_all_method_sel <- ct_frac_all_method_long[ct_frac_all_method_long$method_type == method_name, ]
  
  for(label_var in label_vars){
    ggplot(ct_frac_all_method_sel, aes(x = cell_type, y = value)) + 
      geom_boxplot(alpha = .2) +
      geom_point(size = 0.2) + 
      facet_wrap(~ get(label_var), ncol = 4) +
      labs(title = paste(method_name, 'per', label_var), x = 'cell_type', y = method_name) +
      scale_x_discrete(guide = guide_axis(angle = 90))
    
    
    ggsave(file.path(outp_plot_dir, paste0('label_boxpl_allct_', label_var, '_', method_name, '.png')),
           width = 1500, height = 3000, unit = 'px')
    
    # filtered to immune
    ggplot(ct_frac_all_method_sel[ct_frac_all_method_sel$cell_type %in% ct_names_immune, ], aes(x = cell_type, y = value)) + 
      geom_boxplot(alpha = .2) +
      geom_point(size = 0.2) + 
      facet_wrap(~ get(label_var), ncol = 4) +
      labs(title = paste(method_name, 'per', label_var), x = 'cell_type', y = method_name) +
      scale_x_discrete(guide = guide_axis(angle = 90))
    
    
    ggsave(file.path(outp_plot_dir, paste0('label_boxpl_immune_', label_var, '_', method_name, '.png')),
           width = 1500, height = 3000, unit = 'px')
  }
}


# compare labels ----------------------------------------------------------

labs_all <- ct_frac_all %>%
  select(sample_roi, !!label_vars) %>%
  distinct()

labs_all[labs_all==""]<- "nolabel"
labs_all[is.na(labs_all)]<- "nolabel"


lab_name1 <- 'community_cluster_label_freq0.05'
lab_name2 <- 'network_hub_type_freq0.05'

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

