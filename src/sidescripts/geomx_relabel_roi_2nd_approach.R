library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(patchwork)

# TODO barplots + clustering as for communities

# define variables --------------------------------------------------------
batch_name <- 'batch123'
proj_dir <<- '~/Documents/phd/st'
data_dir <<- '~/Documents/phd/st/data/geomx/batch123/' # batch1 2 and 3
anno_path <<- file.path(data_dir, 'metadata', 'dcc_metadata_batch123_no_tls_cleaned.xlsx') #batch1 and 2 and 3
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808') # batch123

source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))

out_dir <- file.path(output_dir, 'deconvolution', 'relabel-roi-deconv')
dir.create(out_dir, recursive = T)

# outputs from geomx_roi_hubs_integration.R
# TODO first part of that script handling deconv data may be moved here
ct_frac_deconv_path <- file.path(output_dir, 'cycif_integration', 'b123_ct_frac_deconv_bcells.csv')
ct_frac_deconv_roi_path <- file.path(output_dir, 'cycif_integration', 'b123_ct_frac_deconv_roi_bcells.csv')

# just to check - only available for b3tls for now
ct_frac_all_roi_path <- file.path(output_dir, 'cycif_integration', paste0('batch3tls_ct_frac_all_roi_bcells.csv'))


# define names ------------------------------------------------------------

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

# load data ---------------------------------------------------------------

ct_frac_deconv <- as.data.frame(fread(ct_frac_deconv_path))
ct_frac_deconv_roi <- as.data.frame(fread(ct_frac_deconv_roi_path))

ct_frac_all_roi <- fread(ct_frac_all_roi_path)


# plot ct fractions density distributions ---------------------------------



# clustered hmaps with ct fractions ---------------------------------------

deconv_names <- c('sd', 'bp')
ct_names_to_cluster_list <- list(all_immune = ct_names_immune, 
                                 cells4 = c(ct_names_myeloids, ct_names_lymphoids), 
                                 cells3 = c(ct_names_lymphoids, 'myeloids'))

# basic hmap with raw ct frequencies
for(deconv_name in deconv_names){
  for(list_n in 1:length(ct_names_to_cluster_list)){
    ct_names_to_cluster <- unname(unlist(ct_names_to_cluster_list[list_n]))
    out_name <- names(ct_names_to_cluster_list)[list_n]
    
    ct_frac_roi_wide <- ct_frac_deconv_roi[, c('sample_roi', 'cell_type', paste0('ct_frac_', deconv_name))] %>%
      pivot_wider(names_from = cell_type, values_from = !!paste0('ct_frac_', deconv_name)) %>%
      select(sample_roi, !!ct_names_to_cluster) %>%
      column_to_rownames(var = 'sample_roi')
    
    
    # make hmap
    png(filename=file.path(out_dir, paste0('hmap_roi_', out_name, '_', deconv_name, '.png')), 
        width=10, height=6,units="in",res=1000)
    
    ind_heat <- Heatmap(as.matrix(ct_frac_roi_wide), 
                        cluster_columns = T, 
                        cluster_rows= T,
                        show_row_names = F, 
                        show_column_names = TRUE,
                        show_heatmap_legend = TRUE)
    
    
    draw(ind_heat, heatmap_legend_side = "right")
    
    dev.off()
  }
}

# hmap with ct fractions per immune
for(deconv_name in deconv_names){
  for(list_n in 1:length(ct_names_to_cluster_list)){
    ct_names_to_cluster <- unname(unlist(ct_names_to_cluster_list[list_n]))
    out_name <- names(ct_names_to_cluster_list)[list_n]

    ct_frac_roi_wide_immune <- ct_frac_deconv_roi[, c('sample_roi', 'cell_type', paste0('ct_frac_', deconv_name))] %>%
      pivot_wider(names_from = cell_type, values_from = !!paste0('ct_frac_', deconv_name)) %>%
      select(sample_roi, !!ct_names_immune, myeloids, lymphoids, immune) %>%
      mutate(across(ct_names_immune, ~./immune)) %>% # ct fraction as fraction of immune cells
      select(sample_roi, !!ct_names_to_cluster) %>%
      column_to_rownames(var = 'sample_roi')

    # make hmap
    png(filename=file.path(out_dir, paste0('hmap_roi_', out_name, '_', deconv_name, '_immunefrac.png')), 
        width=10, height=6,units="in",res=1000)
    
    ind_heat <- Heatmap(as.matrix(ct_frac_roi_wide_immune), 
                        cluster_columns = T, 
                        cluster_rows= T,
                        show_row_names = F, 
                        show_column_names = TRUE,
                        show_heatmap_legend = TRUE)
    
    
    draw(ind_heat, heatmap_legend_side = "right")
    
    dev.off()
  }
}


# stacked barplots with ct fractions + clustering -------------------------

deconv_name <- 'sd'

# transform to wide + calculate fractions of immune cells
ct_frac_roi_wide_immune <- ct_frac_deconv_roi[, c('sample_roi', 'cell_type', paste0('ct_frac_', deconv_name))] %>%
  pivot_wider(names_from = cell_type, values_from = !!paste0('ct_frac_', deconv_name)) %>%
  select(sample_roi, !!ct_names_immune, myeloids, lymphoids, immune) %>%
  mutate(across(ct_names_immune, ~./immune)) %>%
  select(sample_roi, !!ct_names_immune) %>%
  column_to_rownames(var = 'sample_roi')

# clustering
hclust_avg <- hclust(dist(ct_frac_roi_wide_immune), method = "complete")
plot(hclust_avg, hang = -1, cex = 0.3)
hclust_avg_cut <- cutree(hclust_avg, h = 0.7)
hclust_avg_cut_df <- as.data.frame(hclust_avg_cut) %>%
  rownames_to_column(var = "sample_roi") %>%
  rename(cluster = hclust_avg_cut)

table(hclust_avg_cut_df$cluster)

# transform to long, merge with clustering
ct_frac_roi_long_immune <- ct_frac_roi_wide_immune %>%
  rownames_to_column(var = "sample_roi") %>%
  pivot_longer(-sample_roi, names_to = "cell_type", values_to = "fraction_of_immune") %>%
  left_join(hclust_avg_cut_df) %>%
  arrange(cluster, sample_roi, cell_type, fraction_of_immune)
  

# TODO correctly arrange bars per cluster + add labels for cluster
# stacked bar plot
ggplot(ct_frac_roi_long_immune, aes(x = sample_roi, y = fraction_of_immune, fill = cell_type)) +
  geom_bar(stat = "identity") +
  labs(title = "ROI clustered with cell type fractions of immune", x = "Samples", y = "Counts") +
  theme_minimal() +
  #scale_x_discrete(labels=ct_frac_roi_long_immune$cluster)
  theme(axis.text.x=element_blank()) +
  facet_grid(~ cluster)

plot_list <- lapply(unique(ct_frac_roi_long_immune$cluster), function(cluster_name) {
  
    cluster_name <- as.character(cluster_name)
    dt_cluster <- ct_frac_roi_long_immune[ct_frac_roi_long_immune$cluster == cluster_name, ]
    
    ggplot(dt_cluster, aes(x = sample_roi, y = fraction_of_immune, fill = cell_type)) +
      geom_bar(stat = "identity") +
      labs(title = paste("cluster ", cluster_name),
           x = "ROIs", 
           y = "ct freq") +
      theme(axis.text.x=element_blank())
    
    ggsave(file.path(out_dir, paste0('barplot_cluster_',cluster_name, '_', deconv_name, '_immunefrac.png')))
  })

# Step 4: Combine all plots into one image
combined_plot <- wrap_plots(plot_list) + plot_layout(ncol = 2)

ggsave(file.path(out_dir, paste0('barplot_clustered_', deconv_name, '_immunefrac.png')))

#TODO calculate mean per each cluster and compare across clusters to make labs

#######################################################################
#######################################################################

# # check ct fractions distribution
# for(ct_name in cell_types_immune){
#   roi_ct_long_ct <- roi_ct_frac_immune_long[roi_ct_frac_immune_long$cell_type == ct_name, ]
#   #roi_ct_frac_immune_long_ct <- roi_ct_frac_immune_long[roi_ct_frac_immune_long$cell_type == ct_name, ]
#   #roi_ct_frac_long_ct <- roi_ct_frac_long[roi_ct_frac_long$cell_type == ct_name, ]
#   #aoi_ct_frac_long_ct <- aoi_ct_frac_long[aoi_ct_frac_long$cell_type == ct_name, ]
#   
#   ggplot(data = roi_ct_long_ct) +
#     #geom_density(aes(value)) +
#     geom_histogram(aes(value), bins = 100)
#     ggtitle(ct_name)
#   
#   ggsave(file.path(output_dir, paste0('hist_fraq_immune_', ct_name, '.png')), width = 2000, height = 1000, unit='px')
#   
#   # ggplot(data = aoi_ct_frac_long_ct) +
#   #   geom_density(aes(value, color = Segment)) +
#   #   ggtitle(ct_name)
#   # 
#   # ggsave(file.path(output_dir, paste0('density_fraq_aoi_', ct_name, '.png')), width = 2000, height = 1000, unit='px')
#   # 
# }
# 
# ##############################################
# # cluster ROIs per cell fraction in a hmap
# 
# ct_frac_mtx <- as.matrix(column_to_rownames(roi_ct_frac_immune[, !(names(roi_ct_frac_immune) %in% c('total_cell_nr', 'immune_cell_nr'))],
#                                             'sample_roi'))
# 
# # TODO comment/uncomment for important cells
# ct_frac_mtx <- ct_frac_mtx[, c(cell_types_important)]
# ct_frac_mtx_zscore <- scale(ct_frac_mtx) # zscore by column
# 
# # cluster by hclust
# ct_frac_hclust <- hclust(dist(ct_frac_mtx), method = "average")
# plot(ct_frac_hclust, hang = -1, cex = 0.4)
# ct_frac_hclust_cut <- cutree(ct_frac_hclust, h = 0.3)
# 
# ct_frac_zscore_hclust <- hclust(dist(ct_frac_mtx_zscore), method = "average")
# plot(ct_frac_zscore_hclust, hang = -1, cex = 0.4)
# ct_frac_zscore_hclust_cut <- cutree(ct_frac_zscore_hclust, h = 2)
# ct_frac_zscore_hclust_cut_k <- cutree(ct_frac_zscore_hclust, k = 11)
# 
# #########
# # TODO annotations on hmaps are wrong - only match zscores
# # make hmaps
# ha = HeatmapAnnotation(
#   #ct_label = anno_simple(roi_ct_frac$Annotation_cell),
#   hclust_h = anno_simple(as.character(unname(ct_frac_zscore_hclust_cut))),
#   hclust_k = anno_simple(as.character(unname(ct_frac_zscore_hclust_cut_k))),
#   which = "row", show_legend = TRUE)
# 
# # hmap for ct fraq
# png(filename=file.path(output_dir, paste0('hmap_roi_fraq_immune_important_ct.png')), 
#     width=10, height=6,units="in",res=2000)
# 
# ind_heat <- Heatmap(ct_frac_mtx, cluster_columns = F, cluster_rows= ct_frac_hclust,
#                     show_row_names = TRUE, show_column_names = TRUE,
#                     left_annotation = ha, show_heatmap_legend = TRUE)
# 
# 
# draw(ind_heat, annotation_legend_side = "right", heatmap_legend_side = "right")
# dev.off()
# 
# # hmap for zscore
# png(filename=file.path(output_dir, paste0('hmap_roi_zscore_fraq_immune_important_ct.png')), 
#     width=10, height=6,units="in",res=2000)
# 
# ind_heat <- Heatmap(ct_frac_mtx_zscore, cluster_columns = F, cluster_rows= ct_frac_zscore_hclust,
#                     show_row_names = TRUE, show_column_names = TRUE,
#                     left_annotation = ha, show_heatmap_legend = TRUE)
# 
# 
# draw(ind_heat, annotation_legend_side = "right", heatmap_legend_side = "right")
# dev.off()

# TODO the same for geomx_segment
# TODO compare with deconvoluted fractions
# TODO compare with our labels
# TODO add fractions (from all cells) for cell types (eg macro fraq + dc frac and then: check distrib, label highest ones)



