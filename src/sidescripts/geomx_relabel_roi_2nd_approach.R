library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(patchwork)
library(stringi)
library(PCAtools)
library(uwot)
library(NMF)
library(fpc)

# define variables --------------------------------------------------------
batch_name <- 'batch123'
proj_dir <<- '~/Documents/phd/st'
data_dir <<- '~/Documents/phd/st/data/geomx/batch123/' # batch1 2 and 3
anno_path <<- file.path(data_dir, 'metadata', 'dcc_metadata_batch123_no_tls_cleaned.csv') #batch1 and 2 and 3
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808') # batch123

source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))

out_dir <- file.path(output_dir, 'deconvolution', 'relabel-roi-deconv-dimred')
dir.create(out_dir, recursive = T)

# outputs from geomx_roi_hubs_integration.R
# TODO first part of that script handling deconv data may be moved here
ct_frac_deconv_path <- file.path(output_dir, 'cycif_integration', 'b123_ct_frac_deconv_bcells.csv')
ct_frac_deconv_roi_path <- file.path(output_dir, 'cycif_integration', 'b123_ct_frac_deconv_roi_bcells.csv')

# just to check - only available for b3tls for now
ct_frac_all_roi_path <- file.path(output_dir, 'cycif_integration', paste0('batch3tls_ct_frac_all_roi_bcells.csv'))

# thresholds for cluster labelling high+mid
hithr <- 0.25
midthr <- 0.15

# define names ------------------------------------------------------------

meta_names <- c('dcc_filename', 'Sample', 'NACT_status', 'Annotation_cell', 'Roi_geomx', 
                'Segment_geomx',  'Segment', 'tCycIF_preselection_initial_label') 

meta_names_per_roi <- c('Sample', 'Segment_geomx', 'NACT_status', 'Annotation_cell')

ct_names_all <- c("tumor", "Bcells", "Tcells_CD4", "Tcells_other", "Tcells_CD8", 
                  "Fibroblasts_Mesothelial", "Macrophages_Monocytes", "Mast_cells",
                  "NKcells", "Endothelial_cells", "DCs")

# main immune cells from deconv - also counted in cycif phenotyping
ct_names_immune <- c("Tcells_CD4", "Tcells_CD8", "DCs", "Macrophages_Monocytes", "Bcells", "NKcells") # no Bcells in basic phenotyping eg b1b2 
ct_names_myeloids <- c("DCs", "Macrophages_Monocytes")
ct_names_lymphoids <- c("Tcells_CD4", "Tcells_CD8")

# additional cells from deconv not counted in phenotyping and should be treated as 'other'
ct_names_other <- c("Tcells_other", "Mast_cells") # "Bcells" goes here when basic phenotyping b1b2

# vector to clean cluster labels - just short names for immune cells 
ct_names_toclean <- c(ct_names_immune, 'myeloids')
ct_names_cleaned <- c('CD4', 'CD8', 'CD11', 'Iba1', 'CD20', 'NK', 'Mye')

# load data ---------------------------------------------------------------

ct_frac_deconv <- as.data.frame(fread(ct_frac_deconv_path))
ct_frac_deconv_roi <- as.data.frame(fread(ct_frac_deconv_roi_path))

ct_frac_all_roi <- fread(ct_frac_all_roi_path)

metadt <- fread(anno_path, select = meta_names) %>%
  mutate(sample_roi = paste0(Sample, '_', Roi_geomx)) %>%
  filter(sample_roi %in% ct_frac_deconv_roi$sample_roi)

metadt_roi <- metadt %>%
  select(sample_roi, !!meta_names_per_roi) %>%
  distinct()
  
# clustered hmaps with ct fractions ---------------------------------------

deconv_names <- c('sd', 'bp')
ct_names_to_cluster_list <- list(all_immune = ct_names_immune, 
                                 cells4 = c(ct_names_myeloids, ct_names_lymphoids), 
                                 cells3 = c(ct_names_lymphoids, 'myeloids'),
                                 cells3b = c(ct_names_lymphoids, 'myeloids', 'Bcells'))

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

# calculate fractions of immune_others and play with nr of labels (NK, Bcells to others, DC+macro etc)

#TODO to rmv - just for testing
deconv_name <- 'sd'
hclust_cut <- 2 # 1,5
mye_add <- F
list_n <- 1


deconv_names <- c('sd', 'bp')

#TODO change it to more universal code to use with different names combination 
ct_names_to_cluster_list <- list(mye_lymph_b = c(ct_names_lymphoids, ct_names_myeloids, 'Bcells'), 
                                 mye_lymph = c(ct_names_lymphoids, ct_names_myeloids))

clust_list <- list()

for(deconv_name in deconv_names){
  for(list_n in 1:length(ct_names_to_cluster_list)){
    for(mye_add in c(T, F)){
      for(hclust_cut in c(1.5, 2)){
        
        ct_names_to_cluster <- unname(unlist(ct_names_to_cluster_list[list_n]))
        out_name <- names(ct_names_to_cluster_list)[list_n]
        out_name <- ifelse(mye_add, paste0(out_name, '_myemerged'), out_name)
        out_name <- paste0(deconv_name, '_', out_name, '_hcut', as.character(hclust_cut))

        print(out_name)
        
        # transform to wide + calculate fractions of immune_other cells
        ct_frac_roi_wide_immune <- ct_frac_deconv_roi[, c('sample_roi', 'cell_type', paste0('ct_frac_', deconv_name))] %>%
          pivot_wider(names_from = cell_type, values_from = !!paste0('ct_frac_', deconv_name)) %>%
          select(sample_roi, !!ct_names_immune, !!ct_names_other, immune_other) %>%
          mutate(across(c(ct_names_immune, ct_names_other), ~./immune_other)) %>%
          select(-immune_other) %>%
          column_to_rownames(var = 'sample_roi')
        
        # sum fractions of other cells
        ct_frac_roi_wide_immune <- ct_frac_roi_wide_immune %>%
          mutate(other_immune = rowSums(across(setdiff(colnames(ct_frac_roi_wide_immune), ct_names_to_cluster)))) %>%
          select(!!ct_names_to_cluster, other_immune)
        
        # if needed, sum fractions of myeloids
        if(mye_add){
          ct_frac_roi_wide_immune <- ct_frac_roi_wide_immune %>%
            mutate(myeloids = rowSums(across(ct_names_myeloids))) %>%
            select(-!!ct_names_myeloids)
          
          ct_names_to_cluster <- c(setdiff(ct_names_to_cluster, ct_names_myeloids), 'myeloids')
        }
        
        # clustering
        hclust_avg <- hclust(dist(ct_frac_roi_wide_immune), method = "ward.D2")
        plot(hclust_avg, hang = -1, cex = 0.3)
        hclust_avg_cut <- cutree(hclust_avg, h = hclust_cut) # 1.5 = 9 clusters, 2 = 5 clusters
        hclust_avg_cut_df <- as.data.frame(hclust_avg_cut) %>%
          rownames_to_column(var = "sample_roi") %>%
          rename(cluster = hclust_avg_cut)
        
        print(table(hclust_avg_cut_df$cluster))
        
        # transform to long, merge with clustering
        ct_frac_roi_long_immune_clust <- ct_frac_roi_wide_immune %>%
          rownames_to_column(var = "sample_roi") %>%
          pivot_longer(-sample_roi, names_to = "cell_type", values_to = "fraction_of_immune") %>%
          left_join(hclust_avg_cut_df) %>%
          arrange(cluster, sample_roi, cell_type, fraction_of_immune)
        
        # calculate mean per each cluster and compare across clusters to make labs
        ct_frac_clust_mean <- ct_frac_roi_long_immune_clust %>%
          group_by(cluster, cell_type) %>%
          summarise(mean_ct_frac_of_immune = mean(fraction_of_immune)) %>%
          ungroup()
        
        # create cluster labels based on hi and mid frequency thresholds
        ct_frac_clust_mean_wide <- ct_frac_clust_mean %>%
          pivot_wider(names_from = cell_type, values_from = mean_ct_frac_of_immune)

        ct_frac_clust_mean_wide$cluster_label <- apply(ct_frac_clust_mean_wide, 1, function(clust){
          
          lab <- sapply(ct_names_to_cluster, function(ct_name){
            if(round(as.numeric(clust[[ct_name]]), 2) >= hithr){
              return(paste0('hi', ct_name))
            } else if(round(as.numeric(clust[[ct_name]]), 2) >= midthr){
              return(paste0('mid', ct_name))
            } else{
              return('')
            }
          })
          
          lab <- paste(lab, collapse = '_')
          return(lab)
        })

        # clean labels
        ct_frac_clust_mean_wide$cluster_label <- stri_replace_all_regex(ct_frac_clust_mean_wide$cluster_label,
                               pattern=c(ct_names_toclean, '^_*|_*$', '_+'),
                               replacement=c(ct_names_cleaned, '', '_'),
                               vectorize=FALSE)

        # merge with ctfrac dfs
        ct_frac_roi_long_immune_clust <- left_join(ct_frac_roi_long_immune_clust, ct_frac_clust_mean_wide)
        ct_frac_roi_long_immune_clust$cluster_label_toplot <- paste0(ct_frac_roi_long_immune_clust$cluster, '_', ct_frac_roi_long_immune_clust$cluster_label)
        ct_frac_clust_mean <- left_join(ct_frac_clust_mean, ct_frac_clust_mean_wide)
        ct_frac_clust_mean$cluster_label_toplot <- paste0(ct_frac_clust_mean$cluster, '_', ct_frac_clust_mean$cluster_label)
        
        # stacked bar plot with all samples per cluster
        ggplot(ct_frac_roi_long_immune_clust, aes(x = sample_roi, y = fraction_of_immune, fill = cell_type)) +
          geom_bar(stat = "identity") +
          labs(title = paste0("ROI clustered with ", deconv_name, " cell type fractions of immune"), x = "ROIs", y = "ct_fraction") +
          theme_minimal() +
          theme(axis.text.x=element_blank()) +
          facet_wrap(~ cluster_label_toplot, scales = "free", ncol = 1)
        
        ggsave(file.path(out_dir, paste0('barplot_clust_immunefrac_', out_name,  '.png')))
        
        # stacked barplot for mean ct fraction per cluster
        ggplot(ct_frac_clust_mean, aes(x = cluster_label_toplot, y = mean_ct_frac_of_immune, fill = cell_type)) +
          geom_bar(stat = "identity") +
          labs(title = "mean ct fraction of immune per cluster", x = "ROI clusters", y = "mean ct fraction") +
          theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 6))
        
        ggsave(file.path(out_dir, paste0('barplot_mean_clust_immunefrac_', out_name, '.png')))
        
        # histograms for mean cluster ct fractins
        ggplot(data = ct_frac_clust_mean) +
          geom_histogram(aes(mean_ct_frac_of_immune, fill = cell_type)) +
          ylim(0, 4) +
          facet_wrap(~ cell_type)
        
        ggsave(file.path(out_dir, paste('hist_mean_clust_immunefrac_', out_name, '.png')))
        
        # return clusters for comparison
        clust_df <- ct_frac_roi_long_immune_clust %>%
          select(sample_roi, cluster, cluster_label) %>%
          distinct()
        colnames(clust_df) <- c('sample_roi', out_name, paste0(out_name, '_label'))
        clust_list <- append(clust_list, list(clust_df))
      }
    }
  }
}

clust_df_all <- do.call(cbind, clust_list) %>%
  column_to_rownames(var = 'sample_roi') %>%
  select(-contains('sample_roi')) %>%
  rownames_to_column('sample_roi')

fwrite(clust_df_all, file.path(out_dir, paste0('roi_ctfreq_clusters_df_','mid', as.character(midthr), '_hi', as.character(hithr), '_labs.csv')))


######################################################################
######################################################################
# Totally different approach:
# use different clustering methods and make clusters: Hclust, GMM, DBScan
# use dimentionality reduction methods: PCA, tSNE, UMAP, NMF
# visualise dim red in 2D and color by cluster

# clustering
deconv_name <- 'sd'
ct_names_to_cluster <- c(ct_names_lymphoids, ct_names_myeloids, 'Bcells')
mye_add <- F
out_name <- ifelse(mye_add, 'sd_mye_lymph_b_myemerged', 'sd_mye_lymph_b')
    
print(out_name)
    
# transform to wide + calculate fractions of immune_other cells
ct_frac_roi_wide_immune <- ct_frac_deconv_roi[, c('sample_roi', 'cell_type', paste0('ct_frac_', deconv_name))] %>%
  pivot_wider(names_from = cell_type, values_from = !!paste0('ct_frac_', deconv_name)) %>%
  select(sample_roi, !!ct_names_immune, !!ct_names_other, immune_other) %>%
  mutate(across(c(ct_names_immune, ct_names_other), ~./immune_other)) %>%
  select(-immune_other) %>%
  column_to_rownames(var = 'sample_roi')
    
# sum fractions of other cells
ct_frac_roi_wide_immune <- ct_frac_roi_wide_immune %>%
  mutate(other_immune = rowSums(across(setdiff(colnames(ct_frac_roi_wide_immune), ct_names_to_cluster)))) %>%
  select(!!ct_names_to_cluster, other_immune)
    
# if needed, sum fractions of myeloids
if(mye_add){
  ct_frac_roi_wide_immune <- ct_frac_roi_wide_immune %>%
    mutate(myeloids = rowSums(across(ct_names_myeloids))) %>%
    select(-!!ct_names_myeloids)
  
  ct_names_to_cluster <- c(setdiff(ct_names_to_cluster, ct_names_myeloids), 'myeloids')
}


################################################################
###############################################################
make_and_plot_dimreduction <- function(input_mtx, clusters_df, output_path, dimred_method = c('UMAP', 'PCA'),  id_name = 'sample_roi'){
  
  if(dimred_method == 'PCA'){
    dimred_df <- pca(mat = t(input_mtx))$rotated[, 1:2]
  } else if(dimred_method == 'UMAP'){
    dimred_df <- as.data.frame(umap(input_mtx, n_neighbors = 30))
  }
  
  dimred_df <- dimred_df %>%
    rownames_to_column(var = id_name) %>%
    left_join(clusters_df) 
  
  colnames(dimred_df) <- c(id_name, paste0(dimred_method, '1'), paste0(dimred_method, '2'), 'cluster')
  dimred_df$cluster <- as.character(dimred_df$cluster)
  
  ggplot(dimred_df, aes(x = get(paste0(dimred_method, '1')), y = get(paste0(dimred_method, '2')), color = cluster)) +
    geom_point() +
    labs(title = paste0(dimred_method, " across clusters")) +
    xlab(paste0(dimred_method, '1')) +
    ylab(paste0(dimred_method, '2')) +
    theme_minimal()
  
  ggsave(file.path(output_path))
  
}

####################################
####################################
install.packages(c("cluster", "factoextra"))  # Uncomment if not installed
library(cluster)
library(factoextra)
# clustering: hclust, DBscan GMM
ct_frac_mtx <- as.matrix(ct_frac_roi_wide_immune)

# DBScan
set.seed(220)
dbscan_out <- dbscan(ct_frac_roi_wide_immune, eps = 0.1, MinPts = 10, scale = F)
dbscan_df <- DataFrame(sample_roi = rownames(ct_frac_roi_wide_immune), cluster = dbscan_out$cluster)
table(dbscan_out$cluster)
plot(dbscan_out, ct_frac_roi_wide_immune, main = "DBScan")


##############################################
# iterate through different params for hclust
hclust_cuts = c(1, 1.5, 2, 2.5)

hclust_res <- lapply(hclust_cuts, function(hclust_cut){
  clust_name <- paste0('hclust_cut', as.character(hclust_cut))
  print(clust_name)
  
  hclust_avg <- hclust(dist(ct_frac_mtx), method = "ward.D2")
  hclust_avg_cut <- cutree(hclust_avg, h = hclust_cut) 
  cluster_df <- as.data.frame(hclust_avg_cut) %>%
    rownames_to_column(var = "sample_roi") %>%
    rename(cluster = hclust_avg_cut)
  
  print(table(cluster_df$cluster))
  
  # make dimreduction plots
  make_and_plot_dimreduction(ct_frac_mtx, cluster_df, dimred_method = 'PCA',
                             output_path = file.path(out_dir, paste0('scatter_', out_name,'_', clust_name, '_PCA.png')))
  
  make_and_plot_dimreduction(ct_frac_mtx, cluster_df, dimred_method = 'UMAP',
                             output_path = file.path(out_dir, paste0('scatter_', out_name,'_', clust_name, '_UMAP.png')))
  
  # calculate wcss (inertia)
  wcss <- sum(sapply(unique(cluster_df$cluster), function(cluster) {
    cluster_points <- ct_frac_mtx[cluster_df$cluster == cluster, ]
    return(sum(rowSums((scale(cluster_points) ^ 2))))  # WCSS for each cluster
  }))
  
  # calculate silhouettes
  silhouette_values <- silhouette(cluster_df$cluster, dist(ct_frac_mtx))
  avg_silhouette <- mean(silhouette_values[, 3])
  
  outp_df <- data.frame(method = clust_name, clusters_nr = length(unique(cluster_df$cluster)),
                        inertia = wcss, avg_silh = avg_silhouette)
  
  return(outp_df)
})

hclust_res <- do.call(rbind, hclust_res)


# do the elbow plot with inertia +avg silhouette
ggplot(data=hclust_res, aes(x=method)) +
  geom_line(aes(y=inertia, group=1), color = 'blue')+
  geom_point(aes(y=inertia), color = 'blue') + 
  labs(title = 'Hclust inertia vs cut param')

ggsave(file.path(out_dir, paste0(out_name,'_hclust_inertia.png')))

ggplot(data=hclust_res, aes(x=method)) +
  geom_line(aes(y=avg_silh, group=1), color = 'red')+
  geom_point(aes(y=avg_silh), color = 'red') +
  labs(title = 'Hclust avg silhouette vs cut param')

ggsave(file.path(out_dir, paste0(out_name,'_hclust_avg_silhouette.png')))

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



