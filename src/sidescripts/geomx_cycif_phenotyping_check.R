library(data.table)
library(plyr)
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(ggplot2)
library(tibble)
library(GeomxTools)
library(readxl)
library(tools)
library(tidyr)
library(ggpmisc)
library(raster)

# define paths ------------------------------------------------------------

um_to_pix_ratio <- 1 #0.325 for b2, 0.65 for b3TLS 1 if already in pix

# from master script
proj_dir <<- '~/Documents/phd/st'
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808') # batch123
eyemt_pdrive_dir <- "/home/ad/P-drive/h30492/farkkilab2/9_EyeMT"

metadt_path <- file.path(output_dir, 'metadata_full_SENSITIVE.csv')
geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS') 

deconv_ct_count_path <- file.path(output_dir, 'deconvolution', 'ct_frac_deconv_roi_mid_lvl_ct_updated_4mainimmune.csv')

# phenotyped cells dirs
phenotyped_cells_dirs_paths <- c(file.path(eyemt_pdrive_dir, 'Data/cycif/batch1_adjacent_slides/phenotyped_cells/tribus'),
                                 file.path(eyemt_pdrive_dir, 'Data/cycif/batch2_adjacent_slides/phenotyped_cells/tribus'),
                                 file.path(eyemt_pdrive_dir, 'Data/cycif/batch3_adjacent_slides/phenotyped_cells/tribus'))


out_path_hubs_cells_inroi <- file.path(output_dir, "cycif_integration","batch123_tribus") #  "batch3tls_tribus_and_manualgating"
outp_plot_dir <- file.path(out_path_hubs_cells_inroi, 'plots_newpolygons')
dir.create(outp_plot_dir, recursive = T)


# define cell types -------------------------------------------------------

ct_names_immune <- c("Tcells_CD4", "Tcells_CD8", "DCs", "Macrophages_Monocytes") # all immune ct from phenotyping
ct_names_myeloids <- c("Macrophages_Monocytes", "DCs")
ct_names_lymphoids <- c("Tcells_CD4", "Tcells_CD8")
# additional cells from deconv not counted in phenotyping and should be treated as 'other'
ct_names_other <- c("Tcells_other", "Mast_cells", "NKcells", "Bcells") # "Bcells" goes here when basic phenotyping b1b2 and 'consensus_label_clean'

meta_names_per_roi <- c('Sample', 'Segment_geomx', 'Annotation_cell', "roi_cluster_label_gmm", "roi_cluster_label_hclust")

# load metadata -----------------------------------------------------------

geomx_dcc <- colnames(readRDS(geomx_norm_batch_eff_rm_path))
metadt <- as.data.frame(fread(metadt_path))
metadt <- metadt[metadt$dcc_filename %in% geomx_dcc, ]
rownames(metadt) <- NULL

# fix samplenames in batch1
metadt$Sample_fixed <- ifelse(metadt$main_batch_nr == 1, 
                              mapvalues(metadt$Sample, 
                                        from = c('S015_post', 'S015_pre', 'S027_post', 'S027_pre', 'S032_post', 'S032_pre',
                                                 'S053_post', 'S057_post', 'S065_post', 'S072_post', 'S073_post', 'S076_post',
                                                 'S084_post','S084_pre', 'S139_post', 'S139_pre'),
                                        to = c('S015_iOme', 'S015_pPer', 'S027_iOme', 'S027_pOme', 'S032_iOval', 'S032_pOme',
                                               'S053_iOme', 'S057_iOme', 'S065_iOme', 'S072_iOme', 'S073_iOme', 'S076_iOme',
                                               'S084_iOme','S084_pAdn', 'S139_iOme', 'S139_pPer')),
                              metadt$Sample)


# iterate through samples and select cells within ROIs --------------------

pheno_dir <- phenotyped_cells_dirs_paths[1]
sample_name <- 'S053_iOme'
  
# iterate through batches dirs
pheno_cells_inroi_all <- lapply(phenotyped_cells_dirs_paths, function(pheno_dir){
  batchnr <- as.numeric(substr(gsub('.*batch', '', pheno_dir), 1, 1))
  batch_pheno_all_files <- list.files(pheno_dir, pattern = "tribus_annotated", full.names = T)
  
  metadt_batch <- metadt[metadt$main_batch_nr == batchnr, ]
  
  #iterate through samples in batch
  lapply(unique(metadt_batch$Sample_fixed), function(sample_name){
    
    print(sample_name)
    pheno_cells_sample_path <- batch_pheno_all_files[grepl(sample_name, batch_pheno_all_files)]
    print(pheno_cells_sample_path)
    
    if(length(pheno_cells_sample_path) == 0){
      return(NULL)
    } else{
      pheno_cells_sample <- fread(pheno_cells_sample_path)
      
      # convert to pix
      pheno_cells_sample$X_centroid_px <- pheno_cells_sample$X_centroid / um_to_pix_ratio
      pheno_cells_sample$Y_centroid_px <- pheno_cells_sample$Y_centroid / um_to_pix_ratio

      # count cells within all ROIs in the sample - now in metadata
      roi_coords_sample <- metadt_batch[metadt_batch$Sample_fixed == sample_name, c("Sample","Sample_fixed", "main_batch_nr", "sample_roi", "Roi_geomx",
                                                                                    "roi_c1_X_cycif", "roi_c1_Y_cycif", "roi_c2_X_cycif", "roi_c2_Y_cycif",
                                                                                    "roi_c3_X_cycif", "roi_c3_Y_cycif", "roi_c4_X_cycif", "roi_c4_Y_cycif")]
      
      roi_coords_sample <- distinct(roi_coords_sample)
      
      cells_in_roi_all <- apply(roi_coords_sample, 1, function(row){
        
        #######################################
        # find cells within range
        # coordinates are not longer rectangles, they're a bit rotated - the cells are found inside longer edges of rectangle
        # cells_in_roi <- dplyr::filter(pheno_cells_sample,
        #                               as.numeric(X_centroid_px) >= min(as.numeric(row[['roi_c1_X_cycif']]), as.numeric(row[['roi_c4_X_cycif']])) &
        #                                 as.numeric(X_centroid_px) <= max(as.numeric(row[['roi_c2_X_cycif']]), as.numeric(row[['roi_c3_X_cycif']])) &
        #                                 as.numeric(Y_centroid_px) >= min(as.numeric(row[['roi_c4_Y_cycif']]), as.numeric(row[['roi_c3_Y_cycif']])) &
        #                                 as.numeric(Y_centroid_px) <= max(as.numeric(row[['roi_c1_Y_cycif']]), as.numeric(row[['roi_c2_Y_cycif']])))
        # 
        # 
        # cells_in_roi <- cbind(cells_in_roi, as.data.frame(lapply(row, rep, nrow(cells_in_roi))))
        #######################################
        # new finding cells implementation rewritten from python shapely package to sp R package
        poly_x <- c(row[['roi_c1_X_cycif']], row[['roi_c2_X_cycif']], row[['roi_c3_X_cycif']], row[['roi_c4_X_cycif']])
        poly_y <- c(row[['roi_c1_Y_cycif']], row[['roi_c2_Y_cycif']], row[['roi_c3_Y_cycif']], row[['roi_c4_Y_cycif']])

        cells_in_roi <- dplyr::filter(pheno_cells_sample, sp::point.in.polygon(point.x = as.numeric(X_centroid_px),
                                                                                point.y = as.numeric(Y_centroid_px), 
                                                                                pol.x = poly_x, pol.y = poly_y) != 0)
        
        cells_in_roi <- cbind(cells_in_roi, as.data.frame(lapply(row, rep, nrow(cells_in_roi))))
        #######################################
        return(cells_in_roi)
      })
      
      cells_in_roi_all <- do.call(rbind, cells_in_roi_all)
      print(nrow(cells_in_roi_all))
      table(cells_in_roi_all$Roi_geomx)
      return(cells_in_roi_all)
    }
  })
})

pheno_cells_inroi_all <- unlist(pheno_cells_inroi_all, recursive = FALSE)
pheno_cells_inroi_all_df <- do.call("rbind.fill", pheno_cells_inroi_all)

# clean labels ------------------------------------------------------------

pheno_cells_inroi_all_df$cell_type <- ifelse(grepl('undefined|other', pheno_cells_inroi_all_df$final_label), 'other', pheno_cells_inroi_all_df$final_label)
pheno_cells_inroi_all_df$cell_type <- gsub('^Tumor_|^Stroma_|^Immune_', '', pheno_cells_inroi_all_df$cell_type)
pheno_cells_inroi_all_df$cell_type <- mapvalues(pheno_cells_inroi_all_df$cell_type, 
                                                        from = c("Tumor", "Dcs", "CD8_Tcells", "Stroma", "Macrophages", "CD4_Tcells"),
                                                        to = c("tumor", "DCs", "Tcells_CD8", "stroma", "Macrophages_Monocytes", "Tcells_CD4"))

fwrite(pheno_cells_inroi_all_df, file.path(out_path_hubs_cells_inroi, 'batch123_cells_in_roi_new_polygons.csv'))

#########################################
# to check b3TLS 
# pheno_cells_inroi_all_df <- fread('/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2808/cycif_integration/batch3tls_hubs_cells_inroi_combined_myeloids_min20cells_nomix_dt171715_ct10_dt300.csv')
# pheno_cells_inroi_all_df$cell_type <- ifelse(grepl('undefined|other|Undefined', pheno_cells_inroi_all_df$consensus_label_clean), 'other', pheno_cells_inroi_all_df$consensus_label_clean)
# pheno_cells_inroi_all_df$cell_type <- mapvalues(pheno_cells_inroi_all_df$cell_type, 
#                                                 from = c("Tumor", "Dcs", "CD8_Tcells", "Stroma", "Macrophages", "CD4_Tcells"),
#                                                 to = c("tumor", "DCs", "Tcells_CD8", "stroma", "Macrophages_Monocytes", "Tcells_CD4"))

##########################################

# count cells per roi -----------------------------------------------------

ct_frac_cycif_roi <- as.data.frame(dcast(pheno_cells_inroi_all_df[, c('sample_roi', 'cell_type')], sample_roi ~ cell_type))
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

fwrite(ct_frac_cycif_long_roi, file.path(out_path_hubs_cells_inroi,  'batch123_ct_frac_cycif_roi_new_polygons.csv'))


# merge with cells counted from deconv ------------------------------------

ct_frac_deconv_long_roi <- fread(deconv_ct_count_path)

ct_frac_long_roi <- left_join(ct_frac_cycif_long_roi, ct_frac_deconv_long_roi, by = c('sample_roi', 'cell_type'))
ct_frac_long_roi <- left_join(ct_frac_long_roi, distinct(metadt[, c('sample_roi', meta_names_per_roi)]))

# remove all ROIs with total nr of cells < 100
ct_frac_long_roi <- ct_frac_long_roi[ct_frac_long_roi$total_cell_nr_cycif > 100, ]

# make plots with ct comparisons ------------------------------------------

# scatterplot with geomx vs cycif total cell count
cell_count_roi <- dplyr::select(ct_frac_long_roi, sample_roi, Segment_geomx, total_cell_nr_cycif, total_cell_nr_geomx) %>%
  distinct()

cell_count_scatter <- ggplot(data = cell_count_roi, aes(x = total_cell_nr_cycif, y = total_cell_nr_geomx)) +
  geom_point(aes(color = Segment_geomx)) +
  ggtitle('total cell count per ROI cycif vs geomx') +
  geom_smooth(method='lm', formula= y~x) +
  stat_correlation(method = 'pearson')

ggsave(file.path(outp_plot_dir, paste0('total_cellnr_cycif_vs_geomx.png')),
       width = 2000, height = 2000, unit = 'px')

#######################################################
# compare ct fractions between deconv sd/bp and cycif phenotyping
label_vars_forplots <- c("roi_cluster_label_gmm")

#TODO just to keep the var name
ct_frac_all <- ct_frac_long_roi

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
      
      ggsave(file.path(outp_plot_dir, paste0('scatter_', ct_name, '_', vals[1], '_', vals[2], '_', comp_type, '_', color_var, '.png')),
             width = 2000, height = 2000, unit = 'px')
      
      ct_scatter_by_sample <- ggplot(data = ct_frac_ct, aes(x = get(paste0(comp_type, '_', vals[1])), y = get(paste0(comp_type, '_', vals[2])))) +
        geom_point(aes(color = get(color_var), shape = Segment_geomx)) +
        geom_smooth(method='lm', formula= y~x) +
        stat_correlation(method = 'pearson', output.type = 'text') +
        labs(title = paste(ct_name, comp_type, vals[1], 'vs', vals[2]),
             x = paste0(comp_type, '_', vals[1]), y = paste0(comp_type, '_', vals[2]), color = color_var) +
        facet_wrap(~ Sample, ncol = 10)
      
      ggsave(file.path(outp_plot_dir, paste0('scatter_persample_', ct_name, '_', vals[1], '_', vals[2], '_', comp_type, '_', color_var, '.png')),
             width = 4000, height = 2000, unit = 'px')
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
