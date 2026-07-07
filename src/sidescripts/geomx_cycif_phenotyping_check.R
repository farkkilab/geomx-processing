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

um_to_pix_ratio <- 1 #0.325 for b2, 0.65 for b3TLS, 1 if already in pix
batchname <- 'batch123'

# from master script
proj_dir <<- '~/Documents/phd/st'
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808') # batch123
eyemt_pdrive_dir <- "/home/ad/P-drive/h30492/farkkilab2/9_EyeMT"

metadt_path <- file.path(proj_dir, 'data', 'geomx', 'metadata_full_SENSITIVE.csv')
geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS') 

deconv_ct_count_path <- file.path(output_dir, 'deconvolution', 'ct_frac_deconv_roi_mid_lvl_ct_updated_4mainimmune.csv')

# phenotyped cells dirs
# phenotyped_cells_dirs_paths <- c(file.path(eyemt_pdrive_dir, 'Data/cycif/batch1_adjacent_slides/phenotyped_cells/tribus'),
#                                  file.path(eyemt_pdrive_dir, 'Data/cycif/batch2_adjacent_slides/phenotyped_cells/tribus'),
#                                  file.path(eyemt_pdrive_dir, 'Data/cycif/batch3_adjacent_slides/phenotyped_cells/tribus'))

phenotyped_cells_path <- file.path(eyemt_pdrive_dir, 'Data/cycif/single_cell_datasets', 'phenotypes_original_samplenames_adj.csv')
phenotyped_cells_refined_outpath <- file.path(eyemt_pdrive_dir, 'Data/cycif/single_cell_datasets', 'phenotypes_original_samplenames_adj_consensus_refined.csv')

out_path_hubs_cells_inroi <- file.path(output_dir, "cycif_integration", paste0(batchname, "_adjusted_phenotypes")) #  "batch3tls_tribus_and_manualgating"
outp_plot_dir <- file.path(out_path_hubs_cells_inroi, 'plots')
dir.create(outp_plot_dir, recursive = T)


# define cell types -------------------------------------------------------

ct_names_immune <- c("Tcells_CD4", "Tcells_CD8", "DCs", "Macrophages_Monocytes") # all immune ct from phenotyping
ct_names_myeloids <- c("Macrophages_Monocytes", "DCs")
ct_names_lymphoids <- c("Tcells_CD4", "Tcells_CD8")
# additional cells from deconv not counted in phenotyping and should be treated as 'other'
ct_names_other <- c("Tcells_other", "Mast_cells", "NKcells", "Bcells") # "Bcells" goes here when basic phenotyping b1b2 and 'consensus_label_clean'

meta_names_per_roi <- c('Sample', 'Segment_geomx', 'Annotation_cell', "roi_cluster_label_gmm", "roi_cluster_label_hclust")

b1_samplename_unique <- c('S015_post', 'S015_pre', 'S027_post', 'S027_pre', 'S032_post', 'S032_pre',
                          'S053_post', 'S057_post', 'S065_post', 'S072_post', 'S073_post', 'S076_post',
                          'S084_post','S084_pre', 'S139_post', 'S139_pre')
b1_samplename_orig <- c('S015_iOme', 'S015_pPer', 'S027_iOme', 'S027_pOme', 'S032_iOval', 'S032_pOme',
                        'S053_iOme', 'S057_iOme', 'S065_iOme', 'S072_iOme', 'S073_iOme', 'S076_iOme',
                        'S084_iOme','S084_pAdn', 'S139_iOme', 'S139_pPer')

# load metadata -----------------------------------------------------------

geomx_dcc <- colnames(readRDS(geomx_norm_batch_eff_rm_path))
metadt_all <- as.data.frame(fread(metadt_path))
metadt <- metadt_all[metadt_all$dcc_filename %in% geomx_dcc, ]
rownames(metadt) <- NULL

# fix samplenames in batch1
metadt$Sample_fixed <- ifelse(metadt$main_batch_nr == 1, 
                              mapvalues(metadt$Sample, 
                                        from = b1_samplename_unique,
                                        to = b1_samplename_orig),
                              metadt$Sample)

meta_weird_coords <- metadt[metadt$roi_c1_X_cycif >= metadt$roi_c2_X_cycif | metadt$roi_c4_X_cycif >= metadt$roi_c3_X_cycif |
                            metadt$roi_c1_Y_cycif <= metadt$roi_c4_Y_cycif | metadt$roi_c2_Y_cycif <= metadt$roi_c3_Y_cycif, ]


# load final tiered results and fix sample names --------------------------

# pheno_cells_all <- fread(phenotyped_cells_path)
# 
# unique(pheno_cells_all$imageid)
# pheno_cells_all$Sample <- pheno_cells_all$imageid
# 
# # in batch1 change to _pre _post
# pheno_cells_all$Sample <- mapvalues(pheno_cells_all$Sample,
#                                     from = b1_samplename_orig,
#                                     to = b1_samplename_unique)
# 
# # in batch3 remove additional pre/suffixes
# pheno_cells_all$Sample <- gsub('^9_', '', pheno_cells_all$Sample)
# pheno_cells_all$Sample <- gsub('_[0-9]{1}_[0-9]{1}_.*', '', pheno_cells_all$Sample)
# 
# unique(pheno_cells_all$Sample)
# 
# # add Sample orig and main batch nr
# pheno_cells_all <- left_join(pheno_cells_all, distinct(metadt_all[metadt_all$Sample != '', c('Sample', 'main_batch_nr')])) # to add batchnr also for samples removed from geomx during qc
# pheno_cells_all$Sample_orig <- pheno_cells_all$Sample
# pheno_cells_all$Sample_orig <- ifelse(pheno_cells_all$main_batch_nr == 1,
#                               mapvalues(pheno_cells_all$Sample,
#                                         from = b1_samplename_unique,
#                                         to = b1_samplename_orig),
#                               pheno_cells_all$Sample)
# 
# # manually add S188 sample which was not in geomx
# pheno_cells_all$Sample_orig <- ifelse(pheno_cells_all$imageid == 'S188_iOme', 'S188_iOme', pheno_cells_all$Sample_orig)
# pheno_cells_all$main_batch_nr <- ifelse(pheno_cells_all$imageid == 'S188_iOme', '2', pheno_cells_all$main_batch_nr)
# 
# 
# kk <- distinct(pheno_cells_all[, c('imageid', 'Sample', 'Sample_orig', 'main_batch_nr')])
# 
# fwrite(pheno_cells_all, file.path(eyemt_pdrive_dir, 'Data/cycif/single_cell_datasets', 'phenotypes_original_samplenames_adj.csv'))


# load phenotyped cells and redo consensus label logic --------------------

pheno_cells_all <- fread(phenotyped_cells_path)

pheno_cells_all$aSMA_Vim <- ifelse(pheno_cells_all$aSMA == TRUE | pheno_cells_all$Vimentin == TRUE, TRUE, FALSE)

pheno_cells_all$consensus_label_refined <- pheno_cells_all$final_label_tribus

# If tumor with CD8 relabel to CD8
pheno_cells_all$consensus_label_refined <- ifelse(pheno_cells_all$consensus_label_refined == "Tumor_Tumor" & 
                                       pheno_cells_all$CD8a == TRUE, 
                                              'Tumor_Immune_CD8_Tcells', pheno_cells_all$consensus_label_refined)

# move to stroma/tumor if tum/str marker present
pheno_cells_all$consensus_label_refined <- ifelse(((pheno_cells_all$consensus_label_refined == "Stroma_Immune_CD4_Tcells") & 
                                                 (pheno_cells_all$CD4 == FALSE) & (pheno_cells_all$aSMA_Vim == TRUE)), 
                                              'Stroma_Stroma', pheno_cells_all$consensus_label_refined)
pheno_cells_all$consensus_label_refined <- ifelse(((pheno_cells_all$consensus_label_refined == "Tumor_Immune_CD4_Tcells") & 
                                                 (pheno_cells_all$CD4 == FALSE) & (pheno_cells_all$PanCK == TRUE)), 
                                              'Tumor_Tumor', pheno_cells_all$consensus_label_refined)

pheno_cells_all$consensus_label_refined <- ifelse(((pheno_cells_all$consensus_label_refined == "Stroma_Immune_CD8_Tcells") & 
                                                 (pheno_cells_all$CD8a == FALSE) & (pheno_cells_all$aSMA_Vim == TRUE)), 
                                              'Stroma_Stroma', pheno_cells_all$consensus_label_refined)
pheno_cells_all$consensus_label_refined <- ifelse(((pheno_cells_all$consensus_label_refined == "Tumor_Immune_CD8_Tcells") & 
                                                 (pheno_cells_all$CD8a == FALSE) & (pheno_cells_all$PanCK == TRUE)), 
                                              'Tumor_Tumor', pheno_cells_all$consensus_label_refined)

pheno_cells_all$consensus_label_refined <- ifelse(((pheno_cells_all$consensus_label_refined == "Stroma_Immune_Dcs") & 
                                                 (pheno_cells_all$CD11c == FALSE) & (pheno_cells_all$Iba1 == FALSE) & (pheno_cells_all$aSMA_Vim == TRUE)), 
                                              'Stroma_Stroma', pheno_cells_all$consensus_label_refined)
pheno_cells_all$consensus_label_refined <- ifelse(((pheno_cells_all$consensus_label_refined == "Tumor_Immune_Dcs") & 
                                                 (pheno_cells_all$CD11c == FALSE) & (pheno_cells_all$Iba1 == FALSE) & (pheno_cells_all$PanCK == TRUE)), 
                                              'Tumor_Tumor', pheno_cells_all$consensus_label_refined)

pheno_cells_all$consensus_label_refined <- ifelse(((pheno_cells_all$consensus_label_refined == "Stroma_Immune_Macrophages") & 
                                                 (pheno_cells_all$Iba1 == FALSE) & (pheno_cells_all$CD11c == FALSE) & (pheno_cells_all$aSMA_Vim == TRUE)), 
                                              'Stroma_Stroma', pheno_cells_all$consensus_label_refined)
pheno_cells_all$consensus_label_refined <- ifelse(((pheno_cells_all$consensus_label_refined == "Tumor_Immune_Macrophages") & 
                                                 (pheno_cells_all$Iba1 == FALSE) & (pheno_cells_all$CD11c == FALSE) & (pheno_cells_all$PanCK == TRUE)), 
                                              'Tumor_Tumor', pheno_cells_all$consensus_label_refined)

pheno_cells_all$consensus_label_clean_refined <- gsub('^Tumor_|^Stroma_|^Immune_|^Tumor_Immune_|^Stroma_Immune_', '', 
                                                        pheno_cells_all$consensus_label_refined)

fwrite(pheno_cells_all, phenotyped_cells_refined_outpath)


# bruteforce phenotyping based only on manual gates -----------------------
# BRUTEFORCE

# pheno_cells_all$aSMA_Vim <- ifelse(pheno_cells_all$aSMA == TRUE | pheno_cells_all$Vimentin == TRUE, TRUE, FALSE)
# 
# pheno_cells_all$cell_type_bruteforce <- "other"
# 
# pheno_cells_all$cell_type_bruteforce <- ifelse((pheno_cells_all$PanCK == TRUE),
#                                                         'tumor', pheno_cells_all$cell_type_bruteforce)
# 
# pheno_cells_all$cell_type_bruteforce <- ifelse((pheno_cells_all$aSMA_Vim == TRUE),
#                                                         'stroma', pheno_cells_all$cell_type_bruteforce)
# 
# pheno_cells_all$cell_type_bruteforce <- ifelse((pheno_cells_all$CD4 == TRUE & 
#                                                            pheno_cells_all$CD8a == FALSE & 
#                                                            pheno_cells_all$CD11c == FALSE &
#                                                            pheno_cells_all$Iba1 == FALSE &
#                                                            pheno_cells_all$PanCK == FALSE &
#                                                            pheno_cells_all$aSMA_Vim == FALSE),
#                                                         'Tcells_CD4', pheno_cells_all$cell_type_bruteforce)
# 
# pheno_cells_all$cell_type_bruteforce <- ifelse(((pheno_cells_all$CD11c == TRUE &
#                                                             pheno_cells_all$Iba1 == TRUE) | 
#                                                            (pheno_cells_all$CD11c == FALSE &
#                                                               pheno_cells_all$Iba1 == TRUE)),
#                                                         'Macrophages_Monocytes', pheno_cells_all$cell_type_bruteforce)
# 
# pheno_cells_all$cell_type_bruteforce <- ifelse((pheno_cells_all$CD11c == TRUE &
#                                                            pheno_cells_all$Iba1 == FALSE),
#                                                         'DCs', pheno_cells_all$cell_type_bruteforce)
# 
# pheno_cells_all$cell_type_bruteforce <- ifelse((pheno_cells_all$CD8a == TRUE),
#                                                         'Tcells_CD8', pheno_cells_all$cell_type_bruteforce)
# 
# 
# fwrite(pheno_cells_all, phenotyped_cells_refined_outpath)

# iterate through samples and select cells within ROIs --------------------

pheno_cells_inroi_all <- lapply(unique(metadt$Sample), function(sample_name){
  print('###################################')
  print(sample_name)
  
  pheno_cells_sample <- pheno_cells_all[pheno_cells_all$Sample == sample_name, ]
  
  # convert to pix
  pheno_cells_sample$X_centroid_px <- pheno_cells_sample$X_centroid / um_to_pix_ratio
  pheno_cells_sample$Y_centroid_px <- pheno_cells_sample$Y_centroid / um_to_pix_ratio
  
  # count cells within all ROIs in the sample - now in metadata
  roi_coords_sample <- metadt[metadt$Sample == sample_name, c("Sample","Sample_fixed", "main_batch_nr", "sample_roi", "Roi_geomx",
                                                                                "roi_c1_X_cycif", "roi_c1_Y_cycif", "roi_c2_X_cycif", "roi_c2_Y_cycif",
                                                                                "roi_c3_X_cycif", "roi_c3_Y_cycif", "roi_c4_X_cycif", "roi_c4_Y_cycif")]
  
  roi_coords_sample <- distinct(roi_coords_sample)
  
  # iterate thorugh rois and count cells
  cells_in_roi_all <- apply(roi_coords_sample, 1, function(row){
    
    #new finding cells implementation rewritten from python shapely package to sp R package
    poly_x <- c(row[['roi_c1_X_cycif']], row[['roi_c2_X_cycif']], row[['roi_c3_X_cycif']], row[['roi_c4_X_cycif']])
    poly_y <- c(row[['roi_c1_Y_cycif']], row[['roi_c2_Y_cycif']], row[['roi_c3_Y_cycif']], row[['roi_c4_Y_cycif']])
    
    cells_in_roi <- dplyr::filter(pheno_cells_sample, sp::point.in.polygon(point.x = as.numeric(X_centroid_px),
                                                                           point.y = as.numeric(Y_centroid_px),
                                                                           pol.x = poly_x, pol.y = poly_y) != 0)
    
    cells_in_roi <- cbind(cells_in_roi, as.data.frame(lapply(row, rep, nrow(cells_in_roi))))
    return(cells_in_roi)
  })
  
  cells_in_roi_all <- do.call(rbind, cells_in_roi_all)
  print(nrow(cells_in_roi_all))
  table(cells_in_roi_all$Roi_geomx)

  return(cells_in_roi_all)
})

pheno_cells_inroi_all_df <- do.call("rbind.fill", pheno_cells_inroi_all)


# adjust labels to deconv ctnames
pheno_cells_inroi_all_df$cell_type <- pheno_cells_inroi_all_df$consensus_label_clean_refined
pheno_cells_inroi_all_df$cell_type <- mapvalues(pheno_cells_inroi_all_df$cell_type,
                                                       from = c("Tumor", "Dcs", "CD8_Tcells", "Stroma", "Macrophages", "CD4_Tcells", "Undefined"),
                                                       to = c("tumor", "DCs", "Tcells_CD8", "stroma", "Macrophages_Monocytes", "Tcells_CD4", "other"))


fwrite(pheno_cells_inroi_all_df, file.path(out_path_hubs_cells_inroi, paste0(batchname, '_cells_in_roi.csv')))

labs_count <- group_by(pheno_cells_inroi_all_df, cell_type, final_label_tribus, CD4, CD8a, CD11c, Iba1, PanCK, aSMA_Vim) %>%
  summarise(n = n())

#########################################
# to check b3TLS 
# pheno_cells_inroi_all_df <- fread('/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2808/cycif_integration/batch3tls_hubs_cells_inroi_combined_myeloids_min20cells_nomix_dt171715_ct10_dt300_new_polygons.csv')
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

fwrite(ct_frac_cycif_long_roi, file.path(out_path_hubs_cells_inroi,  paste0(batchname, '_ct_frac_cycif_roi.csv')))


# merge with cells counted from deconv ------------------------------------

ct_frac_deconv_long_roi <- fread(deconv_ct_count_path)

ct_frac_long_roi <- left_join(ct_frac_cycif_long_roi, ct_frac_deconv_long_roi, by = c('sample_roi', 'cell_type'))
ct_frac_long_roi <- left_join(ct_frac_long_roi, distinct(metadt[, c('sample_roi', meta_names_per_roi)]))

# remove all ROIs with total nr of cells < 100
ct_frac_long_roi <- ct_frac_long_roi[ct_frac_long_roi$total_cell_nr_cycif > 100, ]

ct_frac_long_roi$sd_cycif_diff <- abs(ct_frac_long_roi$ct_frac_sd - ct_frac_long_roi$ct_frac_cycif)
ct_frac_long_roi$bp_cycif_diff <- abs(ct_frac_long_roi$ct_frac_bp - ct_frac_long_roi$ct_frac_cycif)
ct_frac_long_roi$bp_sd_diff <- abs(ct_frac_long_roi$ct_frac_bp - ct_frac_long_roi$ct_frac_sd)

ct_frac_long_roi$ct_frac_bpsd <- (ct_frac_long_roi$ct_frac_bp + ct_frac_long_roi$ct_frac_sd)/2

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
for(value_comb in c('bp_sd', 'bp_cycif', 'sd_cycif', 'bpsd_cycif')){
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
