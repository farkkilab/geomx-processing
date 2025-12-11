# integrates geomx and cycif roi coordinates
# integrates cycif cell counts
# integrates hubs 

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

# define paths ------------------------------------------------------------
# from master script
batch <- 'batch123'
proj_dir <<- '~/Documents/phd/st'
data_dir <<- '~/Documents/phd/st/data/geomx/batch123/' # batch1 2 and 3
anno_path <<- file.path(data_dir, 'metadata', 'dcc_metadata_batch123_no_tls_cleaned.xlsx') #batch1 and 2 and 3
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808') # batch123
geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS') 

# coordinates, cell count and hubs pathways for cycif images
eyemt_pdrive_dir <- "/home/ad/P-drive/h30492/farkkilab2/9_EyeMT/"
roi_coords_dir <- file.path(eyemt_pdrive_dir, "Data/geomx/batch2/rois_from_tcycif/cyciF_batch2_ROIs_arrays")
cycif_cell_count_dir <- file.path(eyemt_pdrive_dir, "Data/cycif/batch2_adjacent_slides/phenotyped_cells/tribus/stardist/final_labels_after_NK_gating")
hubs_dir <- file.path(eyemt_pdrive_dir, "/Data_analysis/spatial_analysis/SPACEstat/batch2_interaction_hubs")

# output dirs and paths 
dir.create(file.path(output_dir, "cycif_integration"))
out_path_coords <- file.path(output_dir, "cycif_integration", "batch2_cycif_coordinates.csv")
out_path_cell_count <- file.path(output_dir, "cycif_integration", "batch2_cycif_cell_count_per_roi_stardist.csv")

# load geomx, merge with cleaned metadata ---------------------------------
# TODO run once again in 1811 with already cleaned metadata and just load meta from geomx
meta_geomx <- pData(readRDS(geomx_norm_batch_eff_rm_path))
meta_cleaned <- read_excel(anno_path)

metadt <- left_join(meta_geomx[, c('dcc_filename', 'Slide_Name')], meta_cleaned) # join ensuring order

rm(meta_geomx)
rm(meta_cleaned)

# subset to batch2
metadt <- metadt[metadt$main_batch_nr == 2, ]

# clean and calculate cycif coordinates -----------------------------------

#patient_names <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(list.files(roi_coords_dir)))
sample_names <- file_path_sans_ext(basename(list.files(roi_coords_dir)))

roi_coords_all <- lapply(sample_names, function(sample_name){
  
  # clean roi coords df
  roi_coords <- as.data.frame(fread(file.path(roi_coords_dir, paste0(sample_name, '.csv')), drop = c(6, 11)))
  colnames(roi_coords) <- c('roi_name', 'c1_X', 'c2_X', 'c3_X', 'c4_X', 'c1_Y', 'c2_Y', 'c3_Y', 'c4_Y')
  roi_coords$c1_X <- gsub('\\[|\\]', '', roi_coords$c1_X)
  roi_coords$c1_Y <- gsub('\\[|\\]', '', roi_coords$c1_Y)
  roi_coords <- mutate_at(roi_coords, vars(matches("c[0-9]")), function(x){gsub(" ", "", x)})
  roi_coords <- mutate_at(roi_coords, vars(matches("c[0-9]")), as.numeric)
  roi_coords$roi_name <- as.character(roi_coords$roi_name)
  
  #calculate additional dimentions
  roi_coords <-  apply(roi_coords, 1, function(row){
    row$roi_width <- round(as.numeric(max(as.numeric(row[['c2_X']]), as.numeric(row[['c3_X']])) - min(as.numeric(row[['c1_X']]), as.numeric(row[['c4_X']]))),3)
    row$roi_height <- round(as.numeric(max(as.numeric(row[['c1_Y']]), as.numeric(row[['c2_Y']])) - min(as.numeric(row[['c3_Y']]), as.numeric(row[['c4_Y']]))),3)
    row$roi_center_X <- round(as.numeric(min(as.numeric(row[['c1_X']]), as.numeric(row[['c4_X']])) + (row$roi_width * 0.5)),3)
    row$roi_center_Y <- round(as.numeric(min(as.numeric(row[['c4_Y']]), as.numeric(row[['c3_Y']])) + (row$roi_height * 0.5)),3)
    return(as.data.frame(row))
  })
  
  roi_coords <- as.data.frame(do.call(rbind, roi_coords))
  roi_coords <- mutate_at(roi_coords, vars(matches("c[0-9]")), as.numeric)
  roi_coords$Sample <- sample_name
  
  # TODO save per sample - not sure if needed
  return(roi_coords)
})

roi_coords_all <- do.call(rbind, roi_coords_all)
roi_coords_all <- roi_coords_all[, c(ncol(roi_coords_all), 1:(ncol(roi_coords_all)-1))]
fwrite(roi_coords_all, out_path_coords)

# subset cells to the ones in ROIs ----------------------------------------

cells_in_roi_all_sample <- lapply(unique(roi_coords_all$Sample), function(sample_name){
  
  print(sample_name)
  
  cell_count <- fread(file.path(cycif_cell_count_dir, sample_name,  paste0(sample_name, '_updated.csv')))
  
  # count cells within all ROIs in the sample
  roi_coords_sample <- roi_coords_all[roi_coords_all$Sample == sample_name, ]
  
  
  cells_in_roi_all <- apply(roi_coords_sample, 1, function(row){
    
    # find cells within range
    # coordinates are not longer rectangles, they're a bit rotated
    # the cells are found inside longer edges of rectangle
    print(row[["roi_name"]])

    cells_in_roi <- dplyr::filter(cell_count,
                           as.numeric(X_centroid) >= min(as.numeric(row[['c1_X']]), as.numeric(row[['c4_X']])) &
                             as.numeric(X_centroid) <= max(as.numeric(row[['c2_X']]), as.numeric(row[['c3_X']])) &
                             as.numeric(Y_centroid) >= min(as.numeric(row[['c4_Y']]), as.numeric(row[['c3_Y']])) &
                             as.numeric(Y_centroid) <= max(as.numeric(row[['c1_Y']]), as.numeric(row[['c2_Y']])))
  
    print(nrow(cells_in_roi))
    
    cells_in_roi <- cbind(cells_in_roi, as.data.frame(lapply(row, rep, nrow(cells_in_roi))))
    
    return(cells_in_roi)
  })
  
  cells_in_roi_all <- do.call(rbind, cells_in_roi_all)
  
  rm(cell_count)
  return(cells_in_roi_all)
})

cells_in_roi_all_sample <- do.call(rbind, cells_in_roi_all_sample)




# integrate ROIs and hubs -------------------------------------------------
# hubs_inter <- fread(file.path(hubs_dir, "eyemt_batch2_interactions_dt15171517_ct15.csv"))
# hubs_cells <- fread(file.path(hubs_dir, "eyemt_batch2_cells_dt15171517_ct15.csv"))
hubs_inter <- fread(file.path(hubs_dir, "eyemt_batch2_interactions_dt20222022_ct15_mt20.csv"))
hubs_cells <- fread(file.path(hubs_dir, "eyemt_batch2_cells_dt20222022_ct15_mt20.csv"))
hubs_cells$interaction_hub_ids <- gsub("['", "", hubs_cells$interaction_hub_ids, fixed = T)
hubs_cells$interaction_hub_ids <- gsub("']", "", hubs_cells$interaction_hub_ids, fixed = T)
hubs_cells$interaction_hub_ids <- gsub("',.*", "", hubs_cells$interaction_hub_ids, fixed = F)
# TODO it might be changed back to pix in the original files
hubs_cells$X_centroid_px <- hubs_cells$X_centroid / 0.325
hubs_cells$Y_centroid_px <- hubs_cells$Y_centroid / 0.325

hubs_cells <- left_join(hubs_cells, hubs_inter[, c('hub_id', 'hub_type')], by = c('interaction_hub_ids' = 'hub_id'))

length(unique(hubs_cells$interaction_hub_ids))
length(unique(hubs_inter$hub_id))
length(intersect(hubs_cells$interaction_hub_ids, hubs_inter$hub_id))
setdiff(hubs_cells$interaction_hub_ids, hubs_inter$hub_id)
setdiff(hubs_inter$hub_id, hubs_cells$interaction_hub_ids)
# "S100_iOme_hub_328" "S131_iOme_hub_50"  "S083_iOme2_hub_86"

hubs_cells_in_roi_all <- lapply(unique(roi_coords_all$Sample), function(sample_name){
  print(sample_name)
  
  hubs_cells_sample <- hubs_cells[hubs_cells$imageid == sample_name, ]
  
  # count cells within all ROIs in the sample
  roi_coords_sample <- roi_coords_all[roi_coords_all$Sample == sample_name, ]
  
  
  cells_in_roi_all <- apply(roi_coords_sample, 1, function(row){
    
    # find cells within range
    # coordinates are not longer rectangles, they're a bit rotated
    # the cells are found inside longer edges of rectangle
    print(row[["roi_name"]])
    
    cells_in_roi <- dplyr::filter(hubs_cells_sample,
                                  as.numeric(X_centroid_px) >= min(as.numeric(row[['c1_X']]), as.numeric(row[['c4_X']])) &
                                    as.numeric(X_centroid_px) <= max(as.numeric(row[['c2_X']]), as.numeric(row[['c3_X']])) &
                                    as.numeric(Y_centroid_px) >= min(as.numeric(row[['c4_Y']]), as.numeric(row[['c3_Y']])) &
                                    as.numeric(Y_centroid_px) <= max(as.numeric(row[['c1_Y']]), as.numeric(row[['c2_Y']])))
    
    print(nrow(cells_in_roi))
    
    cells_in_roi <- cbind(cells_in_roi, as.data.frame(lapply(row, rep, nrow(cells_in_roi))))
    
    return(cells_in_roi)
  })
  
  cells_in_roi_all <- do.call(rbind, cells_in_roi_all)
  
  rm(cell_count)
  return(cells_in_roi_all)
})

hubs_cells_in_roi_all <- do.call(rbind, hubs_cells_in_roi_all)
hubs_cells_in_roi_all$interaction_hub_ids <- ifelse(hubs_cells_in_roi_all$interaction_hub_ids == '', 'notinhub', 
                                                    hubs_cells_in_roi_all$interaction_hub_ids)
hubs_cells_in_roi_all$sample_roi <- paste0(hubs_cells_in_roi_all$Sample, '_', as.character(hubs_cells_in_roi_all$roi_name))


table(hubs_cells_in_roi_all$Sample, hubs_cells_in_roi_all$roi_name)
######################################################################
######################################################################
# explore hubs vs ROI

hubs_cells_in_roi_inhub <- hubs_cells_in_roi_all[hubs_cells_in_roi_all$interaction_hub_ids != 'notinhub', ]

table(hubs_cells_in_roi_all$Sample, hubs_cells_in_roi_all$roi_name)
table(hubs_cells_in_roi_inhub$Sample, hubs_cells_in_roi_inhub$roi_name)

# how many rois does not have any hubs?
hubs_unique <- distinct(hubs_cells_in_roi_all, sample_roi, interaction_hub_ids, hub_type)
length(unique(hubs_unique$sample_roi))

# nr of unique hubs in each ROI
hubs_unique_roi <- distinct(hubs_cells_in_roi_inhub, sample_roi, interaction_hub_ids, hub_type)
sort(table(hubs_unique_roi$sample_roi))
length(unique(hubs_unique_roi$sample_roi))

sort(table(hubs_unique_roi$hub_type))

#######################################################
# for sample
sname <- 'S098_iOme'
# % of cells in each roi which belongs to the hub
hubs_cells_sample <- hubs_cells_in_roi_all[hubs_cells_in_roi_all$imageid == sname, ]
cells_sample <- hubs_cells[hubs_cells$imageid == sname, ]

roi_coords_sample <- roi_coords_all[roi_coords_all$Sample == sname, ]

apply(roi_coords_sample,1, function(row){
  print(row[['roi_name']])
  print(as.numeric(row[['c1_X']])*0.325)
  print(as.numeric(row[['c1_Y']])*0.325)
  print(as.numeric(row[['c2_X']])*0.325)
  print(as.numeric(row[['c2_Y']])*0.325)
  print(as.numeric(row[['c3_X']])*0.325)
  print(as.numeric(row[['c3_Y']])*0.325)
  print(as.numeric(row[['c4_X']])*0.325)
  print(as.numeric(row[['c4_Y']])*0.325)
})
######################################################################
######################################################################
# clustering and relabeling based on the cell counts in roi

# cleaning labels
cells_in_roi_all_pt$sample_roi <- paste0(cells_in_roi_all_pt$Sample, '_', cells_in_roi_all_pt$roi_name)

cells_in_roi_all_pt$final_label <-  transmute(cells_in_roi_all_pt, final_label = 
                                                plyr::mapvalues(final_label, c('undefined_Intumor_CD4Tcells',
                                                         'undefined_Intumor_CD8Tcells',
                                                         'Stroma', 'Tumor',
                                                         'undefined_Instroma',
                                                         'undefined_Intumor'), 
                                          c('Intumor_CD4Tcells',
                                            'Intumor_CD8Tcells',
                                            'Instroma_Fibroblasts',
                                            'Intumor_Tumor',
                                            'Instroma_undefined',
                                            'Intumor_undefined')))

cells_in_roi_all_pt$final_label <- ifelse(cells_in_roi_all_pt$final_label == 'NK', 
                                          paste0(cells_in_roi_all_pt$Global, '_NK'), 
                                          cells_in_roi_all_pt$final_label)

# dropping undefined_Global_NK bcs there is just 3 of them
cells_in_roi_all_pt$final_label <- ifelse(cells_in_roi_all_pt$final_label %in% c('undefined_Global_NK', 'undefined_Global'), 
                                          'Global_undefined', cells_in_roi_all_pt$final_label)

cells_in_roi_all_pt$Segment <- gsub('_.*', '', cells_in_roi_all_pt$final_label)
cells_in_roi_all_pt$Segment <- gsub('In', '', cells_in_roi_all_pt$Segment)

cells_in_roi_all_pt$cell_type <- gsub('^.*_', '', cells_in_roi_all_pt$final_label)
cells_in_roi_all_pt$cell_type <-  transmute(cells_in_roi_all_pt, cell_type = 
                                                plyr::mapvalues(cell_type, c('CD11c', 'CD4Tcells',
                                                                               'CD8Tcells', 'Macrophages',
                                                                               'NK', 'Tumor', 'undefined', 'Fibroblasts'), 
                                                                c('DCs', 'Tcells_CD4',
                                                                  'Tcells_CD8', 'Macrophages_Monocytes',
                                                                  'NKs', 'tumor', 'other', 'Fibroblasts')))


fwrite(cells_in_roi_all_pt, output_path)

########################################

cell_types <- c("DCs", "Fibroblasts", "Macrophages_Monocytes", "NKs", "other",
                "Tcells_CD4", "Tcells_CD8", "tumor")
cell_types_important <- c("DCs", "Macrophages_Monocytes", "Tcells_CD4", "Tcells_CD8")
cell_types_immune <- c("DCs", "Macrophages_Monocytes", "NKs", "Tcells_CD4", "Tcells_CD8")

# count nr of cells per ROI and AOI

# count cells per roi
roi_ct <- dcast(cells_in_roi_all_pt, sample_roi ~ cell_type)
roi_ct$total_cell_nr <- rowSums(roi_ct[, cell_types])
roi_ct$immune_cell_nr <- rowSums(roi_ct[, cell_types_immune])

roi_ct_long <- melt(roi_ct[, !(names(roi_ct) %in% c('total_cell_nr', 'immune_cell_nr'))], 
                                        id.vars = c("sample_roi"), variable.name = "cell_type")

roi_ct_frac <- mutate_at(roi_ct, vars(cell_types), funs(. / total_cell_nr))
roi_ct_frac_long <- melt(roi_ct_frac[, !(names(roi_ct_frac) %in% c('total_cell_nr', 'immune_cell_nr'))], 
                         id.vars = c("sample_roi"), variable.name = "cell_type")

roi_ct_frac_immune <- mutate_at(roi_ct[, c("sample_roi", cell_types_immune, "immune_cell_nr")], 
                                vars(cell_types_immune), funs(. / immune_cell_nr))
roi_ct_frac_immune_long <- melt(roi_ct_frac_immune[, !(names(roi_ct_frac_immune) %in% c('immune_cell_nr'))], 
                         id.vars = c("sample_roi"), variable.name = "cell_type")

#TODO merge with our labels

# count cells per aoi
# !!! Global unidentified cells are max 7 in 23 ROIs
# exchanging them to stromal segment to keep the cell count the same
aoi_ct <- cells_in_roi_all_pt
aoi_ct$Segment <- ifelse(aoi_ct$Segment == 'Global', 'stroma', aoi_ct$Segment) # !!!!
aoi_ct <- dcast(aoi_ct, sample_roi+Segment ~ cell_type)
aoi_ct$total_cell_nr <- rowSums(aoi_ct[, 3:ncol(aoi_ct)])

aoi_ct_frac <- mutate_at(aoi_ct, vars(3:10), funs(. / total_cell_nr))
aoi_ct_frac_long <- melt(aoi_ct_frac[, -11], id.vars = c("sample_roi", "Segment"), variable.name = "cell_type")

################################

# TODO use ct number instead of fractions bcs if there are 2 ct abundant, fraction will be lower
# TODO also fraction of immune cells

# check ct fractions distribution
for(ct_name in cell_types_immune){
  roi_ct_long_ct <- roi_ct_frac_immune_long[roi_ct_frac_immune_long$cell_type == ct_name, ]
  #roi_ct_frac_immune_long_ct <- roi_ct_frac_immune_long[roi_ct_frac_immune_long$cell_type == ct_name, ]
  #roi_ct_frac_long_ct <- roi_ct_frac_long[roi_ct_frac_long$cell_type == ct_name, ]
  #aoi_ct_frac_long_ct <- aoi_ct_frac_long[aoi_ct_frac_long$cell_type == ct_name, ]
  
  ggplot(data = roi_ct_long_ct) +
    #geom_density(aes(value)) +
    geom_histogram(aes(value), bins = 100)
    ggtitle(ct_name)
  
  ggsave(file.path(output_dir, paste0('hist_fraq_immune_', ct_name, '.png')), width = 2000, height = 1000, unit='px')
  
  # ggplot(data = aoi_ct_frac_long_ct) +
  #   geom_density(aes(value, color = Segment)) +
  #   ggtitle(ct_name)
  # 
  # ggsave(file.path(output_dir, paste0('density_fraq_aoi_', ct_name, '.png')), width = 2000, height = 1000, unit='px')
  # 
}

##############################################
# cluster ROIs per cell fraction in a hmap

ct_frac_mtx <- as.matrix(column_to_rownames(roi_ct_frac_immune[, !(names(roi_ct_frac_immune) %in% c('total_cell_nr', 'immune_cell_nr'))],
                                            'sample_roi'))

# TODO comment/uncomment for important cells
ct_frac_mtx <- ct_frac_mtx[, c(cell_types_important)]
ct_frac_mtx_zscore <- scale(ct_frac_mtx) # zscore by column

# cluster by hclust
ct_frac_hclust <- hclust(dist(ct_frac_mtx), method = "average")
plot(ct_frac_hclust, hang = -1, cex = 0.4)
ct_frac_hclust_cut <- cutree(ct_frac_hclust, h = 0.3)

ct_frac_zscore_hclust <- hclust(dist(ct_frac_mtx_zscore), method = "average")
plot(ct_frac_zscore_hclust, hang = -1, cex = 0.4)
ct_frac_zscore_hclust_cut <- cutree(ct_frac_zscore_hclust, h = 2)
ct_frac_zscore_hclust_cut_k <- cutree(ct_frac_zscore_hclust, k = 11)

#########
# TODO annotations on hmaps are wrong - only match zscores
# make hmaps
ha = HeatmapAnnotation(
  #ct_label = anno_simple(roi_ct_frac$Annotation_cell),
  hclust_h = anno_simple(as.character(unname(ct_frac_zscore_hclust_cut))),
  hclust_k = anno_simple(as.character(unname(ct_frac_zscore_hclust_cut_k))),
  which = "row", show_legend = TRUE)

# hmap for ct fraq
png(filename=file.path(output_dir, paste0('hmap_roi_fraq_immune_important_ct.png')), 
    width=10, height=6,units="in",res=2000)

ind_heat <- Heatmap(ct_frac_mtx, cluster_columns = F, cluster_rows= ct_frac_hclust,
                    show_row_names = TRUE, show_column_names = TRUE,
                    left_annotation = ha, show_heatmap_legend = TRUE)


draw(ind_heat, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

# hmap for zscore
png(filename=file.path(output_dir, paste0('hmap_roi_zscore_fraq_immune_important_ct.png')), 
    width=10, height=6,units="in",res=2000)

ind_heat <- Heatmap(ct_frac_mtx_zscore, cluster_columns = F, cluster_rows= ct_frac_zscore_hclust,
                    show_row_names = TRUE, show_column_names = TRUE,
                    left_annotation = ha, show_heatmap_legend = TRUE)


draw(ind_heat, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

# TODO the same for geomx_segment
# TODO compare with deconvoluted fractions
# TODO compare with our labels
# TODO add fractions (from all cells) for cell types (eg macro fraq + dc frac and then: check distrib, label highest ones)
