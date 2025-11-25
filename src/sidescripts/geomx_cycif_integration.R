library(data.table)
library(plyr)
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(ggplot2)
library(tibble)

roi_coords_dir <- "/home/ad/P-drive/h30492/farkkilab2/9_EyeMT/CyciF_GeoMx_alignment/CyciF_batch2_ROIs"
cycif_cell_count_dir <- "/home/ad/P-drive/h30492/farkkilab2/9_EyeMT/9_EyeMT_Cycif/batch2_adjacent_slides/subsetting_stardist"

output_dir <- "~/Documents/phd/st/geomx-processing/results/batch2-1903/cycif_cell_count"
output_path <- "~/Documents/phd/st/geomx-processing/results/batch2-1903/cycif_cell_count/batch2_cycif_cell_count_per_roi_stardist.csv"

patient_names <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(list.files(roi_coords_dir)))


cells_in_roi_all_pt <- lapply(patient_names, function(pt_name){
  
  print(pt_name)
  
  cell_count <- fread(file.path(cycif_cell_count_dir, pt_name,  paste0(pt_name, '_updated.csv')))
  
  # clean roi coords df
  roi_coords <- fread(file.path(roi_coords_dir, paste0(pt_name, '.csv')), drop = c(6, 11))
  colnames(roi_coords) <- c('roi_name', 'c1_X', 'c2_X', 'c3_X', 'c4_X', 'c1_Y', 'c2_Y', 'c3_Y', 'c4_Y')
  roi_coords$c1_X <- gsub('\\[|\\]', '', roi_coords$c1_X)
  roi_coords$c1_Y <- gsub('\\[|\\]', '', roi_coords$c1_Y)
  roi_coords <- mutate_at(roi_coords, vars(matches("c[0-9]")), function(x){gsub(" ", "", x)})
  roi_coords <- mutate_at(roi_coords, vars(matches("c[0-9]")), as.numeric)
  roi_coords$roi_name <- as.character(roi_coords$roi_name)
  
  # count cells within all ROIs in the sample
  
  cells_in_roi_all <- apply(roi_coords, 1, function(row){
    
    # find cells within range
    # coordinates are not longer rectangles, they're a bit rotated
    # the cells are found inside longer edges of rectangle
    print(row[["roi_name"]])
    # print(min(row[['c1_X']], row[['c4_X']]))
    # print(max(row[['c2_X']], row[['c3_X']]))
    # print(min(row[['c4_Y']], row[['c3_Y']]))
    # print(max(row[['c1_Y']], row[['c2_Y']]))
    row$roi_width <- as.numeric(max(as.numeric(row[['c2_X']]), as.numeric(row[['c3_X']])) - min(as.numeric(row[['c1_X']]), as.numeric(row[['c4_X']])))
    row$roi_height <- as.numeric(max(as.numeric(row[['c1_Y']]), as.numeric(row[['c2_Y']])) - min(as.numeric(row[['c3_Y']]), as.numeric(row[['c4_Y']])))
    row$roi_center_X <- as.numeric(min(as.numeric(row[['c1_X']]), as.numeric(row[['c4_X']])) + (row$roi_width * 0.5))
    row$roi_center_Y <- as.numeric(min(as.numeric(row[['c4_Y']]), as.numeric(row[['c3_Y']])) + (row$roi_height * 0.5))
    
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

cells_in_roi_all_pt <- do.call(rbind, cells_in_roi_all_pt)

########################

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
