# integrates cycif cell counts + hubs labels by ROI coordinates

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

# define paths ------------------------------------------------------------
# from master script
proj_dir <<- '~/Documents/phd/st'
data_dir <<- '~/Documents/phd/st/data/geomx/batch123/' # batch1 2 and 3
anno_path <<- file.path(data_dir, 'metadata', 'dcc_metadata_batch123_no_tls_cleaned.xlsx') #batch1 and 2 and 3
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808') # batch123
eyemt_pdrive_dir <- "/home/ad/P-drive/h30492/farkkilab2/9_EyeMT"
geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS') 


# cell count and hubs for cycif images pathways 
#cycif_cell_count_dir <- file.path(eyemt_pdrive_dir, "Data/cycif/batch2_adjacent_slides/phenotyped_cells/tribus/stardist/final_labels_after_NK_gating")
#hubs_dir <- file.path(eyemt_pdrive_dir, "/Data_analysis/spatial_analysis/SPACEstat/batch2_interaction_hubs")
#hubs_comm_dir <- file.path(eyemt_pdrive_dir, "/Data_analysis/spatial_analysis/SPACEstat/batch2_communities")
hubs_comm_dir <- file.path(eyemt_pdrive_dir, "Data_analysis/spatial_analysis/SPACEstat/batch3_communities")

# hubs files
dt_hubparam <- 'dt1515151717' # version of distances metrics used for hubs 
ct_hubparam <- 'ct10'
dt_comm_hubparam <- 'dt300'
#hubs_cells_path <- file.path(hubs_dir, paste0("eyemt_batch2_cells_", dt_hubparam, "_", ct_hubparam,  ".csv")) # hubs per cell
hubs_cells_path <- file.path(hubs_comm_dir, paste0("eyemt_batch3_cells_", dt_hubparam, "_", ct_hubparam, "_", dt_comm_hubparam,  "_cmp.csv")) # hubs per cell
#hubs_inter_path <- file.path(hubs_dir, paste0("eyemt_batch2_interactions_", dt_hubparam, "_", ct_hubparam,  ".csv")) # interaction hubs
#hubs_comm_path <- file.path(hubs_comm_dir, paste0("eyemt_batch2_communities_with_annotations_", ct_hubparam, "_", dt_comm_hubparam,  "_", res_hubparam, "_leiden_cmp.csv"))
hubs_comm_path <- file.path(hubs_comm_dir, paste0("eyemt_batch3_network_communities_consensus_", dt_hubparam, "_", ct_hubparam, "_",  dt_comm_hubparam,  "_cmp.csv"))

um_to_pix_ratio <- 0.65 #0.325 for b2, 0.65 for b3TLS

# output dirs and paths 
dir.create(file.path(output_dir, "cycif_integration"))
#out_path_cell_count <- file.path(output_dir, "cycif_integration", "batch2_cycif_cell_count_per_roi_stardist.csv")
# out_path_hubs_cells <- file.path(output_dir, "cycif_integration", 
#                                  paste("batch2_hubs_cells", dt_hubparam, ct_hubparam, dt_comm_hubparam, res_hubparam, ".csv", sep = '_'))
# out_path_hubs_cells_inroi <- file.path(output_dir, "cycif_integration", 
#                                        paste("batch2_hubs_cells_inroi", dt_hubparam, ct_hubparam, dt_comm_hubparam, res_hubparam, ".csv", sep = '_'))
# out_path_hubs_cells <- file.path(output_dir, "cycif_integration", 
#                                  paste("batch3tls_hubs_cells_bcells_", dt_hubparam, dt_comm_hubparam, ".csv", sep = '_'))
out_path_hubs_cells_inroi <- file.path(output_dir, "cycif_integration", 
                                       paste0("batch3tls_hubs_cells_inroi_bcells_", dt_hubparam, '_', ct_hubparam,
                                             '_', dt_comm_hubparam, ".csv"))

####################################
# quickfix
# S084iOme removed during QC bcs of too small GDR
# also S015_pPer (b3) - too small GDR
# also "S139_post" "S139_pre" (b1) super high NTC
b3tls_imageid <- c("S015", "S080", "S081", "S091", "S106", "S112", "S113", "S118", "S120",
                   "S123", "S195", "S225", "S229", "S247", "S309", "S311", "S355", "S378", "S380")

b3tls_samplename <- c("S015_iOme", "S080_iOme2", "S081_iOme", "S091_iOme1", "S106_iOme", "S112_iOme", "S113_iOme", 
                      "S118_iOme",  "S120_iOme", "S123_iOme", "S195_iOme1", "S225_iOme",  "S229_iOme",  "S247_iOme", "S309_iOme", 
                      "S311_iOme", "S355_iOme", "S378_iOme", "S380_iOme") 

# load cleaned metadata ---------------------------------------------------
# TODO run once again in 1811 with already cleaned metadata and just load meta from geomx
meta_geomx <- pData(readRDS(geomx_norm_batch_eff_rm_path))
meta_cleaned <- read_excel(anno_path)

metadt <- left_join(meta_geomx[, c('dcc_filename', 'Slide_Name')], meta_cleaned) # join ensuring order
metadt$sample_roi <- paste0(metadt$Sample, '_', metadt$Roi_geomx)
rm(meta_geomx)
rm(meta_cleaned)

print(colnames(metadt))
# subset to batch2/3
metadt <- metadt[metadt$Sample %in% b3tls_samplename, ]

length(unique(metadt$sample_roi))

# merge all hubs labels per cell ------------------------------------------
hubs_cells <- fread(hubs_cells_path)
hubs_comm <- fread(hubs_comm_path)

double_hubs <- hubs_cells[grepl(',', hubs_cells$interaction_hub_ids), ] # to keep just in case

# in b2 for 294 cells  >1 hub - keep one
# in b3tls for 121 cells  >1 hub - keep one
# in b3tls with bcells for 366 cells >1 hub - keep one
# in b3tls with bcells and ct10 for 856 cells > keep one
# TODO think about this strategy
hubs_cells$interaction_hub_ids <- gsub(",.*", "", hubs_cells$interaction_hub_ids)

# TODO it might be changed back to pix in the original files
hubs_cells$X_centroid_px <- hubs_cells$X_centroid / um_to_pix_ratio
hubs_cells$Y_centroid_px <- hubs_cells$Y_centroid / um_to_pix_ratio

if(!('cluster_label' %in% colnames(hubs_cells))){
  hubs_comm$cluster_label <- paste0('cluster_', hubs_comm$cluster)
}
 
# join with network hubs and communities
hubs_cells <- left_join(hubs_cells, hubs_comm[, c('hubid', 'hub_type', 'hub_size', 'residency', 'community', 'cluster', 'cluster_label')], 
                        by = c('network_hub_id' = 'hubid'))
hubs_cells <-  dplyr::rename(hubs_cells, network_hub_type = hub_type, network_hub_size = hub_size, 
                             community_id = community, community_cluster = cluster, community_cluster_label = cluster_label)

#TODO change imageid to Sample - already handled during phenotyping
hubs_cells$Sample <- mapvalues(hubs_cells$imageid, 
                               from = b3tls_imageid,to = b3tls_samplename)
# filter to shared samples
hubs_cells <- hubs_cells[hubs_cells$Sample %in% b3tls_samplename, ]

# TODO uncomment if needed, a bit to heavy just for checking 
# fwrite(hubs_cells, out_path_hubs_cells)

# filter to cells within ROIs ---------------------------------------------
 
hubs_cells_inroi <- lapply(unique(hubs_cells$Sample), function(sample_name){
  print(sample_name)
  
  hubs_cells_sample <- hubs_cells[hubs_cells$Sample == sample_name, ]
  
  # count cells within all ROIs in the sample - now in metadata
  roi_coords_sample <- metadt[metadt$Sample == sample_name, c("sample_roi", "Roi_geomx", 'Annotation_cell', "roi_c1_X_cycif", "roi_c1_Y_cycif", 
                                                              "roi_c2_X_cycif", "roi_c2_Y_cycif", "roi_c3_X_cycif",
                                                              "roi_c3_Y_cycif", "roi_c4_X_cycif", "roi_c4_Y_cycif")]
  
  cells_in_roi_all <- apply(roi_coords_sample, 1, function(row){
    
    # find cells within range
    # coordinates are not longer rectangles, they're a bit rotated - the cells are found inside longer edges of rectangle
    cells_in_roi <- dplyr::filter(hubs_cells_sample,
                                  as.numeric(X_centroid_px) >= min(as.numeric(row[['roi_c1_X_cycif']]), as.numeric(row[['roi_c4_X_cycif']])) &
                                    as.numeric(X_centroid_px) <= max(as.numeric(row[['roi_c2_X_cycif']]), as.numeric(row[['roi_c3_X_cycif']])) &
                                    as.numeric(Y_centroid_px) >= min(as.numeric(row[['roi_c4_Y_cycif']]), as.numeric(row[['roi_c3_Y_cycif']])) &
                                    as.numeric(Y_centroid_px) <= max(as.numeric(row[['roi_c1_Y_cycif']]), as.numeric(row[['roi_c2_Y_cycif']])))

    cells_in_roi <- cbind(cells_in_roi, as.data.frame(lapply(row, rep, nrow(cells_in_roi))))
    return(cells_in_roi)
  })
  
  cells_in_roi_all <- do.call(rbind, cells_in_roi_all)
  return(cells_in_roi_all)
})

hubs_cells_inroi <- do.call(rbind, hubs_cells_inroi)

# clean labels
hubs_cells_inroi <- mutate(hubs_cells_inroi, across(c(interaction_hub_ids, network_hub_id), ~replace(., . ==  '' , NA)))

fwrite(hubs_cells_inroi, out_path_hubs_cells_inroi)

#########################################################################################
#########################################################################################
# explore hubs vs ROI

# total number of roi:
length(unique(metadt$sample_roi))

# nr of ROI with any cells phenotyped
length(unique(hubs_cells_inroi$sample_roi))

# rois without any cells found
setdiff(metadt$sample_roi, hubs_cells_inroi$sample_roi) 

hubs_cells_inroi_innet <- hubs_cells_inroi[!is.na(hubs_cells_inroi$network_hub_id), ]
hubs_cells_inroi_inhub <- hubs_cells_inroi[!is.na(hubs_cells_inroi$interaction_hub_ids), ]
hubs_cells_inroi_incomm <- hubs_cells_inroi[!is.na(hubs_cells_inroi$community_id), ]

table(hubs_cells_inroi$Sample, hubs_cells_inroi$Roi_geomx)
table(hubs_cells_inroi_innet$Sample, hubs_cells_inroi_innet$roi_name)
table(hubs_cells_inroi_inhub$Sample, hubs_cells_inroi_inhub$roi_name)
table(hubs_cells_inroi_incomm$Sample, hubs_cells_inroi_incomm$Roi_geomx)

# nr of unique interaction hubs in each ROI
hubs_unique_roi <- distinct(hubs_cells_inroi_inhub, sample_roi, interaction_hub_ids)
sort(table(hubs_unique_roi$sample_roi))
length(unique(hubs_unique_roi$sample_roi))
hist(unlist(table(hubs_unique_roi$sample_roi)))

# nr of unique comms in each ROI
comm_unique_roi <- distinct(hubs_cells_inroi_incomm, sample_roi, community_id, community_cluster)
sort(table(comm_unique_roi$sample_roi))
length(unique(comm_unique_roi$sample_roi))
hist(unlist(table(comm_unique_roi$sample_roi)))

# nr of unique networks in each ROI
net_unique_roi <- distinct(hubs_cells_inroi_innet, sample_roi, network_hub_type)
sort(table(net_unique_roi$sample_roi))
length(unique(net_unique_roi$sample_roi))
hist(unlist(table(net_unique_roi$sample_roi)))

# rois with cells but without any communities
roi_wo_comm <- setdiff(hubs_cells_inroi$sample_roi, comm_unique_roi$sample_roi) 

sort(table(hubs_cells_inroi$sample_roi[hubs_cells_inroi$sample_roi %in% roi_wo_comm]))
