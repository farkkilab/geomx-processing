# add AOI signal per ROI --------------------------------------------------

geomx_obj <- readRDS(geomx_qc_path)

# add signal for each Aoi within Roi 

mt <- sData(geomx_obj)
mt <- cbind(paste0(mt$Sample, '_', mt$Roi_geomx), mt)
colnames(mt)[1] <- 'sample_roi'
expr <- exprs(geomx_obj)

# list of vecors with dcc_names per each roi to sum
dcc_roi_list <- lapply(unique(mt$sample_roi), function(roi_name){
  dcc_list <- mt$dcc_filename[mt$sample_roi == roi_name]
  return(roi_name = dcc_list)
})

names(dcc_roi_list) <- unique(mt$sample_roi)

# for each roi, select its dcc columns and sum expr per gene 
roi_sum_dcc_signal <- lapply(names(dcc_roi_list), function(roi_name){
  dcc_names <- dcc_roi_list[[roi_name]]
  roi_dcc_all <- expr[, dcc_roi_list[[roi_name]]]
  
  if(length(dcc_names) > 1){
    roi_dcc_sum <- rowSums(roi_dcc_all)
    return(roi_dcc_sum)
  } else if(length(dcc_names) == 1){
    return(roi_dcc_all)
  } else{
    stop(paste0('no ddc found for ', roi_name))
  }
})

roi_sum_dcc_signal_mtx <- do.call(cbind, roi_sum_dcc_signal)

# ensure dimentions are ok
dim(roi_sum_dcc_signal_mtx)
stopifnot(nrow(expr) == nrow(roi_sum_dcc_signal_mtx))
stopifnot(ncol(roi_sum_dcc_signal) == length(unique(mt$sample_roi)))

colnames(roi_sum_dcc_signal_mtx) <- names(dcc_roi_list)

# adjust matadata
aoi_cols <- c('dcc_filename', 'Segment', 'area', 'nuclei', 'sizefact', 'sizefact_sp', 
              'NegGeoMean_Hs_R_NGS_WTA_v1.0', 'NegGeoSD_Hs_R_NGS_WTA_v1.0', 'LOQ', 'GenesDetected', 'GeneDetectionRate',
              'SampleID', 'Plate_ID', 'Well', "SeqSetId", "Raw", "Trimmed", "Stitched", "Aligned", "umiQ30", "rtsQ30",
              "DeduplicatedReads", "Roi", "Aoi", "NTC_ID", "NTC", "QCFlags", "Trimmed (%)", "Stitched (%)", "Aligned (%)",
              "Saturated (%)", "NegGeoMean")

mt_per_roi <- mt %>% group_by(sample_roi) %>%
  mutate(across(all_of(aoi_cols), ~ paste0(., collapse =';'))) %>%
  ungroup() %>%
  distinct()

# create new geomx object
newassay <- new.env(parent=geomx_obj@assayData)
newassay$exprs <- roi_sum_dcc_signal_mtx

geomx_obj_roibased <- geomx_obj
geomx_obj_roibased@assayData <- newassay
pData(geomx_obj_roibased) <- mt_per_roi[, 1:57]
geomx_obj_roibased@protocolData@data <- mt_per_roi[, 58:ncol(mt_per_roi)]

saveRDS(geomx_obj_roibased, file = geomx_qc_roibased_path)
