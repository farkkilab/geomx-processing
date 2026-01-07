library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(tibble)
library(ComplexHeatmap)
library(GeomxTools)
library(readxl)


# define vars -------------------------------------------------------------

proj_dir <<- '~/Documents/phd/st'
data_dir <<- '~/Documents/phd/st/data/geomx/batch123/' # batch1 2 and 3
anno_path <<- file.path(data_dir, 'metadata', 'dcc_metadata_batch123_no_tls_cleaned.xlsx') #batch1 and 2 and 3
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808') # batch123
geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS') 

bp_cellcounts_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', 'bp_res_mid_lvl_ct_updated_ct_fraction.csv')
sd_cellcounts_path <- file.path(output_dir, 'deconvolution', 'spatial_decon', 'sd_res_mid_lvl_ct_updated_geomxfiltpc_ct_fraction.csv')

cycif_roi_cellnr_path <- file.path(output_dir, "cycif_integration", "batch2_roi_cellnr.csv")
hubs_inroi_path <- file.path(output_dir, "cycif_integration", "batch2_hubs_cells_inroi_dt17191719_ct15.csv")

##############
meta_names <- c('dcc_filename', 'Patient', 'Sample', 'Site', 'NACT_status', 'Annotation_cell', 'Roi', 
                'Segment_geomx',  'Segment') # 'main_batch_nr'

meta_names <- c('dcc_filename', 'Sample', 'Annotation_cell', 'Roi', 
                'Segment_geomx',  'Segment', 'Roi_geomx_original') 

source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))

outp_plot_dir <- file.path(output_dir, 'cycif_integration', 'ct_fractions')
dir.create(outp_plot_dir, recursive = T)

ct_names_all <- c("tumor", "Bcells", "Tcells_CD4", "Tcells_other", "Tcells_CD8", 
                  "Fibroblasts_Mesothelial", "Macrophages_Monocytes", "Mast_cells",
                  "NKcells", "Endothelial_cells", "DCs")

ct_names_stroma <- c("Fibroblasts_Mesothelial", "Endothelial_cells")

ct_names_immune <- c("Bcells", "Tcells_CD4", "Tcells_other", "Tcells_CD8", 
                     "Macrophages_Monocytes", "Mast_cells",
                     "NKcells", "DCs")

# load geomx, merge with cleaned metadata ---------------------------------
# TODO run once again in 1811 with already cleaned metadata and just load meta from geomx
meta_geomx <- pData(readRDS(geomx_norm_batch_eff_rm_path))
meta_cleaned <- read_excel(anno_path)

metadt <- left_join(meta_geomx[, c('dcc_filename', 'Slide_Name')], meta_cleaned) # join ensuring order
metadt$sample_roi <- paste0(metadt$Sample, '_', metadt$Roi)
metadt$sample_roi_orig <- paste0(metadt$Sample, '_', metadt$Roi_geomx_original)

rm(meta_geomx)
rm(meta_cleaned)

# TODO --------------------------------------------------------------------

# TODO compare cell nr (from meta) with cells within roi
# TODO compare sd/bp/pheno cell nr/fractions total/immune
# TODO compare this with hubs labels
# TODO compare with deconv made per roi

# load bind and clean deconv ct fractions ---------------------------------

# max nr of cells = 300 so 0.005 cell fraction is 1 cell/200 cells 1,5 cell/300 cells
min_frac <- 0.005

deconv_list <- list(bp = bp_cellcounts_path, sd = sd_cellcounts_path)

ct_frac_long_both <- lapply(1:length(deconv_list), function(n){
  
  ct_frac_deconv <- fread(deconv_list[[n]], select = c('dcc_filename', ct_names_all))
  ct_frac_deconv <- left_join(ct_frac_deconv, metadt[, meta_names], by = 'dcc_filename')
  
  ct_frac_long <- melt(setDT(ct_frac_deconv), id.vars = meta_names, variable.name = "cell_type")
  
  ct_frac_long$deconv_type <- names(deconv_list[n])
  ct_frac_long$sample_roi <- paste0(ct_frac_long$Sample, '_', ct_frac_long$Roi)
  
  # move unreliable predictions to 0
  ct_frac_long$value_clean <- ifelse(ct_frac_long$value <= min_frac | is.na(ct_frac_long$value), 0, ct_frac_long$value)
  
  ct_frac_long <- dplyr::rename(ct_frac_long, ct_fraction = value, ct_fraction_clean = value_clean)
  
  # add ct number
  ct_frac_long <- left_join(ct_frac_long, metadt[, c('dcc_filename', 'Nuclei')])
  ct_frac_long$ct_number <- round(ct_frac_long$ct_fraction_clean * ct_frac_long$Nuclei)
  
  return(ct_frac_long)
})

ct_frac_long_both <- do.call(rbind, ct_frac_long_both)

# merge per ROI
ct_frac_long_both_roi <- ct_frac_long_both %>%
  group_by(sample_roi, deconv_type, cell_type) %>%
  summarise(Nuclei = sum(Nuclei), ct_number = sum(ct_number), ct_fraction_clean = mean(ct_fraction_clean)) %>%
  ungroup()

# in wide format ----------------------------------------------------------

ct_frac_list <- lapply(1:length(deconv_list), function(n){
  ct_frac <- as.data.frame(fread(deconv_list[[n]], select = c('dcc_filename', ct_names_all)))
  ct_frac <- mutate_at(ct_frac, all_of(ct_names_all), funs(ifelse((is.na(.) | . < 0.005), 0, .)))
  
  ct_frac$stroma <- rowSums(ct_frac[, ct_names_stroma])
  ct_frac$immune <- rowSums(ct_frac[, ct_names_immune])
  
  #calculate cell nr
  ct_nr <- left_join(ct_frac, metadt[, c('dcc_filename', 'Nuclei')])
  ct_nr <- mutate_at(ct_nr, all_of(c(ct_names_all, 'stroma', 'immune')), funs(round(. * Nuclei)))
  ct_nr <- dplyr::rename(ct_nr, total = Nuclei)
  
  colnames(ct_frac) <- c('dcc_filename', paste0(colnames(ct_frac)[2:length(colnames(ct_frac))], '_', names(deconv_list)[n]))
  colnames(ct_nr) <- c('dcc_filename', paste0(colnames(ct_nr)[2:length(colnames(ct_nr))], '_', names(deconv_list)[n]))
  
  return(list(ct_frac = ct_frac, ct_nr = ct_nr))
})

# TODO make it nicer with recursive joining
ct_frac_both <- left_join(ct_frac_list[[1]]$ct_frac, ct_frac_list[[2]]$ct_frac, by = 'dcc_filename')
ct_frac_both <- left_join(ct_frac_both, metadt[, meta_names], by = 'dcc_filename')
ct_frac_both$sample_roi <- paste0(ct_frac_both$Sample, '_', ct_frac_both$Roi)

ct_nr_both <- left_join(ct_frac_list[[1]]$ct_nr, ct_frac_list[[2]]$ct_nr, by = 'dcc_filename')
ct_nr_both <- left_join(ct_nr_both, metadt[, meta_names], by = 'dcc_filename')
ct_nr_both$sample_roi <- paste0(ct_nr_both$Sample, '_', ct_nr_both$Roi)


# merge per ROI

ct_frac_both_roi <- ct_frac_both %>%
  group_by(sample_roi) %>%
  summarise_at(c(paste0(c(ct_names_all, 'stroma', 'immune'), '_bp'), paste0(c(ct_names_all, 'stroma', 'immune'), '_sd')), 
               mean, na.rm = TRUE) %>%
  ungroup()

ct_nr_both_roi <- ct_nr_both %>%
  group_by(sample_roi) %>%
  summarise_at(c(paste0(c(ct_names_all, 'stroma', 'immune', 'total'), '_bp'), paste0(c(ct_names_all, 'stroma', 'immune', 'total'), '_sd')), 
               sum, na.rm = TRUE) %>%
  ungroup()

# calculate cell nr/fractions from phenotyped cells in cycif --------------

cycif_roi_cellnr <- fread(cycif_roi_cellnr_path)

# calculate cell fractions

cycif_roi_cell_fraq <- mutate_all(cycif_roi_cellnr, is.numeric, funs(./total_cell_nr))

cycif_roi_cell_frac <- apply(cycif_roi_cellnr[, -1], 2, function(x){x / cycif_roi_cellnr$total_cell_nr})
# TODO add dcc_filename and remove total cell nr