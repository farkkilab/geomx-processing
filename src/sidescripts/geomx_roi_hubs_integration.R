library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(tibble)
library(ComplexHeatmap)
library(GeomxTools)
library(readxl)

# TODO ensure # "S130_iOme_5" "S130_iOme_6" "S197_iOme_1" - everywhere in metadata - probably removed during QC

# define vars -------------------------------------------------------------

# max nr of cells = 300 so 0.005 cell fraction is 1 cell/200 cells 1,5 cell/300 cells
min_frac <- 0.005

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
source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))

outp_plot_dir <- file.path(output_dir, 'cycif_integration', 'ct_fractions')
dir.create(outp_plot_dir, recursive = T)

output_ct_frac_deconv_path <- file.path(output_dir, 'batch2_ct_frac_deconv.csv')
output_ct_frac_deconv_roi_path <- file.path(output_dir, 'batch2_ct_frac_deconv_roi.csv')
output_ct_frac_cycif_roi_path <- file.path(output_dir, 'batch2_ct_frac_cycif_roi.csv')
output_ct_frac_all_roi_path <- file.path(output_dir, 'batch2_ct_frac_all_roi.csv')

##############
meta_names <- c('dcc_filename', 'Sample', 'Annotation_cell', 'Roi', 
                'Segment_geomx',  'Segment', 'Roi_geomx_original') 

meta_names_per_roi <- c('Sample', 'Segment_geomx', 'Annotation_cell')

ct_names_all <- c("tumor", "Bcells", "Tcells_CD4", "Tcells_other", "Tcells_CD8", 
                  "Fibroblasts_Mesothelial", "Macrophages_Monocytes", "Mast_cells",
                  "NKcells", "Endothelial_cells", "DCs")

# ct from deconv counted as stroma in cycif
ct_names_stroma <- c("Fibroblasts_Mesothelial", "Endothelial_cells")

# main immune cells from deconv - also counted in cycif phenotyping
ct_names_immune <- c("Tcells_CD4", "Tcells_CD8", "Macrophages_Monocytes",
                     "NKcells", "DCs")

# additional cells from deconv not counted in phenotyping and should be treated as 'other'
ct_names_other <- c("Bcells", "Tcells_other", "Mast_cells")

# load geomx, merge with cleaned metadata ---------------------------------
# TODO run once again in 1811 with already cleaned metadata and just load meta from geomx
meta_geomx <- pData(readRDS(geomx_norm_batch_eff_rm_path))
meta_cleaned <- read_excel(anno_path)

metadt <- left_join(meta_geomx[, c('dcc_filename', 'Slide_Name')], meta_cleaned) # join ensuring order
metadt$sample_roi <- paste0(metadt$Sample, '_', metadt$Roi)
metadt$sample_roi_orig <- paste0(metadt$Sample, '_', metadt$Roi_geomx_original)

rm(meta_geomx)
rm(meta_cleaned)

# TODO in b3 ROIs with 2 different names from 2 different scans (16 -> 32)
# have to be adjusted in ROI coordinates + phenotyping to ensure correct ROI

# kk <- distinct(metadt[, c('sample_roi', 'sample_roi_orig')])
# kkk <- kk[kk$sample_roi_orig %in% kk$sample_roi_orig[which(duplicated(kk$sample_roi_orig))], ]

# TODO --------------------------------------------------------------------

# TODO compare cell nr (from meta) with cells within roi
# TODO compare sd/bp/pheno cell nr/fractions
# TODO compare this with hubs labels
# TODO compare with deconv made per roi

# load deconv and transform to long format --------------------------------

deconv_list <- list(bp = bp_cellcounts_path, sd = sd_cellcounts_path)

# loop through sd, bp results, convert to long
ct_frac_deconv_long <- lapply(1:length(deconv_list), function(n){
  ct_frac_deconv <- as.data.frame(fread(deconv_list[[n]], select = c('dcc_filename', ct_names_all)))
  ct_frac_deconv <- mutate_at(ct_frac_deconv, all_of(ct_names_all), funs(ifelse((is.na(.) | . < min_frac), 0, .)))
  
  ct_frac_deconv$stroma <- rowSums(ct_frac_deconv[, ct_names_stroma])
  ct_frac_deconv$immune <- rowSums(ct_frac_deconv[, ct_names_immune])
  ct_frac_deconv$other <- rowSums(ct_frac_deconv[, ct_names_other])
  ct_frac_deconv$immune_other <- rowSums(ct_frac_deconv[, c(ct_names_immune, ct_names_other)])
  
  ct_frac_long <- melt(setDT(ct_frac_deconv), id.vars = 'dcc_filename', variable.name = "cell_type")
  colnames(ct_frac_long)[which(colnames(ct_frac_long) == 'value')] <- paste0('ct_frac_', names(deconv_list[n]))
  
  return(ct_frac_long)
})

ct_frac_deconv_long <- do.call(left_join, ct_frac_deconv_long)

# add metadata
ct_frac_deconv_long <- left_join(ct_frac_deconv_long, metadt[, c(meta_names, 'Nuclei', 'sample_roi', 'sample_roi_orig')],
                               by = 'dcc_filename') %>%
  dplyr::rename(total_cell_nr_geomx = Nuclei)

# add cell nr
ct_frac_deconv_long$ct_nr_bp <- round(ct_frac_deconv_long$ct_frac_bp * ct_frac_deconv_long$total_cell_nr_geomx)
ct_frac_deconv_long$ct_nr_sd <- round(ct_frac_deconv_long$ct_frac_sd * ct_frac_deconv_long$total_cell_nr_geomx)

# merge per ROI
ct_frac_deconv_long_roi <- ct_frac_deconv_long %>%
  group_by(sample_roi, cell_type) %>%
  summarise(total_cell_nr_geomx = sum(total_cell_nr_geomx), ct_nr_bp = sum(ct_nr_bp), ct_nr_sd = sum(ct_nr_sd),
            ct_frac_bp = mean(ct_frac_bp), ct_frac_sd = mean(ct_frac_sd)) %>%
  ungroup()

fwrite(ct_frac_deconv_long, output_ct_frac_deconv_path)
fwrite(ct_frac_deconv_long_roi, output_ct_frac_deconv_roi_path)

# calculate cell nr/fractions from phenotyped cells in cycif --------------

ct_frac_cycif_roi <- as.data.frame(fread(cycif_roi_cellnr_path))
ct_frac_cycif_roi <- dplyr::rename(ct_frac_cycif_roi, immune = immune_cell_nr, immune_other = immune_other_cell_nr, total_cell_nr_cycif = total_cell_nr)

# to eliminate problems for b2. in b3 with non-fixed roi names will throw error
metadt_filt <- metadt[metadt$sample_roi_orig %in% ct_frac_cycif_roi$sample_roi, ]

# transform to long
ct_frac_cycif_long_roi <- melt(setDT(ct_frac_cycif_roi), id.vars = 'sample_roi', variable.name = "cell_type", value.name = "ct_nr_cycif")

# add total cell nr as column, remove total and calculate ct frac
ct_frac_cycif_long_roi <- left_join(ct_frac_cycif_long_roi, ct_frac_cycif_roi[, c('sample_roi', 'total_cell_nr_cycif')]) %>%
  filter(cell_type != 'total_cell_nr_cycif') %>%
  mutate(ct_frac_cycif = ct_nr_cycif / total_cell_nr_cycif) %>%
  dplyr::rename(sample_roi_orig = sample_roi)

# join with metadata sample_roi 
ct_frac_cycif_long_roi <- ct_frac_cycif_long_roi %>%
  filter(sample_roi_orig %in% metadt$sample_roi_orig) %>% # rmv roi not in metadata eg removed during qc
  left_join(distinct(metadt_filt[, c('sample_roi', 'sample_roi_orig')]))

fwrite(ct_frac_cycif_long_roi, output_ct_frac_cycif_roi_path)


# adjust hub dataframe ----------------------------------------------------

hubs_inroi <- as.data.frame(fread(hubs_inroi_path)) %>%
  dplyr::select(sample_roi, interaction_hub_type, network_hub_type) %>%
  dplyr::rename(sample_roi_orig = sample_roi) %>%
  dplyr::filter(interaction_hub_type != '' | network_hub_type != '') %>%
  distinct()

# merge network labels in roi
hubs_inroi <- hubs_inroi %>%
  group_by(sample_roi_orig) %>%
  mutate(network_hub_type = paste0(unique(network_hub_type), collapse = "|")) %>%
  distinct() %>%
  mutate(interaction_hub_type = paste0(unique(interaction_hub_type), collapse = "|")) %>%
  distinct()

hubs_inroi$interaction_hub_type <- gsub("^\\||\\|$", "", hubs_inroi$interaction_hub_type)
hubs_inroi$interaction_hub_type <- gsub("||", "|", hubs_inroi$interaction_hub_type, fixed = T)

hubs_inroi$sample_roi_orig <- gsub(" ", "", hubs_inroi$sample_roi_orig)

# merge all information together ------------------------------------------

# immune_other = all immune + 'other'

# join cycif and deconv ct fraction tables, add metadata and hubs info, reorder
ct_frac_all <- left_join(ct_frac_cycif_long_roi, ct_frac_deconv_long_roi, by = c('sample_roi', 'cell_type')) %>%
  left_join(distinct(metadt_filt[, c(meta_names_per_roi, 'sample_roi')]), by = 'sample_roi') %>%
  left_join(hubs_inroi, by = 'sample_roi_orig') %>%
  mutate(network_hub_type = ifelse(network_hub_type == '', NA, network_hub_type),
         interaction_hub_type = ifelse(interaction_hub_type == '', NA, interaction_hub_type)) %>%
  select(sample_roi, sample_roi_orig, !!meta_names_per_roi, interaction_hub_type, network_hub_type, cell_type,
         total_cell_nr_cycif, total_cell_nr_geomx, ct_nr_cycif, ct_nr_bp, ct_nr_sd, ct_frac_cycif, ct_frac_bp, ct_frac_sd)

fwrite(ct_frac_all, output_ct_frac_all_roi_path)
