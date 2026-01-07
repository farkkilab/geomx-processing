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
hubs_inroi_path <- file.path(output_dir, "cycif_integration", "batch2_hubs_cells_inroi_dt17191719_ct15.csv")

##############
meta_names <- c('dcc_filename', 'Patient', 'Sample', 'Site', 'NACT_status', 'Annotation_cell', 'Roi', 
                'Segment_geomx',  'Segment') # 'main_batch_nr'

meta_names <- c('dcc_filename', 'Sample', 'Annotation_cell', 'Roi', 
                'Segment_geomx',  'Segment') # 'main_batch_nr'

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

rm(meta_geomx)
rm(meta_cleaned)

# TODO --------------------------------------------------------------------

# TODO same sample_roi between meta and hubs
# TODO calculate ct fractions sd+bp per ROI
# TODO compare cell nr (from meta) with cells within roi
# TODO calculate phenotyped cell fractions
# TODO compare sd/bp/pheno cell nr/fractions total/immune
# TODO compare this with hubs labels
# TODO compare with deconv made per roi

# load bind and clean deconv ct fractions ---------------------------------

# max nr of cells = 300 so 0.005 cell fraction is 1 cell/200 cells 1,5 cell/300 cells
min_frac <- 0.005

ct_frac_bp <- fread(bp_cellcounts_path, select = c(meta_names, ct_names_all))
ct_frac_sd <- fread(sd_cellcounts_path, select = c(meta_names, ct_names_all))

deconv_list <- list(bp = ct_frac_bp, sd = ct_frac_sd)

# TODO loop
# long
n <- 1
ct_frac_long <- melt(setDT(deconv_list[[n]]), id.vars = meta_names, variable.name = "cell_type")

ct_frac_long$deconv_type <- names(deconv_list[n])
ct_frac_long$sample_roi <- paste0(ct_frac_long$Sample, '_', ct_frac_long$Roi)

# move unreliable predictions to 0
ct_frac_long$value_clean <- ifelse(ct_frac_long$value <= min_frac | is.na(ct_frac_long$value), 0, ct_frac_long$value)

ct_frac_long <- dplyr::rename(ct_frac_long, ct_fraction = value, ct_fraction_clean = value_clean)

# add ct number
ct_frac_long <- left_join(ct_frac_long, metadt[, c('dcc_filename', 'Nuclei')])
ct_frac_long$ct_number <- round(ct_frac_long$ct_fraction_clean * ct_frac_long$Nuclei)

# in wide format ----------------------------------------------------------

deconv_path <- bp_cellcounts_path 

ct_frac <- as.data.frame(fread(deconv_path, select = c('dcc_filename', ct_names_all)))
ct_frac <- mutate_at(ct_frac, all_of(ct_names_all), funs(ifelse((is.na(.) | . < 0.005), 0, .)))

ct_frac$stroma <- rowSums(ct_frac[, ct_names_stroma])
ct_frac$immune <- rowSums(ct_frac[, ct_names_immune])

# add immune-other, to long, calculate cell nr, merge sd/bp/pheno, merge per roi

ct_frac_long <- melt(setDT(ct_frac), id.vars = "dcc_filename", variable.name = "cell_type")

ct_nr <- left_join(ct_frac, metadt[, c('dcc_filename', 'Nuclei')])
ct_nr <- mutate_at(ct_nr, all_of(c(ct_names_all, 'stroma', 'immune')), funs(round(. * Nuclei)))
ct_nr <- dplyr::rename(ct_nr, total = Nuclei)
