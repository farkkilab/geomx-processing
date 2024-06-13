# TODO check if all packages are needed
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(SpatialDecon)
library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(reshape2)
library(Seurat)
library(tibble)
library(BayesPrism)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)

# define variables --------------------------------------------------------

data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'
output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2'

input_rds_path <- file.path(output_dir, 'geomx_qc_norm.RDS')

deconv_path_mid <- file.path(output_dir, 'deconvolution', 'bp/bp_res_mid_lvl_ct_45.RDS')
deconv_path_ct <- file.path(output_dir, 'deconvolution', 'bp/bp_res_cell_type_45.RDS')

# make dirs and source functions ------------------------------------------

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')


# load deconvolution results ----------------------------------------------

geomx_obj <- readRDS(input_rds_path)
meta_names <- c('dcc_filename', 'Patient', 'Segment', 'Sample', 'NACT status', 'Annotation_cell', 'Roi')

cd8_ct <- 'Tcells'
macro_ct <- 'Macrophages'

deconv_res_mid <- readRDS(deconv_path_mid)
deconv_res_ct <- readRDS(deconv_path_ct)


frac_mid <- get.fraction (bp= deconv_res_mid,
                               which.theta="final",
                               state.or.type="type")

frac_mid <- rownames_to_column(as.data.frame(frac_mid), 'dcc_filename')
frac_mid <- left_join(frac_mid, sData(geomx_obj)[, c(meta_names)])

#TODO doubleneg only if doubleneg not lost in QC etc 
#TODO calculate how many relabeled per sample

# relabel rois per sample -------------------------------------------------

frac_sample_all <- lapply(unique(frac_mid$Sample), function(sample_name){
  frac_sample <- frac_mid[frac_mid$Sample == sample_name,]
  
  frac_sample_roi <- group_by(frac_sample, Roi) %>%
    summarise(macro_sum = sum(Macrophages)/2, cd8_sum = sum(Tcells)/2, nseg = n()) %>%
    filter(nseg == 2) # remove rois where one segment was removed due to qc 
  
  min_macro <- frac_sample_roi$Roi[which.min(frac_sample_roi$macro_sum)]
  min_tcell <- frac_sample_roi$Roi[which.min(frac_sample_roi$cd8_sum)]
  
  double_neg <- min_macro # !!!!!!! here double neg may contain a bit more cd8 than cd8+, but the difference is not big
  frac_sample_roi <- frac_sample_roi[frac_sample_roi$Roi != double_neg, ]
  
  # get min once again after removing doubleneg
  min_macro <- frac_sample_roi$Roi[which.min(frac_sample_roi$macro_sum)]
  min_tcell <- frac_sample_roi$Roi[which.min(frac_sample_roi$cd8_sum)]
  
  # relabel
  frac_sample$Annotation_cell_relabeled <- ifelse(frac_sample$Roi == double_neg, 'negCD8_negIBA1', 
                                                  ifelse(frac_sample$Roi == min_macro, 'posCD8_negIBA1',
                                                         ifelse(frac_sample$Roi == min_tcell, 'negCD8_posIBA1', 'posCD8_posIBA1')))
  
  return(frac_sample)
})

frac_sample_all <- do.call(rbind, frac_sample_all)

# check how many ROIs relabeled

relabeled_rois <- frac_sample_all[frac_sample_all$Annotation_cell != frac_sample_all$Annotation_cell_relabeled, ]

# save relabeled df

fwrite(frac_sample_all, file.path(output_dir, 'deconvolution', 'bp_mid_lvl_ct_relabeled_roi.csv'))
