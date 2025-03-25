library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(GeomxTools)
library(BayesPrism)

source('~/Documents/phd/st/geomx-processing/src/geomx_utils.R')
source('~/Documents/phd/st/st-processing/src/visium_utils.R') # move get treatment hmap to geomx_utils


output_dir <- '~/Documents/phd/st/geomx-processing/results/nact/gsva/caf/'

##################
gsea_res_all_path <- '~/Documents/phd/st/geomx-processing/results/nact/gsva/caf/ssgsea_notnorm_all_caf.csv'
gsea_res_fib_path <- '~/Documents/phd/st/geomx-processing/results/nact/gsva/caf/ssgsea_notnorm_deconv_Fibroblasts_mid_lvl_ct_caf.csv'
deconv_path <- '~/Documents/phd/st/geomx-processing/results/nact/deconvolution/bp/bp_res_mid_lvl_ct_45.RDS'


geomx_obj <- '~/Documents/phd/st/geomx-processing/results/nact/geomx_qc_norm.RDS'

################
# deconv results

deconv <- readRDS(deconv_path)

ct_frac <- get.fraction (bp=deconv,
                         which.theta="final",
                         state.or.type="type")

ct_frac <- as.data.frame(ct_frac)
ct_frac$dcc_filename <- rownames(ct_frac)

fwrite(ct_frac, file.path(output_dir, 'ct_fractions.csv'))

#################
gsea_res_all <- fread(gsea_res_all_path)
gsea_res_fib <- fread(gsea_res_fib_path)


meta <- sData(readRDS(geomx_obj))
meta <- meta[, c('dcc_filename', 'Patient', 'Sample', 'Roi', 'Aoi')]
meta$ROI <- paste(meta$Sample, meta$Roi, sep = '_')

gsea_res_all <- left_join(gsea_res_all, meta[, c('dcc_filename', 'ROI')])
gsea_res_fib <- left_join(gsea_res_fib, meta[, c('dcc_filename', 'ROI')])
# 
# fwrite(gsea_res_all, gsea_res_all_path)
# fwrite(gsea_res_fib, gsea_res_fib_path)




###
# make plots
gsva_df <- gsea_res_fib
gsva_df <- filter(gsva_df, Segment == 'stroma')
gsva_df$pathway <- gsub('_symbol', '', gsva_df$pathway)
gsva_df_post <- filter(gsva_df, `NACT status` == 'post')

####################

pathway_boxplot(gsva_df, 'pathway', 'gsva_score', 'Annotation_cell', c('Segment'), 'ssgsea scores',
                file.path(output_dir, 'box_ssgsea_anno.pdf'))

pathway_boxplot(gsva_df, 'pathway', 'gsva_score', 'NACT status', c('Segment'), 'ssgsea scores',
                file.path(output_dir, 'box_ssgsea_nact_all.pdf'))

pathway_boxplot(gsva_df, 'pathway', 'gsva_score', 'NACT status', c('Annotation_cell', 'Segment'), 'ssgsea scores',
                file.path(output_dir, 'box_ssgsea_nact_peranno.pdf'))

pathway_boxplot(gsva_df_post, 'pathway', 'gsva_score', 'PFS', c('Segment'), 'ssgsea scores',
                file.path(output_dir, 'box_ssgsea_pfs_all.pdf'))

pathway_boxplot(gsva_df_post, 'pathway', 'gsva_score', 'PFS', c('Annotation_cell', 'Segment'), 'ssgsea scores',
                file.path(output_dir, 'box_ssgsea_pfs_peranno.pdf'))
