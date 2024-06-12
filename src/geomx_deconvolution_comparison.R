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

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'deconvolution'), showWarnings = T, recursive = T)

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')

# load deconvolution results ----------------------------------------------

#TODO compare different versions of decon run + cell_types vs mid_lvl_ct

geomx_obj <- readRDS(input_rds_path)
meta_names <- c('dcc_filename', 'Patient', 'Segment', 'Sample', 'NACT status', 'Annotation_cell')

bp_res_files <- list.files(file.path(output_dir, 'deconvolution', 'bp'), full.names = T, pattern = 'RDS')
sd_res_files <- list.files(file.path(output_dir, 'deconvolution', 'sd'), full.names = T, pattern = 'RDS')

#TODO iterate through all files
res_path <- c(bp_res_files, sd_res_files)[4]
####################################
# get deconvolution dirnames etc
deconv_type <- gsub('_.*', '', basename(res_path))
ct_type <- ifelse(grepl('cell_type', basename(res_path)), 'cell_type', 
                  ifelse(grepl('mid_lvl_ct', basename(res_path)), 'mid_lvl_ct', NA))

out_dirname <- gsub(paste0(deconv_type, '_res_|\\.RDS'), '', basename(res_path))
outp_plot_dir <- file.path(output_dir,'deconvolution', deconv_type, out_dirname)

dir.create(outp_plot_dir, showWarnings = T, recursive = T)
dir.create(file.path(outp_plot_dir, 'per_sample'), showWarnings = T, recursive = T)

#############################
# read data
ct_res <- readRDS(res_path)

if(deconv_type == 'bp'){
  ct_frac <- get.fraction (bp=ct_res,
                           which.theta="final",
                           state.or.type="type")
} else if(deconv_type == 'sd'){
  ct_frac <- pData(ct_res)[, 'prop_of_all']
}

ct_frac <- rownames_to_column(as.data.frame(ct_frac), 'dcc_filename')
ct_frac <- left_join(ct_frac, sData(geomx_obj)[, meta_names],
                     by = 'dcc_filename')

#######################################
# # TODO mask ct_frac results if cv > 0.2-0.5 (0.1 thr for bulk, 0.5 for Visium, GeoMx should be in the middle)
# # histogram suggests 0.5 as thr
# ct_frac_cv <- bprism_res@posterior.theta_f@theta.cv
# 
# # extract posterior mean of cell type-specific gene expression count matrix Z  
# # TODO normalise it!
# #Clustering bulk samples by theta or Z (Z can be normalized by vst(round(t(Z.tumor))), 
# # using the vst function from the DESeq2 package.)
# tumor_gene_exp <- get.exp (bp=bprism_res,
#                            state.or.type="type",
#                            cell.name="tumor")

###############################
# check immune and total frac of cells
cells_stroma <- c('Fibroblasts', 'Endothelial cells')
ct_frac$stroma <- rowSums(ct_frac[, cells_stroma])

if(ct_type == 'cell_type'){
  cells_immune <- c('Regulatory T cells', 'Memory B cells', 'CD16- NK cells',
                    'Tem/Trm cytotoxic T cells', 'Tcm/Naive helper T cells', 'Macrophages',
                    'Mast cells', 'Migratory DCs', 'Plasma cells', 'ILC', 'pDC',
                    'CD16+ NK cells', 'NK cells', 'Naive B cells', 'DC1', 'Classical monocytes')
  if('Type 17 helper T cells' %in% colnames(ct_frac)){cells_immune <- c(cells_immune, 'Type 17 helper T cells')} #TODO check if this or bp
  ct_frac$macro_mono <- ct_frac$Macrophages + ct_frac$`Classical monocytes`
  ct_frac$immune <- rowSums(ct_frac[, cells_immune])
  ct_frac$tot <- ct_frac$tumor + ct_frac$stroma + ct_frac$immune
  
  cd8_ct <- 'Tem/Trm cytotoxic T cells'
} else if(ct_type == 'mid_lvl_ct'){
  cells_immune <- c("Tcells", "Bcells", "NKcells", "Macrophages", "Mast cells", "DCs")
  ct_frac$immune <- rowSums(ct_frac[, cells_immune])
  ct_frac$tot <- ct_frac$tumor + ct_frac$stroma + ct_frac$immune + ct_frac$other
  cd8_ct <- 'Tcells'
}

macro_ct <- 'Macrophages'

#############################
# compare basic ct in segment

sapply(c('tumor', 'stroma', 'immune'), function(c){
  ggplot(data = ct_frac, aes(x = Segment, y = get(c))) +
    geom_violin() +
    geom_point(position= position_jitter(),
               size= 0.2, alpha = 0.6) +
    ylab(c) +
    ggtitle(c)
  
  ggsave(file.path(outp_plot_dir, paste0('segment_count_', c, '.pdf')),
         width = 1500, height = 2000, unit = 'px')
})

#############################
# compare macro vs Tcells 

sapply(c('tumor', 'stroma'), function(seg){
  ct_frac_seg <- ct_frac[ct_frac$Segment == seg, ]
  
  ggplot(data = ct_frac_seg, aes(x = get(macro_ct), y = get(cd8_ct), shape = Annotation_cell, color = Sample)) +
    geom_point(alpha = 0.5, size = 2) +
    xlab(macro_ct) +
    ylab(cd8_ct)
  
  ggsave(file.path(outp_plot_dir, paste0('macro_cd8_', seg, '_persample.pdf')),
         width = 1500, height = 2000, unit = 'px')
  
  ggplot(data = ct_frac_seg, aes(x = get(macro_ct), y = get(cd8_ct), color = Annotation_cell)) +
    geom_point(alpha = 0.5, size = 2) +
    xlab(macro_ct) +
    ylab(cd8_ct)
  
  ggsave(file.path(outp_plot_dir, paste0('macro_cd8_', seg, '_peranno.pdf')),
         width = 1500, height = 2000, unit = 'px')

  if(ct_type == 'cell_type'){
    # monocytes
    ggplot(data = ct_frac_seg, aes(x = `Classical monocytes`, y = get(cd8_ct), color = Annotation_cell)) +
      geom_point(alpha = 0.5, size = 2)+
      ylab(cd8_ct)
    
    ggsave(file.path(outp_plot_dir, paste0('mono_cd8_', seg, '_peranno.pdf')),
           width = 1500, height = 2000, unit = 'px')
    
    # mono + macro
    ggplot(data = ct_frac_seg, aes(x = macro_mono, y = get(cd8_ct), color = Annotation_cell)) +
      geom_point(alpha = 0.5, size = 2)+
      ylab(cd8_ct)
    
    ggsave(file.path(outp_plot_dir, paste0('macro_mono_cd8_', seg, '_peranno.pdf')),
           width = 1500, height = 2000, unit = 'px')
  }
  
})
#############################
# macro vs tcells per sample

sapply(unique(ct_frac$Sample), function(s){
  ct_frac_sample <- ct_frac[ct_frac$Sample == s, ]
  
  ggplot(data = ct_frac_sample, aes(x = get(macro_ct), y = get(cd8_ct), color = Annotation_cell, shape = Segment)) +
    geom_point(alpha = 0.5) +
    xlab(macro_ct) +
    ylab(cd8_ct) +
    ggtitle(s)
  
  ggsave(file.path(outp_plot_dir, 'per_sample', paste0(s, '_macro_cd8.pdf')),
         width = 1500, height = 1000, unit = 'px')
  
  if(ct_type == 'cell_type'){
    ggplot(data = ct_frac_sample, aes(x = macro_mono, y = get(cd8_ct), color = Annotation_cell, shape = Segment)) +
      geom_point(alpha = 0.5) +
      ylab(cd8_ct) +
      ggtitle(s)
    
    ggsave(file.path(outp_plot_dir, 'per_sample',  paste0(s, '_macro_mono_cd8.pdf')),
           width = 1500, height = 1000, unit = 'px')
  }
  
})

###########################################
# make hmaps 

col_fun = colorRamp2(c(0, 1), c("white", "red"))

png(filename = file.path(outp_plot_dir, 'hmap_all.png'), width=1000, height=750)
Heatmap(t(as.matrix(ct_frac[, c(cells_immune, cells_stroma, 'tumor')])), col = col_fun) %v%
  HeatmapAnnotation(segment = ct_frac$Segment, 
                    col = list(segment = c("stroma" = "green", "tumor" = "blue")))
dev.off()

sapply(c('stroma', 'tumor'), function(seg){
  ct_frac_seg <- ct_frac[ct_frac$Segment == seg, ]
  
  png(filename = file.path(outp_plot_dir, paste0('hmap_', seg, '_all.png')), width=1000, height=750)
  print(Heatmap(t(as.matrix(ct_frac_seg[, c(cells_immune, cells_stroma, 'tumor')])), col = col_fun) %v%
    HeatmapAnnotation(roi_type = ct_frac_seg$Annotation_cell,
                      sample = ct_frac_seg$Sample,
                      nact_status = ct_frac_seg$`NACT status`))
  dev.off()
  
  png(filename = file.path(outp_plot_dir, paste0('hmap_', seg, '_immune.png')), width=1000, height=750)
  print(Heatmap(t(as.matrix(ct_frac_seg[, c(cells_immune)])), col = col_fun) %v%
    HeatmapAnnotation(roi_type = ct_frac_seg$Annotation_cell,
                      sample = ct_frac_seg$Sample,
                      nact_status = ct_frac_seg$`NACT status`))
  dev.off()
})

#############################################
# stacked barplots with cell composition

# wide to long
ct_frac$roi_name <- paste0(ct_frac$Annotation_cell, '_', ct_frac$Sample)
ct_frac_long <- melt(ct_frac[, c(meta_names, 'roi_name', 'tumor', cells_stroma, cells_immune)],
                             id.vars = c(meta_names, 'roi_name'))

sapply(c('stroma', 'tumor'), function(seg){
  ct_frac_long_seg <- ct_frac_long[ct_frac_long$Segment == seg, ]
  ct_frac_long_seg <- arrange(ct_frac_long_seg, Annotation_cell)
  
  ggplot(data = ct_frac_long_seg, aes(x = dcc_filename, y = value, fill = variable)) +
    geom_bar(position="fill", stat="identity") +
    ggtitle(seg) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4))
  
  ggsave(file.path(outp_plot_dir, paste0('barplot_', seg, '.png')),,
         width = 1500, height = 1000, unit = 'px')
  
  ct_frac_long_seg_immune <- filter(ct_frac_long_seg, variable %in% cells_immune)
  
  ggplot(data = ct_frac_long_seg_immune, aes(x = dcc_filename, y = value, fill = variable)) +
    geom_bar(position="fill", stat="identity") +
    ggtitle(seg) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4))
  
  ggsave(file.path(outp_plot_dir, paste0('barplot_', seg, '_immune.png')),,
         width = 1500, height = 1000, unit = 'px')
})

##########
# paired dotplot compare 2 deconvolution types
bp_path <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2/deconvolution/bp/bp_res_mid_lvl_ct_45.RDS'
sd_path <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2/deconvolution/sd/sd_res_mid_lvl_ct_45_filt_scrna_geomx.RDS'

macro_ct <- 'Macrophages'
cd8_ct <- 'Tcells'

bp_mid_lvl_ct <- get.fraction (bp= readRDS(bp_path),
                               which.theta="final",
                               state.or.type="type")

sd_mid_lvl_ct <- pData(readRDS(sd_path))[, 'prop_of_all']

identical(rownames(bp_mid_lvl_ct), rownames(sd_mid_lvl_ct))

#bp_mid_lvl_ct <- bp_mid_lvl_ct[, c(macro_ct, cd8_ct, 'tumor')]
colnames(bp_mid_lvl_ct) <- paste0('bp_', colnames(bp_mid_lvl_ct))
bp_mid_lvl_ct <- rownames_to_column(as.data.frame(bp_mid_lvl_ct), 'dcc_filename')

#sd_mid_lvl_ct <- sd_mid_lvl_ct[, c(macro_ct, cd8_ct, 'tumor')]
colnames(sd_mid_lvl_ct) <- paste0('sd_', colnames(sd_mid_lvl_ct))
sd_mid_lvl_ct <- rownames_to_column(as.data.frame(sd_mid_lvl_ct), 'dcc_filename')

deconv_both <- left_join(bp_mid_lvl_ct, sd_mid_lvl_ct)

deconv_long <- melt(deconv_both, id.vars = c('dcc_filename'))
deconv_long$cell_type <- gsub('.*_', '', deconv_long$variable)
deconv_long$deconv_type <- gsub('_.*', '', deconv_long$variable)

# make a paired plot
ggplot(deconv_long, aes(x = deconv_type, y = value)) + 
  geom_boxplot(aes(fill = deconv_type), alpha = .2) +
  geom_line(aes(group = dcc_filename)) + 
  geom_point(size = 1) + 
  facet_wrap(~ cell_type)

ggsave(file.path(output_dir, 'deconvolution', paste0('bp_sd_mid_lvl_ct_deconv_comparison_all.png')),,
       width = 1500, height = 1000, unit = 'px')

# separately tumor-stroma
geomx_meta <- sData(readRDS(input_rds_path))[, c('dcc_filename', 'Segment')]
deconv_long <- left_join(deconv_long, geomx_meta)

sapply(c('tumor', 'stroma'), function(seg){
  ggplot(deconv_long[deconv_long$Segment == seg, ], aes(x = deconv_type, y = value)) +
    geom_boxplot(aes(fill = deconv_type), alpha = .2) +
    geom_line(aes(group = dcc_filename)) +
    geom_point(size = 1) +
    facet_wrap(~ cell_type) +
    ggtitle(seg)

  ggsave(file.path(output_dir, 'deconvolution', paste0('bp_sd_mid_lvl_ct_deconv_comparison_all_', seg, '.png')),,
         width = 1500, height = 1000, unit = 'px')
})


###########################################################################
###########################################################################
#############################################################################
# adjust format from adjacent cycif cell counts
#TODO make paired dotplot comparing bp and sd

cycif_counts_path <- '/home/ad/P-drive/h345/afarkkilab/Data/9-EyeMT/Data_integration/Quantified_cells/cell_count_per_AOI.csv'

cycif_counts <- fread(cycif_counts_path)

cycif_counts <- cycif_counts %>%
  mutate(all_intumor_nr = rowSums(dplyr::select(., starts_with("Intumor")))) %>%
  mutate(all_instroma_nr = rowSums(dplyr::select(., starts_with("Instroma")))) %>%
  filter(Note != 'ROTATE')

roi_meta <- sData(readRDS(input_rds_path))[, c('dcc_filename', 'Sample', 'Segment', 'Roi')] %>%
  mutate(Roi = as.numeric(gsub('^0*', '', Roi)))

cycif_counts <- left_join(cycif_counts, roi_meta, by = c('Sample', 'Roi')) %>%
  mutate(cells_nr_aoi = ifelse(Segment == 'tumor', rowSums(dplyr::select(., starts_with("Intumor"))),
         ifelse(Segment == 'stroma', rowSums(dplyr::select(., starts_with("Instroma"))), NA))) %>%
  mutate(frac_aoi_tumor = ifelse(Segment == 'tumor', Intumor_Tumor/cells_nr_aoi, ifelse(Segment == 'stroma',
                                                                         Instroma_Tumor/cells_nr_aoi, NA))) %>%
  mutate(frac_aoi_cd8 = ifelse(Segment == 'tumor', Intumor_CD8Tcells/cells_nr_aoi, ifelse(Segment == 'stroma',
                                                                                          Instroma_CD8Tcells/cells_nr_aoi, NA))) %>%
  mutate(frac_aoi_macro = ifelse(Segment == 'tumor', Intumor_Macrophages/cells_nr_aoi, ifelse(Segment == 'stroma',
                                                                                              Instroma_Macrophages/cells_nr_aoi, NA)))

# add counts from deconvolution
bp_mid_lvl_ct <- get.fraction (bp= readRDS('/media/iganiemi/T7-iga/st/geomx-processing/results/nact2/deconvolution/bp/bp_res_mid_lvl_ct_45.RDS'),
                               which.theta="final",
                               state.or.type="type")

#colnames(bp_mid_lvl_ct) <- paste0('deconv_', colnames(bp_mid_lvl_ct))
bp_mid_lvl_ct <- rownames_to_column(as.data.frame(bp_mid_lvl_ct), 'dcc_filename')
bp_mid_lvl_ct <- bp_mid_lvl_ct[, c('dcc_filename', 'tumor', 'Tcells', 'Macrophages')]
colnames(bp_mid_lvl_ct) <- c('dcc_filename', 'deconv_tumor', 'deconv_cd8', 'deconv_macro')
# bp_cell_type <- get.fraction (bp= readRDS('/media/iganiemi/T7-iga/st/geomx-processing/results/nact2/deconvolution/bp/bp_res_cell_type_45.RDS'),
#                                which.theta="final",
#                                state.or.type="type")

cycif_counts <- left_join(cycif_counts, bp_mid_lvl_ct)

# make it long
cycif_counts_long <- cycif_counts[, c('dcc_filename', 'Sample', 'Roi', 'Segment', 'cells_nr_aoi',
                                      "frac_aoi_tumor", "frac_aoi_cd8", "frac_aoi_macro",
                                      "deconv_tumor", "deconv_cd8", "deconv_macro")]
cycif_counts_long <- melt(cycif_counts_long, id.vars = c('dcc_filename', 'Sample', 'Roi', 'Segment', 'cells_nr_aoi'))
cycif_counts_long$cell_type <- gsub('.*_', '', cycif_counts_long$variable)
cycif_counts_long$count_type <- gsub('_.*', '', cycif_counts_long$variable)

# make a paired plot
ggplot(cycif_counts_long, aes(x = count_type, y = value)) + 
  geom_boxplot(aes(fill = count_type), alpha = .2) +
  geom_line(aes(group = dcc_filename)) + 
  geom_point(size = 2) + 
  facet_wrap(~ cell_type)

ggsave(file.path(output_dir, 'deconvolution', paste0('cycif_count_deconv_comparison_all.png')),,
       width = 1500, height = 1000, unit = 'px')

# only S084_post where are cells in both
ggplot(cycif_counts_long[cycif_counts_long$Sample == 'S084_post', ], aes(x = count_type, y = value)) + 
  geom_boxplot(aes(fill = count_type), alpha = .2) +
  geom_line(aes(group = dcc_filename)) + 
  geom_point(size = 2) + 
  facet_wrap(~ cell_type)+
  ggtitle('S084_post')

ggsave(file.path(output_dir, 'deconvolution', paste0('cycif_count_deconv_comparison_S084_post.png')),,
       width = 1500, height = 1000, unit = 'px')

#############################################
#############################################
# tumor/stroma separately
sapply(c('tumor', 'stroma'), function(seg){
  ggplot(cycif_counts_long[cycif_counts_long$Segment == seg, ], aes(x = count_type, y = value)) + 
    geom_boxplot(aes(fill = count_type), alpha = .2) +
    geom_line(aes(group = dcc_filename)) + 
    geom_point(size = 2) + 
    facet_wrap(~ cell_type) +
    ggtitle(seg)
  
  ggsave(file.path(output_dir, 'deconvolution', paste0('cycif_count_deconv_comparison_all_', seg, '.png')),,
         width = 1500, height = 1000, unit = 'px')
  
  # only S084_post where are cells in both
  ggplot(cycif_counts_long[(cycif_counts_long$Sample == 'S084_post' & cycif_counts_long$Segment == seg), ], 
         aes(x = count_type, y = value)) + 
    geom_boxplot(aes(fill = count_type), alpha = .2) +
    geom_line(aes(group = dcc_filename)) + 
    geom_point(size = 2) + 
    facet_wrap(~ cell_type)+
    ggtitle(paste0(seg, 'S084_post'))
  
  ggsave(file.path(output_dir, 'deconvolution', paste0('cycif_count_deconv_comparison_S084_post_', seg, '.png')),,
         width = 1500, height = 1000, unit = 'px')
})

#############################################
#############################################
# old code
# 
# heatmap(sweep(decon_res_ext@experimentData@other$SpatialDeconMatrix, 1, apply(decon_res_ext@experimentData@other$SpatialDeconMatrix, 1, max), "/"),
#         labRow = NA, margins = c(10, 5))
# 
# res_cols <- c("beta", "p", "t", "se", "prop_of_all", "prop_of_nontumor", "sigmas")
# 
# res_custom <- pData(decon_res_custom)[, res_cols]
# 
# decon <- decon_res_custom
# outp_subdir <- 'custom_cell_types_1x' # 'custom_mid_lvl_ct_5x', 'ext', 'basic' oc_scrnaseq_ref_cell_type_1x_sct
# 
# 
# dt <- sData(decon)[, c('dcc_filename', 'Patient', 'Segment', 'Sample','Nuclei', 'NACT status', 'Annotation_cell')]
# dt$sample_name <- rownames(dt)
# dt <- cbind(dt, res_mtx)
# 
# dt_stroma <- dt[dt$Segment == 'stroma', ]
# dt_tumor <- dt[dt$Segment == 'tumor', ]
# 
# 
# 
# dir.create(file.path(output_dir, 'deconvolution', outp_subdir,  'per_sample'), showWarnings = T, recursive = T)
# 
# macro_name <- 'Macrophages' # 'macrophages' in softTME
# cd8_name <- 'Tem/Trm cytotoxic T cells' # 'Tcells' in mid lvl ct 'CD8.Tcells' in softTME
# 
# res_type <- 'beta' # c('beta', 'cell_counts', 'prop_all')
# 
# #######################
# 
# pData(decon)$sample_name <-  paste0(sData(decon)[['Sample']], '_', 
#                                     sData(decon)[['Roi']], '_', 
#                                     sData(decon)[['Segment']], '_',
#                                     sData(decon)[['Annotation_cell']])
# 
# # mein got, the default functions from spatialdecon are so bad
# # cell counts
# o = hclust(dist(t(decon$cell.counts$cell.counts)))$order
# layout(mat = (matrix(c(1, 2), 1)), widths = c(7, 3))
# pdf(file=file.path(output_dir, 'deconvolution',outp_subdir, 'cell_counts_custom.pdf'), width=1500, height=1000)
# TIL_barplot(t(decon$cell.counts$cell.counts[, o]), draw_legend = TRUE, 
#             cex.names = 0.5)
# dev.off()
# 
# # proportions of cells
# temp = replace(decon$prop_of_nontumor, is.na(decon$prop_of_nontumor), 0)
# o = hclust(dist(temp[decon$Segment == "stroma",]))$order
# pdf(file=file.path(output_dir, 'deconvolution',outp_subdir, 'cell_prop_custom.pdf'), width=1500, height=1000)
# TIL_barplot(t(decon$prop_of_nontumor[decon$Segment == "stroma",])[, o], 
#             draw_legend = TRUE, cex.names = 0.5)
# dev.off()
# 
# 
# 
