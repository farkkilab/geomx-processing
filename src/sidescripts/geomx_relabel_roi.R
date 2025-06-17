# relabel ROI based on deconvolution results

# TODO check if all packages are needed
# library(NanoStringNCTools)
# library(GeomxTools)
# library(GeoMxWorkflows)
# library(SpatialDecon)
# library(reshape2)
# library(Seurat)
# library(tibble)
# library(BayesPrism)
# library(biomaRt)
# library(ComplexHeatmap)
# library(circlize)

library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(tibble)
library(ComplexHeatmap)

# distribution of fraction is continuous - hard to define + and -
# fraq from bp and sd are linearly correlated for most cell types 
# instead of labels we can investigate how certain gsea scores from deconv data correlates with nr of different ct
# eg IFNg GSEA in Tcells vs nr of Macrophages

# some separation may be done based on clustering of averages sd + bp zscores through interesting labels
# per roi and each geomx segment separately
# validation - comparison with ct counts per roi from cycif
# for now clustering is not so clear - will be better to use counts from cycif + ssgsea scores for validation?

# define variables --------------------------------------------------------

batch <- 'batch12'

proj_dir <<- '~/Documents/phd/st'

if(batch == 'batch1'){
  output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch1-1903') # batch1
} else if(batch == 'batch2'){
  output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch2-1903') # batch2
} else if(batch == 'batch12'){
  output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch12-1205-no-counts-shift2') # batch12
} else{
  stop('wrong batch nr')
}

#geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS')

bp_cellcounts_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', 'bp_res_mid_lvl_ct_ct_fraction.csv')
sd_cellcounts_path <- file.path(output_dir, 'deconvolution', 'spatial_decon', 'sd_res_mid_lvl_ct_geomxfiltpc_ct_fraction.csv')


meta_names <- c('dcc_filename', 'Patient', 'Sample', 'Site', 'NACT_status', 'Annotation_cell', 'Roi', 
                'Segment_geomx',  'Segment') # 'main_batch_nr'


ct_names <- c('Bcells', 'DCs', 'Endothelial cells', 'Fibroblasts', 'Macrophages', 'NKcells', 'Tcells', 'tumor')
cells_immune <- c('Bcells', 'DCs', 'Macrophages', 'NKcells', 'Tcells')

label_ct_names <- c('Tcells', 'Macrophages', 'DCs')

source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))

outp_plot_dir <- file.path(output_dir, 'sanity_check', 'ct_fractions')
dir.create(outp_plot_dir, recursive = T)

# make labels match
unite_anno <- function(anno){
  anno_united <- gsub('posCD8_posIBA1', 'CD8_Iba1', anno)
  anno_united <- gsub('negCD8_posIBA1', 'Iba1', anno_united)
  anno_united <- gsub('posCD8_negIBA1', 'CD8', anno_united)
  anno_united <- gsub('negCD8_negIBA1', 'noimmune', anno_united)
}

# load bind and clean deconv ct fractions ---------------------------------

# max nr of cells = 300 so 0.005 cell fraction is 1 cell/200 cells 1,5 cell/300 cells
min_frac <- 0.005

ct_frac_bp <- fread(bp_cellcounts_path)
ct_frac_bp$Annotation_cell_united <- unite_anno(ct_frac_bp$Annotation_cell)
ct_frac_sd <- fread(sd_cellcounts_path)
ct_frac_sd$Annotation_cell_united <- unite_anno(ct_frac_sd$Annotation_cell)

deconv_list <- list(bp = ct_frac_bp, sd = ct_frac_sd)

# transform to long and bind
ct_frac_long_all <- lapply(1:length(deconv_list), function(n){
  ct_frac_long <- melt(setDT(deconv_list[[n]]), id.vars = c(meta_names, 'Annotation_cell_united'), variable.name = "cell_type")
  
  ct_frac_long$deconv_type <- names(deconv_list[n])
  ct_frac_long$Sample_Roi <- paste0(ct_frac_long$Sample, '_', ct_frac_long$Roi)
  
  # move unreliable predictions to 0
  ct_frac_long$value_clean <- ifelse(ct_frac_long$value <= min_frac | is.na(ct_frac_long$value), 0, ct_frac_long$value)
  ct_frac_long$cell_type <- gsub('Endothelial.cells', 'Endothelial cells', ct_frac_long$cell_type)
  ct_frac_long <- ct_frac_long[ct_frac_long$cell_type %in% ct_names, ] # subset to interesting ct
  
  # add fractions per ROI
  ct_frac_long <- group_by(ct_frac_long, Sample_Roi, cell_type) %>%
    mutate(value_roi_mean = sum(value_clean)/n())
  
  return(ct_frac_long)
})

ct_frac_long <- do.call(rbind, ct_frac_long_all)

# per roi
ct_frac_long_roi <- distinct(ct_frac_long, Sample_Roi, NACT_status, Segment_geomx, Annotation_cell, Annotation_cell_united, deconv_type, cell_type, value_roi_mean)

# check fraction distribution across aois and rois ------------------------

for(ct_name in ct_names){
  ct_frac_long_ct <- ct_frac_long[ct_frac_long$cell_type == ct_name, ]
  ct_frac_long_roi_ct <- ct_frac_long_roi[ct_frac_long_roi$cell_type == ct_name, ]
  
  ggplot(data = ct_frac_long_ct) +
    geom_density(aes(value_clean, color = deconv_type)) +
    facet_wrap(~Segment) + 
    ggtitle(ct_name)
  
  ggsave(file.path(outp_plot_dir, paste0('density_fraq_', ct_name, '.png')), width = 2000, height = 1000, unit='px')
  
  ggplot(data = ct_frac_long_roi_ct) +
    geom_density(aes(value_roi_mean, color = deconv_type)) +
    ggtitle(ct_name)
  
  ggsave(file.path(outp_plot_dir, paste0('density_fraq_roi_', ct_name, '.png')), width = 2000, height = 1000, unit='px')
}


# check fractions across current labels -----------------------------------

ct_frac_long_roi$label_T <- ifelse(grepl('CD8|CD4', ct_frac_long_roi$Annotation_cell_united), 'T', '')
ct_frac_long_roi$label_M <- ifelse(grepl('Iba1', ct_frac_long_roi$Annotation_cell_united), 'M', '')
ct_frac_long_roi$label_D <- ifelse(grepl('CD11', ct_frac_long_roi$Annotation_cell_united) & !(grepl('Iba1', ct_frac_long_roi$Annotation_cell_united)), 'D', '')
ct_frac_long_roi$label_TMD <- paste0(ct_frac_long_roi$label_T, ct_frac_long_roi$label_M, ct_frac_long_roi$label_D)

table(ct_frac_long_roi[, c('Annotation_cell_united', 'label_TMD')])

for(ct_name in ct_names){
  ct_frac_long_roi_ct <- ct_frac_long_roi[ct_frac_long_roi$cell_type == ct_name, ]

  ggplot(data = ct_frac_long_roi_ct) +
    geom_density(aes(value_roi_mean, color = Annotation_cell_united)) +
    facet_wrap(~deconv_type) + 
    ggtitle(ct_name)
  
  ggsave(file.path(outp_plot_dir, paste0('anno_density_fraq_roi_', ct_name, '.png')), width = 2000, height = 1000, unit='px')
}

for(ct_name in label_ct_names){
  ct_frac_long_roi_ct <- ct_frac_long_roi[ct_frac_long_roi$cell_type == ct_name, ]
  
  ggplot(data = ct_frac_long_roi_ct) +
    geom_density(aes(value_roi_mean, color = label_TMD)) +
    facet_wrap(~deconv_type) + 
    ggtitle(ct_name)
  
  ggsave(file.path(outp_plot_dir, paste0('anno_tmd_density_fraq_roi_', ct_name, '.png')), width = 2000, height = 1000, unit='px')
}

# for T/M/D separately
ggplot(data = ct_frac_long_roi[ct_frac_long_roi$cell_type == 'Tcells', ]) +
  geom_density(aes(value_roi_mean, color = label_T)) +
  facet_wrap(~deconv_type)

ggsave(file.path(outp_plot_dir, paste0('anno_t_density_fraq_roi.png')), width = 2000, height = 1000, unit='px')

ggplot(data = ct_frac_long_roi[ct_frac_long_roi$cell_type == 'Macrophages', ]) +
  geom_density(aes(value_roi_mean, color = label_M)) +
  facet_wrap(~deconv_type)

ggsave(file.path(outp_plot_dir, paste0('anno_m_density_fraq_roi.png')), width = 2000, height = 1000, unit='px')

ggplot(data = ct_frac_long_roi[ct_frac_long_roi$cell_type == 'DCs', ]) +
  geom_density(aes(value_roi_mean, color = label_D)) +
  facet_wrap(~deconv_type)

ggsave(file.path(outp_plot_dir, paste0('anno_d_density_fraq_roi.png')), width = 2000, height = 1000, unit='px')



# check scatters ----------------------------------------------------------

# bp vs sd per anno
for(ct_name in label_ct_names){
  ct_frac_long_roi_ct <- ct_frac_long_roi[ct_frac_long_roi$cell_type == ct_name, ]
  
  # separate cols for sd and bp vals
  ct_frac_long_roi_ct_deconv <- dcast(setDT(ct_frac_long_roi_ct), 
                                      Sample_Roi+Annotation_cell_united+label_TMD ~ deconv_type, 
                                      value.var = "value_roi_mean")
  
  ggplot(data = ct_frac_long_roi_ct_deconv) +
    geom_point(aes(x = bp, y = sd, color = Annotation_cell_united))
  
  ggsave(file.path(outp_plot_dir, paste0('scatter_label_', ct_name, '.png')), width = 2000, height = 1000, unit='px')
  
  ggplot(data = ct_frac_long_roi_ct_deconv) +
    geom_point(aes(x = bp, y = sd, color = label_TMD))
  
  ggsave(file.path(outp_plot_dir, paste0('scatter_label_tmd_', ct_name, '.png')), width = 2000, height = 1000, unit='px')
}


# make hmap with bp and sd values -----------------------------------------

# separate cols for sd and bp vals
# make clustered hmaps for roi

# per ROI
ct_frac_wide_roi_deconv <- dcast(setDT(ct_frac_long_roi), 
                                    Sample_Roi+Annotation_cell_united+label_TMD+Segment_geomx ~ deconv_type+cell_type, 
                                    value.var = "value_roi_mean")


segm <- 'stroma'
for(segm in unique(ct_frac_wide_roi_deconv$Segment_geomx)){
  
  print(segm)
  
  ct_frac_wide_roi_deconv_segm <- filter(ct_frac_wide_roi_deconv, Segment_geomx == segm)
  
  
  ct_frac_mtx <- as.matrix(column_to_rownames(ct_frac_wide_roi_deconv_segm[, c('Sample_Roi', 'bp_Tcells', 'sd_Tcells', 
                                                                          'bp_Macrophages', 'sd_Macrophages', 'bp_DCs', 'sd_DCs')], 'Sample_Roi'))
  
  ct_frac_mtx_zscore <- scale(ct_frac_mtx) # zscore by column
  
  # cluster by hclust
  ct_frac_zscore_hclust <- hclust(dist(ct_frac_mtx_zscore), method = "average")
  plot(ct_frac_zscore_hclust, hang = -1, cex = 0.4)
  ct_frac_zscore_hclust_cut <- cutree(ct_frac_zscore_hclust, h = 2)
  
  # average zscores per ct
  ct_frac_mtx_zscore_avg <- mutate(as.data.frame(ct_frac_mtx_zscore),
                                   Tcells = (bp_Tcells+sd_Tcells)/2) %>%
    mutate(Macrophages = (bp_Macrophages+sd_Macrophages)/2) %>%
    mutate(DCs = (bp_DCs+sd_DCs)/2) %>%
    select(Tcells, Macrophages, DCs)
  
  # cluster by hclust
  ct_frac_zscore_hclust_avg <- hclust(dist(ct_frac_mtx_zscore_avg), method = "average")
  plot(ct_frac_zscore_hclust_avg, hang = -1, cex = 0.4)
  ct_frac_zscore_hclust_cut_avg <- cutree(ct_frac_zscore_hclust_avg, h = 2)
  
  
  #########
  # make hmaps
  ha = HeatmapAnnotation(
    segment = anno_simple(ct_frac_wide_roi_deconv_segm$Segment_geomx),
    ct_label = anno_simple(ct_frac_wide_roi_deconv_segm$Annotation_cell_united),
    hclust = anno_simple(as.character(unname(ct_frac_zscore_hclust_cut))),
    hclust_avg = anno_simple(as.character(unname(ct_frac_zscore_hclust_cut_avg))),
    which = "row", show_legend = TRUE)
  
  # hmap for sd/bp separately
  png(filename=file.path(outp_plot_dir, paste0('hmap_roi_zscore_', segm, '.png')), 
      width=10, height=6,units="in",res=2000)
  
  ind_heat <- Heatmap(ct_frac_mtx_zscore, cluster_columns = F, cluster_rows= ct_frac_zscore_hclust,
                      show_row_names = TRUE, show_column_names = TRUE,
                      left_annotation = ha, show_heatmap_legend = TRUE)
  
  
  draw(ind_heat, annotation_legend_side = "right", heatmap_legend_side = "right")
  dev.off()
  
  # hmap for averaged sd+bp zscores
  png(filename=file.path(outp_plot_dir, paste0('hmap_roi_zscore_avg_', segm, '.png')), 
      width=10, height=6,units="in",res=2000)
  
  ind_heat <- Heatmap(ct_frac_mtx_zscore_avg, cluster_columns = F, cluster_rows= ct_frac_zscore_hclust_avg,
                      show_row_names = TRUE, show_column_names = TRUE,
                      left_annotation = ha, show_heatmap_legend = TRUE)
  
  
  draw(ind_heat, annotation_legend_side = "right", heatmap_legend_side = "right")
  dev.off()

}


#######################3
# per AOI

ct_frac_wide_deconv <- dcast(setDT(ct_frac_long), 
                             dcc_filename+Segment+Sample_Roi+Annotation_cell_united ~ deconv_type+cell_type, 
                             value.var = "value_clean")

# TODO iterate
#ct_frac_wide_deconv <- filter(ct_frac_wide_deconv, Segment == 'tumor')

# make clustered hmaps for roi
ct_frac_mtx <- as.matrix(column_to_rownames(ct_frac_wide_deconv[, c('dcc_filename', 'bp_Tcells', 'sd_Tcells', 
                                                                        'bp_Macrophages', 'sd_Macrophages', 'bp_DCs', 'sd_DCs')], 'dcc_filename'))

ct_frac_mtx_zscore <- scale(ct_frac_mtx) # by column

# average zscores per ct
ct_frac_mtx_zscore_avg <- mutate(as.data.frame(ct_frac_mtx_zscore),
                                 Tcells = (bp_Tcells+sd_Tcells)/2) %>%
  mutate(Macrophages = (bp_Macrophages+sd_Macrophages)/2) %>%
  mutate(DCs = (bp_DCs+sd_DCs)/2) %>%
  select(Tcells, Macrophages, DCs)

ha = HeatmapAnnotation(
  segment = anno_simple(ct_frac_wide_deconv$Segment),
  ct_label = anno_simple(ct_frac_wide_deconv$Annotation_cell_united),
  which = "row", show_legend = TRUE)


png(filename=file.path(outp_plot_dir, paste0('hmap_aoi_zscore_avg.png')), 
    width=8, height=6,units="in",res=2000)

ind_heat <- Heatmap(ct_frac_mtx_zscore_avg, cluster_columns = F, cluster_rows= T,
                    show_row_names = FALSE, show_column_names = TRUE,
                    left_annotation = ha, show_heatmap_legend = TRUE)


draw(ind_heat)
dev.off()


#####################################################3
# find stromal segment with the highest predicted tumor frac by bp 
# and the biggest difference between bp and sd

# stroma 
ct_frac_wide_deconv_str <- ct_frac_wide_deconv[ct_frac_wide_deconv$Segment == 'stroma', c('dcc_filename', 'Sample_Roi', 'bp_tumor', 'sd_tumor')]
ct_frac_wide_deconv_str$bp_minus_sd_tumor <- ct_frac_wide_deconv_str$bp_tumor - ct_frac_wide_deconv_str$sd_tumor
ct_frac_wide_deconv_str <- arrange(ct_frac_wide_deconv_str, -bp_tumor, -bp_minus_sd_tumor)

fwrite(ct_frac_wide_deconv_str, file.path(outp_plot_dir, paste0('stromal_aoi_with_high_tumor_frac_to_check.csv')))

##########################################
##########################################
# old code from b1
# relabel rois per sample -------------------------------------------------

frac_sample_all <- lapply(unique(frac_mid$Sample), function(sample_name){
  frac_sample <- frac_mid_bp[frac_mid$Sample == sample_name,]
  
  frac_sample_roi <- group_by(frac_sample, Roi) %>%
    dplyr::summarise(macro_sum = sum(Macrophages)/2, cd8_sum = sum(Tcells)/2, nseg = n()) %>%
    dplyr::filter(nseg == 2) # remove rois where one segment was removed due to qc 
  
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

fwrite(frac_sample_all, file.path(output_dir, 'deconvolution', 'sd_mid_lvl_ct_relabeled_roi.csv'))
