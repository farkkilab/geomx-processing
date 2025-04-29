library(dplyr)
library(data.table)
library(ggplot2)
library(GeomxTools)
library(biomaRt)
library(ggpubr)
library(viridis)
library(ggpmisc)
library(reshape2)

# TODO pfs and os with quartiles for every batch 
# heatmaps per annotation
# stacked barplots per annotation

# SD and SD_BG are identical
###########################

proj_dir <<- '~/Documents/phd/st'

#output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch2-1903') # batch2
#output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch1-1903') # batch1
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch12-1004') # batch12

#clinical_dt <<- file.path(proj_dir, 'data/geomx/9_eyemt_patient_clinical_data.csv')

geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS')

ct_names <- c('Bcells', 'DCs', 'Endothelial cells', 'Fibroblasts', 'Macrophages', 'NKcells', 'Tcells', 'tumor')
ct_markers_path <- file.path(proj_dir, 'geomx-processing', 'data', 'signatures', 'ct_markers.csv')

ct_gsea_all_path <- file.path(output_dir, 'pathway_analysis', 'gsea', 'ssgsea_norm_harmony_batch_corr_all_custom_ct_markers.csv.csv')


bp_cellcounts_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', 'bp_res_mid_lvl_ct_ct_fraction.csv')
sd_bg_cellcounts_path <- file.path(output_dir, 'deconvolution', 'spatial_decon', 'sd_res_bg_mid_lvl_ct_geomxfilt_ct_fraction.csv')

# set up metadata variables names -----------------------------------------

aoi_id <<- 'dcc_filename'
roi_id <<- 'Roi'

main_batch_var <- 'main_batch_nr'
batch_var <<- 'batch_nr'


aoi_segment_var <<- "Segment"
main_roi_label <<- "Annotation_cell" 
main_experimental_condition <<- 'NACT_status'
sample_name <<- 'Sample'

other_vars_bio <<- c("Segment_geomx", "Patient", "Site") # 'PFS_months', 'PFS'
other_vars_tech <<- c('Slide_Name', "batch_nr_sample_collection")

# load util functions and create dirs -------------------------------------

source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))

dir.create(file.path(output_dir, 'sanity_check'), recursive = T, showWarnings = F)

#################################################################
################################################################
# cell markers: http://117.50.127.228/CellMarker/CellMarkerBrowse.jsp
# human, ovary


geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)

sort(colnames(sData(geomx_obj)))
dim(geomx_obj@assayData$exprs)

deconv_limma <- readRDS(deconv_limma_path)
deconv_harmony <- readRDS(deconv_harm_path)

# nr of reads per nuclei --------------------------------------------------

raw_sum <- as.data.frame(colSums(geomx_obj@assayData$exprs))
norm_sum <- as.data.frame(colSums(geomx_obj@assayData$deseq2_norm))
limma_sum <- as.data.frame(colSums(geomx_obj@assayData$limma_batch_corr))
harm_sum <- as.data.frame(colSums(geomx_obj@assayData$harmony_batch_corr))

nuclei_sum <- cbind(sData(geomx_obj)[, c('dcc_filename', 'nuclei', 'Segment', 'main_batch_nr')], raw_sum, norm_sum, limma_sum, harm_sum)
colnames(nuclei_sum) <- c('dcc_filename', 'nuclei', 'Segment','main_batch_nr', 'raw', 'norm', 'limma', 'harmony')

for(coln in c('raw', 'norm', 'limma', 'harmony')){
  ggplot(data = nuclei_sum) +
    geom_point(aes(x= nuclei, y = get(coln), shape = Segment, color = as.factor(main_batch_nr))) +
    ylab(coln) +
    ggtitle(paste(coln, 'reads nr per nuclei nr'))
  
  ggsave(file.path(output_dir, 'sanity_check', paste0(coln, '_reads_nr_per_nuclei.png')), width = 2000, height = 1000, unit = 'px')
  
}

# stromal/tumor markers ---------------------------------------------------

ct_markers <- fread(ct_markers_path)
all_markers <-   as.list(ct_markers)
all_markers <- lapply(all_markers, function(l){l[l !=""]})

# not in exprs:
# CCL4 IGLC2 IGLC3 KLRC1 TPSB2 XCL2
# SMA ? ACTA2

ct_gsea_all <- fread(ct_gsea_all_path)
ct_gsea_all$pathway <- ifelse(ct_gsea_all$pathway == 'Nkcells', 'NKcells', ct_gsea_all$pathway)
ct_gsea_all$pathway <- ifelse(ct_gsea_all$pathway == 'CD8_Tcells', 'Tcells_CD8', ct_gsea_all$pathway)


# tum/stromal markers in tum/stromal AOIs
ct_boxpl <- pathway_boxplot(ct_gsea_all, 'pathway', 'ssgsea_score', 'Segment', NULL, 
                            'cell type markers ssgsea score per segment',
                            file.path(output_dir, 'sanity_check', 'ct_markers_per_segment.pdf'))

# macro + Tcell paths activity should be higher in ++ (or single+)
ct_gsea_imm <- filter(ct_gsea_all, pathway %in% c('Tcells', 'Macrophages', 'CD8_Tcells'))

ct_boxpl_anno <- pathway_boxplot(ct_gsea_imm,'pathway', 'ssgsea_score', 'Annotation_cell', facet_var = 'Segment', 
                            'cell type markers ssgsea score per annotation',
                            file.path(output_dir, 'sanity_check', 'ct_markers_per_annotation.pdf'),
                            manual_colours = viridis(14))

# DCs were removed bcs only 3 genes - find sth else to check
# T-cells work as expected (especially batch1, both new labels as well)
# macrophages - not big changes in the markers pathway activity


# deconvolution sanity check ----------------------------------------------

###########################
# ct markers activity should be higher in given ct
ct_gsea_deconv <- list.files(file.path(output_dir, 'pathway_analysis', 'gsea'), pattern = 'ct_markers.*csv', full.names = T)
ct_gsea_deconv <- ct_gsea_deconv[-1]

ct_gsea_deconv <- lapply(ct_gsea_deconv, fread)
ct_gsea_deconv <- do.call(rbind, ct_gsea_deconv)
ct_gsea_deconv$expr_signal <- gsub('deconv_', '', ct_gsea_deconv$expr_signal)

ct_boxpl_deconv <- pathway_boxplot(ct_gsea_deconv,'pathway', 'ssgsea_score', 'expr_signal', facet_var = 'Segment', 
                                 'cell type markers ssgsea score in each deconvoluted data',
                                 file.path(output_dir, 'sanity_check', 'ct_markers_in_deconv.pdf'),
                                 manual_colours = viridis(14))

########################
# ct markers activity in full signal correlated with nr of cells sd-bp
cell_fraq <- list(bp = fread(bp_cellcounts_path), sd = fread(sd_cellcounts_path))
colnames(cell_fraq$sd) <- gsub('.', ' ', colnames(cell_fraq$sd), fixed = T)


sapply(1:length(cell_fraq), function(x){
  print(x)
  cell_fraq_res <- as.data.frame(cell_fraq[[x]])
  ct_gsea_all_fraq <- left_join(ct_gsea_all, cell_fraq_res[, c('dcc_filename', ct_names)])
  
  for(ct_name in ct_names[!ct_names == 'DCs']){
    print(ct_name)
    ct_gsea_all_fraq_ct <- ct_gsea_all_fraq[grepl(ct_name, ct_gsea_all_fraq$pathway), ]
    
    ct_scatter <- ggplot(data = ct_gsea_all_fraq_ct, aes(x = ssgsea_score, y = get(ct_name), color = Segment)) +
      geom_point(aes(shape = Segment)) + 
      ggtitle(paste0(ct_name, ' ssgsea marker activity vs cell count')) + 
      xlab(paste0(ct_name, ' markers ssgsea score')) +
      ylab(paste0(ct_name, ' cell count')) +
      geom_smooth(method='lm', formula= y~x) +
      stat_poly_eq(use_label(c("R2", "p"))) +
      facet_wrap(~pathway, scales = "fixed", dir="v") 
    
    ggsave(file.path(output_dir,'sanity_check', paste0(names(cell_fraq)[x], '_', ct_name, '_ssgsea_vs_cell_count.png')),
           width = 2000, height = 2000, unit = 'px')
  }
})

####################################
# pairwise bp-sd comparison

cell_fraq_bp <- as.data.frame(cell_fraq$bp)
cell_fraq_bp <- cell_fraq_bp[, c('dcc_filename','Segment', 'Annotation_cell', ct_names)]
colnames(cell_fraq_bp) <- c('dcc_filename', 'Segment', 'Annotation_cell', paste0(ct_names, '_bp'))       

cell_fraq_sd <- as.data.frame(cell_fraq$sd)
cell_fraq_sd <- cell_fraq_sd[, c('dcc_filename', ct_names)]
colnames(cell_fraq_sd) <- c('dcc_filename', paste0(ct_names, '_sd'))   

cell_fraq_both <- left_join(cell_fraq_bp, cell_fraq_sd)
cell_fraq_both_long <- melt(cell_fraq_both, id.vars = c('dcc_filename', 'Segment', 'Annotation_cell'),
                            variable.name = 'cell_type', value.name = 'fraction')

cell_fraq_both_long$deconv_type <- ifelse(grepl('bp', cell_fraq_both_long$cell_type), 'bp', 'sd')
cell_fraq_both_long$cell_type <- gsub('_bp|_sd', '', cell_fraq_both_long$cell_type)

# make a paired plot
for(ct in ct_names){
  
  cell_fraq_both_long_ct <- cell_fraq_both_long[cell_fraq_both_long$cell_type == ct, ]
  
  ggplot(cell_fraq_both_long_ct, aes(x = deconv_type, y = fraction)) + 
    geom_boxplot(aes(fill = deconv_type), alpha = .2) +
    geom_line(aes(group = dcc_filename), size = 0.2, alpha = 0.8) + 
    geom_point(size = 0.2) + 
    ggtitle(paste0(ct, ' fractions bp vs sd')) +
    facet_wrap(~ Segment)
  
  ggsave(file.path(output_dir, 'sanity_check', paste0('deconv_comparison_', ct, '_bp_sd.png')),
         width = 1500, height = 1000, unit = 'px')
}

# for all
# ggplot(cell_fraq_both_long, aes(x = deconv_type, y = fraction)) + 
#   geom_boxplot(aes(fill = deconv_type), alpha = .2) +
#   geom_line(aes(group = dcc_filename), size = 0.2, alpha = 0.8) + 
#   geom_point(size = 0.2) + 
#   ggtitle('ct fractions bp vs sd') +
#   facet_wrap(~ cell_type)
# 
# ggsave(file.path(output_dir, 'sanity_check', paste0('deconv_comparison_all_bp_sd.png')),
#        width = 1500, height = 1000, unit = 'px')

#####################
# compare macro vs Tcells - scatter per annotation cell (same per sample?)
for(deconv_type in c('bp', 'sd')){
  macro_tcells_scatter <- ggplot(data = cell_fraq_both, aes(x = get(paste0('Macrophages_', deconv_type)),
                                                            y = get(paste0('Tcells_', deconv_type)), 
                                                            color = Annotation_cell)) +
    geom_point() + 
    ggtitle(paste0(deconv_type, 'Macrophages vs Tcells ct fractions')) + 
    xlab('Macrophages fraq') +
    ylab('Tcells fraq') +
    facet_wrap(~ Segment)
    # geom_smooth(method='lm', formula= y~x) +
    # stat_poly_eq(use_label(c("R2", "p")))
  
  ggsave(file.path(output_dir,'sanity_check', paste0('macro_vs_tcell_', deconv_type, '.png')),
         width = 2000, height = 2000, unit = 'px')
}

