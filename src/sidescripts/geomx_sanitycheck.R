library(plyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(GeomxTools)
library(biomaRt)
library(ggpubr)
library(viridis)
library(ggpmisc)
library(reshape2)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

# PROBLEM: BP predicts a lot more tumor in stromal AOIs, 
# while sd predicts a lot more Tcells 
# maybe binning all Tcells together created highly variable cluster which confuses the algorithm

# careful for genes with variance 0 removed from deconvolution

# TODO pfs and os with quartiles for every batch in clinical table

# TODO histogram x=proportion of genes non-o y = sample count for raw + norm + batch corrected data + deconv 

# SD and SD_BG are identical
###########################

proj_dir <<- '~/Documents/phd/st'

#output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch2-1903') # batch2
#output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch1-1903') # batch1
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch12-1004') # batch12
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch12-1205-no-counts-shift2') # batch12

#clinical_dt <<- file.path(proj_dir, 'data/geomx/9_eyemt_patient_clinical_data.csv')

geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS')

ct_names <- c('Bcells', 'DCs', 'Endothelial cells', 'Fibroblasts', 'Macrophages', 'NKcells', 'Tcells', 'tumor')
cells_immune <- c('Bcells', 'DCs', 'Macrophages', 'NKcells', 'Tcells')

ct_markers_path <- file.path(proj_dir, 'geomx-processing', 'data', 'signatures', 'ct_markers.csv')

ct_gsea_all_path <- file.path(output_dir, 'pathway_analysis', 'gsea', 'ssgsea_norm_harmony_batch_corr_all_custom_ct_markers.csv.csv')


bp_cellcounts_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', 'bp_res_mid_lvl_ct_ct_fraction.csv')
sd_cellcounts_path <- file.path(output_dir, 'deconvolution', 'spatial_decon', 'sd_res_bg_mid_lvl_ct_geomxfilt_ct_fraction.csv')

deconv_raw_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', 'bp_res_mid_lvl_ct.RDS')
deconv_harmony_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', 'bp_res_mid_lvl_ct_expr_mtx_cleaned_vst_harmony_batch_corr.RDS')
deconv_limma_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', 'bp_res_mid_lvl_ct_expr_mtx_cleaned_vst_limma_batch_corr_main_batch_nrbatch_nr_cov_no.RDS')
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
deconv_harmony <- readRDS(deconv_harmony_path)

expr_list <- list(raw = geomx_obj@assayData$exprs, deseq_norm = geomx_obj@assayData$deseq2_norm,
                  q3_norm = geomx_obj@assayData$q3_norm, vst = geomx_obj@assayData$deseq2_vst, 
                  limma_batch_corr = geomx_obj@assayData$limma_batch_corr,
                  harmony_batch_corr = geomx_obj@assayData$harmony_batch_corr)

# proportion of 0 reads ---------------------------------------------------

sapply(names(expr_list), function(x){
  expr <- expr_list[[x]]
  
  # limma and harmony batch effect corrected data are in a log scale so 0 --> 1
  #null_nr <- ifelse(x %in% c('limma_batch_corr', 'harmony_batch_corr'), 1, 0)
  
  prop_0 <- apply(expr, 2, function(x){
    length(which(x == 0))/length(x)
  })
  
  prop_0 <- cbind(sData(geomx_obj)[, c('dcc_filename', 'Segment', 'main_batch_nr')], prop_0)
  
  ggplot(data = prop_0) +
    geom_histogram(aes(x= prop_0, fill = Segment), alpha = 0.5, bins = 200) +
    #facet_wrap(~Segment, scales = "fixed", dir="v") +
    ylab('AOI counts') +
    ggtitle(paste('proportions of 0 reads in ', x, ' expr mtx'), 
            subtitle = paste0('prop of 0 in whole expr mtx: ', round(length(which(expr == 0))/length(expr), 5)))
  
  ggsave(file.path(output_dir, 'sanity_check', paste0('prop_of_0_', x, '.png')), 
         width = 2000, height = 1000, unit = 'px')
})

# nr of reads per nuclei --------------------------------------------------

reads_sum <- lapply(expr_list, colSums)
reads_sum <- as.data.frame(do.call(cbind, reads_sum))
reads_sum <- cbind(sData(geomx_obj)[, c('dcc_filename', 'nuclei', 'Segment', 'main_batch_nr')], reads_sum)
reads_sum$raw <- round(reads_sum$raw) # reverts weird behaviour of colsums to return floats for sum of integers

for(coln in names(expr_list)){
  ggplot(data = reads_sum) +
    geom_point(aes(x= nuclei, y = get(coln), shape = Segment, color = as.factor(main_batch_nr))) +
    ylab(coln) +
    ggtitle(paste(coln, 'reads nr per nuclei nr'))
  
  ggsave(file.path(output_dir, 'sanity_check', paste0('reads_nr_per_nuclei_', coln, '.png')), width = 2000, height = 1000, unit = 'px')
  
}


# gene coverage -----------------------------------------------------------

# what is the % of samples covered by each gene?
loqdt <- sData(geomx_obj)[, c('LOQ')]

# for all data
fdt <- fData(geomx_obj)

# plot gene detection rate
ggplot(data = fdt) +
  geom_histogram(aes(x = fdt$DetectionRate), bins = 100) +
  #geom_density(aes(x = fdt$DetectionRate)) +
  xlim(0, 1) +
  ggtitle('gene detection rate')

ggsave(file.path(output_dir, 'sanity_check', paste0('gene_detection_rate_all.png')))


######################
######################
# for deconvolution
# this doesn't make too much sense
# bprism_res <- readRDS(deconv_raw_path)
# 
# gene_dt_ct_all <- lapply(ct_names, function(ct_name){
#   print(ct_name)
#   
#   # mask unreliable bp results
#   cell_frac_cv <- as.data.frame(bprism_res@posterior.theta_f@theta.cv)
#   cell_to_rm <- rownames(cell_frac_cv)[cell_frac_cv[[ct_name]] > 0.2]
#   
#   deconv_ct <- BayesPrism::get.exp(bp=bprism_res,
#                                    state.or.type="type",
#                                    cell.name=ct_name)
#   
#   deconv_ct_cleaned <- t(deconv_ct[!(rownames(deconv_ct) %in% cell_to_rm), ])
#   
#   print(dim(deconv_ct_cleaned))
#   
#   # genes detected > 0
#   gene_dt_0 <- as.data.frame(rowSums(replace(deconv_ct_cleaned, deconv_ct_cleaned != 0, 1)))
#   colnames(gene_dt_0) <- 'freq'
#   gene_dt_0$cov <- gene_dt_0$freq / ncol(deconv_ct_cleaned)
#   colnames(gene_dt_0) <- paste0(ct_name, '_0_', colnames(gene_dt_0))
#   
#   ggplot(data = gene_dt_0) +
#     geom_histogram(aes(x = get(paste0(ct_name, '_0_cov'))), bins = 100) +
#     xlim(0, 1.1) +
#     ylim(0, 1500) +
#     xlab('gene detection rate') +
#     ggtitle(paste0('gene detection rate (above 0) for ', ct_name))
#   
#   ggsave(file.path(output_dir, 'sanity_check', paste0('gene_detection_rate_deconv_', ct_name, '_above0.png')))
#   
#   # genes detected > LOQ
#   deconv_ct_aboveloq <- lapply(colnames(deconv_ct_cleaned), function(dcc){
#     loq <- loqdt[dcc, ]
#     dcc_above_loq <- replace(deconv_ct_cleaned[, dcc], deconv_ct_cleaned[, dcc] < loq, 0)
#     
#     return(dcc_above_loq)
#   })
#   
#   deconv_ct_aboveloq <- do.call(cbind, deconv_ct_aboveloq)
#   colnames(deconv_ct_aboveloq) <- colnames(deconv_ct_cleaned)
#   
#   gene_dt_aboveloq <- as.data.frame(rowSums(replace(deconv_ct_aboveloq, deconv_ct_aboveloq != 0, 1)))
#   colnames(gene_dt_aboveloq) <- 'freq'
#   gene_dt_aboveloq$cov <- gene_dt_aboveloq$freq / ncol(deconv_ct_aboveloq)
#   colnames(gene_dt_aboveloq) <- paste0(ct_name, '_loq_', colnames(gene_dt_aboveloq))
#   
#   ggplot(data = gene_dt_aboveloq) +
#     geom_histogram(aes(x = get(paste0(ct_name, '_loq_cov'))), bins = 100) +
#     xlim(0, 1.1) +
#     ylim(0, 1500) +
#     xlab('gene detection rate') +
#     ggtitle(paste0('gene detection rate (above LOQ for AOI) for ', ct_name))
#   
#   ggsave(file.path(output_dir, 'sanity_check', paste0('gene_detection_rate_deconv_', ct_name, '_aboveLOQ.png')))
#   
#   stopifnot(identical(rownames(gene_dt_0), rownames(gene_dt_aboveloq)))
#   
#   gene_dt_ct <- cbind(gene_dt_0, gene_dt_aboveloq)
#   
#   return(gene_dt_ct)
# })
# 
# gene_dt_ct_all <- do.call(cbind, gene_dt_ct_all)
# fwrite(gene_dt_ct_all, file.path(output_dir, 'sanity_check', 'gene_coverage_deconvolution.csv'))


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
  cell_fraq_res$stroma <- cell_fraq_res$Fibroblasts + cell_fraq_res$`Endothelial cells`
  
  ct_gsea_all_fraq <- left_join(ct_gsea_all, cell_fraq_res[, c('dcc_filename', c(ct_names, 'stroma'))])
  
  for(ct_name in c(ct_names[!ct_names == 'DCs'], 'stroma')){
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
cell_fraq_bp$other <- cell_fraq_bp$`Mast cells` + cell_fraq_bp$other
cell_fraq_bp <- cell_fraq_bp[, c('dcc_filename','Segment', 'Annotation_cell','other', ct_names)]
colnames(cell_fraq_bp) <- c('dcc_filename', 'Segment', 'Annotation_cell','other_bp', paste0(ct_names, '_bp'))       

cell_fraq_sd <- as.data.frame(cell_fraq$sd)
cell_fraq_sd$other <- cell_fraq_sd$`Mast cells` + cell_fraq_bp$other
cell_fraq_sd <- cell_fraq_sd[, c('dcc_filename','other', ct_names)]
colnames(cell_fraq_sd) <- c('dcc_filename','other_sd', paste0(ct_names, '_sd'))   

cell_fraq_both <- left_join(cell_fraq_bp, cell_fraq_sd)

cell_fraq_both$stroma_bp <- cell_fraq_both$Fibroblasts_bp + cell_fraq_both$`Endothelial cells_bp`
cell_fraq_both$stroma_sd <- cell_fraq_both$Fibroblasts_sd + cell_fraq_both$`Endothelial cells_sd`

cell_fraq_both_long <- melt(cell_fraq_both, id.vars = c('dcc_filename', 'Segment', 'Annotation_cell'),
                            variable.name = 'cell_type', value.name = 'fraction')

cell_fraq_both_long$deconv_type <- ifelse(grepl('bp', cell_fraq_both_long$cell_type), 'bp', 'sd')
cell_fraq_both_long$cell_type <- gsub('_bp|_sd', '', cell_fraq_both_long$cell_type)


# make a paired plot
for(ct in c(ct_names, 'stroma')){
  
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

#########################
# make scatters bp vs sd
for(ct_name in c(ct_names, 'stroma')){
  bp_sd_scatter <- ggplot(data = cell_fraq_both, aes(x = get(paste0(ct_name, '_bp')),
                                                            y = get(paste0(ct_name, '_sd')))) +
    geom_point(aes(color = Annotation_cell)) + 
    ggtitle(paste0(ct_name, ' fractions bp vs sd'),
            subtitle = paste('overall pearson cor: ', round(stats::cor(cell_fraq_both[[paste0(ct_name, '_bp')]],
                                                        cell_fraq_both[[paste0(ct_name, '_sd')]], 
                                                        use = "complete.obs"), 2))) + 
    xlab('bp fraq') +
    ylab('sd fraq') +
    facet_wrap(~ Segment) +
    geom_smooth(method='lm', formula= y~x) +
    # stat_poly_eq(use_label(c("R2", "p"))) +
    stat_correlation(method = 'pearson')
  
  ggsave(file.path(output_dir,'sanity_check', paste0('deconv_comparison_scatter_', ct_name, '.png')),
         width = 2000, height = 2000, unit = 'px')
}

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


# heatmaps and stacked barplots -------------------------------------------

# TODO not really visible - after making better lables - rerun for each Annotation_cell in facet
sapply(c('bp', 'sd'), function(deconv_type){
  sapply(c('stroma', 'tumor'), function(seg){
    
    cell_fraq_both_long_deconv_seg <- cell_fraq_both_long[cell_fraq_both_long$Segment == seg &
                                                            cell_fraq_both_long$deconv_type == deconv_type, ]

    cell_fraq_both_long_deconv_seg <- arrange(cell_fraq_both_long_deconv_seg, Annotation_cell) %>%
      rowid_to_column()
    cell_fraq_both_long_deconv_seg <- cell_fraq_both_long_deconv_seg[!cell_fraq_both_long_deconv_seg$cell_type == 'stroma', ]
    
    anno_color_mapping <- setNames(viridis(length(unique(cell_fraq_both_long_deconv_seg$Annotation_cell))),
                                   unique(cell_fraq_both_long_deconv_seg$Annotation_cell))
    
    ggplot(data = cell_fraq_both_long_deconv_seg, aes(x = reorder(dcc_filename, rowid), 
                                                      y = fraction, fill = cell_type)) +
      geom_bar(position="fill", stat="identity") +
      ggtitle(seg) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 3,
                                       color = anno_color_mapping[cell_fraq_both_long_deconv_seg$Annotation_cell])) +
      scale_x_discrete(labels = cell_fraq_both_long_deconv_seg$Annotation_cell) 
    
    ggsave(file.path(output_dir, 'sanity_check', paste0('barplot_',deconv_type, '_', seg, '.png')),
           width = 1500, height = 1000, unit = 'px')

    cell_fraq_both_long_deconv_seg_imm <- filter(cell_fraq_both_long_deconv_seg, cell_type %in% cells_immune)

    ggplot(data = cell_fraq_both_long_deconv_seg_imm, aes(x = reorder(dcc_filename, rowid), y = fraction, fill = cell_type)) +
      geom_bar(position="fill", stat="identity") +
      ggtitle(seg) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 3,
                                       color = anno_color_mapping[cell_fraq_both_long_deconv_seg$Annotation_cell])) +
      scale_x_discrete(labels = cell_fraq_both_long_deconv_seg$Annotation_cell)

    ggsave(file.path(output_dir, 'sanity_check', paste0('barplot_',deconv_type, '_', seg, '_immune.png')),
           width = 1500, height = 1000, unit = 'px')
  })
})

#########################33
# heatmaps

col_fun = viridis(100)


sapply(c('bp', 'sd'), function(deconv_type){
  
  png(filename = file.path(output_dir, 'sanity_check', paste0('hmap_',deconv_type, '_all.png')), width=1000, height=750)
  Heatmap(t(as.matrix(cell_fraq_both[, paste0(c(ct_names, 'other'), '_', deconv_type)])), col = col_fun) %v%
    HeatmapAnnotation(segment = cell_fraq_both$Segment, 
                      col = list(segment = c("stroma" = "green", "tumor" = "blue")))
  dev.off()
  
  sapply(c('stroma', 'tumor'), function(seg){
    ct_frac_seg <- cell_fraq_both[cell_fraq_both$Segment == seg, ]

    png(filename = file.path(output_dir, 'sanity_check', paste0('hmap_',deconv_type, '_', seg, '.png')), width=1000, height=750)
    print(Heatmap(t(as.matrix(ct_frac_seg[, paste0(c(ct_names, 'other'), '_', deconv_type)])), col = col_fun) %v%
            HeatmapAnnotation(roi_type = ct_frac_seg$Annotation_cell))
    dev.off()

    png(filename = file.path(output_dir, 'sanity_check', paste0('hmap_',deconv_type, '_', seg, '_immune.png')), width=1000, height=750)
    print(Heatmap(t(as.matrix(ct_frac_seg[, paste0(cells_immune, '_', deconv_type)])), col = col_fun) %v%
            HeatmapAnnotation(roi_type = ct_frac_seg$Annotation_cell))
    dev.off()
  })
})

######################333
# for sisana
demo <- fread('/home/iganiemi/Documents/phd/st/gene-regulatory-networks/sisana/sisana/example_input/BRCA_TCGA_20_LumA_LumB_samps_5000_genes_exp.tsv')

geomx_harm <- data.frame(geomx_obj@assayData$harmony_batch_corr)
geomx_harm <- rownames_to_column(geomx_harm, var = 'Target')
write_tsv(geomx_harm, '/home/iganiemi/Documents/phd/st/gene-regulatory-networks/sisana/sisana/geomx_input/geomx_batch12_harmony_corr_expr_mtx2.tsv',
          col_names = T)

meta_segment <- sData(geomx_obj)[, c('dcc_filename', 'Segment')]
meta_segment$dcc_filename <- gsub('\\-', '\\.', meta_segment$dcc_filename)
rownames(meta_segment) <- NULL
colnames(meta_segment) <- NULL

write_csv(meta_segment, '/home/iganiemi/Documents/phd/st/gene-regulatory-networks/sisana/sisana/geomx_input/geomx_meta_segment.csv', 
          col_names = F)
