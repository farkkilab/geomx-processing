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
library(tidyverse)
library(ggpubr)
library(ggpmisc)
library(rstatix)
library(ggcorrplot)

# define variables --------------------------------------------------------

data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'
output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2'

input_rds_path <- file.path(output_dir, 'geomx_qc_norm.RDS')

bp_path <- file.path(output_dir, 'deconvolution', 'bp', 'bp_res_mid_lvl_ct_45.RDS')
sd_path <- file.path(output_dir, 'deconvolution', 'sd', 'sd_res_mid_lvl_ct_nofilt.rds')

deconv_type <- 'bp'
res_path <- bp_path

lm_gsva_all <- fread('/media/iganiemi/T7-iga/st/geomx-processing/results/nact2/gsva/lm_gsva_deconv_all_cell_types/all/bp/lm_gsva_all.csv')
lm_gsva_pre <- fread('/media/iganiemi/T7-iga/st/geomx-processing/results/nact2/gsva/lm_gsva_deconv_all_cell_types/all/bp/lm_gsva_pre.csv')
lm_gsva_post <- fread('/media/iganiemi/T7-iga/st/geomx-processing/results/nact2/gsva/lm_gsva_deconv_all_cell_types/all/bp/lm_gsva_post.csv')

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'cell-types'), showWarnings = T, recursive = T)

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')

# load deconvolution results ----------------------------------------------

geomx_obj <- readRDS(input_rds_path)
meta_names <- c('dcc_filename', 'Patient', 'Segment', 'Sample', 'NACT status', 'Annotation_cell', 'PFS')
meta_less_names <- c('dcc_filename', 'Segment', 'Sample', 'NACT_status')
# cell_types <- c("tumor", "Tcells", "Bcells", "Fibroblasts", "NKcells", "Macrophages",
#                 "DCs", "Endothelial cells", "Mast cells")

cell_types <- c("Tcells", "Bcells",  "NKcells", "Macrophages", "DCs", "Mast cells")

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

ct_frac <- rename(ct_frac, NACT_status = `NACT status`)

# ct_frac$Sample <- as.factor(ct_frac$Sample)
# ct_frac$NACT_status <- as.factor(ct_frac$NACT_status)
# ct_frac$Segment <- as.factor(ct_frac$NACT_status)


#ct_frac <- ct_frac[ct_frac$Annotation_cell == 'posCD8_posIBA1',]
#ct_frac_s <- filter(ct_frac, as.character(Segment) == 'stroma')
ct_frac_s <- ct_frac[ct_frac$Segment == 'stroma',]
ct_frac_t <- ct_frac[ct_frac$Segment == 'tumor',]

ct_frac_post <- ct_frac[ct_frac$NACT_status == 'post',]

ct_frac_long <- melt(ct_frac[, c(meta_less_names, cell_types)], id.vars = meta_less_names)
colnames(ct_frac_long) <- c(meta_less_names, 'cell_type', 'fraq')
ct_frac_long$Segment <- factor(ct_frac_long$Segment, levels = c('stroma', 'tumor'))
ct_frac_long$NACT_status <- factor(ct_frac_long$NACT_status, levels = c('pre', 'post'))

# explore cell nrs --------------------------------------------------------

############################
# ct in segment and nact status
ggplot(data = ct_frac_long, aes(x = cell_type, y = fraq, color = Segment)) +
  geom_boxplot() +
  ylab('fraction') +
  xlab('cell type') +
  #ylim(0, 0.3) +
  geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = TRUE, size = 0.2, label.size = 2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text=element_text(size=10)) +
  facet_wrap(~NACT_status)

ggsave(file.path(output_dir, 'cell-types', paste0('fin_box_segm.png')),
       width = 1500, height = 1000, unit = 'px')

ggplot(data = ct_frac_long, aes(x = cell_type, y = fraq, color = NACT_status)) +
  geom_boxplot() +
  ylab('fraction') +
  xlab('cell type') +
  #ylim(0, 0.3) +
  geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = TRUE, size = 0.2, label.size = 2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text=element_text(size=10)) +
  facet_wrap(~Segment)

ggsave(file.path(output_dir, 'cell-types', paste0('fin_box_nact.png')),
       width = 1500, height = 1000, unit = 'px')

##########################################




####################################3
#Macro vs Tcells
ggplot(data = ct_frac, aes(x = Macrophages, y = Tcells, shape = Segment, color = NACT_status)) +
  geom_point(alpha = 0.5, size = 2)

ggplot(data = ct_frac_post, aes(x = Macrophages, y = Tcells, shape = Segment, color = PFS)) +
  geom_point(alpha = 0.5, size = 2)

#########################################
# boxpl NACT status
sapply(cell_types, function(ct){
  ggplot(data = ct_frac_s, aes(x = NACT_status, y = get(ct))) +
    geom_boxplot() +
    geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = TRUE, size = 0.2, label.size = 2.8) +
    ylab(ct)
  
  ggsave(file.path(output_dir, 'cell-types', paste0('box_nact_stroma_', ct, '.png')),
         width = 1500, height = 1000, unit = 'px')
  
  ggplot(data = ct_frac_t, aes(x = NACT_status, y = get(ct))) +
    geom_boxplot() +
    geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = TRUE, size = 0.2, label.size = 2.8) +
    ylab(ct)
  
  ggsave(file.path(output_dir, 'cell-types', paste0('box_nact_tumor_', ct, '.png')),
         width = 1500, height = 1000, unit = 'px')
})

#########################################
# boxpl Segment facet NACT status
sapply(cell_types, function(ct){
  ggplot(data = ct_frac, aes(x = Segment, y = get(ct))) +
    geom_boxplot() +
    geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = TRUE, size = 0.2, label.size = 2.8) +
    ylab(ct) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_wrap(~NACT_status)
  
  ggsave(file.path(output_dir, 'cell-types', paste0('box_segm_', ct, '.png')),
         width = 1500, height = 1000, unit = 'px')

})


################################################
# boxpl ROI type facet NACT status
sapply(cell_types, function(ct){
  ggplot(data = ct_frac_s, aes(x = Annotation_cell, y = get(ct))) +
    geom_boxplot() +
    geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = TRUE, size = 0.2, label.size = 2.8) +
    ylab(ct) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_wrap(~NACT_status)
  
  ggsave(file.path(output_dir, 'cell-types', paste0('box_anno_stroma_', ct, '.png')),
         width = 1500, height = 1000, unit = 'px')
  
  ggplot(data = ct_frac_t, aes(x = Annotation_cell, y = get(ct))) +
    geom_boxplot() +
    geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = TRUE, size = 0.2, label.size = 2.8) +
    ylab(ct)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_wrap(~NACT_status)
  
  ggsave(file.path(output_dir, 'cell-types', paste0('box_anno_tumor_', ct, '.png')),
         width = 1500, height = 1000, unit = 'px')
})

#########################################
# ct vs ct lm

sapply(c('Macrophages', 'Tcells', 'Bcells', 'NKcells', 'DCs'), function(main_ct){
  sapply(cell_types, function(ct){
    ggplot(data = ct_frac_s, aes(x = get(main_ct), y = get(ct), color = NACT_status)) +
      geom_point(alpha = 0.5, size = 2) +
      stat_cor(aes(x = get(main_ct), y =get(ct), color = NACT_status), method="pearson", label.y.npc="center", label.x.npc = "left") +
      stat_poly_line(aes(x = get(main_ct), y =get(ct), color = NACT_status)) +
      stat_poly_eq(aes(x = get(main_ct), y =get(ct), color = NACT_status)) +
      ylab(ct)+
      xlab(main_ct)
    
    ggsave(file.path(output_dir, 'cell-types', paste0('scatter_', main_ct, '_stroma_', ct, '.png')),
           width = 1500, height = 1000, unit = 'px')
    
    # model = lm(DCs ~ Macrophages:Sample, data = ct_frac_t)
    # summary(model)
    # model.predict <- cbind(ct_frac_t, predict(model, interval = 'confidence'))
    
    ggplot(data = ct_frac_t, aes(x = get(main_ct), y = get(ct), color = NACT_status)) +
      geom_point() +
      stat_cor(aes(x = get(main_ct), y =get(ct), color = NACT_status), method="pearson", label.y.npc="center", label.x.npc = "left") +
      stat_poly_line(aes(x = get(main_ct), y =get(ct)), method='lm') +
      stat_poly_eq(aes(x = get(main_ct), y =get(ct), color = NACT_status)) +
      #geom_line(aes(x = get(main_ct), y=fit)) +
      ylab(ct) +
      xlab(main_ct)
    
    ggsave(file.path(output_dir, 'cell-types', paste0('scatter_', main_ct, '_tumor_', ct, '.png')),
           width = 1500, height = 1000, unit = 'px')
  })
})


#########################################
# corrplots

# stroma, tumor, pre, post
# all, doublepos


sapply(c(TRUE, FALSE), function(doublepos){
  sapply(c('stroma_pre', 'stroma_post', 'tumor_pre', 'tumor_post'), function(dt_name){
    sname <- unlist(strsplit(dt_name, '_'))[1]
    nname <- unlist(strsplit(dt_name, '_'))[2]
    
    ct_df <- ct_frac[ct_frac$Segment == sname & ct_frac$NACT_status == nname, ]
    
    if(doublepos){
      ct_df <- ct_df[ct_df$Annotation_cell == 'posCD8_posIBA1',]
      dposname <- '_doublepos'
    } else{
      dposname <- 'all'
    }

    ct_df <- ct_df[, c(cell_types)]
    
    ct_mat <- cor_mat(ct_df)
    ct_mat <- column_to_rownames(ct_mat, 'rowname')
    ct_pmat <- cor_pmat(ct_df)

    ct_mat[ct_pmat >= 0.05] <- 0 #rmv insignificant correlations
    
    col_fun <- colorRamp2(c(-1, 0, 1), hcl_palette = 'RdBu', reverse = TRUE)
    
    ct_heat <- Heatmap(as.matrix(ct_mat),
                              height = unit(6, "cm") , width = unit(6, "cm"),border="white",
                              rect_gp = gpar(col = "white", lwd = 2), name=paste(sname, nname),
                              cluster_columns = F, cluster_rows= F,
                              column_title_gp = gpar(fontsize = 14),
                              show_heatmap_legend = T, col = col_fun,
                              column_names_gp = gpar(fontsize = 12),
                              row_names_gp = gpar(fontsize = 12),
                              column_names_rot = 45,
                              column_title = paste(sname, nname),
                              heatmap_legend_param = list(title = 'correlation'),
                              na_col = 'white', cell_fun = function(j, i, x, y, width, height, fill) {
                                grid.text(sprintf("%.1f", ct_mat[i, j]), x, y, gp = gpar(fontsize = 5))
                              }) # cor nr in the tiles
    
    #plot(ct_heat)
    
    png(file=file.path(output_dir, 'cell-types', paste0('fin_corrplot_', dt_name, dposname, '.png')),
        width = 1000, height = 1000, units = 'px')
  
    draw(ct_heat, heatmap_legend_side = "left")
    dev.off()
    # 
    # corrplot <- ggcorrplot(ct_mat, hc.order = FALSE, outline.color = "white", p.mat = ct_pmat,
    #                        title = dt_name, tl.cex = 10, pch.cex = 10, lab = T)
    # 
    # pdf(file=file.path(output_dir, 'cell-types', paste0('corrplot_', dt_name, dposname, '.pdf')),
    #     width=5, height=6) 
    # 
    # plot(corrplot)
    # dev.off()
  })
})

################################################################
# corrplots together

dt_name <- 'stroma_pre'

ct_hmap_list <- lapply(c('stroma_pre', 'stroma_post', 'tumor_pre', 'tumor_post'), function(dt_name){
    doublepos <- F
    sname <- unlist(strsplit(dt_name, '_'))[1]
    nname <- unlist(strsplit(dt_name, '_'))[2]
    
    ct_df <- ct_frac[ct_frac$Segment == sname & ct_frac$NACT_status == nname, ]
    
    if(doublepos){
      ct_df <- ct_df[ct_df$Annotation_cell == 'posCD8_posIBA1',]
      dposname <- '_doublepos'
    } else{
      dposname <- 'all'
    }
    
    ct_df <- ct_df[, c(cell_types)]
    
    ct_mat <- cor_mat(ct_df)
    ct_mat <- column_to_rownames(ct_mat, 'rowname')
    ct_pmat <- cor_pmat(ct_df)
    ct_pmat <- column_to_rownames(ct_pmat, 'rowname')
    
    ct_mat[ct_pmat >= 0.05] <- 0 #rmv insignificant correlations
    
    col_fun <- colorRamp2(c(-1, 0, 1), hcl_palette = 'RdBu', reverse = TRUE)
    
    ct_heat <- Heatmap(as.matrix(ct_mat),
                       height = unit(6, "cm") , width = unit(6, "cm"),border="white",
                       rect_gp = gpar(col = "white", lwd = 2), name=paste(sname, nname),
                       cluster_columns = F, cluster_rows= F,
                       column_title_gp = gpar(fontsize = 14),
                       show_heatmap_legend = T, col = col_fun,
                       column_names_gp = gpar(fontsize = 12),
                       row_names_gp = gpar(fontsize = 12),
                       column_names_rot = 45,
                       column_title = paste(sname, nname),
                       heatmap_legend_param = list(title = 'correlation'),
                       na_col = 'white', cell_fun = function(j, i, x, y, width, height, fill) {
                         grid.text(sprintf("%.1f", ct_mat[i, j]), x, y, gp = gpar(fontsize = 5))
                       }) # cor nr in the tiles
    
    
    # png(file=file.path(output_dir, 'cell-types', paste0('fin_corrplot_', dt_name, dposname, '.png')),
    #     width = 1000, height = 1000, units = 'px')
    # 
    # draw(ct_heat, heatmap_legend_side = "left")
    # dev.off()
    return(ct_heat)
  })

#seg_var_heatmap_list <- unlist(seg_var_heatmap_list, recursive = F)

# make HeatmapList object from all heatmaps in a list
all_hmaps = NULL

for(i in seq_along(ct_hmap_list)){
  all_hmaps = all_hmaps + ct_hmap_list[[i]]
}

all_hmap_list <- ct_hmap_list[[1]] %h% ct_hmap_list[[2]] %v% ct_hmap_list[[3]] %h% ct_hmap_list[[4]]

#TODO plot in a grid (?)
# adjust length and height depending on nr of plots
pdf(file=file.path(output_dir, 'cell-types', paste0('fin_corrplot_ct.pdf')),
    width=14,
    height=7)

draw(all_hmaps, ht_gap = unit(1, "cm"), row_km = 2,
     column_title_gp = gpar(fontsize = 15), heatmap_legend_side = "left")

dev.off()

# lm with gsva pathways

rsq_thr <- 0.3

lm_all_filt <- filter(lm_gsva_all, lm_rsq >= rsq_thr, cor_pe_pval <= 0.05)
lm_pre_filt <- filter(lm_gsva_pre, lm_rsq >= rsq_thr, cor_pe_pval <= 0.05)
lm_post_filt <- filter(lm_gsva_post, lm_rsq >= rsq_thr, cor_pe_pval <= 0.05)


