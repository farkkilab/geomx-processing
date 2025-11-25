library(data.table)
library(dplyr)
library(ggplot2)
library(plotly)
library(ggpmisc)

ct_frac <- fread('/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2808/deconvolution/bayes_prism/bp_res_mid_lvl_ct_updated_ct_fraction.csv')
ct_frac_kay <- fread('/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-1811/deconvolution/bayes_prism/bp_res_mid_lvl_ct_updated_ct_fraction.csv')

out_dir <- '/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-1811/deconvolution/bp_vah_vs_haut_reference'
dir.create(out_dir)
# cells_immune <- c("Tcells_CD4","Tcells_CD8", "Tcells_Treg","Tcells_other", "NK", "ILC",
#                   "Plasma_cells", "B_cells", "DC",  "Macrophages",  "Mast_cells")

cells_immune <- c("Tcells_CD8","Tcells_CD4", "Tcells_other", "Bcells", 
                    "Macrophages_Monocytes", "DCs", "Mast_cells", "NKcells")

ct_frac$stroma <- ct_frac$Fibroblasts_Mesothelial + ct_frac$Endothelial_cells
ct_frac$immune <- rowSums(as.data.frame(ct_frac)[, cells_immune], na.rm = T)

cell_frac_long <- melt(ct_frac, id.vars = c('dcc_filename', "Roi", "Segment", 
                                                        "Sample", "NACT_status", "Annotation_cell", "Segment_geomx", 
                                                        "Patient","Site"),
                            variable.name = 'cell_type', value.name = 'fraction')

ct_stroma <- filter(ct_frac, Segment == 'stroma') %>%
  arrange(desc(tumor)) %>%
  select(dcc_filename, tumor, Roi, Sample, Segment_geomx)


#############

for(seg in unique(ct_frac$Segment)){
  
}

ct_frac_str <- ct_frac[ct_frac$Segment == 'stroma', ]
ct_frac_tum <- ct_frac[ct_frac$Segment == 'tumor', ]

ggplot(data = ct_frac, aes(x = tumor, y = stroma)) +
  geom_point(aes(color = Segment_geomx)) + 
  ggtitle(paste0('tumor vs stroma ct fractions')) + 
  xlab('tumor ct fraction') +
  ylab('stroma ct fraction') +
  facet_wrap(~Segment, scales = "fixed", dir="v") 

###
# 3d
plot_ly(x=ct_frac$tumor, y=ct_frac$stroma, z=ct_frac$immune, type="scatter3d", mode="markers", color=ct_frac$Segment, size = 0.5)

plot_ly(x=ct_frac_str$tumor, y=ct_frac_str$stroma, z=ct_frac_str$immune, type="scatter3d", mode="markers", color=ct_frac_str$Segment_geomx, size = 0.5)
plot_ly(x=ct_frac_tum$tumor, y=ct_frac_tum$stroma, z=ct_frac_tum$immune, type="scatter3d", mode="markers", color=ct_frac_tum$Segment_geomx, size = 0.5)


#################3
# annotate stroma segments with highest tumor fractions
ct_frac$Segment_hitumor06 <- ifelse(ct_frac$Segment == 'stroma' & ct_frac$tumor >= 0.6, 'stroma_hitumor', ct_frac$Segment)

top20prc <- as.numeric(quantile(ct_frac$tumor[ct_frac$Segment == 'stroma'], probs = c(0.8))) 

ct_frac$Segment_hitumor20perc <- ifelse(ct_frac$Segment == 'stroma' & ct_frac$tumor >= top20prc, 'stroma_hitumor', ct_frac$Segment)


fwrite(ct_frac[, c('dcc_filename', 'Segment_hitumor06', 'Segment_hitumor20perc')], 
       '/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2808/deconvolution/bayes_prism/bp_hitum_in_stroma.csv')
############################################33
# comparison kay vs erd reference scRNAseq bp ct fractions


ct_frac_erd <- ct_frac[ct_frac$dcc_filename %in% ct_frac_kay$dcc_filename, ]
colnames(ct_frac_erd)[2:12] <- paste0(colnames(ct_frac_erd)[2:12], '_vah')
colnames(ct_frac_erd)
colnames(ct_frac_kay)[2:17] <- paste0(colnames(ct_frac_kay)[2:17], '_haut')
ct_frac_kay <- ct_frac_kay[, 1:17]
colnames(ct_frac_kay)
ct_frac_kay$Fibroblasts_Mesothelial_haut <- ct_frac_kay$Fibroblasts_haut + ct_frac_kay$Mesothelial_haut
ct_frac_kay$Bcells_Plasma_haut <- ct_frac_kay$B_cells_haut + ct_frac_kay$Plasma_cells_haut
ct_frac_kay$Tcells_other_reg_haut <- ct_frac_kay$Tcells_other_haut + ct_frac_kay$Tcells_Treg_haut

ct_frac_both <- left_join(ct_frac_erd, ct_frac_kay)


ct_vah <- c("tumor_vah", "Bcells_vah", "Tcells_CD4_vah", "Tcells_other_vah", "Tcells_CD8_vah", 
            "Fibroblasts_Mesothelial_vah", "Macrophages_Monocytes_vah", "Mast_cells_vah", "NKcells_vah", 
            "Endothelial_cells_vah", "DCs_vah")

ct_haut <- c("tumor_haut", "Bcells_Plasma_haut", "Tcells_CD4_haut", "Tcells_other_reg_haut", "Tcells_CD8_haut",
             "Fibroblasts_Mesothelial_haut", "Macrophages_haut", "Mast_cells_haut", "NK_haut",
             "Endothelial_haut", "DC_haut")

#######
# scatterplots 
for(i in 1:length(ct_vah)){
  ct1 <- ct_vah[i]
  ct2 <- ct_haut[i]
  print(ct1)
  print(ct2)
  
  bp_sd_scatter <- ggplot(data = ct_frac_both, aes(x = get(ct1), y = get(ct2))) +
    geom_point(aes(color = Annotation_cell)) +
    ggtitle(paste0(ct1, ' fractions vah vs haut reference'),
            subtitle = paste('overall pearson cor: ', round(stats::cor(ct_frac_both[[ct1]],
                                                                       ct_frac_both[[ct2]], 
                                                                       use = "complete.obs"), 2))) +
    xlab(ct1) +
    ylab(ct2) +
    facet_wrap(~ Segment) +
    geom_smooth(method='lm', formula= y~x) +
    stat_correlation(method = 'pearson')
  
  ggsave(file.path(out_dir, paste0('deconv_comparison_scatter_', ct1, '.png')),
         width = 2000, height = 2000, unit = 'px')
}

############
# boxplots
# TODO put to long, make names the same..
# for(i in 1:length(ct_vah)){
#   
#   ct1 <- ct_vah[i]
#   ct2 <- ct_haut[i]
#   print(ct1)
#   print(ct2)
#   
#   cell_fraq_both_long_ct <- cell_fraq_both_long[cell_fraq_both_long$cell_type == ct, ]
#   
#   ggplot(cell_fraq_both_long_ct, aes(x = deconv_type, y = fraction)) + 
#     geom_boxplot(aes(fill = deconv_type), alpha = .2) +
#     geom_line(aes(group = dcc_filename), size = 0.2, alpha = 0.8) + 
#     geom_point(size = 0.2) + 
#     ggtitle(paste0(ct, ' fractions bp vs sd')) +
#     facet_wrap(~ Segment)
#   
#   ggsave(file.path(output_dir, 'sanity_check','deconv', paste0('deconv_comparison_', ct, '_bp_sd.png')),
#          width = 1500, height = 1000, unit = 'px')
# }