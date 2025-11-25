# questions

# 1. Are CXCR6+ T cells (CD8 or CD4) found in proximity to CXCL16+ myeloid cells (macrophages or DCs)? 
#   
# 2. Are there differences before chemotherapy and after chemotherapy?
#   
# 3. Are patients that respond better to treatment enriched in CXCL16+ myeloid cells (macrophages or DC) 
# and/or CXCR6+ T cells (CD8 or CD4) compared to patients wth poor responses?
#   
# 4. Are patients that respond better to treatment enriched in CCL5+ T cells or P2X7R+ T cells (CD8 or CD4) 
# compared with patients that do not respond?
#   
# 5. Are patients that respond better to treatment enriched for CCR7+CXCL16+ DCs 
# compared to patients that do not respond? Are they in proximity to CXCR6+ T cells?

# T-cells: CXCR6, CCL5, 'P2RX7'
# macro: CXCL16
# DCs: CXCL16, CCR7

library(data.table)
library(dplyr)
library(tibble)
library(msigdbr)
library(reshape2)
library(ggpubr)
library(ggpmisc)
library(viridis)
library(RColorBrewer)



# set variables -----------------------------------------------------------

res_dir <- '~/Documents/phd/st/geomx-processing/results/batch1-0603/'
out_dir <- '~/Documents/phd/st/geomx-processing/results/batch1-0603/for_cornell'

geomx_path <- file.path(res_dir, "geomx_qc_norm_batch_eff_rm.RDS")
deconv_path <- file.path(res_dir, "deconvolution", "bayes_prism", "bp_res_mid_lvl_ct_expr_mtx_cleaned_vst_harmony_batch_corr.RDS")
deconv_ct_frac_path <- file.path(res_dir, "deconvolution", "bayes_prism", "bp_res_mid_lvl_ct_ct_fraction.csv")
gsea_all_path <- file.path(res_dir, "pathway_analysis", "gsea",  "ssgsea_norm_harmony_batch_corr_all_msigdb.csv")

#TODO gsea_deconv_path

dir.create(out_dir)

source(file.path('~/Documents/phd/st/', 'geomx-processing', 'src', 'geomx_utils.R'))

# set important variables -------------------------------------------------

ct_interest <- c('DCs', 'Tcells', 'Macrophages')
metadt_important <- c('dcc_filename', 'Patient', 'Sample', 'NACT_status', 'PFS_months', 
                      'Segment', 'Annotation_cell', 'Aoi', 'Roi')

genes_interest <- c('CXCR6', 'CXCL16', 'CCL5', 'CCR7', 'P2RX7')

bio_vars_disc <- c('Patient', 'NACT_status', 'Annotation_cell')
bio_vars_cont <- c('PFS_months', ct_interest)

# load files and merge with metadata --------------------------------------

geomx_obj <- readRDS(geomx_path)

# check if all genes are there
all(genes_interest %in% rownames(geomx_obj@assayData$harmony_batch_corr))

metadt <- dplyr::select(sData(geomx_obj), all_of(metadt_important))

# add deconv fractions
deconv_ct_frac <- data.frame(fread(deconv_ct_frac_path))
metadt <- dplyr::left_join(metadt, deconv_ct_frac[, c('dcc_filename', ct_interest)], by = 'dcc_filename')

# add gene expression
genes_expr <- data.frame(t(geomx_obj@assayData$harmony_batch_corr[rownames(geomx_obj@assayData$harmony_batch_corr) 
                                                                  %in% genes_interest, ]))
genes_expr <- rownames_to_column(genes_expr, var = "dcc_filename")
metadt <- dplyr::left_join(metadt, genes_expr, by = 'dcc_filename')


#metadt <- metadt[metadt$Segment == 'stroma', ]


metadt_long <- melt(metadt, id.vars = colnames(metadt)[1: (ncol(metadt) - length(genes_interest))],
                    variable.name = "gene_name", 
                    value.name = "gene_expr")

# add interesting pathways gsea

# check in which pathways there are genes of interest
# msigdb_df <- msigdbr(species = "Homo sapiens")
# msigdb_df <- filter(msigdb_df, gs_cat == 'H' | gs_subcat %in% c('CP:BIOCARTA', 'CP:KEGG','GO:BP'))
# msigdb_df <- filter(msigdb_df, gene_symbol %in% genes_interest)
# 
# all_paths <- unique(msigdb_df$gs_name)


# boxplots with expr across discrete vars ---------------------------------

bio_vars_disc <- c('Patient', 'NACT_status', 'Annotation_cell')



for(color_colname in c('NACT_status', 'Annotation_cell')){
  pathway_boxplot(metadt_long, 'gene_name', 'gene_expr', color_colname, facet_var = 'Segment',
                  plot_title = paste0('genes expr per ', color_colname), 
                  output_path = file.path(out_dir, paste0('boxpl_expr_', color_colname, '_per_segm.pdf')),
                  ymin=0, ymax= max(metadt_long$gene_expr) + 1)
}

pathway_boxplot(metadt_long, 'gene_name', 'gene_expr', 'Annotation_cell', facet_var = c('Segment', 'NACT_status'),
                plot_title = 'genes expr per Annotation_cell', 
                output_path = file.path(out_dir, paste0('boxpl_expr_cell_anno_per_segm_nact.pdf')),
                ymin=0, ymax= max(metadt_long$gene_expr) + 1)


# scatterplots for continuous vars ----------------------------------------

bio_vars_cont <- c('PFS_months', ct_interest)
bio_vars_disc <- c('Patient', 'NACT_status', 'Annotation_cell')

#TODO separately per stroma/tumor segment

# each gene per PFS and ct_fraction

# coloured by c('Patient', 'NACT_status', 'Annotation_cell')

# gname <- genes_interest[1]
# 
# cont_colname <- bio_vars_cont[1]
# color_colname <- bio_vars_disc[2]
# color_colname <- 'NACT_status'
# metadt_long_all <- metadt_long
#metadt_long <- metadt_long_all[metadt_long_all$Segment == 'tumor', ]

for(gname in genes_interest){
  for(cont_colname in bio_vars_cont){
    for(color_colname in c(bio_vars_disc, 'PFS_months')){
      
      print(gname)
      print(cont_colname)
      print(color_colname)
      
      if(color_colname == 'PFS_months'){
        manual_colours <- viridis(length(unique(as.factor(metadt_long$PFS_months))))
        metadt_long$PFS_months <- as.factor(metadt_long$PFS_months)
      } else{
        manual_colours <- brewer.pal(length(unique(as.factor(metadt_long[[color_colname]]))), 'Paired')
      }
      
      if(cont_colname == 'PFS_months'){
        ggplot(data = metadt_long[metadt_long$gene_name == gname,]) +
          geom_point(aes(x = gene_expr, y = get(cont_colname), color = get(color_colname), shape = Segment)) + 
          ggtitle(paste0(gname, ' expression vs ', cont_colname)) +
          xlab(gname) +
          ylab(paste0(cont_colname)) +
          guides(color=guide_legend(title=color_colname)) +
          scale_color_manual(values=manual_colours)
          #geom_smooth(method='lm', formula= y~x) +
          # stat_poly_line() +
          #stat_poly_eq(use_label(c("R2")))
          } else {
          ggplot(data = metadt_long[metadt_long$gene_name == gname,], 
                 aes(x = gene_expr, y = get(cont_colname))) +
            geom_point(aes(color = get(color_colname), shape = Segment)) + 
            ggtitle(paste0(gname, ' expression vs ', cont_colname)) +
            xlab(gname) +
            ylab(paste0(cont_colname)) +
            guides(color=guide_legend(title=color_colname)) +
            scale_color_manual(values=manual_colours) +
            geom_smooth(method='lm', formula= y~x) +
            # stat_poly_line() +
            stat_poly_eq(use_label(c("R2")))
      }

      
      ggsave(file.path(out_dir, paste0('scatter_expr_', cont_colname, '_by_', color_colname, '_', gname, '.png')))
    }
  }
}


# CXCR6 vs CXCL16

for(color_colname in c(bio_vars_disc, 'PFS_months')){

  if(color_colname == 'PFS_months'){
    manual_colours <- viridis(length(unique(as.factor(metadt_long$PFS_months))))
    metadt$PFS_months <- as.factor(metadt$PFS_months)
  } else{
    manual_colours <- brewer.pal(length(unique(as.factor(metadt_long$PFS_months))), 'Paired')
  }
  
  print(gname)
  print(color_colname)

  
  ggplot(data = metadt, aes(x = CXCR6, y = CXCL16)) +
    geom_point(aes(color = as.factor(get(color_colname)), shape = Segment)) + 
    ggtitle('CXCR6 vs CXCL16 expression') +
    xlab('CXCR6') +
    ylab('CXCL16') +
    guides(color=guide_legend(title=color_colname)) +
    scale_color_manual(values=manual_colours) +
    geom_smooth(method='lm', formula= y~x) +
    stat_correlation(method = 'pearson')
    #stat_poly_eq(use_label(c("R2")))
  
  
  ggsave(file.path(out_dir, paste0('scatter_CXCR6_vs_CXCL16_by_', color_colname, '.png')))
}
