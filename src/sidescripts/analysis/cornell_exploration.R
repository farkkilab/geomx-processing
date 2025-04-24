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

# T-cells: CXCR6, CCL5, P2RX7
# macro: CXCL16, TLR9
# DCs: CXCL16, CCR7, TLR9

library(data.table)
library(dplyr)
library(plyr)
library(tibble)
library(msigdbr)
library(reshape2)
library(ggpubr)
library(ggpmisc)
library(viridis)
library(RColorBrewer)
library(GeomxTools)

# TODO rerun for tumor
# TODO rerun for doublepos

# set variables -----------------------------------------------------------

res_dir <- '~/Documents/phd/st/geomx-processing/results/batch1-1903/'
main_out_dir <- '~/Documents/phd/st/geomx-processing/results/batch1-1903/downstream_analysis/for_cornell'

geomx_path <- file.path(res_dir, "geomx_qc_norm_batch_eff_rm.RDS")
deconv_path <- file.path(res_dir, "deconvolution", "bayes_prism", "bp_res_mid_lvl_ct_expr_mtx_cleaned_vst_harmony_batch_corr.RDS")
deconv_ct_frac_path <- file.path(res_dir, "deconvolution", "bayes_prism", "bp_res_mid_lvl_ct_ct_fraction.csv")
gsea_all_path <- file.path(res_dir, "pathway_analysis", "gsea",  "ssgsea_norm_harmony_batch_corr_all_msigdb.csv")

# gsea_deconv_dc_path <- file.path(res_dir, "pathway_analysis", "gsea",  "ssgsea_norm_harmony_batch_corr_deconv_DCs_msigdb.csv")
# gsea_deconv_macro_path <- file.path(res_dir, "pathway_analysis", "gsea",  "ssgsea_norm_harmony_batch_corr_deconv_Macrophages_msigdb.csv")
# gsea_deconv_tcell_path <- file.path(res_dir, "pathway_analysis", "gsea",  "ssgsea_norm_harmony_batch_corr_deconv_Tcells_msigdb.csv")


dir.create(main_out_dir)

source(file.path('~/Documents/phd/st/', 'geomx-processing', 'src', 'geomx_utils.R'))

# set important variables -------------------------------------------------

ct_interest <- c('DCs', 'Tcells', 'Macrophages')
metadt_important <- c('dcc_filename', 'Patient', 'Sample', 'NACT_status', 
                      'Segment', 'Annotation_cell', 'Roi')

genes_interest <- c('CXCR6', 'CXCL16', 'CCL5', 'CCR7', 'P2RX7', 'TLR9')

bio_vars_disc <- c('Patient', 'NACT_status', 'Annotation_cell')
bio_vars_cont <- c(ct_interest)

genes_of_cells <- list(Macrophages = c('CXCL16', 'TLR9'),
                       DCs = c('CXCL16', 'TLR9', 'CCD7'),
                       Tcells = c('CXCR6', 'CCL5', 'P2RX7'))


# load files and merge with metadata for full_data ------------------------

geomx_obj <- readRDS(geomx_path)

# prepare metadata
metadt <- dplyr::select(sData(geomx_obj), all_of(metadt_important))

# add deconv fractions
deconv_ct_frac <- data.frame(fread(deconv_ct_frac_path))
metadt <- dplyr::left_join(metadt, deconv_ct_frac[, c('dcc_filename', ct_interest)], by = 'dcc_filename')

# check if all genes are there
all(genes_interest %in% rownames(geomx_obj@assayData$harmony_batch_corr))

# add gene expression
genes_expr <- data.frame(t(geomx_obj@assayData$harmony_batch_corr[rownames(geomx_obj@assayData$harmony_batch_corr) 
                                                                  %in% genes_interest, ]))
genes_expr <- rownames_to_column(genes_expr, var = "dcc_filename")
genes_expr <- dplyr::left_join(metadt, genes_expr, by = 'dcc_filename')

genes_expr_long <- melt(genes_expr, id.vars = colnames(genes_expr)[1: ncol(metadt)],
                    variable.name = "gene_name", 
                    value.name = "gene_expr")

# add interesting pathways gsea

# check in which pathways there are genes of interest
# msigdb_df <- msigdbr(species = "Homo sapiens")
# msigdb_df <- filter(msigdb_df, gs_cat == 'H' | gs_subcat %in% c('CP:BIOCARTA', 'CP:KEGG','GO:BP'))
# msigdb_df <- filter(msigdb_df, gene_symbol %in% genes_interest)
# 
# all_paths <- unique(msigdb_df$gs_name)


# load files and merge with metadata for deconvoluted data ----------------

deconv_all <- readRDS(deconv_path)

genes_expr_deconv <- lapply(1:length(genes_of_cells), function(x){
  ct_name <- names(genes_of_cells)[x]
  ct_genes <- genes_of_cells[[x]]

  expr_ct <- deconv_all[[ct_name]]
  
  ct_genes_detected <- ct_genes[which(ct_genes %in% rownames(expr_ct))]
  # fix issue if only 1 gene gere
  if(length(ct_genes_detected) == 1){
    expr_ct <- data.frame(expr_ct[ct_genes_detected, ])
  } else if (length(ct_genes_detected) > 1){
    expr_ct <- data.frame(t(expr_ct[ct_genes_detected, ]))
  } else{
    next
  }
  
  colnames(expr_ct) <- paste0(ct_genes_detected, '_', ct_name)
  expr_ct <- rownames_to_column(expr_ct, var = "dcc_filename")
})

# merge all frames

genes_expr_deconv <- join_all(genes_expr_deconv, by='dcc_filename', type='left')

genes_expr_deconv <- dplyr::left_join(metadt, genes_expr_deconv, by = 'dcc_filename')

genes_expr_deconv_long <- melt(genes_expr_deconv, id.vars = colnames(genes_expr_deconv)[1: ncol(metadt)],
                        variable.name = "gene_name", 
                        value.name = "gene_expr")

# boxplots with expr across discrete vars ---------------------------------

dt_type <- 'full_signal' # 'deconvolution', 'full_signal'
seg_type <- 'all_segments'  # 'all_segments' # 'stroma' 'tumor'
anno <- 'double_pos' # 'all', 'double_pos'

for(dt_type in c('full_signal', 'deconvolution')){ 
  for(seg_type in c('all_segments', 'stroma', 'tumor')){
    for(anno in c('all_anno', 'doublepos')){
      
      print(seg_type)
      print(anno)
      
      ###############################
      ##############################3
      out_dir <- file.path(main_out_dir, paste(dt_type, seg_type, anno, sep = '_'))
      dir.create(out_dir)
      
      if(dt_type == 'full_signal'){
        dt_long <- genes_expr_long
        dt_wide <- genes_expr
        xnames <- c('CXCR6','CXCL16')
        ynames <- c('P2RX7', 'CXCL16', 'CXCR6', 'TLR9')
      } else if(dt_type == 'deconvolution'){
        dt_long <- genes_expr_deconv_long
        dt_wide <- genes_expr_deconv
        xnames <- c('CXCR6_Tcells', 'CXCL16_Macrophages')
        ynames <- c('P2RX7_Tcells', 'CXCL16_Macrophages', 'CXCL16_DCs', 'CXCR6_Tcells')
      }
      
      # for(color_colname in c('NACT_status', 'Annotation_cell')){
      #   pathway_boxplot(dt_long, 'gene_name', 'gene_expr', color_colname, facet_var = 'Segment',
      #                   plot_title = paste0('genes expr per ', color_colname), 
      #                   output_path = file.path(out_dir, paste0('boxpl_expr_', color_colname, '_per_segm_',dt_type,  '.pdf')),
      #                   ymin=0, ymax= max(dt_long$gene_expr) + 1)
      # }
      # 
      # pathway_boxplot(dt_long, 'gene_name', 'gene_expr', 'Annotation_cell', facet_var = c('Segment', 'NACT_status'),
      #                 plot_title = 'genes expr per Annotation_cell', 
      #                 output_path = file.path(out_dir, paste0('boxpl_expr_cell_anno_per_segm_nact', dt_type, '.pdf')),
      #                 ymin=0, ymax= max(dt_long$gene_expr) + 1)
      # 
      
      # scatterplots for continuous vars ----------------------------------------
      
      
      if(seg_type == 'stroma'){
        dt_long <- dt_long[dt_long$Segment == 'stroma', ]
        dt_wide <- dt_wide[dt_wide$Segment == 'stroma', ]
      } else if(seg_type == 'tumor'){
        dt_long <- dt_long[dt_long$Segment == 'tumor', ]
        dt_wide <- dt_wide[dt_wide$Segment == 'tumor', ]
      }
      
      if(anno == 'doublepos'){
        dt_long <- dt_long[dt_long$Annotation_cell == "posCD8_posIBA1", ]
        dt_wide <- dt_wide[dt_wide$Annotation_cell == "posCD8_posIBA1", ]
      }
      
      
      
      # gene expr vs cell fraction ----------------------------------------------
      # each gene per PFS and ct_fraction
      # coloured by c('Patient', 'NACT_status', 'Annotation_cell', 'PFS_months')
      
      for(gname in unique(dt_long$gene_name)){
        for(cont_colname in ct_interest){
          for(color_colname in bio_vars_disc){
            
            print(gname)
            print(cont_colname)
            print(color_colname)
            
            # if(color_colname == 'PFS_months'){
            #   manual_colours <- viridis(length(unique(as.factor(dt_long$PFS_months))))
            #   dt_long$PFS_months <- as.factor(dt_long$PFS_months)
            # } else{
            #   manual_colours <- brewer.pal(length(unique(as.factor(dt_long[[color_colname]]))), 'Paired')
            # }
            
            manual_colours <- brewer.pal(12, 'Paired')
            
            dt_gene <- dt_long[dt_long$gene_name == gname,]
            dt_gene <- dt_gene[!is.na(dt_gene$gene_expr) & !is.na(dt_gene[[cont_colname]]), ]
            
            M1 <- lm(get(cont_colname) ~ gene_expr + Patient, data = dt_gene)
            #M2 <- lm(get(cont_colname) ~ gene_expr, data = dt_gene)
            
            #if(color_colname != 'PFS_months'){
            
            ggplot(data = dt_gene, aes(x = gene_expr, y = get(cont_colname))) +
              geom_point(aes(color = get(color_colname), shape = Segment)) + 
              ggtitle(paste0(gname, ' expression vs ', cont_colname),
                      subtitle = paste('lm(y ~ x | Patient) : R2' , round(summary(M1)$r.squared, 2),
                                       'correlation coefficient: ', 
                                       round(cor(dt_gene$gene_expr, dt_gene[[cont_colname]]), 2))) +
              xlab(gname) +
              ylab(paste0(cont_colname)) +
              guides(color=guide_legend(title=color_colname)) +
              scale_color_manual(values=manual_colours) +
              geom_smooth(method='lm', formula= y~x) +
              stat_poly_eq(use_label(c("R2", "p")))
            
            ggsave(file.path(out_dir, paste0('scatter_expr_', cont_colname, '_by_', color_colname, '_', gname,'_', dt_type, '_lm_all.png')),
                   width = 2000, height = 2000, unit = 'px')
            
            #}
            
            #####################
            if(!(color_colname %in% c('Patient', 'PFS_months'))){
              
              ggplot(data = dt_long[dt_long$gene_name == gname,],
                     aes(x = gene_expr, y = get(cont_colname), color = get(color_colname))) +
                geom_point(aes(shape = Segment)) +
                ggtitle(paste0(gname, ' expression vs ', cont_colname)) +
                xlab(gname) +
                ylab(paste0(cont_colname)) +
                guides(color=guide_legend(title=color_colname)) +
                scale_color_manual(values=manual_colours) +
                geom_smooth(method='lm', formula= y~x) +
                # stat_poly_line() +
                stat_poly_eq(use_label(c("R2", "p")))
              
              ggsave(file.path(out_dir, paste0('scatter_expr_', cont_colname, '_by_', color_colname, '_', gname,'_', dt_type, '_lm_per_group.png')),
                     width = 2000, height = 2000, unit = 'px')
            } 
          }
        }
      }
      
      
      # gene vs gene expression -------------------------------------------------
      print(' ^^^^ and now scatters ^^^^')
      for(xname in xnames){
        for(yname in ynames){
          for(color_colname in bio_vars_disc){
            
            # if(color_colname == 'PFS_months'){
            #   manual_colours <- viridis(length(unique(as.factor(dt_long$PFS_months))))
            #   dt_wide$PFS_months <- as.factor(dt_wide$PFS_months)
            # } else{
            #   manual_colours <- brewer.pal(length(unique(as.factor(dt_long$PFS_months))), 'Paired')
            # }
            
            if(xname != yname){
              print(xname)
              print(yname)
              
              manual_colours <- brewer.pal(12, 'Paired')
              
              print(color_colname)
              
              dt_wide_sel <- dt_wide[!is.na(dt_wide[[xname]]) & !is.na(dt_wide[[yname]]), ]
              
              M1 <- lm(get(yname) ~ get(xname) + Patient, data = dt_wide_sel)
              
              ggplot(data = dt_wide_sel, aes(x = get(xname), y = get(yname))) +
                geom_point(aes(color = as.factor(get(color_colname)), shape = Segment)) + 
                ggtitle(paste(xname, 'vs', yname,  'expression'),
                        subtitle = paste('lm(y ~ x | Patient) : R2' , round(summary(M1)$r.squared, 2), 
                                         'correlation coefficient: ', round(cor(dt_wide_sel[[yname]], dt_wide_sel[[xname]]), 2))) +
                xlab(xname) +
                ylab(yname) +
                guides(color=guide_legend(title=color_colname)) +
                scale_color_manual(values=manual_colours) +
                geom_smooth(method='lm', formula= y~x) +
                stat_poly_eq(use_label(c("R2", "p")))
              #stat_correlation(method = 'pearson')
              
              ggsave(file.path(out_dir, paste0('scatter_', xname, '_vs_', yname, '_by_', color_colname, '_', dt_type, '_lm_all.png')),
                     width = 2000, height = 2000, unit = 'px')
              
              if(!color_colname %in% c('Patient', 'PFS_months')){
                ggplot(data = dt_wide_sel, aes(x = get(xname), y = get(yname), color = as.factor(get(color_colname)))) +
                  geom_point(aes(shape = Segment)) + 
                  ggtitle(paste(xname, 'vs', yname,  'expression')) +
                  xlab(xname) +
                  ylab(yname) +
                  guides(color=guide_legend(title=color_colname)) +
                  scale_color_manual(values=manual_colours) +
                  geom_smooth(method='lm', formula= y~x)+
                  stat_poly_eq(use_label(c("R2", "p")))
                
                ggsave(file.path(out_dir, paste0('scatter_', xname, '_vs_', yname, '_by_', color_colname, '_', dt_type, '_lm_per_group.png')),
                       width = 2000, height = 2000, unit = 'px')
              }
            }
            
            #stat_correlation(method = 'pearson')
            #stat_poly_eq(use_label(c("R2")))
          }
        }
      }
      
      #######################################
      #######################################
      
    }
  }
}




# selecting groups for DGE ------------------------------------------------

# full signal
# genes_expr_dge <- genes_expr[genes_expr$Annotation_cell == 'posCD8_posIBA1', ]
# genes_expr_dge$CXCR6_CXCL16_sum <- genes_expr_dge$CXCL16 + genes_expr_dge$CXCR6

geomx_obj2 <- geomx_obj

out_dir_dge <- file.path(main_out_dir, 'dge')
dir.create(out_dir_dge)


dt_type <- 'full_signal'
seg_type <- 'all_segments'

for(dt_type in c('full_signal', 'deconvolution')){ 
  for(seg_type in c('all_segments', 'stroma', 'tumor')){
    
    if(dt_type == 'full_signal'){
      dt_wide <- genes_expr
      cxcr6_name <- 'CXCR6'
      cxcl16_names <- c('CXCL16')
    } else if(dt_type == 'deconvolution'){
      dt_wide <- genes_expr_deconv
      cxcr6_name <- 'CXCR6_Tcells'
      cxcl16_names <- c('CXCL16_Macrophages', 'CXCL16_DCs')
    }
    
    if(seg_type == 'stroma'){
      dt_wide <- dt_wide[dt_wide$Segment == 'stroma', ]
    } else if(seg_type == 'tumor'){
      dt_wide <- dt_wide[dt_wide$Segment == 'tumor', ]
    }

    dt_wide <- dt_wide[dt_wide$Annotation_cell == "posCD8_posIBA1", ]


    for(cxcl16_name in cxcl16_names){
      dt_wide$CXCR6_CXCL16_sum <- dt_wide[[cxcr6_name]] + dt_wide[[cxcl16_name]]
      
      perc_all <- quantile(dt_wide$CXCR6_CXCL16_sum, probs = c(0.25, 0.75), na.rm = T)
      
      dt_wide$CXCR6_CXCL16_sum_label <- ifelse(dt_wide$CXCR6_CXCL16_sum <= as.numeric(perc_all[1]), 'low', 
                                               ifelse(dt_wide$CXCR6_CXCL16_sum >= as.numeric(perc_all[2]), 'high', 'mid'))
      
      ggplot(data = dt_wide, aes(x = get(cxcr6_name), y = get(cxcl16_name), color = CXCR6_CXCL16_sum_label, shape = Segment)) +
        geom_point()
      
      ggsave(file.path(out_dir_dge, paste0(paste(dt_type, seg_type, 'doublepos', cxcr6_name, cxcl16_name, 'quartiles', sep = '_'), '.png')))
      
      # merge with geomx_obj
      dt_dge_lab <- dt_wide[, c('dcc_filename', 'CXCR6_CXCL16_sum_label')]
      colnames(dt_dge_lab) <- c('dcc_filename', paste(dt_type, seg_type, cxcr6_name, cxcl16_name, 'label', sep = '_'))
      
      pData(geomx_obj2) <- left_join(pData(geomx_obj2), dt_dge_lab)
    }
    }}

saveRDS(geomx_obj2, file.path(main_out_dir, 'dge', 'geomx_labels_for_cornell.RDS'))

# # deconv signal
# genes_expr_dge <- genes_expr_deconv[genes_expr_deconv$Annotation_cell == 'posCD8_posIBA1', ]
# genes_expr_dge$CXCR6_CXCL16_sum <- genes_expr_dge$CXCL16_DCs + genes_expr_dge$CXCR6_Tcells
# 
# 
# perc_all <- quantile(genes_expr_dge$CXCR6_CXCL16_sum, probs = c(0.25, 0.75), na.rm = T)
# 
# genes_expr_dge$CXCR6_CXCL16_sum_label <- ifelse(genes_expr_dge$CXCR6_CXCL16_sum <= as.numeric(perc_all[1]), 'low', 
#                                                 ifelse(genes_expr_dge$CXCR6_CXCL16_sum >= as.numeric(perc_all[2]), 'high', 'mid'))
# 
# ggplot(data = genes_expr_dge, aes(x = CXCR6_Tcells, y = CXCL16_DCs, color = CXCR6_CXCL16_sum_label, shape = Segment)) +
#   geom_point()
# 
# ggsave('~/Documents/phd/st/geomx-processing/results/batch1-1903/downstream_analysis/for_cornell/dge/cxcr6_cxcl16macro_deconv.pdf')
# 
# # percentiles for stroma only
# 
# perc_stroma <- quantile(genes_expr_dge$CXCR6_CXCL16_sum[genes_expr_dge$Segment == 'stroma'], probs = c(0.25, 0.75))
# 
# genes_expr_dge$CXCR6_CXCL16_sum_label_stroma <- ifelse(genes_expr_dge$CXCR6_CXCL16_sum <= as.numeric(perc_stroma[1]), 'low', 
#                                                 ifelse(genes_expr_dge$CXCR6_CXCL16_sum >= as.numeric(perc_stroma[2]), 'high', 'mid'))
# 
# genes_expr_dge$CXCR6_CXCL16_sum_label_stroma <- ifelse(genes_expr_dge$Segment == 'stroma', 
#                                                        genes_expr_dge$CXCR6_CXCL16_sum_label_stroma, NA)
# 
# 
# genes_expr_dge <- genes_expr_dge[, c('dcc_filename', 'CXCR6_CXCL16_sum_label', 'CXCR6_CXCL16_sum_label_stroma')]
# 
# # for deconv
# # colnames(genes_expr_dge) <- c('dcc_filename', 'CXCR6_CXCL16macro_deconv_sum_label', 'CXCR6_CXCL16macro_deconv_sum_label_stroma')
# colnames(genes_expr_dge) <- c('dcc_filename', 'CXCR6_CXCL16dc_deconv_sum_label', 'CXCR6_CXCL16dc_deconv_sum_label_stroma')
# 
# 
# # merge with geomx_object
# pData(geomx_obj2) <- left_join(pData(geomx_obj), genes_expr_dge)
# 
# # save RDS object
# 
# saveRDS(geomx_obj, '~/Documents/phd/st/geomx-processing/results/batch1-1903/downstream_analysis/for_cornell/geomx_labels_for_cornell.RDS')
# 

# labels_checkup
# 
# geomx_obj2 <- readRDS('~/Documents/phd/st/geomx-processing/results/batch1-1903/downstream_analysis/for_cornell/geomx_labels_for_cornell.RDS')
# metadt2 <- pData(geomx_obj2)[, c("dcc_filename", "Segment", "Sample", "Annotation_cell", "CXCR6_CXCL16_sum_label"                        
#                                  , "CXCR6_CXCL16_sum_label_stroma", "CXCR6_CXCL16macro_deconv_sum_label"            
#                                  , "CXCR6_CXCL16macro_deconv_sum_label_stroma", "CXCR6_CXCL16dc_deconv_sum_label"               
#                                  , "CXCR6_CXCL16dc_deconv_sum_label_stroma")]
# 
# metadt2 <- left_join(metadt2, genes_expr_dge[, c('dcc_filename', 'CXCR6_CXCL16_sum')])
# metadt2 <- metadt2[metadt2$Annotation_cell == 'posCD8_posIBA1', ]
# 
# 
# 
# min(metadt2$CXCR6_CXCL16_sum[metadt2$CXCR6_CXCL16dc_deconv_sum_label_stroma == 'high'], na.rm = T)
# max(metadt2$CXCR6_CXCL16_sum[metadt2$CXCR6_CXCL16dc_deconv_sum_label_stroma == 'low'], na.rm = T)
