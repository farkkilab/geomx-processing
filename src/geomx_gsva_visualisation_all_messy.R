library(ggplot2)
library(dplyr)
library(data.table)
library(ggpubr)
library(tibble)
library(circlize)
library(ComplexHeatmap)
library(colorspace)
library(ggcorrplot)
library(viridis)
library(RColorBrewer)
library(nichenetr)
library(multinichenetr)

# set up variables --------------------------------------------------------

output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2'

gsva_name <- 'texh_macro_mhc_ifng_myet_forpaper' # or caf or 'texh_macro_mhc' , 'pycr1', 'progeny'
gsva_path <- file.path(output_dir, 'gsva', paste0('gsva_', gsva_name, '.csv')) #, '_forpaper_neggeo_ntc.csv'
progeny_path <- file.path(output_dir, 'progeny', paste0('progeny_perm_neggeo_ntc.csv'))

outp2 <- ifelse(gsva_name == 'progeny', 'progeny', 'gsva')

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')
source('/media/iganiemi/T7-iga/st/st-processing/src/visium_utils.R') # move get treatment hmap to geomx_utils


dir.create(file.path(output_dir, 'corrplots'), showWarnings = T, recursive = T)

# load gsva df ------------------------------------------------------------



# MESSY CODE DOWN THERE
#######################################
#######################################
#######################################
# boxpl for ind genes
ind_genes <- as.data.frame(fread(file.path(output_dir, 'ind_genes/ind_genes_exh_lr_fin.csv')))
ind_genes$gene <- ifelse(ind_genes$gene == 'HAVCR2', 'TIM3', ind_genes$gene)
#ind_genes <- filter(ind_genes, expr < 200) #rmv 3 outliers for plotting
# TODO rerun qc with correct anno instead of manual fix
# median PFS/PFI !!!!
ind_genes$PFS_upd <- ind_genes$PFS
ind_genes$PFS_upd <- ifelse(ind_genes$Patient %in% c('S139'), 'Short', ind_genes$PFS_upd)
ind_genes$PFS_upd <- ifelse(ind_genes$Patient %in% c('S065', 'S027'), 'Long', ind_genes$PFS_upd)

ind_genes$`NACT status` <- factor(ind_genes$`NACT status`, levels = c('pre', 'post'))
ind_genes$PFS <- factor(ind_genes$PFS, levels = c('Short', 'Long'))
ind_genes$PFS_upd <- factor(ind_genes$PFS_upd, levels = c('Short', 'Long'))

ind_genes_post <- filter(ind_genes, `NACT status` == 'post') 

ind_genes$gene <- factor(ind_genes$gene,
                       levels = c('NECTIN2', 'TIGIT', 'CD96', 'CD226', 'HGF', 'CD44', 
                                  'CXCL12', 'CXCR4', 'LAG3', 'PDCD1', 'TIM3', 'CXCL9', 'CXCL10', 'CXCR3'),ordered = TRUE)

pathway_boxplot(ind_genes, 'gene', 'expr', 'NACT status', c('Annotation_cell', 'Segment'), 'gene expression',
                file.path(output_dir, 'ind_genes', paste0('box_ind_genes_full_nact_peranno_fin.pdf')),
                ymin=0, ymax = 750, manual_colours = c('#4169e1ff', '#e8c547ff'))

pathway_boxplot(ind_genes[ind_genes$Annotation_cell == 'posCD8_posIBA1', ], 'gene', 'expr', 'NACT status', 
                c('Annotation_cell', 'Segment'), 'gene expression',
                file.path(output_dir, 'ind_genes', paste0('box_ind_genes_full_nact_peranno_doublepos_fin.pdf')),
                ymin=0, ymax = 750, manual_colours = c('#4169e1ff', '#e8c547ff'))

pathway_boxplot(ind_genes_post, 'gene', 'expr', 'PFS_upd', c('Annotation_cell', 'Segment'), 'gene expression',
                file.path(output_dir, 'ind_genes', paste0('box_ind_genes_full_pfs_peranno_updated.pdf')),
                ymin=0, ymax=210, manual_colours = c('#427aa1ff', '#ebf2faff'))

pathway_boxplot(ind_genes_post[ind_genes_post$Annotation_cell == 'posCD8_posIBA1', ], 'gene', 'expr', 'PFS_upd', 
                c('Annotation_cell', 'Segment'), 'gene expression',
                file.path(output_dir, 'ind_genes', paste0('box_ind_genes_full_pfs_peranno_doublepos_updated.pdf')),
                ymin=0, ymax=210, manual_colours = c('#427aa1ff', '#ebf2faff'))
#######################################
# load gsva dataframe
gsva_df <- as.data.frame(fread(gsva_path))

#TODO  manualfix 
gsva_df$pathway <- ifelse(gsva_df$pathway == 'KEGG_APOPTOSIS', 'KEGG_APOPTOSIS1', gsva_df$pathway)
gsva_df$pathway <- ifelse(gsva_df$pathway == 'TCELL_EXHAUSTION_ZHANG', 'ZHANG_TCELL_EXHAUSTION', gsva_df$pathway)

# TODO maulafix
# median PFS/PFI !!!!
gsva_df$PFS_upd <- gsva_df$PFS
gsva_df$PFS_upd <- ifelse(gsva_df$Patient %in% c('S139'), 'Short', gsva_df$PFS_upd)
gsva_df$PFS_upd <- ifelse(gsva_df$Patient %in% c('S065', 'S027'), 'Long', gsva_df$PFS_upd)


if(gsva_name == 'progeny'){
  gsva_df <- rename(gsva_df, pathway = progeny_path, gsva_score = progeny_score)
}

gsva_df$pathway <- gsub('^[^_]*_', '', gsva_df$pathway)

path_to_rm <- c("APOPTOSIS1", "EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN", "EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP",
                "EXHAUSTED_VS_MEMORY_CD8_TCELL_DN", "EXHAUSTED_VS_MEMORY_CD8_TCELL_UP", "NAIVE_VS_EXHAUSTED_CD8_TCELL_DN",
                "NAIVE_VS_EXHAUSTED_CD8_TCELL_UP", "OVARY_CL13_MONOCYTE_MACROPHAGE")

gsva_df <- dplyr::filter(gsva_df, !(pathway %in% path_to_rm))

path_names <- unique(gsva_df$pathway) 

gsva_df$`NACT status` <- factor(gsva_df$`NACT status`, levels = c('pre', 'post'))
gsva_df$PFS <- factor(gsva_df$PFS, levels = c('Short', 'Long'))
gsva_df$PFS_upd <- factor(gsva_df$PFS_upd, levels = c('Short', 'Long'))

# make it wide
gsva_wide <- dcast(gsva_df, dcc_filename + Segment + Annotation_cell + `NACT status` + PFS + PFS_upd ~ pathway,
                   value.var = 'gsva_score')

#calculate z-score
gsva_wide_zscore <- scale(gsva_wide[, path_names]) 
gsva_wide_zscore <- cbind(gsva_wide[, c('dcc_filename', 'Segment', 'Annotation_cell', 'NACT status', 'PFS', 'PFS_upd')], gsva_wide_zscore)

gsva_df_post <- filter(gsva_df, `NACT status` == 'post')

# make it wide
gsva_wide_post <- dcast(gsva_df_post, dcc_filename + Segment + Annotation_cell + `NACT status` + PFS + PFS_upd ~ pathway,
                   value.var = 'gsva_score')

#calculate z-score
gsva_wide_post_zscore <- scale(gsva_wide_post[, path_names]) 
gsva_wide_post_zscore <- cbind(gsva_wide_post[, c('dcc_filename', 'Segment', 'Annotation_cell', 'NACT status', 'PFS', 'PFS_upd')], gsva_wide_post_zscore)


########################################
# boxplots 
gsva_df$Patient_special <- ifelse(gsva_df$Patient == 'S139', 'S139', 'other')
gsva_df$Patient_special2 <- ifelse(gsva_df$Patient == 'S015', 'S015', 'other')
gsva_df$Patient_special3 <- ifelse(gsva_df$Patient == 'S084', 'S084', 'other')

pathway_boxplot(gsva_df, 'pathway', 'gsva_score', 'Patient_special3', c('Segment'), 'gsva scores',
                file.path(output_dir, outp2, paste0('box_gsva_', gsva_name, '_patient_special3.pdf')))

pathway_boxplot(gsva_df, 'pathway', 'gsva_score', 'Patient_special3', c('Annotation_cell', 'Segment'), 'gsva scores',
                file.path(output_dir, outp2, paste0('box_gsva_', gsva_name, '_patient_special3_peranno.pdf')))

pathway_boxplot(gsva_df, 'pathway', 'gsva_score', 'Annotation_cell', c('Segment'), 'gsva scores',
                file.path(output_dir, outp2, paste0('box_gsva_', gsva_name, '_anno2.pdf')))

pathway_boxplot(gsva_df, 'pathway', 'gsva_score', 'NACT status', c('Segment'), 'gsva scores',
                file.path(output_dir, outp2, paste0('box_gsva_', gsva_name, '_nact_all2.pdf')))

pathway_boxplot(gsva_df, 'pathway', 'gsva_score', 'NACT status', c('Annotation_cell', 'Segment'), 'gsva scores',
                file.path(output_dir, outp2, paste0('box_gsva_', gsva_name, '_nact_peranno2.pdf')))

pathway_boxplot(gsva_df_post, 'pathway', 'gsva_score', 'PFS_upd', c('Segment'), 'gsva scores',
                file.path(output_dir, outp2, paste0('box_gsva_', gsva_name, '_pfs_all.pdf')))

pathway_boxplot(gsva_df_post, 'pathway', 'gsva_score', 'PFS_upd', c('Annotation_cell', 'Segment'), 'gsva scores',
                file.path(output_dir, outp2, paste0('box_gsva_', gsva_name, '_pfs_peranno.pdf')))


#######################################
# make heatmap for values + zscores


for(value_type in c('gsva', 'zscore')){
  
  for(var_name in c('Annotation_cell', 'NACT status', 'PFS_upd')){
    
    if(value_type == 'gsva'){
      gsva <- gsva_wide
      if(var_name == 'PFS_upd'){
        gsva <- gsva[gsva$`NACT status` == 'post',]
      } 
    } else if(value_type == 'zscore'){
      if(var_name == 'PFS_upd'){
        gsva <- gsva_wide_post_zscore
      } else{
        gsva <- gsva_wide_zscore
      }
    }

    # for 1 variable at the time
    gsva_mean <- gsva[, c('Segment', var_name, path_names)] %>%
      #filter(Segment == segment) %>%
      group_by(across(all_of(c('Segment', var_name)))) %>%
      summarise_all(mean, na.rm = TRUE) %>%
      ungroup()
    
    # for 1 variable + Annotation cell
    gsva_mean_peranno <- gsva[, c('Segment', 'Annotation_cell', var_name, path_names)] %>%
      #filter(Segment == segment) %>%
      group_by(across(all_of(c('Segment','Annotation_cell', var_name)))) %>%
      summarise_all(mean, na.rm = TRUE) %>%
      ungroup()
    
    min_val <- min(gsva_mean[, path_names], na.rm = T)
    max_val <- max(gsva_mean[, path_names], na.rm = T)
    min_val_peranno <- min(gsva_mean_peranno[, path_names], na.rm = T)
    max_val_peranno <- max(gsva_mean_peranno[, path_names], na.rm = T)
    
    heat_value_title <- ifelse(gsva_name == 'progeny', 'mean progeny', 'mean gsva')
    heat_value_title <- ifelse(value_type == 'zscore', paste(heat_value_title, 'z-score'), paste(heat_value_title, 'score'))
    
    ###########################
    # make hmaps for 1 variable (split per tumor/stroma)
    var_heatmap_list <- lapply(c('tumor', 'stroma'), function(s){
      heat_seg <- filter(gsva_mean, Segment == s) %>%
        dplyr::select(-Segment) %>%
        column_to_rownames(var_name)
      
      dt_seg_heat <- get_treatment_heatmap(t(as.matrix(heat_seg)), s, min_val, max_val,
                                           heat_value_title, T,
                                           color_scale = 'rb', clust_rows = T)
      
      return(dt_seg_heat)
    })
    
    # make HeatmapList object from all heatmaps in a list
    all_hmaps = NULL
    
    for(i in seq_along(var_heatmap_list)){
      print(i)
      all_hmaps = all_hmaps + var_heatmap_list[[i]]
    }
    
    # adjust length and height depending on nr of plots
    pdf(file=file.path(output_dir, outp2, paste0('heatmap_', var_name, '_', value_type, '.pdf')),
        width=length(var_heatmap_list)*2 + 3,
        height=5) # ,units="in",res=1200 for png
    
    draw(all_hmaps, ht_gap = unit(1, "cm"), 
         column_title = paste0(heat_value_title),
         column_title_gp = gpar(fontsize = 15))
    
    dev.off()
    
    #######################################################
    # make hmaps for annotation + additional variable (split per tumor/stroma)
    
    if(var_name != 'Annotation_cell'){
      seg_var_heatmap_list <- lapply(c('stroma', 'tumor'), function(s){
        heat_seg <- filter(gsva_mean_peranno, Segment == s)
        
        var_heatmap_list <- lapply(unique(unlist(gsva_mean_peranno[, var_name])), function(v){
          v <- as.character(v)
          heat_seg_var <- as.data.frame(heat_seg[heat_seg[,var_name] == v, ])
          rownames(heat_seg_var) <- heat_seg_var$Annotation_cell
          heat_seg_var <- heat_seg_var[, path_names]
          
          dt_seg_heat <- get_treatment_heatmap(t(as.matrix(heat_seg_var)), paste(s, v),
                                               min_val_peranno, max_val_peranno,
                                               heat_value_title, T,
                                               color_scale = 'rb',
                                               clust_rows = T)
          
          return(dt_seg_heat)
        })
        return(var_heatmap_list)
      })
      
      seg_var_heatmap_list <- unlist(seg_var_heatmap_list, recursive = F)
      
      # make HeatmapList object from all heatmaps in a list
      all_hmaps = NULL
      
      for(i in seq_along(seg_var_heatmap_list)){
        all_hmaps = all_hmaps + seg_var_heatmap_list[[i]]
      }
      
      # adjust length and height depending on nr of plots
      pdf(file=file.path(output_dir, outp2, paste0('heatmap_', var_name, '_', value_type, '_peranno.pdf')),
          width=length(seg_var_heatmap_list)*2 + 3,
          height=5) #,units="in",res=1200 for png
      
      draw(all_hmaps, ht_gap = unit(1, "cm"), 
           column_title = paste0(heat_value_title),
           column_title_gp = gpar(fontsize = 15))
      
      dev.off()
    }
  }
}



#####################################
#####################################
#####################################
# dotplots with log2fc
# per anno and between PFS/NACT
nact_paths <- c("TCELL_EXHAUSTION", "PD_1_SIGNALING", "CTLA4_PATHWAY", "REGULATION_OF_T_CELL_APOPTOTIC_PROCESS",
                "CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP", "CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN",
                "JAK_STAT_SIGNALING_PATHWAY", "IL2_PATHWAY")

pfs_paths <- c("CTLA4_INHIBITORY_SIGNALING","CTLA4_PATHWAY",
               "CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP", "CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN",
               "IL2_STAT5_SIGNALING","IL2_PATHWAY", "MTOR_SIGNALING_PATHWAY", "TNFA_SIGNALING_VIA_NFKB",
               "INFLAMMATORY_RESPONSE", "INTERFERON_GAMMA_RESPONSE", 
               "CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION", "MHC_CLASS_II_ANTIGEN_PRESENTATION")


for(var_name in c('NACT status', 'PFS_upd')){
  
  if(var_name %in% c('PFS_upd')){
    gsva_fordot <- gsva_wide[gsva_wide$`NACT status` == 'post',]
    #path_names <- pfs_paths #TODO change only for gsva not for progeny
  } else {
    gsva_fordot <- gsva_wide
    #path_names <- nact_paths
  }
  
  gsva_var1 <- gsva_fordot[gsva_fordot[[var_name]] == unique(gsva_fordot[[var_name]])[1], c('Segment', 'Annotation_cell', var_name, path_names)]
  gsva_var2 <- gsva_fordot[gsva_fordot[[var_name]] == unique(gsva_fordot[[var_name]])[2], c('Segment', 'Annotation_cell', var_name, path_names)]
  
  stats_path <- lapply(path_names, function(path){
    
    stats_segment_anno <- lapply(unique(gsva_fordot$Segment), function(s){
      
      stats_anno <- lapply(unique(gsva_fordot$Annotation_cell), function(a){
        
        test <- tryCatch(expr = wilcox.test(gsva_var1[[path]][gsva_var1$Segment == s & gsva_var1$Annotation_cell == a],
                                            gsva_var2[[path]][gsva_var2$Segment == s & gsva_var2$Annotation_cell == a]), 
                              error = function(e) NA)
        
        log2fc_mean <- log2(mean(gsva_var1[[path]][gsva_var1$Segment == s & gsva_var1$Annotation_cell == a], na.rm = T) /
                              mean(gsva_var2[[path]][gsva_var2$Segment == s & gsva_var2$Annotation_cell == a], na.rm = T))
        
        mean_diff <- mean(gsva_var1[[path]][gsva_var1$Segment == s & gsva_var1$Annotation_cell == a], na.rm = T) -
          mean(gsva_var2[[path]][gsva_var2$Segment == s & gsva_var2$Annotation_cell == a], na.rm = T)
        
        anno_df <- data.frame(pval = as.numeric(test$p.value),
                              wilcox_statistic = as.numeric(test$statistic),
                              log2fc = as.numeric(log2fc_mean),
                              mean_diff = as.numeric(mean_diff),
                              anno = a, 
                              segment = s, 
                              path = path
        ) 
        return(anno_df)
      })
      return(stats_anno)
    })
    
    stats_segment_anno <- unlist(stats_segment_anno, recursive = F)
    return(stats_segment_anno)
  })
  
  stats_path <- unlist(stats_path, recursive = F)
  stats_all_df <- do.call(rbind, stats_path)
  
  
  #fixing very low p-vals
  stats_all_df$pval_for_plot <- ifelse(stats_all_df$pval < 1e-10, 1e-10, stats_all_df$pval)
  
  #ensure order
  stats_all_df$path <-factor(stats_all_df$path, levels = path_names)
  
  #stats_sign_zone_df$sign <- gsub(sign_colname, '', stats_sign_zone_df$sign)
  
  stats_signif_df <- filter(stats_all_df, pval <= 0.05)
  #TODO change for this in ggplot if needed
  
  library(colorspace)
  
  plot <- ggplot(stats_all_df, aes(x = anno, y = path)) +
    geom_point(aes(color = mean_diff, size = pval_for_plot)) +
    theme_classic() +
    xlab(NULL) +
    ylab(NULL) +
    scale_size_area(
      "pval",
      trans = "log10",
      #max_size = ifelse(length(stats_sign_zone_noinf) < 10, 2.5, 2),
      max_size = 2,
      breaks = c(1e-10, 1e-5, 1e-1, 0.05),
      limits = c(1e-10, 0.05)
    ) +
    theme(
      axis.text = element_text(size = rel(0.5)),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = rel(1)),
      strip.placement = "outside",
      strip.background = element_blank(),
      plot.title = element_text(size = 10, face = "bold"),
      aspect.ratio = 1,
      panel.border = element_rect(colour = "black", size = 1.5, fill = NA)
    ) +
    scale_color_continuous_divergingx(palette = 'RdBu', mid = 0, rev = T) + 
    #scale_colour_brewer(palette = 'RdBu')
    # scale_color_gradient2(low="blue", mid="white", high="red", 
    #                       midpoint=0, limits=c(min(stats_sign_zone_df$log2fc, na.rm = T),
    #                                            max(stats_sign_zone_df$log2fc, na.rm = T))) +
    ggtitle(paste(unique(gsva_fordot[[var_name]])[1], 'vs', unique(gsva_fordot[[var_name]])[2])) +
    facet_wrap(~segment, scales = "fixed", dir="h")
  
  
  plot(plot)
  ggsave(file.path(output_dir, outp2, paste0('dotplot_', var_name, '_allpaths.pdf')),
         width=2000, height = 1000, unit='px', device='pdf')
  
}


# make corrplots between all pathways in gsva -----------------------------

# input_df_wide <- gsva_ind_genes_wide

make_corplot <- function(input_df_wide, path_names, cortitle){
  
  rownames(input_df_wide) <- input_df_wide$dcc_filename
  input_df_forcor <- dplyr::select(input_df_wide, all_of(path_names))
  
  input_corr <- round(cor(input_df_forcor, method = 'spearman'), 2)
  
  # filter rows with cor >/<0.5
  thr <- 0.6
  rows_sig <- apply(data.frame(input_corr), 1, function(r) any((r >= thr & r < 1) | (r <= -thr )))
  input_corr_sig <- input_corr[which(rows_sig), which(rows_sig)]
  input_corr_sig[(input_corr_sig > -thr) & (input_corr_sig < thr)] <- 0
  
  
  input_df_forcor_sig <- input_df_forcor[, colnames(input_corr_sig)]
  
  input_pmat <- cor_pmat(input_df_forcor)
  input_pmat_sig <- cor_pmat(input_df_forcor_sig)
  
  corrplot <- ggcorrplot(input_corr_sig, hc.order = FALSE, outline.color = "white", p.mat = input_pmat_sig,
                        title = cortitle, tl.cex = 3, pch.cex = 2) # , lab = TRUE, lab_size = 1)
  
  pdf(file=file.path(output_dir, 'corrplots', paste0('corrplot_signif_06_', cortitle, '.pdf')),
      width=5, height=6) 
  
  plot(corrplot)
  
  dev.off()
  
  return(list(input_corr, input_pmat)) # !!! without filtering
}

####################
# load and update progeny

progeny_df <- as.data.frame(fread(progeny_path))

# TODO maulafix
progeny_df$PFS_upd <- progeny_df$PFS
progeny_df$PFS_upd <- ifelse(progeny_df$Patient %in% c('S139'), 'Short', progeny_df$PFS_upd)
progeny_df$PFS_upd <- ifelse(progeny_df$Patient %in% c('S065', 'S027'), 'Long', progeny_df$PFS_upd)

# make it wide
progeny_wide <- dcast(progeny_df, dcc_filename + Segment + Annotation_cell + `NACT status` + PFS + PFS_upd ~ progeny_path,
                   value.var = 'progeny_score')
# get ind genes
ind_genes_wide <- dcast(ind_genes, dcc_filename + Segment + Annotation_cell + `NACT status` + PFS + PFS_upd ~ gene,
                        value.var = 'expr')

progeny_wide_post <- filter(progeny_wide, `NACT status` == 'post')
# combine gsva with progeny and ind genes

#gsva_ind_genes_wide <- full_join(gsva_wide, ind_genes_wide)
#gsva_ind_genes_wide <- full_join(gsva_ind_genes_wide, progeny_wide)
gsva_ind_genes_wide <- full_join(gsva_wide, progeny_wide) # without individual genes

# make corrplots
path_names <- c(unique(gsva_df$pathway), unique(progeny_df$progeny_path)) #, unique(ind_genes$gene))

#selected paths
path_names <- c('TIGIT', 'NECTIN2', 'TIM3', 
                "CTLA4_PATHWAY", "CTLA4_INHIBITORY_SIGNALING", "TCELL_EXHAUSTION", "PD_1_SIGNALING",
                "JAK_STAT_SIGNALING_PATHWAY", "IL2_STAT5_SIGNALING", "IL2_PATHWAY", "INTERFERON_GAMMA_RESPONSE",
                "TNFA_SIGNALING_VIA_NFKB")

make_corplot(gsva_ind_genes_wide, path_names, 'all')
make_corplot(gsva_ind_genes_wide[gsva_ind_genes_wide$`NACT status` == 'post', ], path_names, 'pre')
make_corplot(gsva_ind_genes_wide[gsva_ind_genes_wide$`NACT status` == 'post', ], path_names, 'post')

gsva_stroma_doublepos <- gsva_ind_genes_wide[gsva_ind_genes_wide$Segment == 'stroma' & gsva_ind_genes_wide$Annotation_cell == 'posCD8_posIBA1', ]
gsva_tumor_doublepos <- gsva_ind_genes_wide[gsva_ind_genes_wide$Segment == 'tumor' & gsva_ind_genes_wide$Annotation_cell == 'posCD8_posIBA1', ]

gsva_stroma_doubleneg <- gsva_ind_genes_wide[gsva_ind_genes_wide$Segment == 'stroma' & gsva_ind_genes_wide$Annotation_cell == 'negCD8_negIBA1', ]
gsva_tumor_doubleneg <- gsva_ind_genes_wide[gsva_ind_genes_wide$Segment == 'tumor' & gsva_ind_genes_wide$Annotation_cell == 'negCD8_negIBA1', ]

make_corplot(gsva_stroma_doublepos, path_names, 'stroma_CD8+IBA1+')
make_corplot(gsva_tumor_doublepos, path_names, 'tumor_CD8+IBA1+')

# pre and post
make_corplot(gsva_stroma_doublepos[gsva_stroma_doublepos$`NACT status` == 'pre', ], path_names, 'stroma_CD8+IBA1+_pre')
make_corplot(gsva_tumor_doublepos[gsva_tumor_doublepos$`NACT status` == 'pre', ], path_names, 'tumor_CD8+IBA1+_pre')

make_corplot(gsva_stroma_doublepos[gsva_stroma_doublepos$`NACT status` == 'post', ], path_names, 'stroma_CD8+IBA1+_post')
make_corplot(gsva_tumor_doublepos[gsva_tumor_doublepos$`NACT status` == 'post', ], path_names, 'tumor_CD8+IBA1+_post')

# long and short
make_corplot(gsva_stroma_doublepos[gsva_stroma_doublepos$`NACT status` == 'post' & gsva_stroma_doublepos$PFS_upd == 'Long', ], path_names, 'stroma_CD8+IBA1+_post_long')
make_corplot(gsva_tumor_doublepos[gsva_tumor_doublepos$`NACT status` == 'post' & gsva_tumor_doublepos$PFS_upd == 'Long', ], path_names, 'tumor_CD8+IBA1+_post_long')

make_corplot(gsva_stroma_doublepos[gsva_stroma_doublepos$`NACT status` == 'post' & gsva_stroma_doublepos$PFS_upd == 'Short', ], path_names, 'stroma_CD8+IBA1+_post_short')
make_corplot(gsva_tumor_doublepos[gsva_tumor_doublepos$`NACT status` == 'post' & gsva_tumor_doublepos$PFS_upd == 'Short', ], path_names, 'tumor_CD8+IBA1+_post_short')

########################
# doubleneg
# pre and post
make_corplot(gsva_stroma_doubleneg[gsva_stroma_doubleneg$`NACT status` == 'pre', ], path_names, 'stroma_doubleneg_pre')
make_corplot(gsva_tumor_doubleneg[gsva_tumor_doubleneg$`NACT status` == 'pre', ], path_names, 'tumor_doubleneg_pre')

make_corplot(gsva_stroma_doubleneg[gsva_stroma_doubleneg$`NACT status` == 'post', ], path_names, 'stroma_doubleneg_post')
make_corplot(gsva_tumor_doubleneg[gsva_tumor_doubleneg$`NACT status` == 'post', ], path_names, 'tumor_doubleneg_post')

# long and short
# make_corplot(gsva_stroma_doubleneg[gsva_stroma_doubleneg$`NACT status` == 'post' & gsva_stroma_doubleneg$PFS_upd == 'Long', ], path_names, 'stroma_doubleneg_post_long')
# make_corplot(gsva_tumor_doubleneg[gsva_tumor_doubleneg$`NACT status` == 'post' & gsva_tumor_doubleneg$PFS_upd == 'Long', ], path_names, 'tumor_doubleneg_post_long')
# 
# make_corplot(gsva_stroma_doubleneg[gsva_stroma_doubleneg$`NACT status` == 'post' & gsva_stroma_doubleneg$PFS_upd == 'Short', ], path_names, 'stroma_doubleneg_post_short')
# make_corplot(gsva_tumor_doubleneg[gsva_tumor_doubleneg$`NACT status` == 'post' & gsva_tumor_doubleneg$PFS_upd == 'Short', ], path_names, 'tumor_doubleneg_post_short')

#############################
path_groups <- data.frame(path = path_names, group = NA)

il2_jak_stat <- c("IL2_STAT5_SIGNALING", "JAK_STAT_SIGNALING_PATHWAY", "IL2_PATHWAY", "JAK.STAT")
texh <- c("CTLA4_PATHWAY", "CTLA4_INHIBITORY_SIGNALING", "PD_1_SIGNALING", 
          "TCELL_EXHAUSTION", "REGULATION_OF_T_CELL_APOPTOTIC_PROCESS")
# "PDCD1", "TIM3", "TIGIT", "CD96", "NECTIN2", "CXCR3", "CXCL9", "IL2RG", "IL2RB" 
tnf_nfkb <- c("TNFA_SIGNALING_VIA_NFKB", "NFkB", "TNFa")
ifng <- c("INTERFERON_GAMMA_RESPONSE")
mtor <- c("MTOR_SIGNALING_PATHWAY", "PI3K")
mapk <- c("MAPK")
tgfb_wnt <- c("TGFb", "WNT")
m1 <- c("MACROPHAGE_M1_VS_M2_UP", "CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP")
m2 <- c("MACROPHAGE_M1_VS_M2_DN", "CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN")
mhc <- c("MHC_PATHWAY", "CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION", "MHC_CLASS_II_ANTIGEN_PRESENTATION")
emt <- c("EPITHELIAL_MESENCHYMAL_TRANSITION")
other <- c("APOPTOSIS", "HYPOXIA", "INFLAMMATORY_RESPONSE",
           "MACROPHAGE_MARKERS", "Androgen", "EGFR", "Estrogen", "Hypoxia", "p53", "Trail", "VEGF",
           "PDCD1", "TIM3", "TIGIT", "CD96", "NECTIN2", "CXCR3", "CXCL9", "IL2RG", "IL2RB")

path_groups <- mutate(path_groups, group = ifelse(path %in% il2_jak_stat, 'il2_jak_stat', group))
path_groups <- mutate(path_groups, group = ifelse(path %in% texh, 'texh', group))
path_groups <- mutate(path_groups, group = ifelse(path %in% tnf_nfkb, 'tnf_nfkb', group))
path_groups <- mutate(path_groups, group = ifelse(path %in% ifng, 'ifng', group))
path_groups <- mutate(path_groups, group = ifelse(path %in% mtor, 'mtor', group))
path_groups <- mutate(path_groups, group = ifelse(path %in% mapk, 'mapk', group))
path_groups <- mutate(path_groups, group = ifelse(path %in% tgfb_wnt, 'tgfb_wnt', group))
path_groups <- mutate(path_groups, group = ifelse(path %in% m1, 'm1', group))
path_groups <- mutate(path_groups, group = ifelse(path %in% m2, 'm2', group))
path_groups <- mutate(path_groups, group = ifelse(path %in% mhc, 'mhc', group))
path_groups <- mutate(path_groups, group = ifelse(path %in% emt, 'emt', group))
path_groups <- mutate(path_groups, group = ifelse(path %in% other, 'other', group))


##############################
get_cor_counts <- function(input_df_wide){
  cor_obj <- make_corplot(input_df_wide, path_names, 'test-cor')
  cor_obj  <- full_join(melt(cor_obj[[1]], value.name = 'cor'), melt(cor_obj[[2]], value.name = 'pval'))
  
  # filter with 0.6 it ignores neg correlation, but there's no in this dataset
  cor_obj <- filter(cor_obj, cor > 0.6 & pval <= 0.05) %>% 
    filter(Var1 != Var2) %>%
    left_join(path_groups, by = c('Var1' = 'path')) %>%
    rename(group_var1 = group) %>%
    left_join(path_groups, by = c('Var2' = 'path')) %>%
    rename(group_var2 = group) 
  
  cor_obj_sel <- distinct(cor_obj, cor, pval, .keep_all = T) %>% # hakierskie, pvals are unique xd
    dplyr::filter(group_var1 != 'other' & group_var2 != 'other')
  
  cor_counts <- as.data.frame(table(cor_obj_sel[, c("group_var1", "group_var2")]))
  
  # count grp1-grp2 + grp2-grp1
  cor_counts$group_var12 <- apply(cor_counts[c('group_var1', 'group_var2')], 1, function(x){paste(sort(x), collapse = '/')})
  cor_counts$Sum <- with(cor_counts, ave(Freq, group_var12, FUN = sum))
  
  cor_counts <- distinct(cor_counts, group_var12, Sum)
  
  return(cor_counts)
}

# pre and post

str_pre <- get_cor_counts(gsva_stroma_doublepos[gsva_stroma_doublepos$`NACT status` == 'pre', ])
tum_pre <- get_cor_counts(gsva_tumor_doublepos[gsva_tumor_doublepos$`NACT status` == 'pre', ])

str_post <- get_cor_counts(gsva_stroma_doublepos[gsva_stroma_doublepos$`NACT status` == 'post', ])
tum_post <- get_cor_counts(gsva_tumor_doublepos[gsva_tumor_doublepos$`NACT status` == 'post', ])


# prepost_comb <- full_join(str_pre, str_post, by = c('group_var1', 'group_var2')) %>%
#   rename(str_pre = Freq.x, str_post = Freq.y) %>%
#   full_join(tum_pre, by = c('group_var1', 'group_var2')) %>%
#   rename(tum_pre = Freq) %>%
#   full_join(tum_post, by = c('group_var1', 'group_var2')) %>%
#   rename(tum_post = Freq) %>%
#   mutate_at(c('str_pre', 'str_post', 'tum_pre', 'tum_post'), ~coalesce(.,0)) %>%
#   filter((str_pre + str_post + tum_pre + tum_post) > 0) %>%
#   mutate(group_var12 = paste0(group_var1, '/', group_var2))

prepost_comb <- full_join(str_pre, str_post, by = c('group_var12')) %>%
  rename(str_pre = Sum.x, str_post = Sum.y) %>%
  full_join(tum_pre, by = c('group_var12')) %>%
  rename(tum_pre = Sum) %>%
  full_join(tum_post, by = c('group_var12')) %>%
  rename(tum_post = Sum) %>%
  mutate_at(c('str_pre', 'str_post', 'tum_pre', 'tum_post'), ~coalesce(.,0)) %>%
  filter((str_pre + str_post + tum_pre + tum_post) > 0) 

# short and long

str_short <- get_cor_counts(gsva_stroma_doublepos[gsva_stroma_doublepos$`NACT status` == 'post' & gsva_stroma_doublepos$PFS_upd == 'Short', ])
tum_short <- get_cor_counts(gsva_tumor_doublepos[gsva_tumor_doublepos$`NACT status` == 'post' & gsva_tumor_doublepos$PFS_upd == 'Short', ])

str_long <- get_cor_counts(gsva_stroma_doublepos[gsva_stroma_doublepos$`NACT status` == 'post' & gsva_stroma_doublepos$PFS_upd == 'Long', ])
tum_long <- get_cor_counts(gsva_tumor_doublepos[gsva_tumor_doublepos$`NACT status` == 'post' & gsva_tumor_doublepos$PFS_upd == 'Long', ])

shortlong_comb <- full_join(str_short, str_long, by = c('group_var12')) %>%
  rename(str_short = Sum.x, str_long = Sum.y) %>%
  full_join(tum_short, by = c('group_var12')) %>%
  rename(tum_short = Sum) %>%
  full_join(tum_long, by = c('group_var12')) %>%
  rename(tum_long = Sum)  %>%
  mutate_at(c('str_short', 'str_long', 'tum_short', 'tum_long'), ~coalesce(.,0)) %>%
  filter((str_short + str_long + tum_short + tum_long) > 0)


fwrite(prepost_comb, file.path(output_dir, 'corrplots', 'corsum_prepost.csv'))
fwrite(shortlong_comb, file.path(output_dir, 'corrplots', 'corsum_shortlong.csv'))

#################
# plots

prepost_comb_long <- melt(prepost_comb)
shortlong_comb_long <- melt(shortlong_comb)

ggplot(prepost_comb_long, aes(x = group_var12, y = value, fill = variable)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5))

ggsave(file.path(output_dir, 'corrplots', 'corsum_prepost.png'), width = 2000, height = 2000, units = 'px')

ggplot(shortlong_comb_long, aes(x = group_var12, y = value, fill = variable)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))

ggsave(file.path(output_dir, 'corrplots', 'corsum_shortlong.png'), width = 2000, height = 2000, units = 'px')

######################
# circle plots

group_name <- 'texh'

groups_nr_df <- as.data.frame(table(path_groups$group))
colnames(groups_nr_df) <- c('group', 'group_nr')

grid_col <- viridis(nrow(groups_nr_df))
names(grid_col) <- as.character(groups_nr_df$group)

comb_df <- prepost_comb
comb_name <- 'prepost'

#comb_list <- list('prepost' = prepost_comb, 'shortlong' = shortlong_comb)

# for 1 group
comb_group_df <- comb_df[grepl(group_name, comb_df$group_var12), ]

comb_group_df <- mutate(comb_group_df, group = gsub(paste0('/', group_name, '|', group_name, '/'), '', group_var12)) %>%
  left_join(groups_nr_df) %>%
  mutate(group_comb_nr = group_nr * as.numeric(groups_nr_df$group_nr[groups_nr_df$group == group_name])) %>%
  select(-c("group_var12", "group_nr"))

comb_group_df_long <- melt(comb_group_df, id.vars = c('group', 'group_comb_nr'), value.name = 'cor_nr') 
comb_group_df_long$main_group <- group_name

# for all groups
comb_df_all <- melt(comb_df, id.vars = c('group_var12'), value.name = 'cor_nr')
comb_df_all <- cbind(comb_df_all, colsplit(string=comb_df_all$group_var12, pattern="/", names=c("group1", "group2")))

for(aoi in colnames(comb_df)[-1]){
  print(aoi)
  #aoi <- 'tum_long'
  aoi_group <- comb_group_df_long[comb_group_df_long$variable == aoi, c('group','main_group', 'cor_nr')]
  aoi_group_rep <- as.data.frame(lapply(aoi_group, rep, aoi_group$cor_nr))
  
  pdf(file.path(output_dir, 'corrplots', paste0('circle_', group_name, '_', comb_name, '_', aoi, '.pdf')))
  par()
  circos.par()
  chordDiagram(aoi_group, grid.col = grid_col)
  circos.clear()
  dev.off()
  
  aoi_all <- comb_df_all[comb_df_all$variable == aoi, c('group1','group2', 'cor_nr')]
  aoi_all_rep <- as.data.frame(lapply(aoi_all, rep, aoi_all$cor_nr))
  
  pdf(file.path(output_dir, 'corrplots', paste0('circle_all_', comb_name, '_', aoi,'.pdf')))
  par()
  circos.par()
  chordDiagram(aoi_all, grid.col = grid_col)
  circos.clear()
  dev.off()
}



#################################
#################################
#################################
# final circle corrplot versions
gsva_stroma_doublepos <- gsva_ind_genes_wide[gsva_ind_genes_wide$Segment == 'stroma' & gsva_ind_genes_wide$Annotation_cell == 'posCD8_posIBA1', ]
gsva_tumor_doublepos <- gsva_ind_genes_wide[gsva_ind_genes_wide$Segment == 'tumor' & gsva_ind_genes_wide$Annotation_cell == 'posCD8_posIBA1', ]

groups_nr_df <- as.data.frame(table(path_groups$group))
colnames(groups_nr_df) <- c('group', 'group_nr')

grid_col <- viridis(nrow(groups_nr_df))
names(grid_col) <- as.character(groups_nr_df$group)

#
group_name <- 'texh'
# 
input_df_list <- list('stroma_pre' = gsva_stroma_doublepos[gsva_stroma_doublepos$`NACT status` == 'pre', ],
                      'stroma_post' = gsva_stroma_doublepos[gsva_stroma_doublepos$`NACT status` == 'post', ],
                      'tumor_pre' = gsva_tumor_doublepos[gsva_tumor_doublepos$`NACT status` == 'pre', ],
                      'tumor_post' = gsva_tumor_doublepos[gsva_tumor_doublepos$`NACT status` == 'post', ],
                      'stroma_short' = gsva_stroma_doublepos[gsva_stroma_doublepos$`NACT status` == 'post' & gsva_stroma_doublepos$PFS_upd == 'Short', ],
                      'stroma_long' = gsva_stroma_doublepos[gsva_stroma_doublepos$`NACT status` == 'post' & gsva_stroma_doublepos$PFS_upd == 'Long', ],
                      'tumor_short' = gsva_tumor_doublepos[gsva_tumor_doublepos$`NACT status` == 'post' & gsva_tumor_doublepos$PFS_upd == 'Short', ],
                      'tumor_long' = gsva_tumor_doublepos[gsva_tumor_doublepos$`NACT status` == 'post' & gsva_tumor_doublepos$PFS_upd == 'Long', ])

cor_thr <- 0.6

sapply(1:length(input_df_list), function(x){
  input_df <- input_df_list[[x]]
  aoi <- names(input_df_list)[x]
  
  print(aoi)
  
  cor_obj <- make_corplot(input_df, path_names, 'test-cor')
  
  cor_obj  <- full_join(melt(cor_obj[[1]], value.name = 'cor'), melt(cor_obj[[2]], value.name = 'pval')) %>%
    mutate(cor = ifelse(pval > 0.05, 0, cor)) %>% # rmv non-significant values
    filter(Var1 != Var2) %>%
    left_join(path_groups, by = c('Var1' = 'path')) %>%
    rename(group_var1 = group) %>%
    left_join(path_groups, by = c('Var2' = 'path')) %>%
    rename(group_var2 = group) %>%
    #distinct(cor, pval, .keep_all = T) %>% # hakierskie rmv duplicated part of the mtx  ver _half of plots
    filter(group_var1 != 'other' & group_var2 != 'other') #rmv group 'other'
  
  cor_all <- select(cor_obj, c(group_var1, group_var2, cor))
  
  # all vars
  pdf(file.path(output_dir, 'corrplots', 'circle-clean-fin', paste0('circle_all_', aoi, '.pdf')))
  par(cex = 1.5, mar = c(0, 0, 0, 0))
  circos.par()
  chordDiagram(cor_all, link.visible = cor_all$cor >= cor_thr, grid.col = grid_col, symmetric = F,
               annotationTrack = c("name", "grid"))
  circos.clear()
  dev.off()
  
  # only with chosen var
  cor_group <- select(cor_obj, c(group_var1, group_var2, cor)) %>%
    #filter(group_var2 == group_name)
    filter(group_var1 == group_name | group_var2 == group_name)
  
  # only categories with any value > 0.6
  cor_group_clean <- cor_group %>%
    group_by(group_var1, group_var2) %>%
    mutate(group_max = max(cor)) %>%
    ungroup() %>%
    filter(group_max >= cor_thr) %>%
    select(-group_max)
    # filter(group_var2 == group_name) %>%
    # arrange(group_var2, group_var1, cor)
  
  
  pdf(file.path(output_dir, 'corrplots', 'circle-clean-fin', paste0('circle_', group_name, '_', aoi, '.pdf')))
  par(cex = 2, mar = c(0, 0, 0, 0))
  circos.par()
  chordDiagram(cor_group_clean, link.visible = cor_group_clean$cor >= cor_thr, grid.col = grid_col, 
               annotationTrack = c("name", "grid"))
  circos.clear()
  dev.off()
  
})

##############################################
##############################################
# circleplots from NicheNetr Package
#input_df_wide <- gsva_stroma_doublepos[gsva_stroma_doublepos$`NACT status` == 'pre', ]
input_df_wide <- gsva_ind_genes_wide[gsva_ind_genes_wide$`NACT status` == 'post' &
                                       gsva_ind_genes_wide$Annotation_cell == 'posCD8_posIBA1' &
                                       gsva_ind_genes_wide$PFS_upd == 'Long' &
                                       gsva_ind_genes_wide$Segment == 'tumor', ]

cor_obj <- make_corplot(input_df_wide, path_names, 'test-cor')
cor_obj  <- full_join(melt(cor_obj[[1]], value.name = 'cor'), melt(cor_obj[[2]], value.name = 'pval'))

# filter with 0.6 it ignores neg correlation, but there's no significant in this dataset
cor_obj <- filter(cor_obj, cor >= 0.7 & pval <= 0.05) %>% 
  filter(Var1 != Var2) %>%
  left_join(path_groups, by = c('Var1' = 'path')) %>%
  rename(group_var1 = group) %>%
  left_join(path_groups, by = c('Var2' = 'path')) %>%
  rename(group_var2 = group) %>%
  dplyr::filter(group_var1 != 'other' & group_var2 != 'other')

cor_obj_sel <- distinct(cor_obj, cor, pval, .keep_all = T)  # hakierskie, pvals are unique xd - HALF
#cor_obj_sel <- cor_obj

#hacking colnames for nichenetr
colnames(cor_obj_sel) <- c('ligand', 'receptor', 'prioritization_score', 'pval', 'sender', 'receiver')
cor_obj_sel$group <- 'oo'
cor_obj_sel$id <- paste(cor_obj_sel$ligand, cor_obj_sel$receptor, cor_obj_sel$sender, cor_obj_sel$receiver, sep = '_')
length(unique(cor_obj_sel$sender))
length(unique(cor_obj_sel$receiver))
intersect(unique(cor_obj_sel$sender), unique(cor_obj_sel$receiver))

# colors_list <- as.list(viridis(length(unique(c(cor_obj_sel$sender, cor_obj_sel$receiver)))))
# names(colors_list) <- unique(c(cor_obj_sel$sender, cor_obj_sel$receiver))

colors_list <- as.list(viridis(length(unique(path_groups$group))))
names(colors_list) <- unique(path_groups$group)

kk <- make_circos_one_group_iga(cor_obj_sel, colors_list, colors_list)
kk

###########
colors_sender <- colors_list
colors_receiver <- colors_list

prioritized_tbl_oi <- cor_obj_sel


make_circos_one_group_iga = function(prioritized_tbl_oi, colors_sender, colors_receiver){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("circlize")
  
  prioritized_tbl_oi = prioritized_tbl_oi %>% dplyr::ungroup() # if grouped: things will be messed up downstream
  
  # Link each cell type to a color
  grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% names(), color_ligand_type = colors_sender)
  grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% names(), color_receptor_type = colors_receiver)
  
  # Make plot
  
  # deal with duplicated sector names
  # dplyr::rename the ligands so we can have the same ligand in multiple senders (and receptors in multiple receivers)
  # only do it with duplicated ones!
  circos_links = prioritized_tbl_oi %>% dplyr::rename(weight = prioritization_score)
  
  df = circos_links
  
  ligand.uni = unique(df$ligand)
  for (i in 1:length(ligand.uni)) {
    df.i = df[df$ligand == ligand.uni[i], ]
    sender.uni = unique(df.i$sender)
    for (j in 1:length(sender.uni)) {
      df.i.j = df.i[df.i$sender == sender.uni[j], ]
      df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
      df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
    }
  }
  receptor.uni = unique(df$receptor)
  for (i in 1:length(receptor.uni)) {
    df.i = df[df$receptor == receptor.uni[i], ]
    receiver.uni = unique(df.i$receiver)
    for (j in 1:length(receiver.uni)) {
      df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
      df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
      df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
    }
  }
  
  intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
  
  while(length(intersecting_ligands_receptors) > 0){
    df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
    df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
    df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
    df = dplyr::bind_rows(df_unique, df_duplicated)
    intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
  }
  
  circos_links = df
  
  # Link ligands/Receptors to the colors of senders/receivers
  circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
  links_circle = circos_links %>% dplyr::distinct(ligand,receptor, weight)
  ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
  grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
  receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
  grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
  grid_col =c(grid_ligand_color,grid_receptor_color)
  names(grid_col) <- c(names(grid_ligand_color), names(grid_receptor_color))
  grid_col <- grid_col[!duplicated(grid_col)]
  
  # Define order of the ligands and receptors and the gaps
  ligand_order = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
    ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
  }) %>% unlist()
  
  receptor_order = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
    receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(receptor) %>% dplyr::distinct(receptor)
  }) %>% unlist()
  
  order = c(ligand_order,receptor_order)
  
  width_same_cell_same_ligand_type = 0.275
  width_different_cell = 3
  width_ligand_receptor = 9
  width_same_cell_same_receptor_type = 0.275
  
  sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
    sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
    gap = width_different_cell
    return(c(sector,gap))
  }) %>% unlist()
  sender_gaps = sender_gaps[-length(sender_gaps)]
  
  receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
    sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
    gap = width_different_cell
    return(c(sector,gap))
  }) %>% unlist()
  receiver_gaps = receiver_gaps[-length(receiver_gaps)]
  
  gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
  
  # print(length(gaps))
  # print(length(union(circos_links$ligand, circos_links$receptor) %>% unique()))
  if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
    warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
  }
  
  grid_col =c(grid_ligand_color,grid_receptor_color)
  #grid_col <- grid_col[!duplicated(grid_col)]
  
  
  links_circle$weight[links_circle$weight == 0] = 0.01
  circos.clear()
  circos.par(gap.degree = gaps)
  chordDiagram(links_circle,
               directional = 0,
               order=order,
               link.sort = TRUE,
               link.decreasing = TRUE,
               grid.col = grid_col,
               # transparency = transparency,
               diffHeight = 0.0075,
               #direction.type = c("diffHeight", "arrows"),
               link.visible = links_circle$weight > 0.01,
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.175),
               grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
               reduce = 0,
               scale = TRUE)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.5)
  }, bg.border = NA) #
  
  p_circos = recordPlot()
  
  
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  grid_col_all = c(colors_receiver, colors_sender)
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                  type = "grid",
                                  #legend_gp = grid::gpar(fill = colors_receiver[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                  legend_gp = grid::gpar(fill = as.character(colors_receiver[prioritized_tbl_oi$receiver %>% unique() %>% sort()])),
                                  title_position = "topleft",
                                  title = "Receiver")
  ComplexHeatmap::draw(legend, just = c("left", "bottom"))

  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                  type = "grid",
                                  #legend_gp = grid::gpar(fill = colors_sender[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                  legend_gp = grid::gpar(fill = as.character(colors_sender[prioritized_tbl_oi$sender %>% unique() %>% sort()])),
                                  title_position = "topleft",
                                  title = "Sender")
  ComplexHeatmap::draw(legend, just = c("left", "top"))

  p_legend = grDevices::recordPlot()

  p_circos$legend = p_legend
  
  return(p_circos)
}

