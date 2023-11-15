library(ggplot2)
library(dplyr)
library(data.table)
library(ggpubr)
library(tibble)
library(circlize)
library(ComplexHeatmap)
library(colorspace)

#TODO adjust for PROGENY

output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact'

gsva_name <- 'texh_macro_mhc' # or caf or 'texh_macro_mhc' , 'pycr1
gsva_path <- file.path(output_dir, 'gsva', paste0('gsva_', gsva_name, '.csv'))
#gsva_path <- file.path(output_dir, 'progeny', paste0('progeny_perm.csv'))

outp2 <- ifelse(gsva_name == 'progeny', 'progeny', 'gsva')

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')
source('/media/iganiemi/T7-iga/st/st-processing/src/visium_utils.R') # move get treatment hmap to geomx_utils

#######################################
# load gsva dataframe
gsva_df <- as.data.frame(fread(gsva_path))

if(gsva_name == 'progeny'){
  gsva_df <- rename(gsva_df, pathway = progeny_path, gsva_score = progeny_score)
}

path_names <- unique(gsva_df$pathway) 
 

gsva_df_post <- filter(gsva_df, `NACT status` == 'post')

# make it wide
gsva_wide <- dcast(gsva_df, dcc_filename + Segment + Annotation_cell + `NACT status` + PFS ~ pathway,
                   value.var = 'gsva_score')

#calculate z-score
gsva_wide_zscore <- scale(gsva_wide[, path_names]) 
gsva_wide_zscore <- cbind(gsva_wide[, c('dcc_filename', 'Segment', 'Annotation_cell', 'NACT status', 'PFS')], gsva_wide_zscore)


########################################
# boxplots 

pathway_boxplot(gsva_df, 'pathway', 'gsva_score', 'Annotation_cell', c('Segment'), 'gsva scores',
                file.path(output_dir, outp2, paste0('box_gsva_', gsva_name, '_anno2.png')))

pathway_boxplot(gsva_df, 'pathway', 'gsva_score', 'NACT status', c('Segment'), 'gsva scores',
                file.path(output_dir, outp2, paste0('box_gsva_', gsva_name, '_nact_all2.png')))

pathway_boxplot(gsva_df, 'pathway', 'gsva_score', 'NACT status', c('Annotation_cell', 'Segment'), 'gsva scores',
                file.path(output_dir, outp2, paste0('box_gsva_', gsva_name, '_nact_peranno2.png')))

pathway_boxplot(gsva_df_post, 'pathway', 'gsva_score', 'PFS', c('Segment'), 'gsva scores',
                file.path(output_dir, outp2, paste0('box_gsva_', gsva_name, '_pfs_all2.png')))

pathway_boxplot(gsva_df_post, 'pathway', 'gsva_score', 'PFS', c('Annotation_cell', 'Segment'), 'gsva scores',
                file.path(output_dir, outp2, paste0('box_gsva_', gsva_name, '_pfs_peranno2.png')))


#######################################
# pairwise t-test, log2fc (doublepos vs doubleneg)

# make heatmap for values + zscores


for(value_type in c('gsva', 'zscore')){
  
  if(value_type == 'gsva'){
    gsva <- gsva_wide
  } else if(value_type == 'zscore'){
    gsva <- gsva_wide_zscore
  }
  
  for(var_name in c('Annotation_cell', 'NACT status', 'PFS')){
    
    if(var_name == 'PFS'){
      gsva <- gsva[gsva$`NACT status` == 'post',]
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
    
    # #TODO pls make it nicer
    # if(value_type == 'gsva'){
    #   min_val <- -1
    #   max_val <- 1
    #   min_val_peranno <- -1
    #   max_val_peranno <- 1
    #   heat_value_title <- 'mean gsva score'
    # } else if(value_type == 'zscore'){
    # 
    #   heat_value_title <- 'mean gsva z-score'
    # }
    # 
    ###########################
    # make hmaps for 1 variable (split per tumor/stroma)
    var_heatmap_list <- lapply(c('tumor', 'stroma'), function(s){
      heat_seg <- filter(gsva_mean, Segment == s) %>%
        dplyr::select(-Segment) %>%
        column_to_rownames(var_name)
      
      dt_seg_heat <- get_treatment_heatmap(t(as.matrix(heat_seg)), s, min_val, max_val,
                                           heat_value_title, T,
                                           color_scale = 'rb', clust_rows = T)
      
      #return(list(heat = dt_seg_heat, heat_mtx = as.matrix(heat_seg))) #TODO why return mtx?
      return(dt_seg_heat)
    })
    
    # combine hmaps
    #var_heatmap_list <- pmap(var_heatmap_list, list)
    
    
    # make HeatmapList object from all heatmaps in a list
    all_hmaps = NULL
    
    for(i in seq_along(var_heatmap_list)){
      print(i)
      all_hmaps = all_hmaps + var_heatmap_list[[i]]
    }
    
    # adjust length and height depending on nr of plots
    png(filename=file.path(output_dir, outp2, paste0('heatmap_', var_name, '_', value_type, '2.png')),
        width=length(var_heatmap_list)*2 + 3,
        height=5,units="in",res=1200)
    
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
      
      #TODO plot in a grid (?)
      # adjust length and height depending on nr of plots
      png(filename=file.path(output_dir, outp2, paste0('heatmap_', var_name, '_', value_type, '_peranno2.png')),
          width=length(seg_var_heatmap_list)*2 + 3,
          height=5,units="in",res=1200)
      
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


for(var_name in c('NACT status', 'PFS')){
  
  if(var_name == 'PFS'){
    gsva_fordot <- gsva_wide[gsva_wide$`NACT status` == 'post',]
  } else {
    gsva_fordot <- gsva_wide
  }
  
  gsva_var1 <- gsva_fordot[gsva_fordot[[var_name]] == unique(gsva_fordot[[var_name]])[1], c('Segment', 'Annotation_cell', var_name, path_names)]
  gsva_var2 <- gsva_fordot[gsva_fordot[[var_name]] == unique(gsva_fordot[[var_name]])[2], c('Segment', 'Annotation_cell', var_name, path_names)]
  
  stats_path <- lapply(path_names, function(path){
    
    stats_segment_anno <- lapply(unique(gsva_fordot$Segment), function(s){
      
      stats_anno <- lapply(unique(gsva_fordot$Annotation_cell), function(a){
        
        test_pval <- tryCatch(expr = wilcox.test(gsva_var1[[path]][gsva_var1$Segment == s & gsva_var1$Annotation_cell == a],
                                            gsva_var2[[path]][gsva_var2$Segment == s & gsva_var2$Annotation_cell == a])$p.value, 
                              error = function(e) NA)
        
        log2fc_mean <- log2(mean(gsva_var1[[path]][gsva_var1$Segment == s & gsva_var1$Annotation_cell == a], na.rm = T) /
                              mean(gsva_var2[[path]][gsva_var2$Segment == s & gsva_var2$Annotation_cell == a], na.rm = T))
        
        mean_diff <- mean(gsva_var1[[path]][gsva_var1$Segment == s & gsva_var1$Annotation_cell == a], na.rm = T) -
          mean(gsva_var2[[path]][gsva_var2$Segment == s & gsva_var2$Annotation_cell == a], na.rm = T)
        
        anno_df <- data.frame(pval = as.numeric(test_pval),
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
  
  #stats_sign_zone_df$sign <- gsub(sign_colname, '', stats_sign_zone_df$sign)
  
  stats_signif_df <- filter(stats_all_df, pval <= 0.05)
  
  library(colorspace)
  
  plot <- ggplot(stats_signif_df, aes(x = anno, y = path)) +
    geom_point(aes(color = mean_diff, size = pval_for_plot)) +
    theme_classic() +
    xlab(NULL) +
    ylab(NULL) +
    scale_size_area(
      "pval",
      trans = "log10",
      #max_size = ifelse(length(stats_sign_zone_noinf) < 10, 2.5, 2),
      max_size = 1,
      breaks = c(1e-10, 1e-5, 1e-1, 0.05),
      limits = c(1e-10, 0.05)
    ) +
    theme(
      axis.text = element_text(size = rel(0.3)),
      axis.text.x = element_text(angle = 45, hjust = 1),
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
  ggsave(file.path(output_dir, outp2, paste0('dotplot_', var_name, '_meandiff_signif2.png')),
         width=2000, height = 1000, unit='px')
  
}

