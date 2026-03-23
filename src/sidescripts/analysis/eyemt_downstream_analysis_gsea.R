library(data.table)
library(plyr)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(viridis)
library(scales)
library(tidyr)
library(circlize)
library(paletteer)
library(car)
library(ggpubr)

#TODO !! 6 AOIs from batch3TLS labels has missing GSEA scores. why?
# no_gsea <- c("DSP-1001660037683-E-H11.dcc", "DSP-1001660037684-F-A08.dcc", "DSP-1001660037685-G-F12.dcc",
#              "DSP-1001660039811-A-F02.dcc", "DSP-1001660039813-C-A06.dcc", "DSP-1001660039813-C-C01.dcc")

# TODO scale the signal in all expr (non-deconv) by the ct frequency (expr*(1-ctfreq)) ?
# TODO for deconv select only pathways specific for given ct, make 1 big hmap
# TODO as above, hmaps per eg CD8 vs Macro pathways 

# set up variables --------------------------------------------------------

ct_of_interest <- c("tumor", "Macrophages_Monocytes", "Tcells_CD8", "Tcells_CD4", "DCs", "Bcells", "Fibroblasts_Mesothelial")
metadt_cols <- c('dcc_filename','Segment', 'Roi_geomx', 'Segment_geomx', 'Sample', 'Patient', 
                 'NACT_status', 'HRP_status', 'BRCA_status', 'PFS_quartile_b123', 'OS_quartile_b123')


# label_cols <- c('Annotation_cell','network_hub_type', 'community_cluster_label', 'network_hub_type_freq0.05',
#                 'community_cluster_label_freq0.05')

proj_dir <<- '~/Documents/phd/st'
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808')
out_dir <- file.path(output_dir, 'downstream', 'gsea_immune') 

#input files

# metadata
metadata_orig_path <- file.path(proj_dir, 'data/geomx/batch123/metadata/dcc_metadata_batch123_no_tls_cleaned.csv')

# ssgsea_scores dir for msigdb and additional
gsea_out_dir <- file.path(output_dir, 'pathway_analysis', 'gsea')

# names of pathways of interest
sigs_path <- file.path(proj_dir, 'geomx-processing', 'data', 'signatures', 'eyemt_immune', 'immune_pathways_names.csv')

# merged cell fractions and roi labels from communities (from geomx_roi_hubs_integration.R)
#ct_frac_all_roi_path <- file.path(output_dir, 'cycif_integration', paste0(batch_name, '_ct_frac_all_roi_bcells.csv'))

# ct fractions of immune per roi (used for clustering) from geomx_relabel_roi_2nd_approach
ct_frac_deconv_roi_path <- file.path(output_dir, 'deconvolution', 'relabel-roi-deconv-dimred', 'sd_mye_lymph_b_ct_fractions_of_immune.csv')

# ROI clusters based on deconvolution ct fractions from geomx_relabel_roi_2nd_approach
ct_frac_clust_path <- file.path(output_dir, 'deconvolution', 'relabel-roi-deconv-dimred', 'sd_mye_lymph_b_all_clustering_results.csv')
clust_types <- c('clusters_gmm_clustnr_5', 'clusters_hclust_cut2')
# descriptive labels for clusters - IN THIS CASE BOTH METHODS HAS THE SAME CLUSTERS DESCRIPTION
clust_labels <- list(CD8_Macro_domin = 1, mixed_w_CD4 = 2, mixed_w_others = 3, Macro_domin = 4, Bcell_domin = 5)

# source and create output dir
source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))
dir.create(out_dir, recursive = T, showWarnings = F)

# load metadata and merge with labels -------------------------------------

metadt <- fread(metadata_orig_path, select = metadt_cols)
metadt$sample_roi <- paste0(metadt$Sample, '_', metadt$Roi_geomx)

ct_frac_immune <- fread(ct_frac_deconv_roi_path)
colnames(ct_frac_immune)[-1] <- paste0('roiimmunefrac_', colnames(ct_frac_immune)[-1])

metadt <- left_join(metadt, ct_frac_immune)

###########################################
# extract communities roi labels and clean rois wo enough cells
# ct_frac_roi <- fread(ct_frac_all_roi_path)
# roi_labels <- ct_frac_roi %>%
#   select(sample_roi, total_cell_nr_cycif, !!label_cols) %>%
#   distinct() %>%
#   filter(total_cell_nr_cycif >= 100)
# 
# # TODO examine for other batches - for now simple cleaning
# roi_labels$community_cluster_label_freq0.05 <- ifelse(roi_labels$community_cluster_label_freq0.05 == 'CD11_Iba1|Iba1', 'CD11_Iba1',
#                                                       roi_labels$community_cluster_label_freq0.05)
# 
# # TODO here only ROIs with any ct phenotyped remains
# # change to do all ROIs clustering wo labels
# # merge ROI labels with metadata - labels per AOI
# metadt_labels_comm <- left_join(metadt, roi_labels, by = 'sample_roi') %>%
#   filter(sample_roi %in% roi_labels$sample_roi)

##########################################
##########################################
# alternative roi labels from deconv clustering
deconv_ct_frac_clust <- fread(ct_frac_clust_path, select = c('sample_roi', clust_types))
metadt_labels_deconv <- as.data.frame(left_join(metadt, deconv_ct_frac_clust, by = 'sample_roi'))
# TODO may need adjustment
metadt_labels_deconv[[paste0(clust_types[1], '_label')]] <- mapvalues(metadt_labels_deconv[[clust_types[1]]], 
                                                                      from=c(unname(unlist(clust_labels))),
                                                                      to=c(names(clust_labels)))
metadt_labels_deconv[[paste0(clust_types[2], '_label')]] <- mapvalues(metadt_labels_deconv[[clust_types[2]]], 
                                                                      from=c(unname(unlist(clust_labels))),
                                                                      to=c(names(clust_labels)))

# write for dge
# TODO move to deconv labelling
# fwrite(metadt_labels_deconv[, c('dcc_filename', paste0(clust_types, '_label'))], 
#        file.path(output_dir, 'deconvolution', 'relabel-roi-deconv', 'dcc_deconv_clusters.csv'))

# which one to use
metadt_labels <- metadt_labels_deconv

# load and filter gsea signatures -----------------------------------------

# load all gsea results filter to signatures and cells of interest, merge with sigs
sigs <- fread(sigs_path)

gsea_all <- lapply(list.files(gsea_out_dir, pattern = 'csv', full.names = T), function(x){
  gsea <- fread(x, select = c('dcc_filename', 'pathway', 'ssgsea_score', 'expr_signal'))
  gsea <- gsea[gsea$pathway %in% sigs$pathway, ]
})

gsea_all <- do.call(rbind, gsea_all)
gsea_all <- gsea_all[gsea_all$expr_signal %in% c('all', paste0('deconv_', ct_of_interest)), ] # filter to deconv ct
gsea_all <- left_join(gsea_all, sigs)


# create paths and dcc annotations ----------------------------------------

path_annots <- c('cell_type', 'immune_effect', 'additional')
path_annots <- c('cell_type', 'immune_effect')
# dcc_annots <- c('Segment', 'Segment_geomx', 'HRP_status', 'PFS_quartile_b123',
#                 'network_hub_type_freq0.05', 'community_cluster_label_freq0.05')
#dcc_annots <- c('Segment', 'HRP_status', 'PFS_quartile_b123', 'network_hub_type_freq0.05', 'community_cluster_label_freq0.05')
#dcc_annots <- c('Segment', 'Segment_geomx', 'HRP_status', 'PFS_quartile_b123', paste0(clust_types, '_label'))
dcc_annots <- c('Segment', 'Segment_geomx', 'HRP_status', 'NACT_status', 'PFS_quartile_b123', paste0(clust_types[1], '_label'))

# row annotations based on pathways
path_annot <- gsea_all %>%
  select('pathway', !!path_annots) %>%
  distinct() %>%
  column_to_rownames(var="pathway") %>%
  mutate(cell_type=replace(cell_type, cell_type=='Tcells_other', 'Tcells')) %>%
  mutate(cell_type=replace(cell_type, cell_type=='immune_other', 'immune')) %>%
  mutate_all(as.factor)

# column annotations based on metadata
dcc_annot <- metadt_labels %>%
  select('dcc_filename', !!dcc_annots) %>%
  column_to_rownames(var="dcc_filename") %>%
  mutate_all(as.factor)


# hmaps with pathways activity across all ROIs ----------------------------

top_var_nr <- 50 # nr of top variable pathways for clustering, NULL for all pathways

#TODO run for tum/str pre/post separately
expr_type <- unique(gsea_all$expr_signal)[1]
pathways_type <- unique(gsea_all$path_type)[1]

for(expr_type in unique(gsea_all$expr_signal)){
  for(pathways_type in unique(gsea_all$path_type)){
    
    out_name <- paste0(expr_type, '_', pathways_type)
    print(out_name)
    
    # filter to all signal, cell state/process pathways and ROIs with labels
    gsea_sel_paths <- gsea_all[gsea_all$expr_signal == expr_type & gsea_all$path_type == pathways_type, ]
    gsea_sel_paths <- gsea_sel_paths[gsea_sel_paths$dcc_filename %in% metadt_labels$dcc_filename, ]
    
    # create matrix with values for each pathway
    gsea_sel_wide <- spread(gsea_sel_paths[, c('dcc_filename','pathway', 'ssgsea_score')],
                                   key = 'pathway', value = 'ssgsea_score') %>%
      column_to_rownames(var="dcc_filename")
    
    # filter to pathways with top variance
    if(!is.null(top_var_nr )){
      gsea_variance <- sapply(gsea_sel_wide, var) 
      top_var_names <- names(sort(gsea_variance, decreasing = T)[1:top_var_nr])
      gsea_sel_wide <- gsea_sel_wide[, top_var_names]
      
      out_name <- paste0(out_name, '_topvar_', top_var_nr)
    }
    
    gsea_sel_wide <- t(as.matrix(gsea_sel_wide))
    
    ########################################333
    # TODO do sth with stupid annot colors
    # # set up colors
    # pal12 <- paletteer_d("colorBlindness::PairedColor12Steps")
    # pal4 <- c('red', 'blue', 'green', 'yellow')
    # pal10 <- paletteer_d("ggsci::category10_d3")
    # pal50 <- paletteer_d("ggsci::default_igv")
    # pal10 <- c('#1F77B4FF', '#FF7F0EFF','#2CA02CFF','#D62728FF', '#9467BDFF','#8C564BFF',
    #            '#E377C2FF', '#7F7F7FFF', '#BCBD22FF', '#17BECFFF')
    # 
    # annotation_colors <- setNames(sample(pal10, length(unique(dcc_annot$sd_mye_lymph_b_hcut2_label))), unique(dcc_annot$sd_mye_lymph_b_hcut2_label))
    # # 
    # dcc_annot_cols <- list(pal4, pal4, pal4, pal4, pal10)
    # 
    # 
    # dcc_annot_colors <- list()
    # for(annot in colnames(dcc_annot)){
    #   annot_colors <- setNames(pal10[1:length(levels(dcc_annot[[annot]]))], levels(dcc_annot[[annot]]))
    #   dcc_annot_colors <- append(dcc_annot_colors, list(annot_colors))
    # }
    # names(dcc_annot_colors) <- colnames(dcc_annot)
    # 
    # path_annot_cols <- list(pal10, pal4)
    # names(path_annot_cols) <- colnames(path_annot)
    # deconv_label_col <- c(pal10[1:length(unique(dcc_annot$sd_mye_lymph_b_hcut2_label))])
    # names(deconv_label_col) <- unique(dcc_annot$sd_mye_lymph_b_hcut2_label)
    
    ###################################3
    
    # filter annots and ensure ordering
    path_annot_toplot <- path_annot[match(rownames(gsea_sel_wide), rownames(path_annot)),] 
    dcc_annot_toplot <- dcc_annot[match(colnames(gsea_sel_wide), rownames(dcc_annot)),] 
    
    # set up annotations
    dcc_ha = HeatmapAnnotation(df = dcc_annot_toplot, which = 'column', na_col = "grey")
    path_ha = HeatmapAnnotation(df = path_annot_toplot, which = 'row', na_col = "grey")
    
    # define colors
    col_fun <- colorRamp2(c(min(gsea_sel_wide), 0, max(gsea_sel_wide)), c("blue", "white", "red"))
    
    # do the hmap
    png(filename=file.path(out_dir, paste0('hmap_',out_name, '_labs_deconv.png')), width=11, height=7,units="in",res=1000)
    
    gsea_heat <- Heatmap(gsea_sel_wide, name = "ssGSEA scorres of immune signalling pathways", 
                         top_annotation = dcc_ha, 
                         left_annotation = path_ha,
                         heatmap_legend_param = list(
                           legend_direction = "horizontal", 
                           legend_width = unit(2, "in")),
                         show_column_names = F,
                         row_names_gp = gpar(fontsize = 4),
                         row_names_max_width = unit(5, "in")
    )
    
    draw(gsea_heat, 
         heatmap_legend_side="bottom", 
         annotation_legend_side="bottom",
         merge_legend = TRUE)
    dev.off()
  }
}


# hmaps with correlations between pathways --------------------------------
corr_thr <- 0.3

#TODO mostly copypasted - add to previous loop
expr_type <- unique(gsea_all$expr_signal)[5]
pathways_type <- unique(gsea_all$path_type)[1]
seg <- 'stroma'
nact_status <- 'post'

# TODO for deconv signal, filter only to pathways relevant for given ct

for(expr_type in unique(gsea_all$expr_signal)){
  for(pathways_type in unique(gsea_all$path_type)){
    for(seg in c('stroma', 'tumor')){
      for(nact_status in c('pre', 'post')){
        
        out_name <- paste(expr_type, pathways_type, seg, nact_status, sep ='_')
        print(out_name)
        
        # filter to all signal, cell state/process, segment and nact status
        gsea_sel_paths <- gsea_all[gsea_all$expr_signal == expr_type & gsea_all$path_type == pathways_type, ]
        gsea_sel_paths_seg <- gsea_sel_paths[gsea_sel_paths$dcc_filename %in% 
                                               metadt_labels$dcc_filename[metadt_labels$Segment == seg & metadt_labels$NACT_status == nact_status], ]
        
        # create matrix with values for each pathway
        gsea_sel_wide_seg <- spread(gsea_sel_paths_seg[, c('dcc_filename','pathway', 'ssgsea_score')],
                                    key = 'pathway', value = 'ssgsea_score') %>%
          column_to_rownames(var="dcc_filename")
        
        # calculate wide df with immunefractions
        immunefrac <- metadt_labels %>%
          filter(dcc_filename %in% rownames(gsea_sel_wide_seg)) %>%
          dplyr::select('dcc_filename', starts_with('roiimmunefrac')) %>%
          column_to_rownames('dcc_filename')
        
        ####################################
        # correlation between pathways
        gsea_seg_corr <- cor(gsea_sel_wide_seg, method = 'pearson')
        
        # correlation between pathways and immune ct fractions
        gsea_seg_corr_ctfrac <- cor(gsea_sel_wide_seg, immunefrac, method = 'spearman')
        
        if(!is.null(corr_thr)){
          gsea_seg_corr[gsea_seg_corr >= -corr_thr & gsea_seg_corr <= corr_thr] <- 0
          cor_abovethr <- (colSums(gsea_seg_corr, na.rm=T) != 1) 
          gsea_seg_corr <- gsea_seg_corr[cor_abovethr, cor_abovethr]
          
          gsea_seg_corr_ctfrac[gsea_seg_corr_ctfrac >= -corr_thr & gsea_seg_corr_ctfrac <= corr_thr] <- 0
          cor_abovethr_ctfrac <- (rowSums(gsea_seg_corr_ctfrac, na.rm=T) != 0) 
          gsea_seg_corr_ctfrac <- gsea_seg_corr_ctfrac[cor_abovethr_ctfrac, ]
          
          out_name <- paste0(out_name, '_filt', as.character(corr_thr))
        }
        
        gsea_corr_list <- list(corr_path = gsea_seg_corr, corr_path_ctfrac = gsea_seg_corr_ctfrac)
        
        #######################################
        
        for(i in 1:length(gsea_corr_list)){
          
          gsea_corr <- as.matrix(gsea_corr_list[[i]])
          
          if(nrow(gsea_corr) > 0){
            path_annot_toplot <- path_annot[match(rownames(gsea_corr), rownames(path_annot)),] # ensure ordering
            
            # set up annotation
            path_ha = HeatmapAnnotation(df = path_annot_toplot, which = 'row', na_col = "grey")
            
            # define colors
            col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
            
            # do the hmap
            png(filename=file.path(out_dir, paste0('hmap_', names(gsea_corr_list)[i], '_', out_name, '.png')), width=11, height=9,units="in",res=1000)
            
            gsea_corr_heat <- Heatmap(gsea_corr, name = paste0("ssGSEA correlations"), 
                                      left_annotation = path_ha,
                                      heatmap_legend_param = list(
                                        legend_direction = "horizontal", 
                                        legend_width = unit(2, "in")),
                                      show_column_names = T,
                                      row_names_gp = gpar(fontsize = 4),
                                      column_names_gp = gpar(fontsize = 4),
                                      row_names_max_width = unit(5, "in"),
                                      col = col_fun
            )
            
            draw(gsea_corr_heat, 
                 heatmap_legend_side="top", 
                 annotation_legend_side="top",
                 merge_legend = TRUE)
            dev.off()
          }
          
        }
      }
    }
  }
}


# hmaps with mean GSEA score per ctfrac cluster ---------------------------
# TODO hmap with mean gsea scores + cluster
# TODO boxplots with wilcox ? (or maybe parametric since scores are normalised)
# TODO hmap/dotplot with gsea scores + wilcox as stars

top_var_nr <- 20
clust_type <- paste0(clust_types[1], '_label')

# TODO repetition - put all plots to main loop
expr_type <- unique(gsea_all$expr_signal)[1]
pathways_type <- unique(gsea_all$path_type)[1]
seg <- 'stroma'
nact_status <- 'post'

for(expr_type in unique(gsea_all$expr_signal)){
  for(pathways_type in unique(gsea_all$path_type)){
    for(seg in c('stroma', 'tumor')){
      for(nact_status in c('pre', 'post')){
        
        out_name <- paste(expr_type, pathways_type, seg, nact_status, sep ='_')
        print(out_name)
        
        # filter to all signal, cell state/process, segment and nact status
        gsea_sel_paths <- gsea_all[gsea_all$expr_signal == expr_type & gsea_all$path_type == pathways_type, ]
        gsea_sel_paths_seg <- gsea_sel_paths[gsea_sel_paths$dcc_filename %in% 
                                               metadt_labels$dcc_filename[metadt_labels$Segment == seg & metadt_labels$NACT_status == nact_status], ]
        
        gsea_sel_paths_seg <- left_join(gsea_sel_paths_seg, metadt_labels_deconv[, c('dcc_filename', clust_type, 'Sample')])
        
        # calculate mean ssgsea score per cluster
        gsea_sel_clust_mean <- gsea_sel_paths_seg %>%
          group_by(!!!syms(clust_type), pathway) %>%
          summarise(mean_ssgsea_score = mean(ssgsea_score))
        
        # transfrom to wide and create mtx
        gsea_sel_clust_mean_wide <- spread(gsea_sel_clust_mean[, c(clust_type,'pathway', 'mean_ssgsea_score')],
                                           key = 'pathway', value = 'mean_ssgsea_score') %>%
          column_to_rownames(var=clust_type) %>%
          t()
        
        # filter to pathways with top variance
        # if(!is.null(top_var_nr )){
        #   gsea_variance <- sapply(as.data.frame(t(gsea_sel_clust_mean_wide)), var) 
        #   top_var_names <- names(sort(gsea_variance, decreasing = T)[1:top_var_nr])
        #   gsea_sel_clust_mean_wide <- gsea_sel_clust_mean_wide[top_var_names, ] #rows
        #   
        #   out_name <- paste0(out_name, '_topvar_', top_var_nr)
        # }
        
        # filter to pathways with significant ancova score
        path_aov_pval <- lapply(unique(gsea_sel_paths_seg$pathway), function(path_name){
          gsea_clust_path <- gsea_sel_paths_seg[gsea_sel_paths_seg$pathway == path_name, ]
          
          ancova_model <- aov(gsea_clust_path$ssgsea_score ~ factor(gsea_clust_path[[clust_type]]) + factor(gsea_clust_path$Sample))
          ancova_res <- Anova(ancova_model, type="III") 
          aov_pval_sc <- ancova_res["factor(gsea_clust_path[[clust_type]])", "Pr(>F)"]

          path_df <- data.frame(pathway = path_name, aov_pval = aov_pval_sc)
        })
        
        path_aov_pval <- do.call(rbind, path_aov_pval)
        path_aov_pval <- path_aov_pval[path_aov_pval$aov_pval <= 0.05, ]
        
        # filter to pathways with significant differences
        gsea_sel_clust_mean_wide <- gsea_sel_clust_mean_wide[rownames(gsea_sel_clust_mean_wide) %in% path_aov_pval$pathway, ]
        
        # filter annots and ensure ordering
        path_annot_toplot <- path_annot[match(rownames(gsea_sel_clust_mean_wide), rownames(path_annot)),] 
        
        # set up annotations
        path_ha = HeatmapAnnotation(df = path_annot_toplot, which = 'row', na_col = "grey")
        
        # define colors
        # TODO maybe min/max(gsea_all) to have the same scale over all hmaps?
        col_fun <- colorRamp2(c(min(gsea_sel_clust_mean_wide), 0, max(gsea_sel_clust_mean_wide)), c("blue", "white", "red"))
        
        # do the hmap
        png(filename=file.path(out_dir, paste0('hmap_mean_clust_',out_name, '.png')), width=11, height=7,units="in",res=1000)
        
        gsea_heat <- Heatmap(gsea_sel_clust_mean_wide, name = "ssGSEA scores", 
                             left_annotation = path_ha,
                             col = col_fun,
                             cluster_columns = F,
                             heatmap_legend_param = list(
                               legend_direction = "horizontal", 
                               legend_width = unit(2, "in")),
                             row_names_gp = gpar(fontsize = 4),
                             row_names_max_width = unit(5, "in")
        )
        
        draw(gsea_heat, 
             heatmap_legend_side="bottom", 
             annotation_legend_side="bottom",
             merge_legend = TRUE)
        dev.off()
        
      }
    }
  }
}


##########################
# compute anova per each pathway

# TODO above hmap with means, but restricted to signif pathways instead of variance

# TODO repetition - put all plots to main loop
expr_type <- unique(gsea_all$expr_signal)[1]
pathways_type <- unique(gsea_all$path_type)[2]
seg <- 'stroma'
nact_status <- 'post'

out_name <- paste(expr_type, pathways_type, seg, nact_status, sep ='_')
print(out_name)

# filter to all signal, cell state/process, segment and nact status
gsea_sel_paths <- gsea_all[gsea_all$expr_signal == expr_type & gsea_all$path_type == pathways_type, ]
gsea_sel_paths_seg <- gsea_sel_paths[gsea_sel_paths$dcc_filename %in% 
                                       metadt_labels$dcc_filename[metadt_labels$Segment == seg & metadt_labels$NACT_status == nact_status], ]

gsea_sel_paths_clust <- gsea_sel_paths_seg %>%
  left_join(metadt_labels_deconv[, c('dcc_filename', clust_type, 'Sample')])

# do levene Test to check if variance is similar
leveneTest(gsea_sel_paths_clust$ssgsea_score ~ factor(gsea_sel_paths_clust[[clust_type]]))

# iterate through pathways and compute aov pval

path_aov_pval <- lapply(unique(gsea_sel_paths_clust$pathway), function(path_name){
  gsea_clust_path <- gsea_sel_paths_clust[gsea_sel_paths_clust$pathway == path_name, ]
  
  # do anova per cluster 
  # gsea_clust_aov <- aov(gsea_clust_path$ssgsea_score ~ factor(gsea_clust_path[[clust_type]]))
  # summary(gsea_clust_aov)
  # aov_pval_sc <- summary(gsea_clust_aov)[[1]][["Pr(>F)"]][1]
  
  # fit ANCOVA model with Sample as covariate
  # https://www.statology.org/ancova-in-r/
  ancova_model <- aov(gsea_clust_path$ssgsea_score ~ factor(gsea_clust_path[[clust_type]]) + factor(gsea_clust_path$Sample))
  #summary(ancova_model)
  ancova_res <- Anova(ancova_model, type="III") 
  
  # view summary of model
  aov_pval_sc <- ancova_res["factor(gsea_clust_path[[clust_type]])", "Pr(>F)"]
  
  #perform Tukey post-hoc test
  #TukeyHSD(gsea_clust_aov)

  path_df <- data.frame(pathway = path_name, aov_pval = aov_pval_sc)
})

path_aov_pval <- do.call(rbind, path_aov_pval)
path_aov_pval <- path_aov_pval[path_aov_pval$aov_pval <= 0.05, ]

# filter to pathways with significant differences
gsea_sel_paths_clust <- gsea_sel_paths_clust[gsea_sel_paths_clust$pathway %in% path_aov_pval$pathway]

for(path_name in unique(gsea_sel_paths_clust$pathway)){
  print(path_name)
  path_df <- gsea_sel_paths_clust[gsea_sel_paths_clust$pathway == path_name, ]
  
  ggplot(path_df, aes(x = factor(get(clust_type)), y = ssgsea_score, fill = factor(get(clust_type)))) +
    geom_boxplot(color = "black", alpha = 0.7) +
    labs(title = path_name, x = "cluster", y = "ssGSEA score") +
    theme_minimal() +
    theme(legend.position = "top") + 
    geom_pwc(method = "tukey_hsd", label = "p.signif", hide.ns = TRUE, size = 0.2, label.size = 2.8)
    #facet_wrap(~pathway, scales = "fixed", dir="v")
  
  ggsave(file.path(out_dir, paste0('boxpl_clusters_stroma_post_all_signal_', path_name, '.png')))
}



# compare the models
# library(AICcmodavg)
# model.set <- list(mtcars_aov, mtcars_aov2)
# model.names <- c("mtcars_aov", "mtcars_aov2")
# aictab(model.set, modnames = model.names)

# ANCOVA tutorial

