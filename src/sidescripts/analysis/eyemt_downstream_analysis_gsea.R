library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ComplexHeatmap)
library(viridis)
library(scales)
library(tidyr)
library(circlize)
library(paletteer)
library(car)
library(ggpubr)
library(GeomxTools)

#TODO !! 6 AOIs from batch3TLS labels has missing GSEA scores. why?
# no_gsea <- c("DSP-1001660037683-E-H11.dcc", "DSP-1001660037684-F-A08.dcc", "DSP-1001660037685-G-F12.dcc",
#              "DSP-1001660039811-A-F02.dcc", "DSP-1001660039813-C-A06.dcc", "DSP-1001660039813-C-C01.dcc")

# TODO scale the signal in all expr (non-deconv) by the ct frequency (expr*(1-ctfreq)) ? - maybe not the best idea

# TODO clean the code and put everything to 1 loop + functions

# TODO corr vs ct number - make one table and analyse systematically

# set up variables --------------------------------------------------------

#ct_of_interest <- c("tumor", "Macrophages_Monocytes", "Tcells_CD8", "Tcells_CD4", "DCs", "Bcells")
ct_of_interest <- c("Macrophages_Monocytes")
ct_names_immune <- c("Tcells_CD4", "Tcells_CD8", "DCs", "Macrophages_Monocytes", "Bcells", "NKcells")
ct_names_myeloids <- c("DCs", "Macrophages_Monocytes")
ct_names_lymphoids <- c("Tcells_CD4", "Tcells_CD8", "Bcells")

min_frac <- 0.01 # gsea scores computed for dcc with smaller fraction, will be removed

proj_dir <<- '~/Documents/phd/st'
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808')
out_dir <- file.path(output_dir, 'downstream', 'gsea_immune_Macro') 

#input files
geomx_path <- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS')

# metadata compiled in metadata_join.R
metadt_path <- file.path(output_dir, 'metadata_full_SENSITIVE.csv')

# ssgsea_scores dir for msigdb and additional
gsea_out_dir <- file.path(output_dir, 'pathway_analysis', 'gsea')

# names of pathways of interest
sigs_path <- file.path(proj_dir, 'geomx-processing', 'data', 'signatures', 'eyemt_immune', 'immune_pathways_names.csv')
# RDS object with list of pathways and their genes
signs_genes_list_path <- file.path(proj_dir, 'geomx-processing', 'data', 'signatures', 'eyemt_immune', 'sign_list_immune.RDS')

# source and create output dir
source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))
dir.create(out_dir, recursive = T, showWarnings = F)


# load metadata and filter to dcc after qc --------------------------------

geomx_dcc <- colnames(readRDS(geomx_path))
metadt <- as.data.frame(fread(metadt_path))
metadt <- metadt[metadt$dcc_filename %in% geomx_dcc, ]

# load and filter gsea signatures -----------------------------------------

# load signatures and remove the ones which are subset of others
# check which pathways are subsets of other pathways - choose longest
sigs <- fread(sigs_path)
sign_genes_list <- readRDS(signs_genes_list_path)

unique_paths <- unlist(sapply(1:length(sign_genes_list), function(i){
  sign <- sign_genes_list[[i]]
  sign_inters <- sapply(sign_genes_list[-i], function(sign2){
    length(intersect(sign, sign2))
  })
  
  if(max(sign_inters) < length(sign)){
    return(names(sign_genes_list)[i])
  } else{
    return()
  }
}))

sigs <- sigs[sigs$pathway %in% unique_paths, ]

# load all gsea results filter to signatures and cells of interest, merge with sigs
gsea_all <- lapply(list.files(gsea_out_dir, pattern = 'csv', full.names = T), function(x){
  gsea <- fread(x, select = c('dcc_filename', 'pathway', 'ssgsea_score', 'expr_signal'))
  gsea <- gsea[gsea$pathway %in% sigs$pathway, ]
})

gsea_all <- do.call(rbind, gsea_all)
gsea_all <- gsea_all[gsea_all$expr_signal %in% c('all', paste0('deconv_', ct_of_interest)), ] # filter to deconv ct
gsea_all <- left_join(gsea_all, sigs)

# for deconvoluted signal, filter to pathways expressed by given ct
# filter to dcc with high enough fraction of given cell
gsea_deconv_sel <- lapply(ct_of_interest, function(ct){
  
  ct_groups <- c('all', grep(ct, unique(sigs$path_cell_type), value = T))
  if(ct %in% ct_names_immune){ct_groups <- c(ct_groups, 'immune')} 
  if(ct %in% ct_names_lymphoids){ct_groups <- c(ct_groups, 'lymphoids')} 
  if(ct %in% ct_names_myeloids){ct_groups <- c(ct_groups, 'myeloids')}
  
  gsea_ct <- gsea_all %>%
    filter(expr_signal == paste0('deconv_', ct)) %>%
    filter(path_cell_type %in% ct_groups) %>%
    filter(dcc_filename %in% metadt$dcc_filename[metadt[[paste0("ct_frac_sd_aoi_", ct)]] >= min_frac])
  
  return(gsea_ct)
})

gsea_all <- rbind(gsea_all[gsea_all$expr_signal == 'all', ], do.call(rbind, gsea_deconv_sel))

# create paths and dcc annotations ----------------------------------------

path_annots <- c('path_cell_type', 'immune_effect', 'additional')
path_annots <- c('path_cell_type', 'immune_effect')
dcc_annots <- c('Segment', 'Segment_geomx', 'HRP_status', 'NACT_status', 'roi_cluster_label_gmm')

# row annotations based on pathways
path_annot <- gsea_all %>%
  select('pathway', !!path_annots) %>%
  distinct() %>%
  column_to_rownames(var="pathway") %>%
  mutate_all(as.factor)

# column annotations based on metadata
dcc_annot <- metadt %>%
  select('dcc_filename', !!dcc_annots) %>%
  column_to_rownames(var="dcc_filename") %>%
  mutate_all(as.factor)


# hmaps with pathways activity across all ROIs ----------------------------

top_var_nr <- 50 # nr of top variable pathways for clustering, NULL for all pathways

################
# for testing
expr_type <- unique(gsea_all$expr_signal)[1]
pathways_type <- unique(gsea_all$path_type)[1]
seg <- 'stroma'
nact_status <- 'pre'
################

for(expr_type in unique(gsea_all$expr_signal)){
    for(seg in c('both', 'stroma', 'tumor')){
      for(nact_status in c('pre', 'post')){

        # out_name <- ifelse(seg == 'both', paste(expr_type, pathways_type, sep ='_'),
        #                    paste(expr_type, pathways_type, seg, nact_status, sep ='_'))
        
        out_name <- ifelse(seg == 'both', expr_type,
                           paste(expr_type, seg, nact_status, sep ='_'))
        print(out_name)
        
        # filter to expr signal type, cell state/process pathways, segment and nact status
        #gsea_sel_paths <- gsea_all[gsea_all$expr_signal == expr_type & gsea_all$path_type == pathways_type, ]
        gsea_sel_paths <- gsea_all[gsea_all$expr_signal == expr_type, ]
        if(seg != 'both'){
          gsea_sel_paths <- gsea_sel_paths[gsea_sel_paths$dcc_filename %in% 
                                                 metadt$dcc_filename[metadt$Segment == seg & metadt$NACT_status == nact_status], ]
        }
        
        if(length(unique(gsea_sel_paths$dcc_filename)) > 20){
          
          print(paste0(length(unique(gsea_sel_paths$dcc_filename)), ' AOI for comparison'))
          # create matrix with values for each pathway
          gsea_sel_wide <- spread(gsea_sel_paths[, c('dcc_filename','pathway', 'ssgsea_score')],
                                  key = 'pathway', value = 'ssgsea_score') %>%
            column_to_rownames(var="dcc_filename")
          
          # filter to pathways with top variance
          if(!is.null(top_var_nr )){
            if(top_var_nr < ncol(gsea_sel_wide)){
              gsea_variance <- sapply(gsea_sel_wide, var) 
              top_var_names <- names(sort(gsea_variance, decreasing = T)[1:top_var_nr])
              gsea_sel_wide <- gsea_sel_wide[, top_var_names]
              
              out_name <- paste0(out_name, '_topvar_', top_var_nr)
            }
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
          png(filename=file.path(out_dir, paste0('hmap_scores_',out_name, '.png')), width=11, height=7,units="in",res=1000)
          
          gsea_heat <- Heatmap(gsea_sel_wide, name = "ssGSEA scores of immune signalling pathways", 
                               top_annotation = dcc_ha, 
                               left_annotation = path_ha,
                               col = col_fun,
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
    }
}


# hmaps with correlations between pathways and pathways vs ct frac --------

corr_thr <- 0.3

#TODO mostly copypasted - add to previous loop
##################
# testing
expr_type <- unique(gsea_all$expr_signal)[1]
pathways_type <- unique(gsea_all$path_type)[1]
seg <- 'stroma'
nact_status <- 'post'
#################

for(expr_type in unique(gsea_all$expr_signal)){
    for(seg in c('both', 'stroma', 'tumor')){
      for(nact_status in c('pre', 'post')){
        
        #out_name <- paste(expr_type, pathways_type, seg, nact_status, sep ='_')
        out_name <- ifelse(seg == 'both', expr_type,
                           paste(expr_type, seg, nact_status, sep ='_'))
        print(out_name)
        
        # filter to all signal, cell state/process, segment and nact status
        #gsea_sel_paths <- gsea_all[gsea_all$expr_signal == expr_type & gsea_all$path_type == pathways_type, ]
        gsea_sel_paths <- gsea_all[gsea_all$expr_signal == expr_type, ]
        
        if(seg != 'both'){
          gsea_sel_paths <- gsea_sel_paths[gsea_sel_paths$dcc_filename %in% 
                                                 metadt$dcc_filename[metadt$Segment == seg & metadt$NACT_status == nact_status], ]
        }
        
        if(length(unique(gsea_sel_paths$dcc_filename)) > 20){
          
          # create matrix with values for each pathway
          gsea_sel_wide_seg <- spread(gsea_sel_paths[, c('dcc_filename','pathway', 'ssgsea_score')],
                                      key = 'pathway', value = 'ssgsea_score') %>%
            column_to_rownames(var="dcc_filename")
          
          # wide df with immunefractions
          roi_immunefrac <- metadt %>%
            filter(dcc_filename %in% rownames(gsea_sel_wide_seg)) %>%
            dplyr::select('dcc_filename', starts_with('ct_immunefrac_sd_roi')) %>%
            column_to_rownames('dcc_filename')
          
          aoi_ctfrac <- metadt %>%
            filter(dcc_filename %in% rownames(gsea_sel_wide_seg)) %>%
            dplyr::select('dcc_filename', starts_with('ct_frac_sd_aoi')) %>%
            column_to_rownames('dcc_filename')
          
          ####################################
          # correlation between pathways
          gsea_seg_corr <- cor(gsea_sel_wide_seg, method = 'pearson')
          
          # correlation between pathways and immune ct fractions
          gsea_seg_corr_roi_immunefrac <- cor(gsea_sel_wide_seg, roi_immunefrac, method = 'spearman')
          gsea_seg_corr_aoi_ctfrac <- cor(gsea_sel_wide_seg, aoi_ctfrac, method = 'spearman')
          
          gsea_seg_corr_roi_immunefrac[is.na(gsea_seg_corr_roi_immunefrac)] <- 0
          gsea_seg_corr_aoi_ctfrac[is.na(gsea_seg_corr_aoi_ctfrac)] <- 0
          
          if(!is.null(corr_thr)){
            gsea_seg_corr[gsea_seg_corr >= -corr_thr & gsea_seg_corr <= corr_thr] <- 0
            cor_abovethr <- (colSums(gsea_seg_corr, na.rm=T) != 1) # cor btw same pathway will be 1
            gsea_seg_corr <- gsea_seg_corr[cor_abovethr, cor_abovethr]
            
            # TODO rmv duplication
            gsea_seg_corr_roi_immunefrac[gsea_seg_corr_roi_immunefrac >= -corr_thr & gsea_seg_corr_roi_immunefrac <= corr_thr] <- 0
            cor_abovethr_roi_immunefrac_row <- (rowSums(gsea_seg_corr_roi_immunefrac, na.rm=T) != 0) 
            cor_abovethr_roi_immunefrac_col <- (colSums(gsea_seg_corr_roi_immunefrac, na.rm=T) != 0)
            gsea_seg_corr_roi_immunefrac <- as.matrix(gsea_seg_corr_roi_immunefrac[cor_abovethr_roi_immunefrac_row,  cor_abovethr_roi_immunefrac_col])
            
            gsea_seg_corr_aoi_ctfrac[gsea_seg_corr_aoi_ctfrac >= -corr_thr & gsea_seg_corr_aoi_ctfrac <= corr_thr] <- 0
            cor_abovethr_aoi_ctfrac_row <- (rowSums(gsea_seg_corr_aoi_ctfrac, na.rm=T) != 0) 
            cor_abovethr_aoi_ctfrac_col <- (colSums(gsea_seg_corr_aoi_ctfrac, na.rm=T) != 0) 
            gsea_seg_corr_aoi_ctfrac <- as.matrix(gsea_seg_corr_aoi_ctfrac[cor_abovethr_aoi_ctfrac_row, cor_abovethr_aoi_ctfrac_col])
            
            out_name <- paste0(out_name, '_filt', as.character(corr_thr))
          }
          
          gsea_corr_list <- list(corr_path = gsea_seg_corr, corr_path_roi_immunefrac = gsea_seg_corr_roi_immunefrac,
                                 corr_path_aoi_ctfrac = gsea_seg_corr_aoi_ctfrac)
          
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
                                        row_names_gp = gpar(fontsize = ifelse(nrow(gsea_corr) < 20, 6, 4)),
                                        column_names_gp = gpar(fontsize = ifelse(ncol(gsea_corr) < 20, 6, 4)),
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


# heatmaps with correlations between 2 cell types deconv pathways ---------

# TODO for deconv select only pathways specific for given ct, make 1 big hmap (not sure if needed)

################
# testing
seg <- 'tumor'
nact_status <- 'post'
deconv_ct1 <- 'Tcells_CD8'
deconv_ct2 <- 'Bcells'
#################

# corr between specific deconv
corr_thr <- 0.7

#for(pathways_type in unique(gsea_all$path_type)){
for(cells_comb in combn(ct_of_interest, 2, simplify = F)){
  for(seg in c('both', 'stroma', 'tumor')){
    for(nact_status in c('pre', 'post')){
      
      deconv_ct1 <- cells_comb[1]
      deconv_ct2 <- cells_comb[2]
      
      out_name <- ifelse(seg == 'both', paste('deconv', deconv_ct1, 'vs', deconv_ct2, sep ='_'),
                         paste('deconv', deconv_ct1, 'vs', deconv_ct2, seg, nact_status, sep ='_'))
      print(out_name)
      #out_name <- paste('deconv', deconv_ct1, 'vs', deconv_ct2, pathways_type, seg, nact_status, sep ='_')
      
      # filter to cell state/process, segment and nact status
      # gsea_sel_paths <- gsea_all[gsea_all$path_type == pathways_type, ]
      
      if(seg != 'both'){
        gsea_sel_paths <- gsea_all[gsea_all$dcc_filename %in% 
                                     metadt$dcc_filename[metadt$Segment == seg & metadt$NACT_status == nact_status], ]
      } else{
        gsea_sel_paths <- gsea_all
      }
      
      gsea_deconv_ct1 <- gsea_sel_paths[gsea_sel_paths$expr_signal == paste0('deconv_', deconv_ct1), ]
      gsea_deconv_ct2 <- gsea_sel_paths[gsea_sel_paths$expr_signal == paste0('deconv_', deconv_ct2), ]
      
      print(paste0("nr of shared AOIs:", as.character(length(intersect(gsea_deconv_ct1$dcc_filename, gsea_deconv_ct2$dcc_filename)))))
      
      if(length(intersect(gsea_deconv_ct1$dcc_filename, gsea_deconv_ct2$dcc_filename)) > 20){
        
        # create matrix with values for each pathway
        gsea_deconv_wide_ct1 <- spread(gsea_deconv_ct1[, c('dcc_filename','pathway', 'ssgsea_score')],
                                       key = 'pathway', value = 'ssgsea_score') %>%
          column_to_rownames(var="dcc_filename")
        colnames(gsea_deconv_wide_ct1) <- paste0(colnames(gsea_deconv_wide_ct1), '_', deconv_ct1)
        
        gsea_deconv_wide_ct2 <- spread(gsea_deconv_ct2[, c('dcc_filename','pathway', 'ssgsea_score')],
                                       key = 'pathway', value = 'ssgsea_score') %>%
          column_to_rownames(var="dcc_filename")
        colnames(gsea_deconv_wide_ct2) <- paste0(colnames(gsea_deconv_wide_ct2), '_', deconv_ct2)
        
        # filter to dcc common in both mtx and ensure ordering
        gsea_deconv_wide_ct1 <- gsea_deconv_wide_ct1[intersect(rownames(gsea_deconv_wide_ct1), rownames(gsea_deconv_wide_ct2)),]
        gsea_deconv_wide_ct2 <- gsea_deconv_wide_ct2[intersect(rownames(gsea_deconv_wide_ct1), rownames(gsea_deconv_wide_ct2)),]
        
        # correlation between pathways
        gsea_deconv_corr <- cor(gsea_deconv_wide_ct1, gsea_deconv_wide_ct2, method = 'pearson')
        
        print(paste0("max corr: ", as.character(max(gsea_deconv_corr))))
        
        # filter to values above thr
        if(!is.null(corr_thr)){
          gsea_deconv_corr[gsea_deconv_corr >= -corr_thr & gsea_deconv_corr <= corr_thr] <- 0
          cor_abovethr_col <- (colSums(gsea_deconv_corr, na.rm=T) != 0) 
          cor_abovethr_row <- (rowSums(gsea_deconv_corr, na.rm=T) != 0) 
          gsea_deconv_corr <- as.matrix(gsea_deconv_corr[cor_abovethr_row, cor_abovethr_col])
          # restore colnames - will throw error when 1col and 1row
          if(length(which(cor_abovethr_col)) == 1){colnames(gsea_deconv_corr) <- names(which(cor_abovethr_col))} 
          if(length(which(cor_abovethr_row)) == 1){colnames(gsea_deconv_corr) <- names(which(cor_abovethr_row))}
          
          out_name <- paste0(out_name, '_filt', as.character(corr_thr))
        }
        
        if(nrow(gsea_deconv_corr) > 0){
          # do hmap
          # TODO annots for rows + cols paths separately
          # path_annot_toplot <- path_annot[match(rownames(gsea_deconv_corr), rownames(path_annot)),] # ensure ordering
          # # set up annotation
          # path_ha = HeatmapAnnotation(df = path_annot_toplot, which = 'row', na_col = "grey")
          
          # define colors
          col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
          
          # do the hmap
          png(filename=file.path(out_dir, paste0('hmap_corr_path_btwcells_', out_name, '.png')), width=11, height=9,units="in",res=1000)
          
          gsea_corr_heat <- Heatmap(gsea_deconv_corr, name = paste0("ssGSEA correlations ", deconv_ct1, " vs ", deconv_ct2), 
                                    #left_annotation = path_ha,
                                    heatmap_legend_param = list(
                                      legend_direction = "horizontal", 
                                      legend_width = unit(2, "in")),
                                    show_column_names = T,
                                    row_names_gp = gpar(fontsize = 6),
                                    column_names_gp = gpar(fontsize = 6),
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
#}


# hmaps with mean GSEA score per ctfrac cluster ---------------------------
# TODO hmap with mean gsea scores + cluster
# TODO boxplots with tuckey test
# TODO hmap/dotplot with gsea scores + tuckey as stars - but to what since they're many groups?

top_var_nr <- 20
clust_type <- paste0(clust_types[1], '_label')

# TODO repetition - put all plots to main loop
expr_type <- unique(gsea_all$expr_signal)[1]
pathways_type <- unique(gsea_all$path_type)[1]
seg <- 'stroma'
nact_status <- 'post'

for(expr_type in unique(gsea_all$expr_signal)[7:8]){
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
          ancova_pval <- tryCatch({Anova(ancova_model, type="III")["factor(gsea_clust_path[[clust_type]])", "Pr(>F)"]},
                                 error = function(e){1})

          path_df <- data.frame(pathway = path_name, aov_pval = ancova_pval)
        })
        
        path_aov_pval <- do.call(rbind, path_aov_pval)
        path_aov_pval <- path_aov_pval[path_aov_pval$aov_pval <= 0.05, ]
        
        # filter to pathways with significant differences
        gsea_sel_clust_mean_wide <- gsea_sel_clust_mean_wide[rownames(gsea_sel_clust_mean_wide) %in% path_aov_pval$pathway, ]
        
        if(nrow(gsea_sel_clust_mean_wide) > 0){
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
}


##########################
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

################################################
################################################

# clusters distributions vs ct frac and clinical vars ---------------------

# each ct freq in tumor/stroma segment across nact status + freq cluster
#clust_type <- paste0(clust_types[1], '_label')
clust_type <- 'roi_cluster_label_hclust'
clust_type_name <- 'hclust' # for plotting

aoi_ctfrac_long <- metadt %>%
  dplyr::select('dcc_filename', starts_with('ct_frac_sd_aoi')) %>%
  pivot_longer(cols = starts_with('ct_frac_sd_aoi'), names_to = 'cell_type', values_to = 'ct_frac_sd_aoi') %>%
  mutate(cell_type = gsub('ct_frac_sd_aoi_', '', cell_type)) %>%
  left_join(metadt[, c('dcc_filename', 'Segment', 'NACT_status', clust_type)]) %>%
  mutate(segment_nact = paste0(Segment, '_', NACT_status)) %>%
  filter(cell_type %in% !!ct_names_immune)


# boxpl all cells at once, color by segment_nact
ggplot(aoi_ctfrac_long, aes(x = factor(cell_type), y = ct_frac_sd_aoi, fill = factor(segment_nact))) +
  geom_boxplot() +
  labs(title = 'ct freq aross segment and nact status', x = "cell type", y = "ct frac sd in AOI") +
  theme_minimal()

ggsave(file.path(out_dir, paste0('ct_frac_aoi_segment_nact.png')))

# boxpl faceted by cell, color by segment_nact all clusters at once
ggplot(aoi_ctfrac_long, aes(x = factor(get(clust_type)), y = ct_frac_sd_aoi, fill = factor(segment_nact))) +
  geom_boxplot() +
  labs(title = 'ct freq aross segment and nact status', x = "cell type", y = "ct frac sd in AOI") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  facet_wrap(~cell_type, scales = "fixed", dir="v")

ggsave(file.path(out_dir, paste0('ct_frac_aoi_per_cluster_', clust_type_name, '_segment_nact.png')))


###########################################
# ct frac clusters distribution across samples (pre, post, HRD, PFS, OS)

vars_labels <- c('NACT_status', 'HRP_status', 'primary_treatment_response')
vars_cont <- c('TMB', 'ovaHRDscar_score', 'PFS_days', 'OS_days')

# count total nr of clusters in dataset (for ordering)
cluster_labels_count <- metadt %>%
  select(dcc_filename, Sample, !!vars_labels, !!vars_cont, !!clust_type) %>%
  group_by(get(clust_type)) %>%   
  mutate(clust_name_occur_total = n()) %>%
  ungroup() %>%
  as.data.frame()

# all combinations
clust_allcombs <- tidyr::expand(metadt, Sample, get(clust_type))
colnames(clust_allcombs) <- c('Sample', clust_type) #fixing stupid names

# count ROI label frequency per sample (per ROI not AOIs!)
cluster_freqs_per_sample <- metadt[, c('Sample', 'sample_roi', clust_type)] %>%
  distinct() %>%
  group_by(Sample, get(clust_type)) %>%
  mutate(clust_nr_per_sample = n()) %>%
  ungroup() %>%
  group_by(Sample) %>%
  mutate(clust_freq_per_sample = clust_nr_per_sample/n()) %>%
  select(Sample, !!clust_type, clust_nr_per_sample, clust_freq_per_sample) %>%
  distinct() %>%
  full_join(clust_allcombs, by = c('Sample', clust_type)) %>% # join with all combs to get 0
  replace(is.na(.), 0) %>%
  left_join(distinct(metadt[, c('Sample', vars_labels, vars_cont)]))


# stacked barplot for nr of clusters across samples faceted by discrete vars
for(label_var in vars_labels){

  # nrs of AOIs from given ROI cluster 
  ggplot(cluster_labels_count, aes(x = reorder(Sample, clust_name_occur_total), fill = get(clust_type))) +
    geom_bar(stat = "count") +
    labs(title = paste0("nr of AOIs per sample across ", label_var), x = "Sample", y = "AOI nr", fill='ROI cluster type') +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 8)) +
    facet_wrap(~get(label_var), dir="v", scales = "free")
  
  ggsave(file.path(out_dir, paste0('aoi_nr_cluster_', clust_type_name, '_color_', label_var, '.png')))
  
  # frequencies of ROIs clusters - stacked barplots
  ggplot(cluster_freqs_per_sample, aes(x = Sample, y = clust_freq_per_sample, fill = get(clust_type))) +
    geom_bar(stat = "identity") +
    labs(title = paste0("frequencies of ROI clusters per sample across ", label_var), x = "Sample", y = "ROI cluster frequency", fill='ROI cluster type') +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 8)) +
    facet_wrap(~get(label_var), dir="v", scales = "free")
  
  ggsave(file.path(out_dir, paste0('roi_freq_cluster_', clust_type_name, '_color_', label_var, '.png')))
  
  # frequencies of ROIs clusters - boxplots
  ggplot(cluster_freqs_per_sample, aes(x = get(clust_type), y = clust_freq_per_sample, fill = get(label_var))) +
    geom_boxplot() +
    geom_point(position= position_jitterdodge(dodge.width = 1, jitter.width= .3, jitter.height = 0),
               size= 0.5, alpha = 0.6) +
    geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = TRUE, size = 0.2, label.size = 2.8) +
    labs(title = paste0("frequencies of ROI clusters per sample across ", label_var), x = "ROI cluster", y = "ROI cluster frequency", fill= paste0(label_var)) +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 8))
  
  ggsave(file.path(out_dir, paste0('roi_freq_cluster_', clust_type_name, '_boxpl_color_', label_var, '.png')))
}

# scatterplots for continuous vars
for(cont_var in vars_cont){
  
  ggplot(cluster_freqs_per_sample, aes(x = clust_freq_per_sample, y = get(cont_var), color = get(clust_type))) +
    geom_point(size = 3) +
    xlab("ROI cluster frequency in sample") +
    ylab(cont_var) +
    scale_color_discrete(name = clust_type) +
    theme_bw()
  
  ggsave(file.path(out_dir, paste0('roi_freq_cluster_', clust_type_name, '_scatter_color_', cont_var, '.png')))
}

#fwrite(metadt_labels_count, file.path(out_dir, 'metadt_labels.csv'))
