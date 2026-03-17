library(data.table)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(viridis)
library(scales)
library(tidyr)
library(circlize)
library(paletteer)

#TODO !! 6 AOIs from batch3TLS labels has missing GSEA scores. why?
# no_gsea <- c("DSP-1001660037683-E-H11.dcc", "DSP-1001660037684-F-A08.dcc", "DSP-1001660037685-G-F12.dcc",
#              "DSP-1001660039811-A-F02.dcc", "DSP-1001660039813-C-A06.dcc", "DSP-1001660039813-C-C01.dcc")

#TODO rescue 20 ROIs which 1 AOI was removed 
# 'b123_ct_frac_deconv_bcells.csv' vs 'b123_ct_frac_deconv_roi_bcells.csv' from geomx_roi_hubs_integration.R)

# set up variables --------------------------------------------------------
batch_name <- 'batch3tls'

ct_of_interest <- c("tumor", "Macrophages_Monocytes", "Tcells_CD8", "Tcells_CD4", "DCs", "Bcells", "Fibroblasts_Mesothelial")
metadt_cols <- c('dcc_filename','Segment', 'Roi_geomx', 'Segment_geomx', 'Sample', 'Patient', 
                 'NACT_status', 'HRP_status', 'BRCA_status', 'PFS_quartile_b123', 'OS_quartile_b123')
label_cols <- c('Annotation_cell','network_hub_type', 'community_cluster_label', 'network_hub_type_freq0.05',
                'community_cluster_label_freq0.05')

proj_dir <<- '~/Documents/phd/st'
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808')
out_dir <- file.path(output_dir, 'downstream', 'gsea_immune') 

#input files

# metadata
metadata_orig_path <- file.path(proj_dir, 'data/geomx/batch123/metadata/dcc_metadata_batch123_no_tls_cleaned.csv')
# merged cell fractions and roi labels (from geomx_roi_hubs_integration.R)
ct_frac_all_roi_path <- file.path(output_dir, 'cycif_integration', paste0(batch_name, '_ct_frac_all_roi_bcells.csv'))
# ssgsea_scores dir for msigdb and additional
gsea_out_dir <- file.path(output_dir, 'pathway_analysis', 'gsea')
# names of pathways of interest
sigs_path <- file.path(proj_dir, 'geomx-processing', 'data', 'signatures', 'eyemt_immune', 'immune_pathways_names.csv')
# ROI clusters based on ct fractions
ct_frac_clust_path <- file.path(output_dir, 'deconvolution', 'relabel-roi-deconv', 'roi_ctfreq_clusters_df_mid0.15_hi0.25_labs.csv')
clust_types <- c('sd_mye_lymph_b_myemerged_hcut2', 'sd_mye_lymph_b_hcut2', 'bp_mye_lymph_b_myemerged_hcut2', 'bp_mye_lymph_b_hcut2')

# source and create output dir
source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))
dir.create(out_dir, recursive = T, showWarnings = F)


# load metadata and merge with labels -------------------------------------

metadt <- fread(metadata_orig_path, select = metadt_cols)
metadt$sample_roi <- paste0(metadt$Sample, '_', metadt$Roi_geomx)

# extract communities roi labels and clean rois wo enough cells
ct_frac_roi <- fread(ct_frac_all_roi_path)
roi_labels <- ct_frac_roi %>%
  select(sample_roi, total_cell_nr_cycif, !!label_cols) %>%
  distinct() %>%
  filter(total_cell_nr_cycif >= 100)

# TODO examine for other batches - for now simple cleaning
roi_labels$community_cluster_label_freq0.05 <- ifelse(roi_labels$community_cluster_label_freq0.05 == 'CD11_Iba1|Iba1', 'CD11_Iba1',
                                                      roi_labels$community_cluster_label_freq0.05)

# TODO here only ROIs with any ct phenotyped remains
# change to do all ROIs clustering wo labels
# merge ROI labels with metadata - labels per AOI
metadt_labels_comm <- left_join(metadt, roi_labels, by = 'sample_roi') %>%
  filter(sample_roi %in% roi_labels$sample_roi)

##########################################
##########################################
# alternative roi labels from deconv clustering
deconv_ct_frac_clust <- fread(ct_frac_clust_path, select = c('sample_roi', clust_types, paste0(clust_types, '_label')))
metadt_labels_deconv <- as.data.frame(left_join(metadt, deconv_ct_frac_clust, by = 'sample_roi'))

# write for dge
# TODO move to deconv labelling
fwrite(metadt_labels_deconv[, c('dcc_filename', paste0(clust_types, '_label'))], 
       file.path(output_dir, 'deconvolution', 'relabel-roi-deconv', 'dcc_deconv_clusters.csv'))

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



# hmaps with pathways activity across all ROIs ----------------------------

path_annots <- c('cell_type', 'immune_effect', 'additional')
path_annots <- c('cell_type', 'immune_effect')
# dcc_annots <- c('Segment', 'Segment_geomx', 'HRP_status', 'PFS_quartile_b123',
#                 'network_hub_type_freq0.05', 'community_cluster_label_freq0.05')
#dcc_annots <- c('Segment', 'HRP_status', 'PFS_quartile_b123', 'network_hub_type_freq0.05', 'community_cluster_label_freq0.05')
#dcc_annots <- c('Segment', 'Segment_geomx', 'HRP_status', 'PFS_quartile_b123', paste0(clust_types, '_label'))
dcc_annots <- c('Segment', 'Segment_geomx', 'HRP_status', 'NACT_status', 'PFS_quartile_b123', 'sd_mye_lymph_b_hcut2_label')


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
    
    gsea_sel_wide <- t(as.matrix(gsea_sel_wide))
    
    # row annotations based on pathways
    path_annot <- gsea_sel_paths %>%
      select('pathway', !!path_annots) %>%
      distinct() %>%
      column_to_rownames(var="pathway") %>%
      mutate(cell_type=replace(cell_type, cell_type=='Tcells_other', 'Tcells')) %>%
      mutate(cell_type=replace(cell_type, cell_type=='immune_other', 'immune')) %>%
      mutate_all(as.factor)
    
    path_annot <- path_annot[match(rownames(gsea_sel_wide), rownames(path_annot)),] # ensure ordering
    
    # column annotations based on metadata
    dcc_annot <- metadt_labels %>%
      select('dcc_filename', !!dcc_annots) %>%
      filter(dcc_filename %in% gsea_sel_paths$dcc_filename) %>% # filter back bcs some labelled ROIs has no GSEA
      column_to_rownames(var="dcc_filename") %>%
      mutate_all(as.factor)
    
    dcc_annot <- dcc_annot[match(colnames(gsea_sel_wide), rownames(dcc_annot)),] # ensure ordering
    
    ########################################333
    # TODO do sth with stupid colors
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
    
    # set up annotations
    dcc_ha = HeatmapAnnotation(df = dcc_annot, which = 'column', na_col = "grey")
    path_ha = HeatmapAnnotation(df = path_annot, which = 'row', na_col = "grey")
    
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
expr_type <- unique(gsea_all$expr_signal)[1]
pathways_type <- unique(gsea_all$path_type)[1]
seg <- 'stroma'
nact_status <- 'pre'
    
out_name <- paste(expr_type, pathways_type, seg, nact_status, sep ='_')
print(out_name)
    
# filter to all signal, cell state/process pathways and ROIs with labels
gsea_sel_paths <- gsea_all[gsea_all$expr_signal == expr_type & gsea_all$path_type == pathways_type, ]
gsea_sel_paths_seg <- gsea_sel_paths[gsea_sel_paths$dcc_filename %in% 
                                       metadt_labels$dcc_filename[metadt_labels$Segment == seg & metadt_labels$NACT_status == nact_status], ]

# create matrix with values for each pathway
gsea_sel_wide_seg <- spread(gsea_sel_paths_seg[, c('dcc_filename','pathway', 'ssgsea_score')],
                        key = 'pathway', value = 'ssgsea_score') %>%
  column_to_rownames(var="dcc_filename")

gsea_seg_corr <- cor(gsea_sel_wide_seg, method = 'pearson')

# move all cor < thr to 0 and select only pathways with any cors
gsea_seg_corr_filt <- gsea_seg_corr
gsea_seg_corr_filt[gsea_seg_corr_filt >= -corr_thr & gsea_seg_corr_filt <= corr_thr] <- 0

cor_abovethr <- (colSums(gsea_seg_corr_filt, na.rm=T) != 1) # 1 because there always cor =1 to itself 
gsea_seg_corr_filt <- gsea_seg_corr_filt[cor_abovethr, cor_abovethr]


gsea_seg_corr <- gsea_seg_corr_filt
# row annotations based on pathways
path_annot <- gsea_sel_paths_seg %>%
  select('pathway', !!path_annots) %>%
  distinct() %>%
  column_to_rownames(var="pathway") %>%
  mutate(cell_type=replace(cell_type, cell_type=='Tcells_other', 'Tcells')) %>%
  mutate(cell_type=replace(cell_type, cell_type=='immune_other', 'immune')) %>%
  mutate_all(as.factor)

path_annot <- path_annot[match(rownames(gsea_seg_corr), rownames(path_annot)),] # ensure ordering

# set up annotation
path_ha = HeatmapAnnotation(df = path_annot, which = 'row', na_col = "grey")

# do the hmap
png(filename=file.path(out_dir, paste0('hmap_corr_', out_name, '_labs_deconv_filt03.png')), width=11, height=7,units="in",res=1000)

gsea_corr_heat <- Heatmap(gsea_seg_corr, name = paste0("pearson correlation of ssGSEA scores of ", out_name), 
                     left_annotation = path_ha,
                     heatmap_legend_param = list(
                       legend_direction = "horizontal", 
                       legend_width = unit(2, "in")),
                     show_column_names = T,
                     row_names_gp = gpar(fontsize = 4),
                     column_names_gp = gpar(fontsize = 4),
                     row_names_max_width = unit(5, "in")
)

draw(gsea_corr_heat, 
     heatmap_legend_side="top", 
     annotation_legend_side="top",
     merge_legend = TRUE)
dev.off()
