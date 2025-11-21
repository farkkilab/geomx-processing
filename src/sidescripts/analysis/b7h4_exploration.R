library(readxl)
library(data.table)
library(plyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(reshape2)
library(GeomxTools)
library(NanoStringNCTools)

# boxplot across main values - Segment, NACT status, Site, PFS/OS
# paired boxplot pre vs post coloured by PFS/OS
# boxplot for each deconv ct
# TODO scatter expr vs cell_type_fraq for AOI/matched Roi

scrna_anno <- 'mid_lvl_ct_updated'
ct_of_interest <<- c("tumor", "Macrophages_Monocytes", "Tcells_CD8", "Tcells_CD4", "DCs", "Bcells", "Fibroblasts_Mesothelial")

proj_dir <<- '~/Documents/phd/st'
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808')
out_dir <- file.path(output_dir, 'exploration', 'vtcn1') # for plots

geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS') 
norm_type <- 'harmony_batch_corr_q3_norm'

bp_cellcounts_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', paste0('bp_res_', scrna_anno, '_ct_fraction.csv'))
bp_deconv_path <- file.path(output_dir, 'deconvolution', 'bayes_prism', paste0('bp_res_', scrna_anno, '_expr_mtx_cleaned_deseq2_vst_harmony_corr.RDS'))

#TODO loop through many genes
goi <- c('VTCN1')

source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))
dir.create(out_dir, recursive = T, showWarnings = F)

# TODO sth happened to Roi column in metadata - check it from beginning - now just merging
metadata_orig_path <- '/home/iganiemi/Documents/phd/st/data/geomx/batch123/metadata/dcc_metadata_batch123_no_tls.xlsx' 

#########
# names
aoi_id <<- 'dcc_filename'
roi_id <<- 'Roi'

main_batch_var <- 'main_batch_nr'
batch_var <<- 'batch_nr'


aoi_segment_var <<- "Segment"
main_roi_label <<- "Annotation_cell" 
main_experimental_condition <<- 'NACT_status'
sample_name <<- 'Sample'

other_vars_bio <<- c("Segment_geomx", "Patient", "Site") # 'PFS_months', 'PFS' , "tls_status"
other_vars_tech <<- c('Slide_Name')

# tu use for plots
meta_names <- c(aoi_id, roi_id, aoi_segment_var, sample_name, main_experimental_condition, 
                other_vars_bio, 'paired_status', 'PFS_quartile_b123', 'OS_quartile_b123',
                'PFS_median_b123', 'OS_median_b123', 'PFS_quartile_paired', 'OS_quartile_paired',
                'PFS_median_paired', 'OS_median_paired')


########################################
# load data

geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)

deconv_res <- readRDS(bp_deconv_path)
deconv_res <- deconv_res[names(deconv_res) %in% ct_of_interest]

# name fix
geomx_obj@phenoData$Site <- ifelse(geomx_obj@phenoData$Site == 'Peritroneum', 'Peritoneum', geomx_obj@phenoData$Site)

# TODO!! - important - check when the bug was introduced and fix + rerun DGE 
geomx_obj@phenoData$PFS_quartile_b123 <- ifelse(geomx_obj@phenoData$Patient == 'S032', 4,
                                                ifelse(geomx_obj@phenoData$Patient == 'S309', 4, geomx_obj@phenoData$PFS_quartile_b123))
geomx_obj@phenoData$OS_quartile_b123 <- ifelse(geomx_obj@phenoData$Patient == 'S032', 4,
                                                ifelse(geomx_obj@phenoData$Patient == 'S309', 2, geomx_obj@phenoData$OS_quartile_b123))

# add stuff to metadata
# TODO move somewhere

geomx_obj@phenoData$PFS_median_b123 <- ifelse(geomx_obj@phenoData$PFS_quartile_b123 %in% c(1, 2), 1, 2)
geomx_obj@phenoData$OS_median_b123 <- ifelse(geomx_obj@phenoData$OS_quartile_b123 %in% c(1, 2), 1, 2)

geomx_obj@phenoData$OS_quartile_paired <- mapvalues(geomx_obj@phenoData$Patient, 
                                                    from=c("S015", "S027", "S032", "S069", "S084", "S139",
                                                           "S229", "S333"), 
                                                    to=c(4, 3, 4, 2, 1, 1, 3, 2))

geomx_obj@phenoData$OS_quartile_paired <- ifelse(geomx_obj@phenoData$OS_quartile_paired %in% c(1, 2, 3, 4), 
                                                 geomx_obj@phenoData$OS_quartile_paired, 0)

geomx_obj@phenoData$PFS_quartile_paired <- mapvalues(geomx_obj@phenoData$Patient, 
                                                     from=c("S015", "S027", "S032", "S069", "S084", "S139",
                                                            "S229", "S333"), 
                                                     to=c(3, 2, 4, 2, 1, 1, 4, 3))

geomx_obj@phenoData$PFS_quartile_paired <- ifelse(geomx_obj@phenoData$PFS_quartile_paired %in% c(1, 2, 3, 4), 
                                                  geomx_obj@phenoData$PFS_quartile_paired, 0)

geomx_obj@phenoData$PFS_median_paired <- ifelse(geomx_obj@phenoData$PFS_quartile_paired %in% c(1, 2), 1,
                                                ifelse(geomx_obj@phenoData$PFS_quartile_paired %in% c(3, 4), 2, 0))

geomx_obj@phenoData$OS_median_paired <- ifelse(geomx_obj@phenoData$OS_quartile_paired %in% c(1, 2), 1,
                                               ifelse(geomx_obj@phenoData$OS_quartile_paired %in% c(3, 4), 2, 0))

meta_orig <- read_excel(metadata_orig_path)

geomx_obj@phenoData@data <- left_join(geomx_obj@phenoData@data, meta_orig[, c('dcc_filename', 'Roi')])

###################################################
# make gene expr df for all + deconv_data

genes_expr_list <- lapply(c('all', ct_of_interest), function(ct){
  
  if(ct == 'all'){
    expr <- geomx_obj@assayData[[norm_type]]
  } else{
    expr <- deconv_res[[ct]]
  }

  goi_in_dt <- goi[which(goi %in% rownames(expr))]
  
  if(length(goi_in_dt) > 1){
    genes_expr <- data.frame(t(expr[rownames(expr) %in% goi_in_dt, ]))
  } else if(length(goi_in_dt) == 1){
    genes_expr <- data.frame(expr[rownames(expr) %in% goi_in_dt, ])
    colnames(genes_expr) <- goi_in_dt
  } else(
    stop('not enough goi in dataset')
  )
  
  genes_expr <- rownames_to_column(genes_expr, var = aoi_id)
  genes_expr$cell_type <- ct
  return(genes_expr)
})

genes_expr_list <- do.call(rbind, genes_expr_list)

genes_expr_long <- dplyr::left_join(geomx_obj@phenoData@data[meta_names], genes_expr_list, by = aoi_id)

genes_expr_long <- melt(genes_expr_long, id.vars = c(colnames(genes_expr_long)[1: length(meta_names)], 'cell_type'),
                        variable.name = "gene_name", 
                        value.name = "gene_expr")

#############################################################################
# boxplot for expression for main vars
manual_colours <- brewer.pal(12, 'Paired')

for(ct in c('all', 'tumor')){
  
  genes_expr_long_ct <- genes_expr_long[genes_expr_long$cell_type == ct, ]
  
  genes_expr_long_post_ct <- genes_expr_long_ct[genes_expr_long_ct$NACT_status == 'post', ]
  genes_expr_long_pre_ct <- genes_expr_long_ct[genes_expr_long_ct$NACT_status == 'pre', ]
  
  # TODO loop through gene
  for(color_colname in c('Segment')){
    pathway_boxplot(genes_expr_long_ct, 'gene_name', 'gene_expr', color_colname, facet_var = c(),
                    plot_title = paste0('genes expr per ', color_colname),
                    output_path = file.path(out_dir, paste0('boxpl_expr_', color_colname, '_', ct, '.pdf')),
                    ymin=0, ymax= max(genes_expr_long_ct$gene_expr) + 1, manual_colours = manual_colours)
  }
  
  
  for(color_colname in c('NACT_status')){
    pathway_boxplot(genes_expr_long_ct, 'gene_name', 'gene_expr', color_colname, facet_var = 'Segment',
                    plot_title = paste0('genes expr per ', color_colname),
                    output_path = file.path(out_dir, paste0('boxpl_expr_', color_colname, '_', ct, '_per_segm.pdf')),
                    ymin=0, ymax= max(genes_expr_long_ct$gene_expr) + 1, manual_colours = manual_colours)
  }
  
  for(color_colname in c('Site', 'PFS_quartile_b123', 'OS_quartile_b123', 'PFS_median_b123', 'OS_median_b123')){
    pathway_boxplot(genes_expr_long_post_ct, 'gene_name', 'gene_expr', color_colname, facet_var = 'Segment',
                    plot_title = paste0('genes expr per ', color_colname),
                    output_path = file.path(out_dir, paste0('boxpl_expr_', color_colname, '_', ct, '_per_segm_post.pdf')),
                    ymin=0, ymax= max(genes_expr_long_post_ct$gene_expr) + 1, manual_colours = manual_colours)
  }
  
  for(color_colname in c('Site', 'PFS_median_paired', 'OS_median_paired')){
    pathway_boxplot(genes_expr_long_pre_ct, 'gene_name', 'gene_expr', color_colname, facet_var = 'Segment',
                    plot_title = paste0('genes expr per ', color_colname),
                    output_path = file.path(out_dir, paste0('boxpl_expr_', color_colname, '_', ct, '_per_segm_pre.pdf')),
                    ymin=0, ymax= max(genes_expr_long_post_ct$gene_expr) + 1, manual_colours = manual_colours)
  }
  

}


###############################################################
# paired boxplot pre vs post coloured by OS/PFS

for(ct in c('all', 'tumor')){
  
  genes_expr_long_ct <- genes_expr_long[genes_expr_long$cell_type == ct, ]
  
  genes_expr_long_ct_paired <- genes_expr_long_ct[genes_expr_long_ct$paired_status == 'paired', ]
  
  genes_expr_long_ct_paired_mean <- group_by(genes_expr_long_ct_paired, Patient, NACT_status, gene_name, Segment) %>%
    summarise(gexpr_mean = mean(gene_expr)) %>%
    left_join(genes_expr_long_ct_paired, by = c('Patient', 'NACT_status', 'gene_name', 'Segment')) %>%
    distinct(Patient, NACT_status, gene_name, Segment, .keep_all = TRUE)
  
  # TODO loop through gene
  for(color_colname in c('OS_quartile_paired', 'PFS_quartile_paired', 
                         'OS_median_paired', 'PFS_median_paired')){
    ggplot(genes_expr_long_ct_paired_mean, aes(x = NACT_status, y = gexpr_mean, fill = as.factor(get(color_colname)))) + 
      #geom_boxplot(alpha = .2) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
      geom_line(aes(group = Patient), size = 0.2, alpha = 0.8) + 
      ggtitle(paste0('paired gene_expr')) +
      guides(fill=guide_legend(title=color_colname)) +
      facet_wrap(~ Segment)
    
    ggsave(file.path(out_dir, paste0('paired_comparison_', color_colname, '_', ct,  '_.png')),
           width = 1500, height = 1000, unit = 'px')
  }
}


####################################################################
# expression in different cts
genes_expr_long_post <- genes_expr_long[genes_expr_long$NACT_status == 'post', ]

color_colname <- 'cell_type'
pathway_boxplot(genes_expr_long_post, 'gene_name', 'gene_expr', color_colname, facet_var = 'Segment',
                plot_title = paste0('genes expr per ', color_colname),
                output_path = file.path(out_dir, paste0('boxpl_expr_', color_colname, '_per_segm_post.pdf')),
                ymin=0, ymax= max(genes_expr_long_post$gene_expr) + 1, manual_colours = manual_colours)

# version without wilxoc
gene_boxpl <- ggplot(data = genes_expr_long_post, aes(x = gene_name, y = gene_expr, fill = as.factor(get(color_colname)))) +
  geom_boxplot() +
  # geom_point(position= position_jitterdodge(dodge.width = 1, jitter.width= .3, jitter.height = 0),
  #            size= 0.2, alpha = 0.6) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 5)) +
  ggtitle(paste0('genes expr per ', color_colname))+
  xlab('gene_name') +
  ylab('gene_expr') +
  guides(fill=guide_legend(title=color_colname)) +
  scale_fill_manual(values=manual_colours) +
  facet_wrap(~get('Segment'), scales = "fixed", dir="v")


pdf(file= file.path(out_dir, paste0('boxpl_expr_', color_colname, '_per_segm_post_no_wilcox.pdf')), width=8, height=5)
plot(gene_boxpl)
dev.off()

###################################################################
# merge expr df with bp cell counts
genes_expr_long$sample_roi <- paste0(genes_expr_long$Sample, '_', genes_expr_long$Roi)


