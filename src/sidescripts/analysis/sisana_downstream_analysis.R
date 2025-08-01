library(data.table)
library(fgsea)
library(ComplexHeatmap)
library(RColorBrewer)
library(tibble)
library(plyr)
library(dplyr)
library(umap)
library(Rtsne)
library(GeomxTools)
library(preprocessCore)
library(PCAtools)
library(DESeq2)
library(msigdbr)
library(biomaRt)
library(parallel)
library(ComplexHeatmap)

outp_dir <- '~/Documents/phd/st/gene-regulatory-networks/sisana/sisana/geomx_output_deseq2_harmony'

geomx_obj_path <- '~/Documents/phd/st/geomx-processing/results/batch12-1205-no-counts-shift2/geomx_qc_norm_batch_eff_rm.RDS'


#####
dge_path <- '~/Documents/phd/st/geomx-processing/results/batch12-1205-no-counts-shift2/dge/dge_within_slide_Segment_bin_FALSE__/dge_all_dge_within_slide_Segment_bin_FALSE__.csv'
#######

# for checking gse enrichment on my dge
gsea_dge_dir_path <- '~/Documents/phd/st/geomx-processing/results/batch12-1205-no-counts-shift2/dge/dge_within_slide_Segment_bin_FALSE__/gsea_enrichment/kegg/'
gsea_dge_path <- file.path(gsea_dge_dir_path, 'gsea_dge_clust_msigdb__dge_all_fc0.5_nofiltering.csv')

# expression and indegrees difference from sisana
expr_mean_diff_path <- file.path(outp_dir, 'compare_means', 'comparison_mw_between_stroma_tumor_expression.txt')
ind_mean_diff_path <- file.path(outp_dir, 'compare_means', 'comparison_mw_between_stroma_tumor_indegree.txt')

# indegree results
ind_vals_path <- file.path(outp_dir, 'network', 'lioness_indegree.csv')
out_vals_path <- file.path(outp_dir, 'network', 'lioness_outdegree.csv')

# gsea on expr/ind diff 
gsea_dirs <- list.dirs(file.path(outp_dir, 'gsea'), recursive = F)

gsea_clust_df_path <- file.path(outp_dir, 'gsea_from_geomx', 'gsea_indegrees_gobp', 'gsea_sign_fdr0.05_clust.csv')
geneset_filt_edges_path <- file.path(outp_dir, 'gsea_from_geomx', 'lioness_filtered_for_APP_genes.csv')

fdr_thr <- 0.05
fwer_thr <- 0.01
jaccard_hclust_cuts <- c(0.5, 1, 1.2, 1.5)
cutnr <- 1.2

# umap_vars <- c(aoi_segment_var, main_roi_label, main_experimental_condition, sample_name, 
#                main_batch_var, batch_var, other_vars_bio, other_vars_tech)

umap_vars <- c('Segment', 'NACT_status', 'Annotation_cell', 'Sample', 
               'main_batch_nr', 'batch_nr', 'Segment_geomx', 'Patient', 'Site')

source(file.path('~/Documents/phd/st', 'geomx-processing', 'src', 'geomx_utils.R'))

dir.create(file.path(outp_dir, 'clustering'), recursive = T, showWarnings = F)

#############################################################################
#############################################################################
# prepare input dfs for sisana
# output_dir <- '/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2706'
# dir.create(file.path(output_dir, 'sisana', 'geomx_input'), recursive = T)
# geomx_obj <- readRDS(file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS'))
# 
# geomx_harm <- data.frame(geomx_obj@assayData$harmony_batch_corr)
# geomx_harm <- rownames_to_column(geomx_harm, var = 'Target')
# geomx_harm_demo <- geomx_harm[, 1:10]
# write_tsv(geomx_harm, file.path(output_dir, 'sisana', 'geomx_input', 'geomx_batch123_harmony_corr_expr_mtx.tsv'),
#           col_names = T)
# write_tsv(geomx_harm_demo, file.path(output_dir, 'sisana', 'geomx_input', 'geomx_batch123_harmony_corr_expr_mtx_demo.tsv'),
#           col_names = T)

# meta_segment <- sData(geomx_obj)[, c('dcc_filename', 'Segment')]
# meta_segment$dcc_filename <- gsub('\\-', '\\.', meta_segment$dcc_filename)
# rownames(meta_segment) <- NULL
# colnames(meta_segment) <- NULL
# 
# write_csv(meta_segment, file.path(output_dir, 'sisana', 'geomx_input', 'geomx_meta_segment.csv'), 
#           col_names = F)

############################################################################
############################################################################
# cluster gsea pathways for sisana results --------------------------------

for(gsea_dir_path in gsea_dirs){
  
  print(gsea_dir_path)
  
  # load gsea results
  gsea_res <- fread(file.path(gsea_dir_path, 'gseapy.gene_set.prerank.report.csv'))
  
  # filter to significant results
  gsea_sign <- gsea_res[gsea_res$`FDR q-val` <= fdr_thr, ]
  gsea_sign <- gsea_sign[gsea_sign$`FWER p-val` <= fwer_thr]
  
  if(nrow(gsea_sign) > 0){
    
    # cluster pathways by jaccard idx and make heatmap
    hmap_outpath <- file.path(gsea_dir_path, paste0('hmap_fdr', as.character(fdr_thr), '_fwer',as.character(fwer_thr), '.png'))
    gsea_clust <- cluster_gsea_enrichment(gsea_sign, 'Lead_genes', 'Term', hmap_outpath = hmap_outpath, hmap_title = 'stroma vs tumor')
    
    # get best pathway per cluster
    for(clust_cutnr in jaccard_hclust_cuts){
      
      print(clust_cutnr)
      
      clust_colname <- paste0('path_cluster_cut_', gsub('\\.', '', as.character(clust_cutnr)))
      print(length(unique(gsea_clust[[clust_colname]])))
      
      gsea_best <- best_pathway_clust(gsea_clust, clust_colname, nes_colname = 'NES')
      
      fwrite(gsea_best, file.path(gsea_dir_path, paste0('gsea_sign_fdr', as.character(fdr_thr),  '_clust_best_cut_', as.character(clust_cutnr), '.csv')))
      
    }
  }
}

#############################################################################
#############################################################################
# run compare means results through geomx pipeline gsea enrichment --------

geomx_obj <- readRDS(geomx_obj_path)
diff_expr <- fread(expr_mean_diff_path)
diff_ind <- fread(ind_mean_diff_path)

# TODO iterate by expr + indegrees + gobp + kegg
gene_diff_df <- diff_ind

#'CP:KEGG_MEDICUS', 'GO:BP'
compute_hallmark <- F
msigdb_subcat <- c('CP:KEGG_MEDICUS')
min_sign_gene_nr <- 10

outp_gsea_dir <- file.path(outp_dir, 'gsea_from_geomx', 'gsea_indegrees_kegg')
dir.create(outp_gsea_dir, recursive = T)

#####################
# make signatures
sign_list <- prepare_msigdb_sign_list(adjust_synonym = T, geomx_obj = geomx_obj, hal = compute_hallmark, 
                                      db_subcat_list = msigdb_subcat)

sign_list <- sign_list[sapply(sign_list, length) >= min_sign_gene_nr]


# rank and do gsea enrichment
gsea_res <- rank_genes_and_do_gsea_enrichment(gene_diff_df, 'difference_of_means_(tumor-stroma)',
                                              'mw_pvalue', 'Target', sign_list, gsea_scoretype = 'std')

fwrite(gsea_res, file.path(outp_gsea_dir, 'gsea_res_all.csv'))

# filter to significant results and main pathways
gsea_sign <- gsea_res[gsea_res$padj <= fdr_thr, ]
gsea_sign <- gsea_sign[gsea_sign$is_main_pathway == 'yes', ]

# cluster pathways by jaccard idx and make heatmap
hmap_outpath <- file.path(outp_gsea_dir, paste0('hmap_fdr', as.character(fdr_thr), '.png'))

gsea_clust <- cluster_gsea_enrichment(gsea_sign, 'leadingEdge', 'pathway', hmap_outpath = hmap_outpath, hmap_title = 'stroma vs tumor')

fwrite(gsea_clust, file.path(outp_gsea_dir, paste0('gsea_sign_fdr', as.character(fdr_thr),  '_clust.csv')))

# get one path per cluster ------------------------------------------------

for(clust_cutnr in jaccard_hclust_cuts){
  
  print(clust_cutnr)
  
  clust_colname <- paste0('path_cluster_cut_', gsub('\\.', '', as.character(clust_cutnr)))
  print(length(unique(gsea_clust[[clust_colname]])))
  
  gsea_best <- best_pathway_clust(gsea_clust, clust_colname, nes_colname = 'NES')
  
  fwrite(gsea_best, file.path(outp_gsea_dir, paste0('gsea_sign_fdr', as.character(fdr_thr),  '_clust_best_cut_', as.character(clust_cutnr), '.csv')))
  
}


#####################################################
#####################################################
# best clustered pathway for gsea from DGE

gsea_dge <- fread(gsea_dge_path)

for(clust_cutnr in jaccard_hclust_cuts){
  
  print(clust_cutnr)
  
  clust_colname <- paste0('path_cluster_cut_', gsub('\\.', '', as.character(clust_cutnr)))
  print(length(unique(gsea_dge[[clust_colname]])))
  
  gsea_best <- best_pathway_clust(gsea_dge, clust_colname, nes_colname = 'NES')
  
  fwrite(gsea_best, file.path(outp_dir,'gsea', paste0('gsea_from_dge_sign_fdr', as.character(fdr_thr),  '_clust_best_cut_', as.character(clust_cutnr), '_kegg.csv')))
  
}


######################################
#######################################
# choosing interesting gene set for edges exploration

gsea_clust_df <- fread(gsea_clust_df_path)
clust_cutnr <- 1.2 #1.2 - clust 6, # 1.5 clust 4
clust_name <- 6
lead_genes_colname <- 'leadingEdge'
lead_genes_split <- '|'

clust_colname <- paste0('path_cluster_cut_', gsub('\\.', '', as.character(clust_cutnr)))

sort(table(gsea_clust_df[[clust_colname]]))
length(unique(gsea_clust_df[[clust_colname]]))

path_clust <- gsea_clust_df[gsea_clust_df[[clust_colname]] == clust_name, ]

pathclust_genes <- lapply(path_clust[[lead_genes_colname]], function(x){
  genelist <- unlist(strsplit(x, split=lead_genes_split, fixed=T)) # TODO wtf it looks like any split works
})

pathclust_genes <- unique(unlist(pathclust_genes))

# save in txt to use in sisana extract command
fwrite(list(pathclust_genes), file = file.path(outp_dir, 'gsea_from_geomx', 'APP_genes.txt'))

#######################################################
#######################################################
# edges exploration
geomx_obj_meta <- readRDS(geomx_obj_path)@phenoData@data
geomx_obj_meta$dcc_filename <- gsub('\\-', '\\.', geomx_obj_meta$dcc_filename)
top_edge_nr <- 100

# file from sisana extract genes (using lioness.pickle and APP_genes.txt)
geneset_edges <- fread(geneset_filt_edges_path)

# move to long 
geneset_edges_long <- melt(setDT(geneset_edges), id.vars = c('TF', 'Target'), variable.name = "dcc_filename")
geneset_edges_long <- full_join(as.data.frame(geneset_edges_long), geomx_obj_meta[, c('dcc_filename', 'Segment')])
geneset_edges_long$TF_Target <- paste0(geneset_edges_long$TF, '_', geneset_edges_long$Target)

# select top n edges per sample
geneset_edges_long_top <- geneset_edges_long %>% 
  arrange(desc(value)) %>% 
  group_by(dcc_filename) %>% slice_head(n = top_edge_nr)

geneset_edges_long_top_str <- geneset_edges_long_top[geneset_edges_long_top$Segment == 'stroma', ]
geneset_edges_long_top_tum <- geneset_edges_long_top[geneset_edges_long_top$Segment == 'tumor', ]

# divide by segment and back to wide for cytoscape
geneset_edges_top_wide_str <- dcast(setDT(geneset_edges_long_top[geneset_edges_long_top$Segment == 'stroma', ]),
                               TF+Target ~ dcc_filename,
                               value.var = "value")

geneset_edges_top_wide_tum <- dcast(setDT(geneset_edges_long_top[geneset_edges_long_top$Segment == 'tumor', ]),
                                    TF+Target ~ dcc_filename,
                                    value.var = "value")

fwrite(geneset_edges_top_wide_str, file = file.path(outp_dir, 'gsea_from_geomx', 'lioness_APP_genes_top_100_per_sample_str.csv'))
fwrite(geneset_edges_top_wide_tum, file = file.path(outp_dir, 'gsea_from_geomx', 'lioness_APP_genes_top_100_per_sample_tum.csv'))


sort(table(geneset_edges_long_top_str$TF), decreasing = T)[1:50]
sort(table(geneset_edges_long_top_tum$TF), decreasing = T)[1:50]

sort(table(geneset_edges_long_top_str$Target), decreasing = T)[1:50]
sort(table(geneset_edges_long_top_tum$Target), decreasing = T)[1:50]

length(unique(geneset_edges_long_top_str$TF))
length(unique(geneset_edges_long_top_tum$TF))
length(intersect(geneset_edges_long_top_str$TF, geneset_edges_long_top_tum$TF))

length(unique(geneset_edges_long_top$TF_Target))
length(intersect(geneset_edges_long_top_str$TF_Target, geneset_edges_long_top_tum$TF_Target))

####################
# that makes more sense I think
# sum edges by TF and gene per stroma and tumor separately

# for all 
geneset_edges_long_tfsum <- geneset_edges_long %>%
  group_by(Segment, TF) %>%
  summarise(TF_sum_segm = sum(value))

geneset_edges_long_targetsum <- geneset_edges_long %>%
  group_by(Segment, Target) %>%
  summarise(Target_sum_segm = sum(value))

# for top per sample
geneset_edges_long_top_tfsum <- geneset_edges_long_top %>%
  group_by(Segment, TF) %>%
  summarise(TF_sum_segm = sum(value))

geneset_edges_long_top_targetsum <- geneset_edges_long_top %>%
  group_by(Segment, Target) %>%
  summarise(Target_sum_segm = sum(value))


sort(geneset_edges_long_tfsum$TF_sum_segm[geneset_edges_long_tfsum$Segment == 'stroma'], decreasing = T)[1:50]
sort(geneset_edges_long_tfsum$TF_sum_segm[geneset_edges_long_tfsum$Segment == 'tumor'], decreasing = T)[1:50]


# top 100 edges per gene (sum by all samples together, tum/str separately)
geneset_edges_long_edgesum <- geneset_edges_long_top %>%
  group_by(Segment, TF_Target) %>%
  summarise(TF_Target_sum_segm = sum(value)) %>%
  mutate(Target = gsub('.*_', '', TF_Target))


geneset_edges_long_top_edges_per_gene <- geneset_edges_long_edgesum %>%
  arrange(desc(TF_Target_sum_segm)) %>% 
  group_by(Segment, Target) %>%
  slice_head(n = top_edge_nr)

##############################################################################
##############################################################################
##############################################################################
# networks indegree heatmap and clustering --------------------------------

top_var <- 1000

for(res_type in c('expr', 'ind')){
  metadata <- sData(readRDS(geomx_obj_path))
  
  if(res_type == 'expr'){
    ind_vals <- readRDS(geomx_obj_path)@assayData$harmony_batch_corr
  } else if(res_type == 'ind'){
    metadata$dcc_filename <- gsub('\\-', '\\.', metadata$dcc_filename) # change names to match indegrees
    
    ind_vals <- fread(ind_vals_path)
    ind_vals <- as.matrix(column_to_rownames(ind_vals, 'Target'))
  }
  
  identical(metadata$dcc_filename, colnames(ind_vals))
  
  per_gene_variance <- apply(ind_vals, 1, stats::var)
  top_var_ind <- names(sort(per_gene_variance, decreasing = T)[1:top_var])
  
  ind_vals_var <- ind_vals[rownames(ind_vals) %in% top_var_ind, ]
  
  ind_vals_var_zscore <- scale(ind_vals_var) # by column
  
  # do the heatmap
  #######
  # make clustered heatmap
  segment_anno <- sapply(colnames(ind_vals_var_zscore), function(s){
    seg <- metadata$Segment[metadata$dcc_filename == s]
  })
  
  patient_anno <- sapply(colnames(ind_vals_var_zscore), function(s){
    seg <- metadata$Patient[metadata$dcc_filename == s]
  })
  
  ha = HeatmapAnnotation(
    segment = anno_simple(segment_anno, col = c("stroma" = "green", "tumor" = "blue")),
    patient = anno_simple(patient_anno),
    annotation_name_side = "left")
  
  
  png(filename=file.path(outp_dir, paste0('hmap_', res_type, '_clustered', '_topvar', as.character(top_var), '.png')), 
      width=8, height=6,units="in",res=2000)
  
  ind_heat <- Heatmap(ind_vals_var_zscore, cluster_columns = T, cluster_rows= T,
                      show_row_names = FALSE, show_column_names = FALSE, show_row_dend = FALSE,
                      top_annotation = ha)
  
  
  draw(ind_heat)
  dev.off()
}


# networks UMAP -----------------------------------------------------------

# plot
plot_umap_tsne <- function(pheno_data, xcol, ycol, 
                           color_var, shape_var = 'Segment',
                           output_name){
  
  pheno_data[[color_var]] <- as.character(pheno_data[[color_var]])
  
  ggplot(pheno_data,
         aes(x = get(xcol), 
             y = get(ycol), 
             color = get(color_var), shape = get(shape_var))) +
    geom_point(size = 3) +
    scale_color_discrete(name = color_var) + 
    scale_shape_discrete(name = shape_var) + 
    theme_bw()
  
  ggsave(output_name, width = 2000, height = 1500, unit='px', device='png')
}
##############
# indegrees (maybe also outdegrees) to UMAP + visualisation + clustering
# clustering done using PCA from most variable indegrees - 2000 most variable, 50 PCs 

# 1 select 2000 most variable indegrees
# do PCA
# use top50 PCA as an input to UMAP

metadata <- sData(readRDS(geomx_obj_path))

# already after log transformation so no need to log again
ind_vals <- fread(ind_vals_path)
ind_vals <- as.matrix(column_to_rownames(ind_vals, 'Target')) 

# get 2000 most variable indegrees
per_gene_variance <- apply(ind_vals, 1, stats::var)
top_var_ind <- names(sort(per_gene_variance, decreasing = T)[1:2000])

ind_vals_var <- ind_vals[rownames(ind_vals) %in% top_var_ind, ]

# normalise
ind_vals_qnorm <- normalize.quantiles(ind_vals_var, keep.names = T)

# do PCA
ind_vals_var_pca <- pca(ind_vals_qnorm)
ind_vals_var_pca <- t(ind_vals_var_pca$rotated)
ind_vals_var_pca <- ind_vals_var_pca[1:50, ]


#######################

dim(ind_vals)
hist(as.matrix(ind_vals), breaks = 1000)

dim(ind_vals_var)
hist(as.matrix(ind_vals_var), breaks = 1000)

dim(ind_vals_qnorm)
hist(as.matrix(ind_vals_qnorm), breaks = 1000)

dim(ind_vals_var_pca)
hist(as.matrix(ind_vals_var_pca), breaks = 1000)

#########
dt_mtx_list <- list(indegrees = ind_vals, indegrees_qnorm = ind_vals_qnorm, 
                    indegrees_qnorm_var_pca = ind_vals_var_pca)

for(dt_mtx_name in names(dt_mtx_list)){
  
  dt_mtx <- dt_mtx_list[[dt_mtx_name]]
  
  custom_umap <- umap::umap.defaults
  custom_umap$random_state <- 42
  
  # make umap
  umap_out <- umap(t(dt_mtx), config = custom_umap)
  umap_out_res <- umap_out$layout
  colnames(umap_out_res) <- c('UMAP_1', 'UMAP_2')
  
  
  # set the seed for tSNE 
  set.seed(42) 
  
  # make tsne
  # perplexity not be bigger than 3 * perplexity < nrow(X) - 1
  tsne_out <- Rtsne(t(dt_mtx), perplexity = ncol(dt_mtx)*.15)
  tsne_out_res <- tsne_out$Y
  colnames(tsne_out_res) <- c('tSNE_1', 'tSNE_2')
  
  metadata <- cbind(metadata[, c('dcc_filename', umap_vars)], umap_out_res, tsne_out_res)

  for(colvar in umap_vars){
    plot_umap_tsne(metadata, 'UMAP_1', 'UMAP_2', color_var = colvar, 
                   output_name = file.path(outp_dir, 'clustering', paste0('umap_', dt_mtx_name, '_', colvar, '.png')))
    
    plot_umap_tsne(metadata, 'tSNE_1', 'tSNE_2', color_var = colvar, 
                   output_name = file.path(outp_dir, 'clustering', paste0('tsne_', dt_mtx_name, '_', colvar, '.png')))
  }
  
}

#####################################################
#####################################################
# quick checks

gsea_dge <- fread(gsea_dge_path)

gsea_dge_genes_list <- lapply(gsea_dge$leadingEdge, function(x){
  genelist <- unlist(strsplit(x, split='|', fixed=T))
})

names(gsea_dge_genes_list) <- gsea_dge$pathway

######################
# compare on the gene lvl for all leading edge genes

gsea_indg_genes <- unique(unlist(paths_genes_list))
gsea_dge_genes <- unique(unlist(gsea_dge_genes_list))

length(intersect(gsea_indg_genes, gsea_dge_genes))

# not so much - 101 out of 421 are the same

# compare expr dge and indegree means -------------------------------------

diff_expr <- fread(expr_mean_diff_path)
diff_ind <- fread(ind_mean_diff_path)

dge_res <- fread(dge_path)

diff_expr_sign <- diff_expr[diff_expr$FDR <= 0.05, ]
diff_expr_sign <- arrange(diff_expr_sign, -`difference_of_means_(tumor-stroma)`)
dge_sign <- dge_res[dge_res$FDR <= 0.05, ]
dge_sign <- arrange(dge_sign, Estimate)

length(intersect(diff_expr_sign$Target[1:100], dge_sign$Gene[1:100])) # 98/100 from overexpr tumor
length(intersect(diff_expr_sign$Target[(nrow(diff_expr_sign)-100): nrow(diff_expr_sign)], 
                 dge_sign$Gene[(nrow(dge_sign)-100): nrow(dge_sign)])) # 97/100 from overexpr stroma



