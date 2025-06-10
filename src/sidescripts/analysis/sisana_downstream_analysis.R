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

outp_dir <- '~/Documents/phd/st/gene-regulatory-networks/sisana/sisana/geomx_output'

dge_path <- '~/Documents/phd/st/geomx-processing/results/batch12-1205-no-counts-shift2/dge/dge_within_slide_Segment_bin_FALSE__/dge_all_dge_within_slide_Segment_bin_FALSE__.csv'
gsea_dge_dir_path <- '~/Documents/phd/st/geomx-processing/results/batch12-1205-no-counts-shift2/dge/dge_within_slide_Segment_bin_FALSE__/gsea_enrichment/hallmark_gobp'
gsea_dge_path <- file.path(gsea_dge_dir_path, 'gsea_dge_clust_msigdb__dge_deconv_Tcells_fc0.5_nofiltering.csv')

geomx_obj_path <- '~/Documents/phd/st/geomx-processing/results/batch12-1205-no-counts-shift2/geomx_qc_norm_batch_eff_rm.RDS'

dge_path <- '~/Documents/phd/st/geomx-processing/results/batch12-1205-no-counts-shift2/dge/dge_within_slide_Segment_bin_FALSE__/dge_all_dge_within_slide_Segment_bin_FALSE__.csv'
#######
# TODO iterate through all 
gsea_dir_path <- file.path(outp_dir, 'gsea', 'gsea_indegree_kegg')

gsea_res_path <- file.path(gsea_dir_path, 'gseapy.gene_set.prerank.report.csv')

fdr_thr <- 0.05
fwer_thr <- 0.01
jaccard_hclust_cuts <- c(0.5, 1, 1.2, 1.5)

expr_mean_diff_path <- file.path(outp_dir, 'compare_means', 'comparison_mw_between_stroma_tumor_expression.txt')
ind_mean_diff_path <- file.path(outp_dir, 'compare_means', 'comparison_mw_between_stroma_tumor_indegree.txt')
ind_vals_path <- file.path(outp_dir, 'network', 'lioness_indegree.csv')



# umap_vars <- c(aoi_segment_var, main_roi_label, main_experimental_condition, sample_name, 
#                main_batch_var, batch_var, other_vars_bio, other_vars_tech)

umap_vars <- c('Segment', 'NACT_status', 'Annotation_cell', 'Sample', 
               'main_batch_nr', 'batch_nr', 'Segment_geomx', 'Patient', 'Site')

dir.create(file.path(outp_dir, 'clustering'), recursive = T, showWarnings = F)


# load results ------------------------------------------------------------

gsea_res <- fread(gsea_res_path)

gsea_sign <- gsea_res[gsea_res$`FDR q-val` <= fdr_thr, ]
gsea_sign <- gsea_sign[gsea_sign$`FWER p-val` <= fwer_thr]


# collapse to main pathways -----------------------------------------------
# have to be hacked from fgsea
# collapsedPathways <- collapsePathways(gsea_res, sign_list, dge_sub_rank)
# gsea_res$is_main_pathway <- ifelse(gsea_res$pathway %in% collapsedPathways$mainPathways, 'yes', 'no')


# cluster pathways based on jaccard idx -----------------------------------

paths_genes_list <- lapply(gsea_sign$Lead_genes, function(x){
  genelist <- unlist(strsplit(x, split=';', fixed=T))
})

names(paths_genes_list) <- gsea_sign$Term

# calculate jaccard score between each pathway leading gene set
path_jaccard <- lapply(paths_genes_list, function(x){
  p1 <- lapply(paths_genes_list, function(y){
    jacc_idx <- as.numeric(round(length(intersect(x, y)) / length(union(x,y)), digits = 4))
  })
  return(unlist(p1))
})

path_jaccard_mtx <- do.call('cbind', path_jaccard)

# clustering with hclust
path_hclust <- hclust(dist(path_jaccard_mtx), method = "average")
#plot(path_hclust, hang = -1, cex = 0.4)

for(cutnr in jaccard_hclust_cuts){
  # cut the hclust tree at given point
  path_hclust_cut <- cutree(path_hclust, h = cutnr)
  
  # merge with gsea result
  if(identical(gsea_sign$Term, names(path_hclust_cut))){
    gsea_sign[[paste0('path_cluster_cut_', gsub('\\.', '', as.character(cutnr)))]] <- path_hclust_cut
  }
}


heatmap(path_jaccard_mtx)
#######
# make clustered heatmap
nes_anno <- sapply(colnames(path_jaccard_mtx), function(path){
  nes <- ifelse(sign(gsea_sign$NES[gsea_sign$Term == path]) == 1, 'pos', 'neg')
})

ha = HeatmapAnnotation(
  NES = anno_simple(nes_anno, col = c("pos" = "green", "neg" = "blue")),
  annotation_name_side = "left")

png(filename=file.path(gsea_dir_path, paste0('hmap_fdr', as.character(fdr_thr), '_fwer',as.character(fwer_thr), '.png')), 
    width=8, height=6,units="in",res=1000)

condition_heat <- Heatmap(as.matrix(path_jaccard_mtx), border="white",
                          rect_gp = gpar(col = "white", lwd = 2), column_title = 'tumor vs stroma',
                          cluster_columns = T, cluster_rows= T, col = brewer.pal(5, "YlOrRd"),
                          show_heatmap_legend = F, top_annotation = ha,
                          row_names_gp = gpar(fontsize = 6),
                          column_names_gp = gpar(fontsize = 6))

draw(condition_heat)
dev.off()
#######


# get one path per cluster ------------------------------------------------

cutnr <- 1.2
cut_colname <- paste0('path_cluster_cut_', gsub('\\.', '', as.character(cutnr)))

length(unique(gsea_sign$path_cluster_cut_05))
length(unique(gsea_sign$path_cluster_cut_1))
length(unique(gsea_sign$path_cluster_cut_12))
length(unique(gsea_sign$path_cluster_cut_15))


for(sign in c('pos', 'neg')){
  s <- ifelse(sign == 'pos', 1, -1)

  gsea_main <- gsea_sign[sign(gsea_sign$NES) == s, ]
  gsea_main <- arrange(gsea_main, NES)
  
  # keep 1st pathway from cluster
  gsea_main <- gsea_main[!(duplicated(gsea_main[[cut_colname]])), ]
  
  if(nrow(gsea_main) > 0){
    fwrite(gsea_main, file.path(gsea_dir_path, paste0('gsea_clustered_fdr', as.character(fdr_thr), '_fwer',as.character(fwer_thr), '_',
                                                      sign, '_', cut_colname, '.csv')))
  }

}


##################

# check gsea from DGE

cutnr <- 1.2
cut_colname <- paste0('path_cluster_cut_', gsub('\\.', '', as.character(cutnr)))


dge <- fread(dge_path)
gsea_dge <- fread(gsea_dge_path)


for(sign in c('pos', 'neg')){
  s <- ifelse(sign == 'pos', 1, -1)
  
  gsea_dge_main <- gsea_dge[sign(gsea_dge$NES) == s, ]
  gsea_dge_main <- arrange(gsea_dge_main, -NES)
  
  # keep 1st pathway from cluster
  gsea_dge_main <- gsea_dge_main[!(duplicated(gsea_dge_main[[cut_colname]])), ]
  
  if(nrow(gsea_dge_main) > 0){
    # TODO better name
    fwrite(gsea_dge_main, file.path(gsea_dge_dir_path, paste0('gsea_clustered_best_fc05_Tcells', '_', sign, '_', cut_colname, '.csv')))
  }
  
}

# for now only neg bcs there were only neg found from indegrees
gsea_dge_neg <- gsea_dge[sign(gsea_dge$NES) == -1, ]
gsea_dge_neg <- arrange(gsea_dge_neg, NES)

# keep 1st pathway from cluster
gsea_dge_main <- gsea_dge_neg[!(duplicated(gsea_dge_neg[[cut_colname]])), ]


gsea_dge_genes_list <- lapply(gsea_dge_neg$leadingEdge, function(x){
  genelist <- unlist(strsplit(x, split='|', fixed=T))
})

names(gsea_dge_genes_list) <- gsea_dge_neg$pathway

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



# networks indegree heatmap and clustering --------------------------------
top_var <- 1000

metadata <- sData(readRDS(geomx_obj_path))
metadata$dcc_filename <- gsub('\\-', '\\.', metadata$dcc_filename)

ind_vals <- fread(ind_vals_path)
ind_vals <- as.matrix(column_to_rownames(ind_vals, 'Target'))

# for expression (comments dcc names change)
#ind_vals <- readRDS(geomx_obj_path)@assayData$harmony_batch_corr

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


png(filename=file.path(outp_dir, paste0('hmap_expr_clustered', '_topvar', as.character(top_var), '.png')), 
    width=8, height=6,units="in",res=2000)

ind_heat <- Heatmap(ind_vals_var_zscore, cluster_columns = T, cluster_rows= T,
              show_row_names = FALSE, show_column_names = FALSE, show_row_dend = FALSE,
              top_annotation = ha)


draw(ind_heat)
dev.off()
#######


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
