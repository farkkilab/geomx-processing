# unbiased GSEA/ORA on DGE results - to be added to main pipeline !!
#library(org.Hs.eg.db) # for enrichGO
library(fgsea)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# set variables -----------------------------------------------------------

# proj_dir <<- '~/Documents/phd/st'
# output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch12-1205-no-counts-shift2')
# geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS')
# 
# signature_type <<- 'msigdb' # c('msigdb', 'custom')
# custom_sign_path <<- file.path(proj_dir, 'geomx-processing', 'data', 'signatures',
#                                'ct_markers.csv')
# signature_name <<- ifelse(signature_type == 'custom', gsub('.csv', '', basename(custom_sign_path)), '')
# 
# dge_name <- 'dge_within_slide_Segment_bin_FALSE__NACT_status'
# 
# # dge_dir_path <- file.path(output_dir, 'dge', dge_name)
# # dge_df_path <- file.path(dge_dir_path, paste0('dge_', dge_inp_data, '_', dge_name, '.csv'))
# dge_dir_path <- '/home/iganiemi/Documents/phd/st/geomx-processing/results/batch12-1205-no-counts-shift2/dge/dge_within_slide_Segment_bin_FALSE__NACT_status'
# dge_df_path <- '/home/iganiemi/Documents/phd/st/geomx-processing/results/batch12-1205-no-counts-shift2/dge/dge_within_slide_Segment_bin_FALSE__NACT_status/dge_all_dge_within_slide_Segment_bin_FALSE__NACT_status.csv'
# 
# dge_inp_data <- 'dge_all'  # as for dge in deconv 'deconv_Macrophages itp'
# 
# 
# 

################
################
# select thr
fc_thr <- 0.5
pval_thr <- 0.05 # for DEG genes
gsea_padj_thr <- 0.01 # for GSEA results

# jaccard idx hclust cuts 
# hlust dendrogram height cut for clustering pathways by jaccard idx
jaccard_hclust_cuts <- c(0.5, 1, 1.2, 1.5)

adj_synonym <- T # whether or not adjust synonyms genes
# around 300 genes can be rescued this way but ensembl does not always work
# if there are issues, turn it off

# signatures and DEG results with less nr of genes will be removed
min_sign_gene_nr <- 5
compute_hallmark <- T # should GSEA for msigdb hallmark be computed
msigdb_subcat <- c('CP:BIOCARTA', 'CP:KEGG','GO:BP')

source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))

# prepare signatures list -------------------------------------------------

geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path) # only mtx with rownames needed

if(signature_type == 'msigdb'){
  # signatures from all Hallmark + selected CP from msigDB 
  sign_list <- prepare_msigdb_sign_list(adjust_synonym = adj_synonym, geomx_obj = geomx_obj, hal = compute_hallmark, 
                                        db_subcat_list = msigdb_subcat)
  out_name <- 'msigdb'
} else if(signature_type == 'custom'){
  # signatures from custom file
  sign_list <- prepare_custom_sign_list(fread(custom_sign_path), adjust_synonym = F,
                                        geomx_obj = geomx_obj)
  out_name <- paste0('custom_', gsub('//.csv', '', basename(custom_sign_path)))
} else{
  stop("signature_type parameter can only be 'msigb' or 'custom'")
}

sign_list <- sign_list[sapply(sign_list, length) >= min_sign_gene_nr]


# load dge files ----------------------------------------------------------

dge_dir_path <- file.path(output_dir, 'dge', dge_name)
dge_df_list <- list.files(dge_dir_path, pattern = paste0(dge_name, '.csv'), full.names = T)

# loop through all dge results
for(dge_df_path in dge_df_list){
  
  dge_inp_data <- gsub(paste0( '_',dge_name, '.csv'), '', basename(dge_df_path))
  print(paste0('##### ', dge_inp_data, ' #####'))
  
  # filter DGE results ------------------------------------------------------
  
  dge_df <- fread(dge_df_path)
  
  # filter to significant
  dge_df_sig <- filter(dge_df, FDR <= pval_thr & (Estimate >= fc_thr | Estimate <= -fc_thr))
  dge_df_sig$data_group <- ifelse(is.na(dge_df_sig$data_group), 'onegroup', dge_df_sig$data_group)
  
  # needed for ORA
  #dge_df_sig <- left_join(dge_df_sig, gene_entrez_universe, by = c('Gene' = 'external_gene_name'))
  
  
  # GSEA with fgsea ---------------------------------------------------------
  
  gsea_res_all <- lapply(unique(dge_df_sig$Contrast), function(cont){
    lapply(unique(dge_df_sig$data_group), function(dt_group){
      
      # subset to data and contrast group
      dge_sub <- dge_df_sig[dge_df_sig$data_group == dt_group & dge_df_sig$Contrast == cont, ]
      
      # continue if > min_sign_gene_nr genes
      if(nrow(dge_sub) >= min_sign_gene_nr){
        
        # different from DGE calculation cause fdr was for all segments and contrasts. this is appropriate one
        dge_sub$FDR_adj <- p.adjust(dge_sub$`Pr(>|t|)`, method = "fdr") 
        dge_sub$rank_p_fcval <- dge_sub$Estimate*(-log10(dge_sub$`Pr(>|t|)`))
        
        # sort by ranking by fc and pval (because fdr has many ties)
        dge_sub_rank <- dge_sub$rank_p_fcval
        names(dge_sub_rank) <- dge_sub$Gene
        dge_sub_rank <- sort(dge_sub_rank, decreasing = T)
        
        # fix infinite ranks if needed
        # Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
        max_ranking <- max(dge_sub_rank[is.finite(dge_sub_rank)])
        min_ranking <- min(dge_sub_rank[is.finite(dge_sub_rank)])
        dge_sub_rank <- replace(dge_sub_rank, dge_sub_rank > max_ranking, max_ranking * 10)
        dge_sub_rank <- replace(dge_sub_rank, dge_sub_rank < min_ranking, min_ranking * 10)
        dge_sub_rank <- sort(dge_sub_rank, decreasing = TRUE) # sort genes by ranking
        
        # do gsea on ranked dge gene list
        gsea_res <- fgsea(pathways = sign_list, # List of gene sets to check
                          stats = dge_sub_rank,
                          scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                          minSize = 10,
                          maxSize = 500,
                          nproc = 18) # for parallelisation
        
        # padj is stochastic - different runs on the same data give slightly different results
        gsea_res <- arrange(gsea_res, padj) %>%
          filter(padj <= gsea_padj_thr) 
        
        # assign independent pathways, removing redundancies/similar pathways
        collapsedPathways <- collapsePathways(gsea_res, sign_list, dge_sub_rank)
        
        gsea_res$is_main_pathway <- ifelse(gsea_res$pathway %in% collapsedPathways$mainPathways, 'yes', 'no')
        
        gsea_res$data_group <- dt_group
        gsea_res$Contrast <- cont
        
        if(nrow(gsea_res) > 0){
          return(gsea_res)
        } else{
          return()
        }

      } else{
        return()
      }
    })
  })
  
  gsea_res_all <- do.call(rbind, unlist(gsea_res_all, recursive=FALSE))
  
  if(!is.null(gsea_res_all)){
    fwrite(gsea_res_all, file.path(dge_dir_path, paste0('gsea_dge_', signature_type,
                                                        '_', signature_name, '_', dge_inp_data,
                                                        '_fc', as.character(fc_thr),'.csv')))
    
    
    # cluster gsea signatures by jaccard idx ----------------------------------
    # clustering based on jaccard idx - nr of common elements in a set / union of sets
    
    gsea_res_clust_all <- lapply(unique(gsea_res_all$Contrast), function(cont){
      lapply(unique(gsea_res_all$data_group), function(dt_group){
        
        # subset to dt group + contrast
        gsea_subset <- gsea_res_all[gsea_res_all$data_group == dt_group & gsea_res_all$Contrast == cont, ]
        gsea_subset_name <- gsub(' ', '', paste0(dt_group, '_', cont))
        
        if(signature_type == 'msigdb'){
          gsea_subset <- gsea_subset[gsea_subset$is_main_pathway == 'yes', ]
        }
        
        # get leading genes for each pathway
        paths_genes_list <- lapply(gsea_subset$leadingEdge, function(x){
          genelist <- unlist(strsplit(x, split='|', fixed=T))
        })
        
        names(paths_genes_list) <- gsea_subset$pathway
        
        # calculate jaccard score between each pathway leading gene set
        path_jaccard <- lapply(paths_genes_list, function(x){
          p1 <- lapply(paths_genes_list, function(y){
            jacc_idx <- as.numeric(round(length(intersect(x, y)) / length(union(x,y)), digits = 4))
          })
          return(unlist(p1))
        })
        
        path_jaccard_mtx <- do.call('cbind', path_jaccard)
        
        # make clustered heatmap
        nes_anno <- sapply(colnames(path_jaccard_mtx), function(path){
          nes <- ifelse(sign(gsea_subset$NES[gsea_subset$pathway == path]) == 1, 'pos', 'neg')
        })
        
        ha = HeatmapAnnotation(
          NES = anno_simple(nes_anno, col = c("pos" = "green", "neg" = "blue")),
          annotation_name_side = "left")
        
        png(filename=file.path(dge_dir_path, paste0('hmap_',gsea_subset_name, '_', signature_type,
                                                    '_', signature_name, '_', dge_inp_data,
                                                    '_fc', as.character(fc_thr), '.png')), 
            width=8, height=6,units="in",res=1000)
        
        condition_heat <- Heatmap(as.matrix(path_jaccard_mtx), border="white",
                                  rect_gp = gpar(col = "white", lwd = 2), column_title = gsea_subset_name,
                                  cluster_columns = T, cluster_rows= T, col = brewer.pal(5, "YlOrRd"),
                                  show_heatmap_legend = F, top_annotation = ha,
                                  row_names_gp = gpar(fontsize = 6),
                                  column_names_gp = gpar(fontsize = 6))
        
        draw(condition_heat)
        dev.off()
        
        #heatmap(path_jaccard_mtx)
        
        # clustering with hclust
        path_hclust <- hclust(dist(path_jaccard_mtx), method = "average")
        #plot(path_hclust, hang = -1, cex = 0.4)
        
        for(cutnr in jaccard_hclust_cuts){
          # cut the hclust tree at given point
          path_hclust_cut <- cutree(path_hclust, h = cutnr)
          
          # add subset name to cluster name
          path_hclust_cut <- paste0(gsea_subset_name, '_', as.character(path_hclust_cut))
          
          # merge with gsea result
          gsea_subset[[paste0('path_cluster_cut_', gsub('\\.', '', as.character(cutnr)))]] <- path_hclust_cut
          
        }
        
        return(gsea_subset)
        
      })
    })
    
    gsea_res_clust_all <- do.call(rbind, unlist(gsea_res_clust_all, recursive=FALSE))
    
    fwrite(gsea_res_clust_all, file.path(dge_dir_path, paste0('gsea_dge_clust_', signature_type,
                                                              '_', signature_name, '_', dge_inp_data,
                                                              '_fc', as.character(fc_thr), '.csv')))
  } else{
    print('no GSEA enrichment for this DEG list')
  }
  }

# make entrez gene universe -----------------------------------------------


# it just have to be expr mtx  for deconv sice only rownames are taken
# geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)
# 
# # convert to Entrez ID
# ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "https://useast.ensembl.org")
# 
# # bcg genes - all from geomx dataset
# gene_entrez_universe <- getBM(attributes=c('external_gene_name', 'entrezgene_id'),
#                               filters = 'external_gene_name',
#                               values = rownames(geomx_obj),
#                               mart = ensembl)
# 
# #rmv duplicates
# gene_entrez_universe <- gene_entrez_universe[!duplicated(gene_entrez_universe$external_gene_name),]
# 
# 
# # TODO rmv just for testing when ensembl does not work
# fwrite(gene_entrez_universe, file.path(proj_dir, 'geomx-processing', 'data', 'signatures',
#                                        'enterz_universe_all.csv'))
# gene_entrez_universe <- fread(file.path(proj_dir, 'geomx-processing', 'data', 'signatures',
#                                         'enterz_universe_all.csv'))

# ORA ---------------------------------------------------------------------

# ORA on GO
# go_res <- as.data.frame(enrichGO(gene = as.character(unlist(dge_sub$entrezgene_id)),
#                                  ont = "BP",
#                                  OrgDb ="org.Hs.eg.db",
#                                  universe = as.character(unlist(gene_entrez_universe$entrezgene_id)),
#                                  readable=TRUE,
#                                  pvalueCutoff = pval_thr))


# #TODO maybe it should be done separately for upregulated and downregulated?
# # ORA on H+CP
# msigdb_df <- msigdbr(species = "Homo sapiens")
# 
# # TODO loop
# dge_sub_pos <- dge_sub$entrezgene_id[dge_sub$Estimate > 0]
# dge_sub_neg <- dge_sub$entrezgene_id[dge_sub$Estimate < 0]
# 
# 
# msigdb_res <- enricher(
#   gene = as.character(unlist(dge_sub_pos)),
#   pvalueCutoff = pval_thr, # Can choose a FDR cutoff
#   pAdjustMethod = "BH", 
#   universe = as.character(unlist(gene_entrez_universe$entrezgene_id)), 
#   TERM2GENE = dplyr::select(msigdb_df, gs_name, entrez_gene)
# )
# 
# msigdb_res <- data.frame(msigdb_res@result) %>%
#   filter(p.adjust <= pval_thr)


# inspect GSEA results ----------------------------------------------------


# anno_stroma <- fread(file.path("/media/iganiemi/T7-iga/st/geomx-processing/results/nact/dge/gsea/gsea_clust/pos",
#                                "gsea_post_anno_rank_p_fcval_clust_stroma.csv"))
# 
# anno_tumor <- fread(file.path("/media/iganiemi/T7-iga/st/geomx-processing/results/nact/dge/gsea/gsea_clust/pos",
#                               "gsea_post_anno_rank_p_fcval_clust_tumor.csv"))
# 
# pfi_stroma <- fread(file.path("/media/iganiemi/T7-iga/st/geomx-processing/results/nact/dge/gsea/gsea_clust/neg",
#                               "gsea_post_pfi_rank_p_fcval_clust_stroma.csv"))
# 
# pfi_tumor <- fread(file.path("/media/iganiemi/T7-iga/st/geomx-processing/results/nact/dge/gsea/gsea_clust/neg",
#                              "gsea_post_pfi_rank_p_fcval_clust_tumor.csv"))
# 
# gsea_subset_clust <- pfi_stroma
# gsea_subset_clust <- arrange(gsea_subset_clust, desc(NES))
# gsea_subset_clust <- arrange(gsea_subset_clust, NES)
# 
# #gsea_subset_clust <- gsea_subset_clust[grepl('GOBP', gsea_subset_clust$pathway)] 
# 
# #keep 1st pathway from cluster
# gsea_subset_clust_1 <- gsea_subset_clust[!(duplicated(gsea_subset_clust$path_cluster_cut_1)), ]
# gsea_subset_clust_12 <- gsea_subset_clust[!(duplicated(gsea_subset_clust$path_cluster_cut_12)), ]
# gsea_subset_clust_15 <- gsea_subset_clust[!(duplicated(gsea_subset_clust$path_cluster_cut_15)), ]
# 
# 
# # for ++ vs each other ROI type separately
# # analyse which pathways are upregulated in 3 and 2
# sharing <- stack(sapply(unique(gsea_subset_clust$pathway), function(x){
#   gsea_p <- gsea_subset_clust$Contrast[gsea_subset_clust$pathway == x]
#   
#   if(length(gsea_p) == 3){
#     y <- 'shared_3'
#   } else if(length(gsea_p) == 2){
#     if(!("negCD8_negIBA1 - posCD8_posIBA1" %in% gsea_p)){
#       y <- 'shared_2_no_doubleneg'
#     } else if(!("negCD8_posIBA1 - posCD8_posIBA1" %in% gsea_p)){
#       y <- 'shared_2_no_negpos'
#     } else if(!("posCD8_negIBA1 - posCD8_posIBA1" %in% gsea_p)){
#       y <- 'shared_2_no_negpos'
#     }
#   } else{
#     y <- paste0('shared_1_', gsea_p)
#   }
#   return(y)
# }))
# 
# gsea_subset_fin <- gsea_subset_clust_1[gsea_subset_clust_1$pathway %in% sharing$ind[sharing$values == 'shared_3'] &
#                                          grepl('posCD8_posIBA1', gsea_subset_clust_1$Contrast), ]
# 
# gsea_subset_fin <- arrange(gsea_subset_fin, desc(NES))
# 
# #TODO filter for 1st appearance on pathway - fir different comparison it may had different
# # leading edge genes and haven't been clustered before
# 
# fwrite(gsea_subset_fin, file.path("/media/iganiemi/T7-iga/st/geomx-processing/results/nact/dge/gsea/gsea_clust/fin",
#                                   "gsea_post_anno_rank_p_fcval_clust_tumor_pos_pospos_shared_fin.csv"))

############################################################################
###########################################################################


# clustering with dbscan
# dbscan::kNNdistplot(path_jaccard_mtx, k =  5)
# 
# dbscan_res <- dbscan::dbscan(path_jaccard_mtx, eps = 1.5)
# 
# path_clust <- dbscan_res$cluster
# names(path_clust) <- colnames(path_jaccard_mtx)
# 
# cl0 <- names(path_clust[path_clust == 0])
# cl1 <- names(path_clust[path_clust == 1])
# cl2 <- names(path_clust[path_clust == 2])

