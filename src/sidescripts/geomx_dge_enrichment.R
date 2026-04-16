# unbiased GSEA/ORA on DGE results - to be added to main pipeline !!
# devtools::install_github("Evotec-Bioinformatics/evoGO")
library(fgsea)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(evoGO)

# set variables -----------------------------------------------------------

# proj_dir <<- '~/Documents/phd/st'
# output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch12-1205-no-counts-shift2')
# geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS')
# 
# signature_type <<- 'msigdb' # c('msigdb', 'custom')
# custom_sign_path <<- file.path(proj_dir, 'geomx-processing', 'data', 'signatures',
#                                'ct_markers.csv')
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
qval_thr <- 0.2 # for GO
gsea_padj_thr <- 0.05 # for GSEA results

# jaccard idx hclust cuts 
# hlust dendrogram height cut for clustering pathways by jaccard idx
jaccard_hclust_cuts <- c(0.5, 1, 1.2, 1.5)

adj_synonym <- T # whether or not adjust synonyms genes
# around 300 genes can be rescued this way but ensembl does not always work
# if there are issues, turn it off

# signatures and DEG results with less nr of genes will be removed
min_sign_gene_nr <- 10
#compute_hallmark <- T # should GSEA for msigdb hallmark be computed
#msigdb_subcat <- c('CP:BIOCARTA', 'CP:KEGG_MEDICUS','GO:BP')
#msigdb_subcat <- c('GO:BP', 'CP:KEGG_MEDICUS')

scrna_ref_cleaned_path <- file.path(output_dir, 'deconvolution', gsub('.RDS', '_cleaned_for_deconv.RDS', basename(scrna_ref_path)))

source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))

dge_dir_path <- file.path(output_dir, 'dge', dge_name)
dir.create(file.path(dge_dir_path, "gsea_enrichment"))

dge_df_list <- list.files(dge_dir_path, pattern = paste0(dge_name, '.csv'), full.names = T)

# download the latest version of go annot
# TODO only for evoGO, also may not work bc of ensembl
# goAnnotation <- getGOAnnotation("hsapiens")
# goAnnotation <- loadGOAnnotation("hsapiens")

# prepare signatures list -------------------------------------------------

geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path) # only mtx with rownames needed

if(signature_type == 'msigdb'){
  # signatures from all Hallmark + selected CP from msigDB 
  # sign_list <- prepare_msigdb_sign_list(adjust_synonym = adj_synonym, geomx_obj = geomx_obj, hal = compute_hallmark, 
  #                                       db_subcat_list = msigdb_subcat)
  
  sign_list <- prepare_msigdb_sign_list(adjust_synonym = adj_synonym, geomx_obj = geomx_obj, msigdb_subcat = msigdb_subcat)
  out_name <- paste0('msigdb', '_', msigdb_subcat)
} else if(signature_type == 'custom'){
  # signatures from custom file
  sign_list <- prepare_custom_sign_list(fread(custom_sign_path), adjust_synonym = F,
                                        geomx_obj = geomx_obj)
  out_name <- paste0('custom_', gsub('//.csv', '', basename(custom_sign_path)))
} else{
  stop("signature_type parameter can only be 'msigb' or 'custom'")
}

sign_list <- sign_list[sapply(sign_list, length) >= min_sign_gene_nr]

print(paste0(length(sign_list), ' signatures will be used'))

# load dge files ----------------------------------------------------------

# loop through all dge results
for(dge_df_path in dge_df_list){
  
  dge_inp_data <- gsub(paste0( '_',dge_name, '.csv'), '', basename(dge_df_path))
  print(paste0('##### ', dge_inp_data, ' #####'))
  
  dir.create(file.path(dge_dir_path, "gsea_enrichment", dge_inp_data))
  
  # read DGE results --------------------------------------------------------
  
  dge_df <- fread(dge_df_path)
  dge_df$data_group <- ifelse(is.na(dge_df$data_group), 'onegroup', dge_df$data_group) # add to avoid bugs
  
  # GSEA with fgsea ---------------------------------------------------------
  
  gsea_res_all <- lapply(unique(dge_df$Contrast), function(cont){
    lapply(unique(dge_df$data_group), function(dt_group){
      
      # subset to data and contrast group
      dge_sub <- dge_df[dge_df$data_group == dt_group & dge_df$Contrast == cont, ]
      
      # rank by log2fc * -log10(pval)
      # TODO check if collapsepathways would work with custom signatures
      gsea_res <- rank_genes_and_do_gsea_enrichment(dge_sub, 'Estimate', 'Pr(>|t|)', 'Gene', sign_list)
      
      gsea_res$data_group <- dt_group
      gsea_res$Contrast <- cont
      
      if(nrow(gsea_res) > 0){
        return(gsea_res)
      } else{
        return()
      }
    })
  })
  
  gsea_res_all <- do.call(rbind, unlist(gsea_res_all, recursive=FALSE))
  
  if(!is.null(gsea_res_all)){
    fwrite(gsea_res_all, file.path(dge_dir_path, 'gsea_enrichment', dge_inp_data,
                                   paste0('gsea_dge_', out_name, 
                                          '_', dge_inp_data, '_fc', as.character(fc_thr),'_nofiltering.csv')))
    
    
    # cluster gsea signatures by jaccard idx ----------------------------------
    # clustering based on jaccard idx - nr of common elements in a set / union of sets
    
    gsea_res_clust_all <- lapply(unique(gsea_res_all$Contrast), function(cont){
      lapply(unique(gsea_res_all$data_group), function(dt_group){
        
        # subset to dt group + contrast
        gsea_subset <- gsea_res_all[gsea_res_all$data_group == dt_group & gsea_res_all$Contrast == cont, ]
        gsea_subset_name <- gsub(' ', '', paste0(dt_group, '_', cont))
        print(gsea_subset_name)
        
        # filter to padj
        # padj is stochastic - different runs on the same data give slightly different results
        gsea_subset <- arrange(gsea_subset, padj) %>%
          filter(padj <= gsea_padj_thr) 
        
        if(nrow(gsea_subset) > 0){
          if(signature_type == 'msigdb'){
            gsea_subset <- gsea_subset[gsea_subset$is_main_pathway == 'yes', ]
          }
          
          # only clustering more than 1 pathways makes sense
          if(nrow(gsea_subset) > 1){
            hmap_outpath <- file.path(dge_dir_path, 'gsea_enrichment',dge_inp_data,
                                      paste0('hmap_',gsea_subset_name, '_', out_name, '_', dge_inp_data, 
                                             '_fc', as.character(fc_thr), '.png'))
            
            # cluster pathways by jaccard idx and make heatmap
            gsea_subset <- cluster_gsea_enrichment(gsea_subset, 'leadingEdge', 'pathway', 
                                                  hmap_outpath = hmap_outpath, hmap_title = gsea_subset_name, 
                                                  lead_genes_split = '|')
          }
          
          return(gsea_subset)
        } else{
          return()
        }
      })
    })
    
    gsea_res_clust_all <- do.call(rbind.fill, unlist(gsea_res_clust_all, recursive=FALSE))
    
    
    if(!is.null(gsea_res_clust_all)){
    fwrite(gsea_res_clust_all, file.path(dge_dir_path, 'gsea_enrichment', dge_inp_data,
                                         paste0('gsea_dge_clust_', out_name,
                                                '_', dge_inp_data, '_fc', as.character(fc_thr), '.csv')))
    }
    
  } else{
    print('no GSEA enrichment for this DEG list')
  }
}

# make entrez gene universe -----------------------------------------------

# TODO save it during dge script and load now
# geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)
# scrna_ref_obj <- readRDS(scrna_ref_cleaned_path)
# 
# geomx_filt <- remove_low_complex_and_noncoding_genes(geomx_obj, scrna_ref_obj, raw_counts_layer = 'counts')
# 
# # convert to Entrez ID
# ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "https://useast.ensembl.org")
# 
# # bcg genes - all from geomx dataset
# gene_entrez_universe <- getBM(attributes=c('external_gene_name', 'entrezgene_id'),
#                               filters = 'external_gene_name',
#                               values = rownames(geomx_filt),
#                               mart = ensembl)
# 
# #rmv duplicates
# gene_entrez_universe <- gene_entrez_universe[!duplicated(gene_entrez_universe$external_gene_name),]
# 
# rm(geomx_obj)
# rm(scrna_ref_obj)
# rm(geomx_filt)

# TODO rmv just for testing when ensembl does not work
# fwrite(gene_entrez_universe, file.path(proj_dir, 'geomx-processing', 'data', 'signatures',
#                                        'enterz_universe_geomx_pc.csv'))
gene_entrez_universe <- fread(file.path(proj_dir, 'geomx-processing', 'data', 'signatures',
                                        'enterz_universe_geomx_pc.csv'))


# do GO enrichment --------------------------------------------------------

dge_df_path <- dge_df_list[1]
cont <- "Bcell_domin - CD8_Macro_domin"
dt_group <- 'stroma_post'

dir.create(file.path(dge_dir_path, 'go_enrichment'))

# loop through all dge results
for(dge_df_path in dge_df_list){
  
  dge_inp_data <- gsub(paste0( '_',dge_name, '.csv'), '', basename(dge_df_path))
  print(paste0('##### ', dge_inp_data, ' #####'))
  
  dir.create(file.path(dge_dir_path, 'go_enrichment', dge_inp_data))
  
  # read DGE results --------------------------------------------------------
  
  dge_df <- fread(dge_df_path)
  dge_df$data_group <- ifelse(is.na(dge_df$data_group), 'onegroup', dge_df$data_group) # add to avoid bugs
  
  # GO enrichment with clusterprofiler --------------------------------------

  go_res_all <- lapply(unique(dge_df$Contrast), function(cont){
    lapply(unique(dge_df$data_group), function(dt_group){
      
      # subset to data and contrast group
      dge_sub <- dge_df[dge_df$data_group == dt_group & dge_df$Contrast == cont, ]
      
      dge_subset_name <- gsub(' ', '', paste0(dt_group, '_', cont))
      print(dge_subset_name)
      
      # subset to differential genes using set up thresholds
      dge_sub_signif_pos <- dge_sub[dge_sub$Estimate >= fc_thr & dge_sub$`Pr(>|t|)` <= pval_thr, ]
      dge_sub_signif_neg <- dge_sub[dge_sub$Estimate <= -fc_thr & dge_sub$`Pr(>|t|)` <= pval_thr, ]
      
      dge_sub_signif_list <- list(pos = dge_sub_signif_pos, neg = dge_sub_signif_neg)
      
      lapply(1:length(dge_sub_signif_list), function(i){
        dge_sub_signif <- dge_sub_signif_list[[i]]
        dge_signif_name <- names(dge_sub_signif_list)[i]
        
        if(nrow(dge_sub_signif) >= min_sign_gene_nr){
          
          # perform GO enrichment
          go_res_obj <- enrichGO(gene = dge_sub_signif$Gene,
                                 OrgDb = org.Hs.eg.db,
                                 keyType = "SYMBOL",
                                 ont = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = pval_thr,
                                 qvalueCutoff = qval_thr,
                                 readable = T,
                                 universe = gene_entrez_universe$external_gene_name,
                                 minGSSize = min_sign_gene_nr)
          
          go_res <- go_res_obj@result
          go_res <- go_res[go_res$p.adjust <= pval_thr, ]
          
          if(nrow(go_res) > 0){
            
            print(paste0(dge_signif_name, " - nr of significant go terms: ", as.character(nrow(go_res))))
            
            # make visualisation plot
            png(filename=file.path(dge_dir_path, 'go_enrichment', dge_inp_data,
                                   paste0('dotplot_go_', dge_inp_data, '_', dge_subset_name, '_', dge_signif_name,
                                                                         '_fc', as.character(fc_thr), '_pval', as.character(pval_thr), '.png')), 
                width=12, height=6,units="in",res=1000)
            
            plot(dotplot(go_res_obj, x = "GeneRatio", color = "p.adjust", title = paste0("Top 15 of GO Enrichment", dge_signif_name),
                    showCategory = 15, label_format = 80))
            dev.off()
            
            go_res$data_group <- dt_group
            go_res$Contrast <- cont
            go_res$direction <- dge_signif_name
            
            return(go_res)
          } else{
            return()
          }
        } else{
          return()
        }
      })
    })
  })
  
  # combine into 1 df
  go_res_all_flat <- list_flatten(list_flatten(go_res_all)) # flatten to single list of dfs
  go_res_all_flat <- keep(go_res_all_flat, ~ is.data.frame(.x)) # keep only not-null
  go_res_all_fin <- do.call(rbind, go_res_all_flat)
  
  if(!is.null(go_res_all_fin)){
    fwrite(go_res_all_fin, file.path(dge_dir_path, 'go_enrichment',
                                   paste0('go_dge_', dge_inp_data, '_fc', as.character(fc_thr),'_pval', as.character(pval_thr), '.csv')))
  } else{
    print('no GO enrichment for this DEG list')
  }
}


# do KEGG enrichment analysis ---------------------------------------------

dge_df_path <- dge_df_list[1]
cont <- "Bcell_domin - CD8_Macro_domin"
dt_group <- 'stroma_post'

dir.create(file.path(dge_dir_path, 'kegg_enrichment'))

# loop through all dge results
for(dge_df_path in dge_df_list){
  
  dge_inp_data <- gsub(paste0( '_',dge_name, '.csv'), '', basename(dge_df_path))
  print(paste0('##### ', dge_inp_data, ' #####'))
  
  dir.create(file.path(dge_dir_path, 'kegg_enrichment', dge_inp_data))
  
  # read DGE results --------------------------------------------------------
  
  dge_df <- fread(dge_df_path)
  dge_df$data_group <- ifelse(is.na(dge_df$data_group), 'onegroup', dge_df$data_group) # add to avoid bugs
  
  
  # convert to Entrez ID
  ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "https://useast.ensembl.org")
  
  # bcg genes - all from geomx dataset
  dge_entrez <- getBM(attributes=c('external_gene_name', 'entrezgene_id'),
                      filters = 'external_gene_name',
                      values = dge_df$Gene,
                      mart = ensembl)
  
  #rmv duplicates
  dge_entrez <- dge_entrez[!duplicated(dge_entrez$external_gene_name),]
  
  dge_df <- left_join(dge_df, dge_entrez, by = c('Gene' = 'external_gene_name'))
  
  # GO enrichment with clusterprofiler --------------------------------------
  
  kegg_res_all <- lapply(unique(dge_df$Contrast), function(cont){
    lapply(unique(dge_df$data_group), function(dt_group){
      
      # subset to data and contrast group
      dge_sub <- dge_df[dge_df$data_group == dt_group & dge_df$Contrast == cont, ]
      
      dge_subset_name <- gsub(' ', '', paste0(dt_group, '_', cont))
      print(dge_subset_name)
      
      # subset to differential genes using set up thresholds
      dge_sub_signif_pos <- dge_sub[dge_sub$Estimate >= fc_thr & dge_sub$`Pr(>|t|)` <= pval_thr, ]
      dge_sub_signif_neg <- dge_sub[dge_sub$Estimate <= -fc_thr & dge_sub$`Pr(>|t|)` <= pval_thr, ]
      
      dge_sub_signif_list <- list(pos = dge_sub_signif_pos, neg = dge_sub_signif_neg)
      
      lapply(1:length(dge_sub_signif_list), function(i){
        dge_sub_signif <- dge_sub_signif_list[[i]]
        dge_signif_name <- names(dge_sub_signif_list)[i]
        
        if(nrow(dge_sub_signif) >= min_sign_gene_nr){
          # perform KEGG enrichment
          kegg_res_obj <- enrichKEGG(dge_sub_signif$entrezgene_id, 
                                     organism = "hsa", 
                                     keyType = "kegg", 
                                     pvalueCutoff = pval_thr, 
                                     pAdjustMethod = "BH", 
                                     universe = as.character(gene_entrez_universe$entrezgene_id), 
                                     minGSSize = min_sign_gene_nr, 
                                     qvalueCutoff = qval_thr)

          kegg_res <- kegg_res_obj@result
          kegg_res <- kegg_res[kegg_res$p.adjust <= pval_thr, ]
          
          if(nrow(kegg_res) > 0){
            
            print(paste0(dge_signif_name, " - nr of significant kegg terms: ", as.character(nrow(kegg_res))))
            
            # make visualisation plot
            png(filename=file.path(dge_dir_path, 'kegg_enrichment', dge_inp_data,
                                   paste0('dotplot_kegg_', dge_inp_data, '_', dge_subset_name, '_', dge_signif_name,
                                          '_fc', as.character(fc_thr), '_pval', as.character(pval_thr), '.png')), 
                width=12, height=6,units="in",res=1000)
            
            plot(dotplot(kegg_res_obj, x = "GeneRatio", color = "p.adjust", title = paste0("Top 15 of KEGG Enrichment DGE ", dge_signif_name),
                         showCategory = 15, label_format = 80))
            dev.off()
            
            kegg_res$data_group <- dt_group
            kegg_res$Contrast <- cont
            kegg_res$direction <- dge_signif_name
            
            return(kegg_res)
          } else{
            return()
          }
        } else{
          return()
        }
      })
    })
  })
  
  # combine into 1 df
  kegg_res_all_flat <- list_flatten(list_flatten(kegg_res_all)) # flatten to single list of dfs
  kegg_res_all_flat <- keep(kegg_res_all_flat, ~ is.data.frame(.x)) # keep only not-null
  kegg_res_all_fin <- do.call(rbind, kegg_res_all_flat)
  
  if(!is.null(kegg_res_all_fin)){
    fwrite(kegg_res_all_fin, file.path(dge_dir_path, 'kegg_enrichment',
                                     paste0('kegg_dge_', dge_inp_data, '_fc', as.character(fc_thr),'_pval', as.character(pval_thr), '.csv')))
  } else{
    print('no KEGG enrichment for this DEG list')
  }
}



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

