#TODO check if all needed

# WARNING: DGE with mixed model will take around 30G RAM
# best to run in >10 cores
# 

# problems with Matrix package - i has to be lower that 1.7 to work with lmer
#devtools::install_version("Matrix","1.6.4")

library(data.table)
#library(clusterProfiler)
# library(msigdbr)
# library(progeny)
# library(biomaRt)
# library(GSVA)
# library(ggpubr)
# library(topGO)
# library(fgsea)
# library(fpc)
# library(dbscan)


# get variables -----------------------------------------------------------
# data_dir <- '~/Documents/phd/st/data/geomx/geomx_batch2_1124/'
# output_dir <- '~/Documents/phd/st/geomx-processing/results/batch2'
# 


comparison_type <- 'within' 
# 'within' when you compare different ROI types within sample
# between - comparisons between slides

cofounder_name <- 'Sample' # don't change it
# then 'Sample' is added as a cofounder (random intercept in LLM model)

main_var_name <- "Annotation_cell" 
# main_var - main variable to make comparison 

main_var_is_bin <- TRUE 
# if main_var_is_bin is True, main_var_main_val will be compared 
# with all other categories in main_var
# if False - each category in main_var will be compared with every other one
# main_var_main_val is set to NULL

main_var_main_val <- 'CD4_CD8_CD11_Iba1' 
# main_var_main_val value among main_var which needs to be compare against all other vals
# or part of the value eg 'CD8' within values for comparison


dge_categories <- c('Segment', 'NACT status')
# dge_categories - all conditions for which we want to make DGE separately

################
# comparison_type <- 'within'
# cofounder_name <- 'Sample'
# main_var_name <- "Annotation_cell" 
# main_var_is_bin <- TRUE 
# main_var_main_val <- 'CD4_CD8_CD11_Iba1' 
# dge_categories <- c('Segment', 'NACT status')
################
comparison_type <- 'between'
cofounder_name <- 'Sample'
main_var_name <- "NACT status"
main_var_is_bin <- FALSE
dge_categories <- c('Segment', 'Annotation_cell')
###############

norm_type <- 'q3_norm' # either 'q3_norm', 'quant_norm' or 'deseq2_norm'

multicore = TRUE # if Linux or macOS, for Windows multicore = FALSE

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'dge'), showWarnings = T, recursive = T)

#source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')

# load geomx obj from rds -------------------------------------------------

geomx_obj <- readRDS(geomx_norm_path)

# convert normalized counts to log scale
assayDataElement(object = geomx_obj, elt = paste0("log_", norm_type)) <-
  assayDataApply(geomx_obj, 2, FUN = log, base = 2, elt = norm_type)


# DGE with main variable comparison ---------------------------------------

if(main_var_is_bin){
  # make binary vector - either main variable has the desired value or not
  pData(geomx_obj)$main_var <- ifelse(grepl(main_var_main_val, pData(geomx_obj)[, main_var_name]),
                                          main_var_main_val, 'other_roi_type')
} else{
  pData(geomx_obj)$main_var <- pData(geomx_obj)[, main_var_name]
}


# convert test variables to factors
for(col in c(dge_categories, 'main_var')){
  pData(geomx_obj)[[paste0(col, "_factor")]] <- factor(pData(geomx_obj)[[col]])
}

pData(geomx_obj)$cofounder_factor <- factor(pData(geomx_obj)[[cofounder_name]])

# make variable with all dge categories
pData(geomx_obj)$dge_group <- apply(pData(geomx_obj), 1, function(row){
  group <- sapply(dge_categories, function(var){
    paste(row[var])
  })
  group <- paste(group, collapse = '_')
  return(group)
})

# create formula for the LLM model:
# Sample is used as a mixed effect (cofounder)
if(comparison_type == 'within'){
  # within slide analysis - with random slope in LLM
  model_formula <- ~ main_var_factor + (1 + main_var_factor | cofounder_factor) # random slope + random intercept
} else if(comparison_type == 'between'){
  model_formula <- ~ main_var_factor + (1 | cofounder_factor) # random intercept
} else{stop('comparison type can be either "within" or "between"')}


# run LMM:
# formula follows conventions defined by the lme4 package

dge_results <- c()
for(data_group in unique(pData(geomx_obj)[, 'dge_group'])){
  
  print(data_group)
  ind <- geomx_obj@phenoData@data$dge_group == data_group
  
  mixed_result <- tryCatch({
    mixedOutmc <- mixedModelDE(
      geomx_obj[, ind],
      elt = paste0("log_", norm_type),
      modelFormula = model_formula, 
      groupVar = 'main_var_factor',
      nCores = (parallel::detectCores() - 2),
      multiCore = multicore
    )
    mixedOutmc  # Return the result of mixedModelDE
  }, error = function(e) {
    # Return an empty dataframe if an error occurs eg to little ROIs
    data.frame()
    
  })
  
  gc()
  
  if(nrow(mixed_result) > 1){
    # format results as data.frame
    r_test <- do.call(rbind, mixed_result["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests
    
    # use lapply in case you have multiple levels of your test factor to
    # correctly associate gene name with it's row in the results table
    r_test$Gene <-
      unlist(lapply(colnames(mixed_result),
                    rep, nrow(mixed_result["lsmeans", ][[1]])))
    r_test$data_group <- data_group
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
    r_test <- r_test[, c("Gene", "data_group",  "Contrast", "Estimate",
                         "Pr(>|t|)", "FDR")]
    dge_results <- rbind(dge_results, r_test)
  } else{
    print('error while computing dge. probably too little ROI for comparison')
    dge_results <- dge_results
  }
  
}

fwrite(dge_results, file.path(output_dir, 'dge', 
                          paste('dge', comparison_type, 'slide', main_var_name, 
                                 'bin', main_var_is_bin, paste0(dge_categories, collapse = '_'), 
                                sep = '_')))

# enrichment on DGE -------------------------------------------------------

dge_data_dir <- file.path("/media/iganiemi/T7-iga/st/geomx-processing/results/nact/dge")

list.files(file.path(dge_data_dir, 'dge'))

# ROI type doublepos vs all other together pre post separately
dge_df <- fread(file.path(dge_data_dir, "dge",  "dge_annotation_dual_cell_pre_post_separately.csv"))
# doublepost post short vs long pfi
dge_df <- fread(file.path(dge_data_dir, "dge",  "dge_pfs_doublepos_updated.csv"))
# doublepos pre vs post
dge_df <- fread(file.path(dge_data_dir, "dge",  "dge_pre_post_doublepos.csv"))

# select thr
fc_thr <- 1
pval_thr <- 0.05
# filter to post and significant

dge_df_sig <- filter(dge_df, Subset == 'post' & FDR <= pval_thr & (Estimate >= fc_thr | Estimate <= -fc_thr))

table(dge_df$Contrast)

# subset to anno group
dge_sub <- dge_df_sig[dge_df_sig$Contrast == 'negCD8_negIBA1 - posCD8_posIBA1', ]

# convert to Entrez ID
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

# bcg genes - all from geomx dataset
gene_entrez_universe <- getBM(attributes=c('external_gene_name', 'entrezgene_id'),
                    filters = 'external_gene_name',
                    values = rownames(geomx_obj),
                    mart = ensembl)

#rmv duplicates
gene_entrez_universe <- gene_entrez_universe[!duplicated(gene_entrez_universe$external_gene_name),]

dge_sub <- left_join(dge_sub, gene_entrez_universe, by = c('Gene' = 'external_gene_name'))


# filter msigdb -----------------------------------------------------------

# Hallmark +CP
msigdb_df <- msigdbr(species = "Homo sapiens")
msigdb_df <- filter(msigdb_df, gs_cat %in% c("H", "C2", "C5") & !(gs_subcat %in% c("CGP", "GO:CC", "GO:MF", "HPO", "CP")))

msigdb_list <- lapply(unique(msigdb_df$gs_name), function(x){
  ls <- msigdb_df$gene_symbol[msigdb_df$gs_name == x]
})

names(msigdb_list) <- unique(msigdb_df$gs_name)

#TODO add T-cell exhaustion pathway here and rerun for previous DGE 

# ORA ---------------------------------------------------------------------

# ORA on GO
go_res <- as.data.frame(enrichGO(gene = as.character(unlist(dge_sub$entrezgene_id)),
                   ont = "BP",
                   OrgDb ="org.Hs.eg.db",
                   universe = as.character(unlist(gene_entrez_universe$entrezgene_id)),
                   readable=TRUE,
                   pvalueCutoff = pval_thr))


# ORA on H+CP
msigdb_res <- enricher(
  gene = as.character(unlist(dge_sub$entrezgene_id)),
  pvalueCutoff = pval_thr, # Can choose a FDR cutoff
  pAdjustMethod = "BH", 
  universe = as.character(unlist(gene_entrez_universe$entrezgene_id)), 
  TERM2GENE = dplyr::select(msigdb_df, gs_name, entrez_gene)
)

msigdb_res <- data.frame(msigdb_res@result) %>%
  filter(p.adjust <= pval_thr)

# GSEA with fgsea ---------------------------------------------------------

rank_type <- 'rank_p' # c('rank_p', 'rank_p_fcval', 'rank_fdr_fcval')

dge_df <- filter(dge_df, Subset == 'post') #for anno comparison

gsea_all <- lapply(unique(dge_df$Contrast), function(contrast){
  lapply(unique(dge_df$Segment), function(segment){
    
    dge_all_sub <- filter(dge_df, Contrast == contrast &
                            Segment == segment)
    
    # different from DGE calculation cause fdr was for all segments and contrasts. this is appropriate one
    dge_all_sub$FDR2 <- p.adjust(dge_all_sub$`Pr(>|t|)`, method = "fdr") 
    
     dge_all_sub$rank_p <- sign(dge_all_sub$Estimate)*(-log10(dge_all_sub$`Pr(>|t|)`))
     dge_all_sub$rank_p_fcval <- dge_all_sub$Estimate*(-log10(dge_all_sub$`Pr(>|t|)`))
     dge_all_sub$rank_fdr_fcval <- dge_all_sub$Estimate*(-log10(dge_all_sub$FDR2))
    
     # fdr rank alnone makes no sense - too many ties
     # plot(dge_all_sub$rank_p, dge_all_sub$rank_fdr) #linear + plateau - may be used interchangeably
     # plot(dge_all_sub$rank_p, dge_all_sub$rank_p_fcval) #linear + plateau - may be used interchangeably
     # plot(dge_all_sub$rank_fdr, dge_all_sub$rank_fdr_fcval) # linear
    
     dge_sub_rank <- dge_all_sub[[rank_type]]
     names(dge_sub_rank) <- dge_all_sub$Gene
     dge_sub_rank <- sort(dge_sub_rank, decreasing = T)
    
     plot(dge_sub_rank)
    
     # fix infinite ranks if needed
     # Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
     max_ranking <- max(dge_sub_rank[is.finite(dge_sub_rank)])
     min_ranking <- min(dge_sub_rank[is.finite(dge_sub_rank)])
     dge_sub_rank <- replace(dge_sub_rank, dge_sub_rank > max_ranking, max_ranking * 10)
     dge_sub_rank <- replace(dge_sub_rank, dge_sub_rank < min_ranking, min_ranking * 10)
     dge_sub_rank <- sort(dge_sub_rank, decreasing = TRUE) # sort genes by ranking
    
    
     gsea_res <- fgsea(pathways = msigdb_list, # List of gene sets to check
                       stats = dge_sub_rank,
                       scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                       minSize = 10,
                      maxSize = 500,
                       nproc = 18) # for parallelisation
    
    
     gsea_res <- arrange(gsea_res, padj) %>%
       filter(padj <= 0.01)
    
     # Select only independent pathways, removing redundancies/similar pathways
     collapsedPathways <- collapsePathways(gsea_res, msigdb_list, dge_sub_rank)
     mainPathways <- gsea_res[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
     gsea_res_main <- filter(gsea_res, pathway %in% mainPathways)
    
     gsea_res_main$Segment <- segment
     gsea_res_main$Contrast <- contrast

     # plotEnrichment(msigdb_list[[head(gsea_res[order(padj), ], 1)$pathway]],
     #                dge_sub_rank) +
     #   labs(title = head(gsea_res[order(padj), ], 1)$pathway)

    print(paste(contrast, segment, 'pass'))
    return(gsea_res_main)
  })
})


gsea_all <- do.call(rbind, unlist(gsea_all, recursive=FALSE))

fwrite(gsea_all, file.path(dge_data_dir, 'gsea', paste0('gsea_post_anno_dual_', rank_type, '.csv')))


# inspect GSEA results ----------------------------------------------------

# doublepos vs others
gsea_res_path <- "/media/iganiemi/T7-iga/st/geomx-processing/results/nact/dge/gsea/gsea_post_anno_dual_rank_p_fcval.csv"

# doublepos vs each other separately
gsea_res_path <- "/media/iganiemi/T7-iga/st/geomx-processing/results/nact/dge/gsea/gsea_post_anno_rank_p_fcval.csv"

# doublepos : pfi long vs short
gsea_res_path <- "/media/iganiemi/T7-iga/st/geomx-processing/results/nact/dge/gsea/gsea_post_pfi_rank_p_fcval.csv"

gsea_res <- fread(gsea_res_path)

# filter only for doublepositive comparisons
#gsea_res <- gsea_res[grepl('posCD8_posIBA1', gsea_res$Contrast), ]

# NES < 0 upregulated in ++
# NES > 0 downregulated in ++

intersect(gsea_res$pathway[gsea_res$NES > 0], gsea_res$pathway[gsea_res$NES < 0])

# TODO change +/- here
gsea_res <- gsea_res[gsea_res$NES < 0, ]

length(unique(gsea_res$pathway))

gsea_stroma <- gsea_res[gsea_res$Segment == 'stroma', ]
gsea_tumor <- gsea_res[gsea_res$Segment == 'tumor', ]

# clustering based on jaccard idx - nr of common elements in a set / union of sets

# TODO change str/tum here
gsea_subset <- gsea_tumor

length(unique(gsea_subset$pathway))

paths_genes_list <- lapply(gsea_subset$leadingEdge, function(x){
  genelist <- unlist(strsplit(x, split='|', fixed=T))
})

names(paths_genes_list) <- gsea_subset$pathway


path_jaccard <- lapply(paths_genes_list, function(x){
  p1 <- lapply(paths_genes_list, function(y){
    jacc_idx <- as.numeric(round(length(intersect(x, y)) / length(union(x,y)), digits = 4))
  })
  return(unlist(p1))
})

path_jaccard_mtx <- do.call('cbind', path_jaccard)

hist(path_jaccard_mtx, breaks = 100)
heatmap(path_jaccard_mtx)

# clustering with hclust
path_hclust <- hclust(dist(path_jaccard_mtx), method = "average")
plot(path_hclust, hang = -1, cex = 0.2)

path_hclust_cut_1 <- cutree(path_hclust, h = 1)
path_hclust_cut_12 <- cutree(path_hclust, h = 1.2)
path_hclust_cut_15 <- cutree(path_hclust, h = 1.5)

length(unique(unname(path_hclust_cut_1)))
length(unique(unname(path_hclust_cut_12)))
length(unique(unname(path_hclust_cut_15)))
sort(table(path_hclust_cut_1), decreasing = T)

path_hclust_cut_1[path_hclust_cut_1 == '9']

# merge with gsea result
if(identical(gsea_subset$pathway, names(path_hclust_cut_1))){
  gsea_subset$path_cluster_cut_1 <- path_hclust_cut_1
  gsea_subset$path_cluster_cut_12 <- path_hclust_cut_12
  gsea_subset$path_cluster_cut_15 <- path_hclust_cut_15
}

fwrite(gsea_subset, file.path("/media/iganiemi/T7-iga/st/geomx-processing/results/nact/dge/gsea/gsea_clust/neg",
                              "gsea_post_anno_rank_p_fcval_clust_tumor.csv"))


anno_stroma <- fread(file.path("/media/iganiemi/T7-iga/st/geomx-processing/results/nact/dge/gsea/gsea_clust/pos",
                               "gsea_post_anno_rank_p_fcval_clust_stroma.csv"))

anno_tumor <- fread(file.path("/media/iganiemi/T7-iga/st/geomx-processing/results/nact/dge/gsea/gsea_clust/pos",
                               "gsea_post_anno_rank_p_fcval_clust_tumor.csv"))

pfi_stroma <- fread(file.path("/media/iganiemi/T7-iga/st/geomx-processing/results/nact/dge/gsea/gsea_clust/neg",
                               "gsea_post_pfi_rank_p_fcval_clust_stroma.csv"))

pfi_tumor <- fread(file.path("/media/iganiemi/T7-iga/st/geomx-processing/results/nact/dge/gsea/gsea_clust/neg",
                              "gsea_post_pfi_rank_p_fcval_clust_tumor.csv"))

gsea_subset_clust <- pfi_stroma
gsea_subset_clust <- arrange(gsea_subset_clust, desc(NES))
gsea_subset_clust <- arrange(gsea_subset_clust, NES)

#gsea_subset_clust <- gsea_subset_clust[grepl('GOBP', gsea_subset_clust$pathway)] 

#keep 1st pathway from cluster
gsea_subset_clust_1 <- gsea_subset_clust[!(duplicated(gsea_subset_clust$path_cluster_cut_1)), ]
gsea_subset_clust_12 <- gsea_subset_clust[!(duplicated(gsea_subset_clust$path_cluster_cut_12)), ]
gsea_subset_clust_15 <- gsea_subset_clust[!(duplicated(gsea_subset_clust$path_cluster_cut_15)), ]


# for ++ vs each other ROI type separately
# analyse which pathways are upregulated in 3 and 2
sharing <- stack(sapply(unique(gsea_subset_clust$pathway), function(x){
  gsea_p <- gsea_subset_clust$Contrast[gsea_subset_clust$pathway == x]

  if(length(gsea_p) == 3){
    y <- 'shared_3'
  } else if(length(gsea_p) == 2){
    if(!("negCD8_negIBA1 - posCD8_posIBA1" %in% gsea_p)){
      y <- 'shared_2_no_doubleneg'
    } else if(!("negCD8_posIBA1 - posCD8_posIBA1" %in% gsea_p)){
      y <- 'shared_2_no_negpos'
    } else if(!("posCD8_negIBA1 - posCD8_posIBA1" %in% gsea_p)){
      y <- 'shared_2_no_negpos'
    }
  } else{
    y <- paste0('shared_1_', gsea_p)
  }
  return(y)
}))

gsea_subset_fin <- gsea_subset_clust_1[gsea_subset_clust_1$pathway %in% sharing$ind[sharing$values == 'shared_3'] &
                                          grepl('posCD8_posIBA1', gsea_subset_clust_1$Contrast), ]

gsea_subset_fin <- arrange(gsea_subset_fin, desc(NES))

#TODO filter for 1st appearance on pathway - fir different comparison it may had different
# leading edge genes and haven't been clustered before

fwrite(gsea_subset_fin, file.path("/media/iganiemi/T7-iga/st/geomx-processing/results/nact/dge/gsea/gsea_clust/fin",
                                  "gsea_post_anno_rank_p_fcval_clust_tumor_pos_pospos_shared_fin.csv"))

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

