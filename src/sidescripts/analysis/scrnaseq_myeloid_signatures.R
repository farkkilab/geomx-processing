library(Seurat)
library(data.table)
library(plyr)
library(dplyr)
library(GSVA)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(multinichenetr)
library(SingleCellExperiment)
library(tidyr)
library(tibble)

#TODO check why gsea scores are so big
#TODO check in deconvoluted bulk RNAseq


proj_dir <- '~/Documents/phd/st'
outdir <- '~/Documents/phd/st/geomx-processing/results/batch123-2808/downstream/myelonets_signatures_validation/'

#ct_markers_path <- '/home/iganiemi/Documents/phd/st/geomx-processing/data/signatures/ct_markers.csv'
myelonets_sign_path <- file.path(outdir, 'myelonets_signatures_list.RDS')


source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))
########################
# Vaharautio lab scRNAseq ref from Erdogan

dataset_name <- 'GSE266577'
gsea_outpath <- file.path(outdir, paste0('ssgsea_myelonets_signatures_', dataset_name, '.csv'))


mtx_path <- '/home/iganiemi/Documents/phd/st/data/scrna/GSE266577_counts_raw.mtx'
bar_path <- '/home/iganiemi/Documents/phd/st/data/scrna/GSE266577_barcodes.txt'
ft_path <- '/home/iganiemi/Documents/phd/st/data/scrna/GSE266577_seurat_features.txt'

meta_path <- '/home/iganiemi/Documents/phd/st/data/scrna/GSE266577_metadata.txt'

sample_colname <- 'publication_sample_code_final'

sc_ref <- CreateSeuratObject(
  ReadMtx(mtx_path, bar_path, ft_path, feature.column = 1),
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = fread(meta_path),
  min.cells = 10,
  min.features = 200)

dim(sc_ref)

sc_ref <- subset(sc_ref, subset = percent.mt < 10)

dim(sc_ref)
# log normalise
sc_ref <- NormalizeData(object = sc_ref)

# filter to post samples and macrophages
sc_ref_post <- subset(sc_ref, subset = treatment_stage == 'IDS')
sc_ref_post_macro <- subset(sc_ref_post, subset = cell_type %in% c('Macrophages')) # , 
sc_ref_post_macro <- subset(sc_ref_post_macro, subset = publication_patient_code_final != 'S014')

# retrieved from publication
pt_pfi <- data.frame(sample_id = c("S022", "S001", "S002", "S025", "S008", 
                                   "S009", "S010", "S027", "S011", "S012", 
                                   "S028", "S029", "S014", "S015", "S017", 
                                   "S030", "S018", "S019", "S020", "S021", 
                                   "S031", "S032"),
                     PFI_group = c('long', 'short', 'long', 'long', 'long',
                                   'short', 'short', 'short', 'short', 'short',
                                   'long', 'short', 'short', 'long', 'short',
                                   'short', 'short', 'short', 'short', 'short',
                                   'long', 'short'))

meta_postmacro <- sc_ref_post_macro@meta.data %>%
  mutate(sample_id = publication_patient_code_final) %>%
  mutate(celltype_id = cell_type) %>%
  left_join(pt_pfi)

table(meta_postmacro$sample_id, meta_postmacro$PFI_group)

# meta_postmacro <- sc_ref_post_macro@meta.data %>%
#   mutate(sample_id = patient_id) %>%
#   mutate(cell_name = cell) %>%
#   mutate(celltype_id = cell_subtype) %>%
#   left_join(pt_pfi)
# 
# pt_pfi_macro <- pt_pfi[pt_pfi$sample_id %in% meta_postmacro$sample_id]
# 
# table(meta_postmacro$sample_id, meta_postmacro$PFI_group)

#####################################
# Hautaniemi lab scRNAseq ref from Kaiyang

# dataset_name <- 'GSE165897'
# gsea_outpath <- paste0('~/Documents/phd/st/geomx-processing/results/batch123-2808/lr_interactions/multi_niche_netr/ssgsea_myelonets_signatures_', dataset_name, '.csv')
# 
# mtx_path <- '/home/iganiemi/Documents/phd/st/data/scrna/GSE165897_UMIcounts_HGSOC.tsv'
# meta_path <-  '/home/iganiemi/Documents/phd/st/data/scrna/GSE165897_cellInfo_HGSOC.tsv'
# tcell_sub_path <- '/home/iganiemi/Documents/phd/st/data/scrna/GSE165897_Tcell_subtypes.tsv'
# clinical_path <- '/home/iganiemi/Documents/phd/st/data/scrna/GSE165897_clinical_data.csv'
# 
# sc_ref <- CreateSeuratObject(
#   counts = as.matrix(fread(mtx_path), rownames=1),
#   assay = "RNA",
#   names.field = 1,
#   names.delim = "_",
#   meta.data = fread(meta_path),
#   min.cells = 10,
#   min.features = 200)
# 
# dim(sc_ref)
# gc()
# 
# sc_ref <- subset(sc_ref, subset = percent.mt < 7.5) 
# 
# dim(sc_ref)
# 
# # log normalise
# sc_ref <- NormalizeData(object = sc_ref)
# 
# 
# # filter to post samples and macrophages
# sc_ref_post <- subset(sc_ref, subset = treatment_phase == 'post-NACT')
# sc_ref_post_macro <- subset(sc_ref_post, subset = cell_subtype == 'Macrophages')
# sc_ref_post_macro <- subset(sc_ref_post_macro, subset = patient_id != 'EOC1005')
# 
# 
# 
# # retrieved from publication
# pt_pfi <- fread(clinical_path) %>%
#   mutate(PFI_group = ifelse(PFIdays >= 365, 'long', 'short')) %>%
#   mutate(sample_id = PatientID)
# 
# meta_postmacro <- sc_ref_post_macro@meta.data %>%
#   mutate(sample_id = patient_id) %>%
#   mutate(cell_name = cell) %>%
#   mutate(celltype_id = cell_subtype) %>%
#   left_join(pt_pfi)
# 
# pt_pfi_macro <- pt_pfi[pt_pfi$sample_id %in% meta_postmacro$sample_id]
# 
# table(meta_postmacro$sample_id, meta_postmacro$PFI_group)


# -------------------------------------------------------------------------

# filter to sufficiently expressed genes and make a pseudobulk

expr_mtx <- SeuratObject::GetAssayData(object = sc_ref_post_macro,
                                       assay = "RNA",
                                       layer = "data")

# creating single cell experiment object
sce <- SingleCellExperiment(
  assays = list(counts = expr_mtx),
  colData = meta_postmacro
)

# make sure that gene symbols used in the expression data are updated
sce = makenames_SCE(alias_to_symbol_SCE(sce, "human")) 

# make sure names are valid
SummarizedExperiment::colData(sce)$sample_id = SummarizedExperiment::colData(sce)$sample_id %>% make.names()
SummarizedExperiment::colData(sce)$celltype_id = SummarizedExperiment::colData(sce)$celltype_id %>% make.names()


# filter sufficiently expressed genes
frq_list = get_frac_exprs(
  sce = sce,
  sample_id = 'sample_id', celltype_id =  'celltype_id', group_id = 'PFI_group',
  batches = NA,
  min_cells = 10,
  fraction_cutoff = 0.05, min_sample_prop = 0.25)

# keep genes that are expressed by at least one cell type
genes_oi = frq_list$expressed_df %>%
  filter(expressed == TRUE) %>% pull(gene) %>% unique()
sce = sce[genes_oi, ]

dim(sce)

# calculate pseudobulk
celltype_info = get_avg_pb_exprs(
  sce = sce, 
  sample_id = 'sample_id', celltype_id =  'celltype_id', group_id = 'PFI_group',
  batches = NA, 
  min_cells = 10)

pb_df <- celltype_info$avg_df

# TODO works only if 1 cell type
pb_expr_mtx <- pivot_wider(pb_df[, 1:3], names_from = gene, values_from = average_sample) %>%
  column_to_rownames(var = 'sample') %>%
  t() 

# -------------------------------------------------------------------------

# check if signature genes are matching
signatures_list <- readRDS(myelonets_sign_path)

for(sign in signatures_list){
  print(length(sign))
  print(length(intersect(sign, rownames(pb_expr_mtx))))
  print('XXXX')
}

# calculate gsea
gsea_sign <- gsva(ssgseaParam(pb_expr_mtx, signatures_list, minSize = 10, normalize = T))

# adjust df and save
gsea_long <- melt(gsea_sign)
colnames(gsea_long) <- c('signature', 'sample_id', 'ssgsea_score')
gsea_long <- left_join(gsea_long, pt_pfi) 


fwrite(gsea_long, file.path(outdir, paste0('gsea_per_sample_avg_', dataset_name, '_Macro.csv')))

# all sign, per group
ggplot(data = gsea_long, aes(x = signature, y = ssgsea_score, fill = PFI_group)) +
  geom_boxplot() +
  #geom_violin() +
  geom_point(position= position_jitterdodge(dodge.width = 1, jitter.width= .3, jitter.height = 0),
             size= 0.2, alpha = 0.6) +
  stat_summary(fun = "mean", geom = "point", colour = "red", position = position_dodge(0.9), size=0.3) +
  geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = FALSE, size = 0.2, label.size = 2.8) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 5))


ggsave(file.path(outdir, paste0("mye_vs_PFI_", dataset_name, "_avg_Macro.png")))


# sign <- names(signatures_list)[4]
# gsea_long_sign <- gsea_long[gsea_long$signature == sign, ]
# plot(gsea_long_sign$ssgsea_score, gsea_long_sign$PFIdays)
# 
# for(sign in names(signatures_list)){
#   gsea_long_sign <- gsea_long[gsea_long$signature == sign, ]
#   print(sign)
#   print(cor(gsea_long_sign$ssgsea_score, gsea_long_sign$PFIdays, method = 'spearman'))
#   print(cor(gsea_long_sign$ssgsea_score, gsea_long_sign$PFIdays, method = 'pearson'))
#   print('XXX')
# }



################################################################################
################################################################################3
# GSEA in each cell
# expr_mtx_log <- log2(expr_mtx + 1)

# do gsea
# gsea_sign <- gsva(ssgseaParam(expr_mtx, signatures_list, minSize = 10, normalize = T))
# 
# # adjust df and save
# gsea_long <- melt(gsea_sign)
# colnames(gsea_long) <- c('signature', 'cell_name', 'ssgsea_score')
# gsea_long <- left_join(gsea_long, meta_postmacro)
# 
# fwrite(gsea_long, file.path(outdir, paste0('ssgsea_myelonets_signatures_', dataset_name, '_percell.csv')))
# 
# 
# # plot scores over samples
# 
# sign_name <- names(signatures_list)[1]
# 
# gsea_long_sign <- gsea_long[gsea_long$signature == sign_name, ]
# 
# ##########################
# # per patient
# ggplot(data = gsea_long_sign, aes(x = sample_id, y = ssgsea_score)) +
#   #geom_boxplot() +
#   geom_violin() +
#   geom_point(position= position_jitterdodge(dodge.width = 1, jitter.width= .3, jitter.height = 0),
#              size= 0.2, alpha = 0.6) +
#   stat_summary(fun = "mean", geom = "point", colour = "red", position = position_dodge(0.9), size=0.3) +
#   #geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = TRUE, size = 0.2, label.size = 2.8) +
#   theme(axis.text.x = element_text(angle=45, hjust=1, size = 5))
# 
# ##############################
# # all sign, per group
# ggplot(data = gsea_long, aes(x = signature, y = ssgsea_score, fill = PFI_group)) +
#   #geom_boxplot() +
#   geom_violin() +
#   # geom_point(position= position_jitterdodge(dodge.width = 1, jitter.width= .3, jitter.height = 0),
#   #            size= 0.2, alpha = 0.6) +
#   stat_summary(fun = "mean", geom = "point", colour = "red", position = position_dodge(0.9), size=0.3) +
#   geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = TRUE, size = 0.2, label.size = 2.8) +
#   theme(axis.text.x = element_text(angle=45, hjust=1, size = 5))+
#   ylim(min(gsea_long$ssgsea_score), 1.2)
# 
# ggsave(file.path(outdir, paste0("mye_vs_PFI_", dataset_name, "_percell_Macro.png")))

# ########################################
# # downsample to 200 cells per sample and 7 short samples
# # TODO sample per group
# x <- 200
# short_pfi_random8 <- sample(pt_pfi$sample_id[pt_pfi$PFI_group == 'short'], 4)
# 
# gsea_long_downsampled <- gsea_long %>% 
#   group_by(sample_id, signature) %>%
#   sample_n(min(n(), x))
#   #filter(!(sample_id %in% short_pfi_random8))
#   
# 
# table(gsea_long_downsampled$sample_id, gsea_long_downsampled$PFI_group)
# table(gsea_long$sample_id, gsea_long$PFI_group)
# 
# ggplot(data = gsea_long_downsampled, aes(x = signature, y = ssgsea_score, fill = PFI_group)) +
#   #geom_boxplot() +
#   geom_violin() +
#   # geom_point(position= position_jitterdodge(dodge.width = 1, jitter.width= .3, jitter.height = 0),
#   #            size= 0.2, alpha = 0.6) +
#   stat_summary(fun = "mean", geom = "point", colour = "red", position = position_dodge(0.9), size=0.3) +
#   geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = TRUE, size = 0.2, label.size = 2.8) +
#   theme(axis.text.x = element_text(angle=45, hjust=1, size = 5))+
#   ylim(min(gsea_long_downsampled$ssgsea_score), 1.2)
# 
# ggsave(paste0("~/Documents/phd/st/geomx-processing/results/batch123-2808/lr_interactions/multi_niche_netr/mye_vs_PFI_", dataset_name, "_down_200cell5.png"))
# 
# #########################################
# # mean per patient and then per group
# gsea_long_mean <- gsea_long %>%
#   group_by(sample_id, signature) %>%
#   summarise(mean_ssgsea = mean(ssgsea_score)) %>%
#   ungroup() %>%
#   left_join(pt_pfi) 
# 
# ggplot(data = gsea_long_mean, aes(x = signature, y = mean_ssgsea, fill = PFI_group)) +
#   geom_boxplot() +
#   #geom_violin() +
#   geom_point(position= position_jitterdodge(dodge.width = 1, jitter.width= .3, jitter.height = 0),
#              size= 0.2, alpha = 0.6) +
#   stat_summary(fun = "mean", geom = "point", colour = "red", position = position_dodge(0.9), size=0.3) +
#   geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = FALSE, size = 0.2, label.size = 2.8) +
#   theme(axis.text.x = element_text(angle=45, hjust=1, size = 5)) +
#   ylim(min(gsea_long_mean$mean_ssgsea), 0.7)
# 
# ggsave(paste0("~/Documents/phd/st/geomx-processing/results/batch123-2808/lr_interactions/multi_niche_netr/mye_vs_PFI_", dataset_name, "_meanperpt.png"))
# 
# ##############################################################################
# #############################################################################
# #############################################################################
# # check sd vs cycif once and for all