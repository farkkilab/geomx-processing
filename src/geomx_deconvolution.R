# TODO check if all packages are needed
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(SpatialDecon)
library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(reshape2)
library(Seurat)
library(tibble)
library(BayesPrism)
library(biomaRt)

# define variables --------------------------------------------------------

data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'
output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2'

input_rds_path <- file.path(output_dir, 'geomx_qc_norm.RDS')
output_rds_path <- file.path(output_dir, 'geomx_qc_norm_deconv.RDS')

scrna_ref_path <- '/media/iganiemi/T7-iga/st/data/scrna/vaharautio_scrnaseq_dataset_downsampled_for_iga_processed.RDS'
scrna_ref_cleaned_path <- file.path(output_dir,'deconvolution', 'scrna_ref_cleaned.RDS')
#output_scrna_mtx_path <- file.path(output_dir, 'oc_scrna_ref_mtx_for_spatialdecon.RDS')

norm_type <- 'q3_norm'
scrna_anno <- 'mid_lvl_ct' # either 'cell_type' or 'mid_lvl_ct'
ct_nr_thr <- 20 # best 20 or 45 to rmv cell states not abundant enough in scrnaseq

# imp_vars <- c("Segment", "Annotation_cell", "NACT status", "PFS") # vals used for sankey, detection rate plots, 
# main_var <- "Annotation_cell" # legend in sankey, 
# umap_vars <- c(imp_vars, "Patient", "Sample")

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'deconvolution'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'deconvolution', 'spatial_decon', scrna_anno), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'prism'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'deconvolution', 'bayes_prism'), showWarnings = T, recursive = T)

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')

# prepare scrnaseq reference dataset --------------------------------------
#TODO can be moved to some other script
geomx_obj <- readRDS(input_rds_path)
scrna_ref_obj <- readRDS(scrna_ref_path)

###################################
###################################
# repair synonymuous gene names
#TODO export to function
length(rownames(geomx_obj@assayData$exprs))
length(rownames(scrna_ref_obj@assays$RNA@data))
length(intersect(rownames(geomx_obj@assayData$exprs), rownames(scrna_ref_obj@assays$RNA@data)))

geo_non_ex <- setdiff(rownames(geomx_obj@assayData$exprs), rownames(scrna_ref_obj@assays$RNA@data))

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

geo_non_ex_syn <- getBM(attributes = c('external_gene_name', 'external_synonym'),
                        filters = 'external_gene_name',
                        values = geo_non_ex,
                        mart = ensembl)


geo_syn_in_scrna <- filter(geo_non_ex_syn, external_synonym %in% rownames(scrna_ref_obj@assays$RNA@data)) %>%
  distinct(external_gene_name, .keep_all = T) %>% # it'll remove a handful of weird genes with multiple synonyms simultaneously present in scrna, may be ignored
  distinct(external_synonym, .keep_all = T)
  
common_genes <- sapply(rownames(scrna_ref_obj@assays$RNA@data), function(x){
  if(x %in% geo_syn_in_scrna$external_synonym){
    gname <- geo_syn_in_scrna$external_gene_name[geo_syn_in_scrna$external_synonym == x]
  } else{
    gname <- x
  }
  return(gname)
})

#  make a new assay with renamed genes
RNA_common_genes <- scrna_ref_obj@assays$RNA
RNA_common_genes@counts@Dimnames[[1]] <- common_genes
RNA_common_genes@data@Dimnames[[1]] <- common_genes
scrna_ref_obj@assays$RNA_common_genes <- RNA_common_genes

length(intersect(rownames(geomx_obj@assayData$exprs), rownames(scrna_ref_obj@assays$RNA@data)))
length(intersect(rownames(geomx_obj@assayData$exprs), rownames(scrna_ref_obj@assays$RNA_common_genes@data)))

###########################
###########################

scrna_ref_obj@meta.data$cell_type <- ifelse(scrna_ref_obj@meta.data$cell_type == 'Epithelial cells', 
                                            'tumor', scrna_ref_obj@meta.data$cell_type)

# cell states - clustering tumor cells by patient
scrna_ref_obj@meta.data$cell_state <- ifelse(scrna_ref_obj@meta.data$cell_type == 'tumor', 
                                             paste0('tumor_', scrna_ref_obj@meta.data$patient), 
                                             scrna_ref_obj@meta.data$cell_type)




#fix  mid-lvl-ct
scrna_ref_obj@meta.data$mid_lvl_ct <- ifelse(scrna_ref_obj@meta.data$mid_lvl_ct == 'Plasma cells', 
                                             'Bcells', scrna_ref_obj@meta.data$mid_lvl_ct)
scrna_ref_obj@meta.data$mid_lvl_ct <- ifelse(scrna_ref_obj@meta.data$mid_lvl_ct == 'Classical monocytes', 
                                             'Macrophages', scrna_ref_obj@meta.data$mid_lvl_ct)
scrna_ref_obj@meta.data$mid_lvl_ct <- ifelse(scrna_ref_obj@meta.data$mid_lvl_ct == 'Epithelial cells', 
                                             'tumor', scrna_ref_obj@meta.data$mid_lvl_ct)

# TODO check if 'other' cells should be removed or not
# mid_ct_other_cells <- scrna_ref_obj@meta.data$cell_name[scrna_ref_obj@meta.data$mid_lvl_ct == 'other']
# scrna_ref_obj <- scrna_ref_obj[, !colnames(scrna_ref_obj) %in% mid_ct_other_cells]

########################################
########################################
# QC of cell states

# TODO think of changing labels for mast cells, Th17, tumor_H103

# plot.cor.phi (input=t(scrna_ref_obj@assays$RNA@data),
#               input.labels=scrna_ref_obj@meta.data$cell_state,
#               title="cell state correlation",
#               #specify pdf.prefix if need to output to pdf
#               #pdf.prefix="gbm.cor.cs",
#               cexRow=0.6, cexCol=0.6,
#               margins=c(6,6))
# 
# dev.off()
# 
# plot.cor.phi (input=t(scrna_ref_obj@assays$RNA@data),
#               input.labels=scrna_ref_obj@meta.data$mid_lvl_ct,
#               title="cell type correlation",
#               #specify pdf.prefix if need to output to pdf
#               #pdf.prefix="gbm.cor.ct",
#               cexRow=0.5, cexCol=0.5,
# )
# 
# dev.off()

#################################
#################################

# check genes outliers
# TODO important for BayesPrism, check how it affects SpatialDecon 

scrna_stat <- plot.scRNA.outlier(
  input=t(scrna_ref_obj@assays$RNA_common_genes@data), #make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels=scrna_ref_obj@meta.data$cell_type,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE #return the data used for plotting.
  #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)

View(scrna_stat)

geomx_stat <- plot.bulk.outlier(
  bulk.input=t(geomx_obj@assayData$exprs),#make sure the colnames are gene symbol or ENSMEBL ID
  sc.input=t(scrna_ref_obj@assays$RNA_common_genes@data), #make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels=scrna_ref_obj@meta.data$cell_type,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE
  #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
)

View(geomx_stat)

# filter out outlier genes
scrna_filt <- cleanup.genes (input=t(scrna_ref_obj@assays$RNA_common_genes@data),
                             input.type="count.matrix",
                             species="hs", 
                             gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                             exp.cells=5)

dim(t(scrna_ref_obj@assays$RNA_common_genes@data))
dim(scrna_filt)

# geomx doen't have to be filtered since later on they took only intersection of genes
# TODO but maybe it should? check if it improves SpatialDecon
geomx_stat_to_rm <- geomx_stat[ rowSums(geomx_stat[, -c(1,2)]) >= 1, ]
geomx_filtered <- geomx_obj[!(rownames(geomx_obj) %in% geomx_stat_to_rm),  ]

# check expr concordance for different gene types
#plot.bulk.vs.sc (sc.input = scrna_filt, bulk.input = geomx_raw)

# subset to protein coding genes
scrna_filt_pc <-  select.gene.type(scrna_filt, gene.type = "protein_coding")

#  make a new assay with filtered genes
RNA_common_genes_filt_pc <- scrna_ref_obj@assays$RNA_common_genes
RNA_common_genes_filt_pc@counts <- RNA_common_genes_filt_pc@counts[rownames(RNA_common_genes_filt_pc@counts) %in% colnames(scrna_filt_pc),  ]
RNA_common_genes_filt_pc@data <- RNA_common_genes_filt_pc@data[rownames(RNA_common_genes_filt_pc@data) %in% colnames(scrna_filt_pc),  ]
scrna_ref_obj@assays$RNA_common_genes_filt_pc <- RNA_common_genes_filt_pc

# save adjusted scRNAseq file
saveRDS(scrna_ref_obj, file = scrna_ref_cleaned_path)

# TODO takes > 64G of memory, have to be run on linux machine
# subset to signature genes (differentially expressed trough cell types)
# diff_exp_stat <- get.exp.stat(sc.dat=scrna_raw[,colSums(scrna_raw>0)>3],# filter genes to reduce memory use
#                               cell.type.labels=scrna_ref_obj@meta.data$cell_type,
#                               cell.state.labels=scrna_ref_obj@meta.data$cell_state,
#                               pseudo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
#                               cell.count.cutoff=20, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
#                               n.cores=8 #number of threads
# )

# scrna_filt_pc_sig <- select.marker (sc.dat=scrna_filt_pc,
#                                          stat=diff_exp_stat,
#                                          pval.max=0.01,
#                                          lfc.min=0.1)

# dim(scrna_filt_pc_sig)

############################################################################
############################################################################
############################################################################
###############################################################################
###############################################################################
# BayesPrism
# deconvolution by bayesprism ---------------------------------------------
# https://github.com/Danko-Lab/BayesPrism/blob/main/tutorial_deconvolution.html

# load cleaned scrna and geomx
geomx_obj <- readRDS(input_rds_path)
scrna_ref_obj <- readRDS(scrna_ref_cleaned_path)

# remove cells from cell states with < thr cells in ref scrnaseq
ct_freq <- as.data.frame(table(scrna_ref_obj@meta.data$cell_state))
low_ct_cells <- scrna_ref_obj@meta.data$cell_name[scrna_ref_obj@meta.data$cell_state %in% 
                                                    as.character(ct_freq$Var1[ct_freq$Freq < ct_nr_thr])]

scrna_ref_obj <- scrna_ref_obj[, !colnames(scrna_ref_obj) %in% low_ct_cells]

# make a prism object
prism_obj <- new.prism(
  reference=t(scrna_ref_obj@assays$RNA_common_genes_filt_pc@data), 
  mixture=t(geomx_obj@assayData$exprs),
  input.type="count.matrix", 
  cell.type.labels = scrna_ref_obj@meta.data[[scrna_anno]], 
  cell.state.labels = scrna_ref_obj@meta.data$cell_state,
  key="tumor",
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

# run bayesprism
bprism_res <- run.prism(prism = prism_obj, n.cores=12)

# save res
saveRDS(bprism_res, file = file.path(output_dir,'deconvolution', 'bayes_prism', 
                                     paste0('bp_res_', scrna_anno, '_', ct_nr_thr, '.RDS')))

##########################
###########################
# exploration
slotNames(bprism_res)

mean_ct_frac <- get.fraction (bp=bprism_res,
                              which.theta="final",
                              state.or.type="type")

# TODO mask ct_frac results if cv > 0.2-0.5 (0.1 thr for bult, 0.5 for Visium, GeoMx should be in the middle)
# histogram suggests 0.5 as thr
ct_frac_cv <- bprism_res@posterior.theta_f@theta.cv

# extract posterior mean of cell type-specific gene expression count matrix Z  
# TODO normalise it!
#Clustering bulk samples by theta or Z (Z can be normalized by vst(round(t(Z.tumor))), 
# using the vst function from the DESeq2 package.)
tumor_gene_exp <- get.exp (bp=bprism_res,
                           state.or.type="type",
                           cell.name="tumor")

###############################################################################
###############################################################################
############################################################################
# prepare data for SpatialDecon -------------------------------------------
# from
# https://bioconductor.org/packages/release/bioc/vignettes/SpatialDecon/inst/doc/SpatialDecon_vignette_NSCLC.html

#TODO check geomx_obj and geomx_obj_filtered 

featureType(geomx_obj) <- "Target"

sampleNames(geomx_obj) <- sData(geomx_obj)[['dcc_filename']]

# prepare cell profile matrix from reference scRNAseq

# format annotations
scrna_anno_dt <- scrna_ref_obj@meta.data[, c('cell_name', scrna_anno)]
rownames(scrna_anno_dt) <- NULL
colnames(scrna_anno_dt) <- c('cell_name', 'cell_type')

# TODO examine scalingFactor: 1 or 5 or what?
custom_oc_mtx <- create_profile_matrix(mtx = scrna_ref_obj@assays$SCT@data,            # cell x gene count matrix
                                       cellAnnots = scrna_anno_dt,  # cell annotations with cell type and cell name as columns
                                       cellTypeCol = "cell_type",  # column containing cell type
                                       cellNameCol = "cell_name",           # column containing cell ID/name
                                       matrixName = "oc_scrnaseq_ref_cell_type_sct", # name of final profile matrix
                                       outDir = output_dir,                    # path to desired output directory, set to NULL if matrix should not be written
                                       normalize = FALSE,                # Should data be normalized?
                                       minCellNum = 50,                   # minimum number of cells of one type needed to create profile, exclusive
                                       minGenes = 10,                    # minimum number of genes expressed in a cell, exclusive
                                       scalingFactor = 1,                # what should all values be multiplied by for final matrix
                                       discardCellTypes = TRUE)          # should cell types be filtered for types like mitotic, doublet, low quality, unknown, etc.




######################################################
######################################################
# all bells and whistles from the vignette
# prepare data for SpatialDecon -------------------------------------------
# from
# https://bioconductor.org/packages/release/bioc/vignettes/SpatialDecon/inst/doc/SpatialDecon_vignette_NSCLC.html

#TODO check geomx_obj and geomx_obj_filtered 

featureType(geomx_obj) <- "Target"

sampleNames(geomx_obj) <- sData(geomx_obj)[['dcc_filename']]

# get negative probes (aggregated to features already) names
negativeProbefData <- subset(fData(geomx_obj), CodeClass == "Negative")

# estimate bcg for every segment based on neg probes 
# TODO re-check if >1 module
geomx_bg <- derive_GeoMx_background(norm = geomx_obj@assayData[[norm_type]],
                             probepool = fData(geomx_obj)$Module,
                             negnames = negativeProbefData$TargetName)

# load pre-defined TME cell profile matrix
tme_mtx <- download_profile_matrix(species = "Human",
                                   age_group = "Adult", 
                                   matrixname = "ImmuneTumor_safeTME")

data("safeTME")
data("safeTME.matches")

# prepare cell profile matrix from reference scRNAseq ---------------------

scrna_ref_obj <- readRDS(scrna_ref_path)

# format annotations
scrna_anno_dt <- scrna_ref_obj@meta.data[, c('cell_name', scrna_anno)]
rownames(scrna_anno_dt) <- NULL
colnames(scrna_anno_dt) <- c('cell_name', 'cell_type')

# TODO examine scalingFactor: 1 or 5 or what?
custom_oc_mtx <- create_profile_matrix(mtx = scrna_ref_obj@assays$SCT@data,            # cell x gene count matrix
                                    cellAnnots = scrna_anno_dt,  # cell annotations with cell type and cell name as columns
                                    cellTypeCol = "cell_type",  # column containing cell type
                                    cellNameCol = "cell_name",           # column containing cell ID/name
                                    matrixName = "oc_scrnaseq_ref_cell_type_sct", # name of final profile matrix
                                    outDir = output_dir,                    # path to desired output directory, set to NULL if matrix should not be written
                                    normalize = FALSE,                # Should data be normalized?
                                    minCellNum = 50,                   # minimum number of cells of one type needed to create profile, exclusive
                                    minGenes = 10,                    # minimum number of genes expressed in a cell, exclusive
                                    scalingFactor = 1,                # what should all values be multiplied by for final matrix
                                    discardCellTypes = TRUE)          # should cell types be filtered for types like mitotic, doublet, low quality, unknown, etc.



# run basic SpatialDecon --------------------------------------------------
decon_res <-  runspatialdecon(object = geomx_obj,
                      norm_elt = norm_type,
                      raw_elt = "exprs",
                      X = tme_mtx,
                      align_genes = TRUE)

#heatmap(t(decon_res$beta), cexCol = 0.5, cexRow = 0.7, margins = c(10,7))

# run extended SpatialDecon -----------------------------------------------

#give info about tumor
geomx_obj$istumor = geomx_obj$Segment == "tumor"

# TODO nuclei counts from geomx are unreliable - match with info from cycif
# TODO examine if istumor should be used - it's not pure in our case
# TODO examine n_tumor_clusters param with different n

decon_res_ext <- runspatialdecon(object = geomx_obj,
                          norm_elt = norm_type,                # normalized data
                          raw_elt = "exprs",                      # expected background counts for every data point in norm
                          X = safeTME,                            # safeTME matrix, used by default
                          cellmerges = safeTME.matches,           # safeTME.matches object, used by default
                          #cell_counts = geomx_obj$Nuclei,      # nuclei counts, used to estimate total cells
                          #is_pure_tumor = geomx_obj$istumor,   # identities of the Tumor segments/observations
                          n_tumor_clusters = 5)               # how many distinct tumor profiles to append to safeTME

heatmap(sweep(decon_res_ext@experimentData@other$SpatialDeconMatrix, 1, apply(decon_res_ext@experimentData@other$SpatialDeconMatrix, 1, max), "/"),
        labRow = NA, margins = c(10, 5))

# run extended SpatialDecon with custom oc mtx ----------------------------

decon_res_custom <- runspatialdecon(object = geomx_obj,
                                 norm_elt = norm_type,                # normalized data
                                 raw_elt = "exprs",                      # expected background counts for every data point in norm
                                 X = custom_oc_mtx,                            # safeTME matrix, used by default
                                 #cell_counts = geomx_obj$Nuclei,      # nuclei counts, used to estimate total cells
                                 #is_pure_tumor = geomx_obj$istumor,   # identities of the Tumor segments/observations
                                 n_tumor_clusters = 5)               # how many distinct tumor profiles to append to safeTME

# run extended SpatialDecon with custom oc mtx and estimated bg------------

#TODO to use bg geomx obj have to be converted to seurat
decon_res_custom_bg <- runspatialdecon(object = geomx_obj,
                                    bg = geomx_bg,                      # expected background counts for every data point in norm
                                    X = custom_oc_mtx,                            # safeTME matrix, used by default
                                    #cell_counts = geomx_obj$Nuclei,      # nuclei counts, used to estimate total cells
                                    #is_pure_tumor = geomx_obj$istumor,   # identities of the Tumor segments/observations
                                    n_tumor_clusters = 5)               # how many distinct tumor profiles to append to safeTME


# save spatialdecon results -----------------------------------------------

# colnames from output
res_cols <- c("beta", "p", "t", "se", "prop_of_all", "prop_of_nontumor")

res_ext <- pData(decon_res_ext)[, c(res_cols, "sigma")]
saveRDS(res_ext, file = file.path(output_dir, 'deconvolution', 'spatial_decon', scrna_anno, 
                              'spat_dec_res_ext.rds'))

res_custom <- pData(decon_res_custom)[, c(res_cols, "sigmas")]
saveRDS(res_ext, file = file.path(output_dir, 'deconvolution', 'spatial_decon', scrna_anno, 
                                  'spat_dec_res_custom.rds'))


# messy code

# visualise ---------------------------------------------------------------

decon <- decon_res_custom
outp_subdir <- 'custom_cell_types_1x' # 'custom_mid_lvl_ct_5x', 'ext', 'basic' oc_scrnaseq_ref_cell_type_1x_sct

dir.create(file.path(output_dir, 'deconvolution', outp_subdir,  'per_sample'), showWarnings = T, recursive = T)

macro_name <- 'Macrophages' # 'macrophages' in softTME
cd8_name <- 'Tem/Trm cytotoxic T cells' # 'Tcells' in mid lvl ct 'CD8.Tcells' in softTME

res_type <- 'beta' # c('beta', 'cell_counts', 'prop_all')

#######################

pData(decon)$sample_name <-  paste0(sData(decon)[['Sample']], '_', 
                                            sData(decon)[['Roi']], '_', 
                                            sData(decon)[['Segment']], '_',
                                            sData(decon)[['Annotation_cell']])

# mein got, the default functions from spatialdecon are so bad
# cell counts
o = hclust(dist(t(decon$cell.counts$cell.counts)))$order
layout(mat = (matrix(c(1, 2), 1)), widths = c(7, 3))
pdf(file=file.path(output_dir, 'deconvolution',outp_subdir, 'cell_counts_custom.pdf'), width=1500, height=1000)
TIL_barplot(t(decon$cell.counts$cell.counts[, o]), draw_legend = TRUE, 
            cex.names = 0.5)
dev.off()

# proportions of cells
temp = replace(decon$prop_of_nontumor, is.na(decon$prop_of_nontumor), 0)
o = hclust(dist(temp[decon$Segment == "stroma",]))$order
pdf(file=file.path(output_dir, 'deconvolution',outp_subdir, 'cell_prop_custom.pdf'), width=1500, height=1000)
TIL_barplot(t(decon$prop_of_nontumor[decon$Segment == "stroma",])[, o], 
            draw_legend = TRUE, cex.names = 0.5)
dev.off()

# scatter

if(res_type == 'beta'){
  res_mtx <- sData(decon)$beta
} else if(res_type == 'cell_counts'){
  res_mtx <- sData(decon)$cell.counts$cell.counts
} else if(res_type == 'prop_all'){
  res_mtx <- sData(decon)$prop_of_all
}

# colnames(sData(decon_res_ext))
# View(sData(decon_res_ext)[, c(2:7, 22:26, 54, 48:59, 63:76)])
# View(sData(decon_res_ext)[, c(2:7, 22:26, 48:59, 63:76)])

dt <- sData(decon)[, c('dcc_filename', 'Patient', 'Segment', 'Sample','Nuclei', 'NACT status', 'Annotation_cell')]
dt$sample_name <- rownames(dt)
dt <- cbind(dt, res_mtx)

dt_stroma <- dt[dt$Segment == 'stroma', ]
dt_tumor <- dt[dt$Segment == 'tumor', ]

ggplot(data = dt_stroma, aes(x = get(macro_name), y = get(cd8_name), shape = Annotation_cell, color = Sample)) +
  geom_point(alpha = 0.5, size = 2)

ggsave(file.path(output_dir, 'deconvolution',outp_subdir,  paste0(res_type, '_stroma.pdf')),
       width = 1500, height = 1000, unit = 'px')

ggplot(data = dt_tumor, aes(x = get(macro_name), y = get(cd8_name), shape = Annotation_cell, color = Sample)) +
  geom_point(alpha = 0.5, , size = 2)

ggsave(file.path(output_dir, 'deconvolution', outp_subdir,  paste0(res_type, '_tumor.pdf')),
       width = 1500, height = 1000, unit = 'px')

#################
# plot per sample - validation of ROI selection

for(s in unique(dt$Sample)){
  dt_s <- dt[dt$Sample == s, ]
  
  ggplot(data = dt_s, aes(x = get(macro_name), y = get(cd8_name), color = Annotation_cell, shape = Segment)) +
    geom_point(alpha = 0.5) +
    ggtitle(s)
  
  ggsave(file.path(output_dir, 'deconvolution', outp_subdir,  'per_sample', paste0(res_type, '_', s, '.pdf')),
         width = 1500, height = 1000, unit = 'px')
}


##################
# stacked barplots

# wide to long
dt_long <- melt(dt, id.vars = c(c('dcc_filename', 'Patient', 'Segment', 'Sample','Nuclei', 'NACT status', 'Annotation_cell', 'sample_name')))

dt_long <- arrange(dt_long, Annotation_cell)

ggplot(data = dt_long[dt_long$Segment == 'stroma',], aes(x = sample_name, y = value, fill = variable)) +
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6))

ggsave(file.path(output_dir, 'deconvolution', outp_subdir, paste0(res_type,'_barplot_stroma.pdf')),
       width = 1500, height = 1000, unit = 'px')

ggplot(data = dt_long[dt_long$Segment == 'tumor',], aes(x = sample_name, y = value, fill = variable)) +
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6))

ggsave(file.path(output_dir, 'deconvolution', outp_subdir, paste0(res_type,'_barplot_tumor.pdf')),
       width = 1500, height = 1000, unit = 'px')


# prepare data for PRISM --------------------------------------------------

geomx_obj <- readRDS(input_rds_path)
scrna_ref_obj <- readRDS(scrna_ref_path)

geomx_raw <- geomx_obj@assayData$exprs
scrna_raw <- scrna_ref_obj@assays$RNA@data

dim(geomx_raw)
dim(scrna_raw)

#TODO 700 genes from geomx not in scrna but they are synonyms - look for them and repair
rownames(geomx_raw)[which(!(rownames(geomx_raw) %in% rownames(scrna_raw)))]

# find common genes and filter matrices
common_genes <- intersect(rownames(geomx_raw), rownames(scrna_raw))

scrna_raw <- scrna_raw[rownames(scrna_raw) %in% common_genes, ]
geomx_raw <- geomx_raw[rownames(geomx_raw) %in% common_genes, ]

# reorder scrna
reorder_idx <- match(rownames(geomx_raw), rownames(scrna_raw))
scrna_raw <- scrna_raw[reorder_idx, ]  

identical(rownames(geomx_raw), rownames(scrna_raw))

# change format and write
scrna_raw <- as.data.frame(scrna_raw)
geomx_raw <- as.data.frame(geomx_raw)

scrna_raw <- rownames_to_column(scrna_raw, var = 'gene')
geomx_raw <- rownames_to_column(geomx_raw, var = 'gene')

View(scrna_raw[1:10, 1:10])

table(scrna_ref_obj@meta.data$cell_type)

# remove trash cell type 
scrna_raw <- scrna_raw[, !names(scrna_raw) %in% 
                          c(scrna_ref_obj@meta.data$cell_name[scrna_ref_obj@meta.data$cell_type == 'Late erythroid'])]

#TODO fread give some weird but this is terribly slow
write.table(scrna_raw, file = file.path(output_dir, 'prism', 'scrna_raw_shared.tsv'), row.names = F, col.names = T, sep = '\t')
write.table(geomx_raw, file = file.path(output_dir, 'prism', 'geomx_raw_shared.tsv'), row.names = F, col.names = T, sep = '\t')

# prepare binary cell type weights table
cell_types <- scrna_ref_obj@meta.data[, c('cell_name', 'cell_type', 'mid_lvl_ct')]
cell_types <- cell_types[cell_types$cell_type != 'Late erythroid', ]

cell_types$low_lvl_ct <- ifelse(cell_types$mid_lvl_ct == 'Epithelial cells', 'tumor',
                                ifelse(cell_types$mid_lvl_ct %in% c('Endothelial cells', 'Fibroblasts'), 'stroma', 'immune'))


make_bin_weights_mtx <- function(dt, id_name, feature_name, output_path){
  weights <- unique(dt[, c(id_name, feature_name)])
  weights <- dcast(weights, formula = get(id_name) ~ get(feature_name), fun.aggregate = length)
  colnames(weights)[1] <- id_name #fix name
  
  reorder_idx <- match(rownames(dt), weights[[id_name]])
  weights <- weights[reorder_idx, ]  
  
  print(identical(rownames(dt), weights[[id_name]]))
  
  write.table(weights, file = output_path, row.names = F, col.names = T, sep = '\t')
}

make_bin_weights_mtx(cell_types, 'cell_name', 'cell_type', file.path(output_dir, 'prism', 'weights_high_lvl_ct.tsv'))
make_bin_weights_mtx(cell_types, 'cell_name', 'mid_lvl_ct', file.path(output_dir, 'prism', 'weights_mid_lvl_ct.tsv'))
make_bin_weights_mtx(cell_types, 'cell_name', 'low_lvl_ct', file.path(output_dir, 'prism', 'weights_low_lvl_ct.tsv'))

# check if bin weights have the same order as scrna mtx

scrna_mtx <- fread(file.path(output_dir, 'prism', 'scrna_raw_shared.tsv'))
weights_low <- fread(file.path(output_dir, 'prism', 'weights_low_lvl_ct.tsv'))
weights_mid <- fread(file.path(output_dir, 'prism', 'weights_mid_lvl_ct.tsv'))
weights_high <- fread(file.path(output_dir, 'prism', 'weights_high_lvl_ct.tsv'))

identical(colnames(scrna_mtx)[-1], weights_low$cell_name)
identical(colnames(scrna_mtx)[-1], weights_mid$cell_name)
identical(colnames(scrna_mtx)[-1], weights_high$cell_name)

# prism commands 
#TODO change into R code when they fix the package

# prism-gain -H /media/iganiemi/T7-iga/st/geomx-processing/results/nact2/prism/scrna_raw_shared.tsv \
# -G /media/iganiemi/T7-iga/st/geomx-processing/results/nact2/prism/scrna_gains.tsv

# this won't be used since we already have labels
# prism-clust -H /media/iganiemi/T7-iga/st/geomx-processing/results/nact2/prism/scrna_raw_shared.tsv -g /media/iganiemi/T7-iga/st/geomx-processing/results/nact2/prism/scrna_gains.tsv -T /media/iganiemi/T7-iga/st/geomx-processing/results/nact2/prism/sc_tree.tsv

# prism-decom -H /media/iganiemi/T7-iga/st/geomx-processing/results/nact2/prism/geomx_raw_shared.tsv /media/iganiemi/T7-iga/st/geomx-processing/results/nact2/prism/scrna_raw_shared.tsv \
# -g /media/iganiemi/T7-iga/st/geomx-processing/results/nact2/prism/scrna_gains.tsv \
# -w /media/iganiemi/T7-iga/st/geomx-processing/results/nact2/prism/weights_high_lvl_ct.tsv \
# -Z /media/iganiemi/T7-iga/st/geomx-processing/results/nact2/prism/results/geomx_decom_high_lvl_ct.tsv \
# -G /media/iganiemi/T7-iga/st/geomx-processing/results/nact2/prism/results/geomx_gains_high_lvl_ct.tsv \
# -W /media/iganiemi/T7-iga/st/geomx-processing/results/nact2/prism/results/geomx_weights_high_lvl_ct.tsv 

high_deconv <- fread(file.path(output_dir, 'prism','results', 'geomx_decom_mid_lvl_ct.tsv'))
high_gains <- fread(file.path(output_dir, 'prism','results', 'geomx_gains_mid_lvl_ct.tsv'))
high_weights <- fread(file.path(output_dir, 'prism','results', 'geomx_weights_mid_lvl_ct.tsv'))





##############################################################
##############################################################
##############################################################

# compare with ROI type 

ct_frac <- as.data.frame(mean_ct_frac)
ct_frac$stroma <- ct_frac$Fibroblasts + ct_frac$`Endothelial cells`
ct_frac$immune <- ct_frac$`Regulatory T cells` + ct_frac$`Memory B cells` + ct_frac$`CD16- NK cells` + 
  ct_frac$`Tem/Trm cytotoxic T cells` + ct_frac$`Tcm/Naive helper T cells` + ct_frac$Macrophages + ct_frac$`Mast cells` +
  ct_frac$`Migratory DCs` + ct_frac$`Plasma cells` + ct_frac$ILC + ct_frac$pDC + ct_frac$`Type 17 helper T cells` +
  ct_frac$`CD16+ NK cells` + ct_frac$`NK cells` + ct_frac$`Naive B cells` + ct_frac$DC1 + ct_frac$`Classical monocytes`

#ct_frac <- mutate(ct_frac, immune = rowSums(select(ct_frac, -tumor, -stroma, -Fibroblasts, -`Endothelial cells`)))
ct_frac$tot <- ct_frac$tumor + ct_frac$stroma + ct_frac$immune

ct_frac <- rownames_to_column(ct_frac, 'dcc_filename')
ct_frac <- left_join(ct_frac, sData(geomx_obj)[, c('dcc_filename', 'Patient', 'Segment', 'Sample','Nuclei', 'NACT status', 'Annotation_cell')],
                     by = 'dcc_filename')


ggplot(data = ct_frac, aes(x = Segment, y = tumor)) +
  geom_violin() 
  # geom_point(position= position_jitterdodge(dodge.width = 1, jitter.width= .3, jitter.height = 0),
  #            size= 0.2, alpha = 0.6) 


ggplot(data = ct_frac, aes(x = Segment, y = stroma)) +
  geom_violin() 

ggplot(data = ct_frac, aes(x = Segment, y = immune)) +
  geom_violin() 


ct_frac_stroma <- ct_frac[ct_frac$Segment == 'stroma', ]
ct_frac_tumor <- ct_frac[ct_frac$Segment == 'tumor', ]

ggplot(data = ct_frac_stroma, aes(x = Macrophages, y = `Tem/Trm cytotoxic T cells`, shape = Annotation_cell, color = Sample)) +
  geom_point(alpha = 0.5, size = 2)

ggsave(file.path(output_dir, 'bayes-prism', paste0('macro_cd8_stroma.pdf')),
       width = 1500, height = 2000, unit = 'px')

ggplot(data = ct_frac_tumor, aes(x = Macrophages, y = `Tem/Trm cytotoxic T cells`, shape = Annotation_cell, color = Sample)) +
  geom_point(alpha = 0.5, size = 2)

ggsave(file.path(output_dir, 'bayes-prism', paste0('macro_cd8_tumor.pdf')),
       width = 1500, height = 2000, unit = 'px')

###############

ggplot(data = ct_frac_stroma, aes(x = Macrophages, y = `Tem/Trm cytotoxic T cells`, color = Annotation_cell)) +
  geom_point(alpha = 0.5, size = 2)

ggsave(file.path(output_dir, 'bayes-prism', paste0('macro_cd8_stroma2.pdf')),
       width = 1500, height = 2000, unit = 'px')

ggplot(data = ct_frac_tumor, aes(x = Macrophages, y = `Tem/Trm cytotoxic T cells`, color = Annotation_cell)) +
  geom_point(alpha = 0.5, size = 2)

ggsave(file.path(output_dir, 'bayes-prism', paste0('macro_cd8_tumor2.pdf')),
       width = 1500, height = 2000, unit = 'px')

