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

# define variables --------------------------------------------------------

data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'
output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2'

input_rds_path <- file.path(output_dir, 'geomx_qc_norm.RDS')
output_rds_path <- file.path(output_dir, 'geomx_qc_norm_deconv.RDS')

scrna_ref_path <- '/media/iganiemi/T7-iga/st/data/scrna/vaharautio_scrnaseq_dataset_downsampled_for_iga_processed.RDS'
#output_scrna_mtx_path <- file.path(output_dir, 'oc_scrna_ref_mtx_for_spatialdecon.RDS')

norm_type <- 'q3_norm'
scrna_anno <- 'cell_type' # either 'cell_type' or 'mid_lvl_ct'

# imp_vars <- c("Segment", "Annotation_cell", "NACT status", "PFS") # vals used for sankey, detection rate plots, 
# main_var <- "Annotation_cell" # legend in sankey, 
# 
# umap_vars <- c(imp_vars, "Patient", "Sample")

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'deconvolution'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'prism'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'bayes-prism'), showWarnings = T, recursive = T)

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')


# read geomx object -------------------------------------------------------

geomx_obj <- readRDS(input_rds_path)


# prepare data for SpatialDecon -------------------------------------------
# from
# https://bioconductor.org/packages/release/bioc/vignettes/SpatialDecon/inst/doc/SpatialDecon_vignette_NSCLC.html

featureType(geomx_obj) <- "Target"

# sample+ROI+AOI as sample name
sampleNames(geomx_obj) <-  paste0(sData(geomx_obj)[['Sample']], '_', 
                                  sData(geomx_obj)[['Roi']], '_', 
                                  sData(geomx_obj)[['Segment']], '_',
                                  sData(geomx_obj)[['Annotation_cell']])

# get negative probes (aggregated to features already) names
negativeProbefData <- subset(fData(geomx_obj), CodeClass == "Negative")

# estimate bcg for every segment based on neg probes 
# TODO re-check if >1 module
# TODO how to use bg?
bg <- derive_GeoMx_background(norm = geomx_obj@assayData[[norm_type]],
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

#data("mini_singleCell_dataset")

# TODO examine scalingFactor: 1 or 5 or what?
custom_oc_mtx <- create_profile_matrix(mtx = scrna_ref_obj@assays$SCT@data,            # cell x gene count matrix
                                    cellAnnots = scrna_anno_dt,  # cell annotations with cell type and cell name as columns
                                    cellTypeCol = "cell_type",  # column containing cell type
                                    cellNameCol = "cell_name",           # column containing cell ID/name
                                    matrixName = "oc_scrnaseq_ref_cell_type_1x_sct", # name of final profile matrix
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

str(pData(decon_res))


heatmap(t(decon_res$beta), cexCol = 0.5, cexRow = 0.7, margins = c(10,7))

# run extended SpatialDecon -----------------------------------------------

#give info about tumor
geomx_obj$istumor = geomx_obj$Segment == "tumor"

# TODO nuclei counts from geomx are unreliable - match with info from cycif
# TODO examine n_tumor_clusters param with different n

decon_res_ext <- runspatialdecon(object = geomx_obj,
                          norm_elt = norm_type,                # normalized data
                          raw_elt = "exprs",                      # expected background counts for every data point in norm
                          X = safeTME,                            # safeTME matrix, used by default
                          cellmerges = safeTME.matches,           # safeTME.matches object, used by default
                          cell_counts = geomx_obj$Nuclei,      # nuclei counts, used to estimate total cells
                          is_pure_tumor = geomx_obj$istumor,   # identities of the Tumor segments/observations
                          n_tumor_clusters = 5)               # how many distinct tumor profiles to append to safeTME

names(decon_res_ext@assayData)
str(pData(decon_res_ext))

heatmap(sweep(decon_res_ext@experimentData@other$SpatialDeconMatrix, 1, apply(decon_res_ext@experimentData@other$SpatialDeconMatrix, 1, max), "/"),
        labRow = NA, margins = c(10, 5))

# TODO add it at the beginning and rmv from here
# pData(decon_res_ext)$sample_name <-  paste0(sData(decon_res_ext)[['Sample']], '_', 
#                                   sData(decon_res_ext)[['Roi']], '_', 
#                                   sData(decon_res_ext)[['Segment']], '_',
#                                   sData(decon_res_ext)[['Annotation_cell']])



# run extended SpatialDecon with custom oc mtx ----------------------------

decon_res_custom <- runspatialdecon(object = geomx_obj,
                                 norm_elt = norm_type,                # normalized data
                                 raw_elt = "exprs",                      # expected background counts for every data point in norm
                                 X = custom_oc_mtx,                            # safeTME matrix, used by default
                                 cell_counts = geomx_obj$Nuclei,      # nuclei counts, used to estimate total cells
                                 is_pure_tumor = geomx_obj$istumor,   # identities of the Tumor segments/observations
                                 n_tumor_clusters = 5)               # how many distinct tumor profiles to append to safeTME

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

