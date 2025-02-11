library(Biobase, quietly =T)
library(NanoStringNCTools, quietly =T)
library(GeomxTools, quietly =T)
library(GeoDiff, quietly =T)
library(Biobase, quietly =T)
library(DESeq2, quietly =T)
library(SpatialDecon, quietly =T)
#install preprocessCore manually from source
# BiocManager::install("preprocessCore", configure.args = c(preprocessCore = "--disable-threading"), 
# force= TRUE, update=TRUE, type = "source")
library(preprocessCore, quietly =T)
library(umap, quietly =T)
library(Rtsne, quietly =T)

# problems with Matrix package - i has to be lower that 1.7 to work with lmer
# devtools::install_version("Matrix","1.6.4")
library(Matrix)

library(limma)
library(pvca)
library(harmony)

library(ggforce, quietly =T)
library(plyr, quietly =T)
library(dplyr, quietly =T)
library(ggplot2, quietly =T)
library(cowplot, quietly =T)
library(reshape2, quietly =T)
library(data.table, quietly =T)
library(tibble, quietly =T)

library(Seurat, quietly =T)
library(BayesPrism, quietly =T)
library(biomaRt, quietly =T)

# define variables and paths ----------------------------------------------
proj_dir <- '~/Documents/phd/st'
data_dir <- '~/Documents/phd/st/data/geomx/geomx_batch2_1124/'

dcc_path <- dir(file.path(data_dir, "dcc"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
pkc_path <- file.path(data_dir, 'metadata', 'Hs_R_NGS_WTA_v1.0.pkc')

# anno file have to contain sheet named 'Sheet1' and following column names:
# 'Sample_ID', 'Aoi', 'Roi', 'Sample', 'Slide_Name'
# '_' instead of whitespace in all column names!!!
anno_path <- file.path(data_dir, 'metadata', 'dcc_metadata_all_batch2_1124.xlsx')

# path to reference scRNAseq dataset for deconvolution
# have to contain 'cell_type' column name
scrna_ref_path <- file.path(proj_dir, 'data/scrna/vaharautio_scrnaseq_dataset_downsampled_for_iga_processed.RDS')

output_dir <- file.path(proj_dir, 'geomx-processing', 'results', 'batch2')

# load utils functions ----------------------------------------------------

source(file.path(proj_dir, 'st-processing', 'src', 'visium_utils.R')) #TODO add needed functions to geomx_utils
source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))

dir.create(output_dir, recursive = T, showWarnings = F)

# define intermediate output paths ----------------------------------------

geomx_qc_path <<- file.path(output_dir, 'geomx_qc.RDS')
geomx_norm_path <<- file.path(output_dir, 'geomx_qc_norm.RDS')
geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS')
geomx_deconvolution_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm_deconv.RDS')

# start the pipeline ------------------------------------------------------

print('#############')
print('GeoMx pipeline starting :O')
print('#############')

# conditionally run preprocessing -----------------------------------------

run_unless_exists('Preprocessing', geomx_qc_path, 
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_qc.R'))

# conditionally run normalisation -----------------------------------------

run_unless_exists('Normalisation', geomx_norm_path, 
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_normalisation.R'))


# conditonally run batch effect removal -----------------------------------

run_unless_exists('Batch effect removal', geomx_norm_batch_eff_rm_path, 
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_batch_effect_rmv.R'))


# conditionally run differential gene expression --------------------------

#TODO add comparison between 1 and 2nd group (eg 'CD8_IBA1' vs CD8-_IBA1)

comparison_type <<- 'within' 
# 'within' when you compare different ROI types within sample
# between - comparisons between slides
cofounder_name <<- 'Sample' # better don't change
main_var_name <<- 'Annotation_cell' # main variable to make comparison between
main_var_is_bin <<- TRUE # should variable be compared with all others at once (TRUE) or with each other separately
# if FALSE all labels in main_var_name will be compared as they are
main_var_main_val <- 'CD8_.*Iba1' # if main_var_is_bin - TRUE - name of the main value (or regex)
dge_categories <<- c('Segment', 'NACT_status') # categories to divide to when making DGE separately

dge_logs_path <<- file.path(output_dir, 'dge', 
                            paste0('dge_', comparison_type, '_slide_', main_var_name, 
                                   '_bin_', main_var_is_bin, '_', main_var_main_val, '_',
                                   paste0(dge_categories, collapse = '_'), 
                                   '_logs.txt'))

run_unless_exists('Differential Gene Expression', dge_logs_path, 
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_dge.R'))


# conditionally run deconvolution -----------------------------------------

scrna_anno <<- 'mid_lvl_ct' # either 'cell_type' or 'mid_lvl_ct'
# column name of cell type label in scRNAseq metadata

run_unless_exists('Deconvolution', geomx_deconvolution_path, 
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_deconvolution.R'))
