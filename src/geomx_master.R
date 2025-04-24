# main packages for all scripts
library(plyr, quietly =T)
library(dplyr, quietly =T)
library(data.table, quietly =T)
library(tibble, quietly =T)
library(tools, quietly = T)
library(parallel, quietly = T)

library(ggforce, quietly =T)
library(ggplot2, quietly =T)
library(cowplot, quietly =T)
library(ggrepel, quietly =T)
library(reshape2, quietly =T)

# problems with Matrix package - i has to be lower that 1.7 to work with lmer
# devtools::install_version("Matrix","1.6.4")
# dependencies:  c("grDevices", "graphics", "grid", "lattice", "methods",  "stats", "utils")  
library(Matrix)

# script specific packages
library(Biobase, quietly =T)
library(NanoStringNCTools, quietly =T)
library(GeomxTools, quietly =T)
library(GeoDiff, quietly =T)
library(DESeq2, quietly =T)
library(SpatialDecon, quietly =T)
#install preprocessCore manually from source
# BiocManager::install("preprocessCore", configure.args = c(preprocessCore = "--disable-threading"), 
# force= TRUE, update=TRUE, type = "source")
library(preprocessCore, quietly =T)
library(umap, quietly =T)
library(Rtsne, quietly =T)

library(limma)
library(pvca)
library(harmony)

library(Seurat, quietly =T)
library(BayesPrism, quietly =T)
library(biomaRt, quietly =T)
library(msigdbr, quietly =T)

library(GSVA, quietly =T)
library(clusterProfiler, quietly =T)
library(progeny, quietly =T)

# or batch1, batch2, batch12 etc
batch <<- 'demo' # just for running demo data


# define variables and paths ----------------------------------------------

proj_dir <<- '~/Documents/phd/st/geomx-processing'
pdrive_dir <<- '/media/Pdrive/h30492/farkkilab2/9_EyeMT'

data_dir <- file.path(proj_dir, 'demo_data')
output_dir <<- file.path(proj_dir, 'results', 'demo_batch')

# input data
dcc_path <<- dir(file.path(data_dir, "dcc"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
pkc_path <<- file.path(data_dir, 'metadata', 'Hs_R_NGS_WTA_v1.0.pkc')

# anno file have to contain sheet named 'Sheet1' and following column names:
# 'Sample_ID', 'Slide_Name',  'Aoi', 'Roi' and 'Panel'
# and dcc_name of proper NTC in 'NTC_ID' column if theres no 1NTC/batch
anno_path <<- file.path(data_dir, 'metadata', 'dcc_metadata_demobatch.xlsx') 

# path to reference scRNAseq dataset for deconvolution
# have to contain 'cell_type' column name in metadata
scrna_ref_path <<- file.path(pdrive_dir, '9_EyeMT_reference_scRNAseq/vaharautio_scrnaseq_dataset_downsampled_for_iga_processed.RDS')
scrna_ref_path <<- '/home/iganiemi/Documents/phd/st/data/scrna/vaharautio_scrnaseq_dataset_downsampled_for_iga_processed.RDS'

# path to csv file with custom gene signatures
custom_sign_path <<- file.path(data_dir, 'signatures', 'ct_markers.csv')


# set up metadata variables names -----------------------------------------

aoi_id <<- 'dcc_filename'
roi_id <<- 'Roi'

main_batch_var <- 'main_batch_nr'
batch_var <<- 'batch_nr'


aoi_segment_var <<- "Segment"
main_roi_label <<- "Annotation_cell" 
main_experimental_condition <<- NULL # usually 'NACT_status', but NULL for demo data
sample_name <<- 'Sample'

other_vars_bio <<- c("Patient", "Site")
other_vars_tech <<- c('Slide_Name')

# load util functions and create dirs -------------------------------------

source(file.path(proj_dir, 'src', 'geomx_utils.R'))

dir.create(output_dir, recursive = T, showWarnings = F)

# define intermediate output paths ----------------------------------------

geomx_qc_path <<- file.path(output_dir, 'geomx_qc.RDS')
geomx_norm_path <<- file.path(output_dir, 'geomx_qc_norm.RDS')
geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS')

# start the pipeline ------------------------------------------------------

print('#############')
print('GeoMx pipeline starting :O')
print('#############')

# conditionally run preprocessing -----------------------------------------

run_unless_exists('Preprocessing', geomx_qc_path, 
                  file.path(proj_dir, 'src', 'geomx_qc.R'))

# conditionally run normalisation -----------------------------------------

run_unless_exists('Normalisation', geomx_norm_path, 
                  file.path(proj_dir, 'src', 'geomx_normalisation.R'))


# conditonally run batch effect removal -----------------------------------
# !!! check throughfully the 1st PVCA plots - if another variables are responsible for variance 
# primary_batch_var and secondary_batch_var values should be changed
# secondary batch variable has to be INDEPENDENT from the primary_batch_var

# should be the same as in batch effect rm script
# if analysing each batch separately, only batch_var is considered
# if analysing many big batches together, both main_batch_var and batch_var are considered
primary_batch_var <<- ifelse(batch %in% c('batch1', 'batch2', 'batch3'), batch_var, main_batch_var)
if(batch %in% c('batch1', 'batch2', 'batch3')){secondary_batch_var <<- NULL} else{secondary_batch_var <<- batch_var}

run_unless_exists('Batch effect removal', geomx_norm_batch_eff_rm_path, 
                  file.path(proj_dir, 'src', 'geomx_batch_effect_rmv.R'))

# conditionally run deconvolution -----------------------------------------

# column name of cell type label in scRNAseq metadata
scrna_anno <<- 'mid_lvl_ct' # either 'cell_type' or 'mid_lvl_ct'

deconv_logs_path <<- file.path(output_dir,'deconvolution', 
                             paste0('deconv_', scrna_anno, '_logs.txt'))

run_unless_exists('Deconvolution', deconv_logs_path, 
                  file.path(proj_dir, 'src', 'geomx_deconvolution.R'))

# conditionally run pathway analysis --------------------------------------

# cell type column name in reference scRNAseq dataset
scrna_anno <<- 'mid_lvl_ct' # either 'cell_type' or 'mid_lvl_ct'

pathway_inp_data_type <<- c('all', 'bp') # within c('all', 'bp')
# all - full geomx data (not-deconvoluted)
# bp - bayes prism deconvoluted data

ct_of_interest <<- c("tumor", "Tcells", "Bcells", "Fibroblasts", "NKcells", 
                     "Macrophages", "DCs", "Endothelial cells")
# if running for 'bp' (bayes prism deconvolution results) 
# specifies for which cell types GSEA should be computed (as in scrna_anno column in scRNAseq reference ds)
# if ct_of_interest <<- NULL - GSEA will be computed for all cell types

signature_type <<- 'custom' # c('msigdb', 'custom')
# msigdb - on all pathways from msigdb (Hallmark + CP)
# custom - on custom signatures list specified in custom_sign_path

gsea_type <<- 'ssgsea' # 'gsva' or 'ssgsea'

signature_name <<- ifelse(signature_type == 'custom', gsub('.csv', '', basename(custom_sign_path)), '')

gsea_logs_path <<- file.path(output_dir,'pathway_analysis', 'gsea', 
                             paste0(gsea_type, '_', signature_type, signature_name, '_logs.txt'))


run_unless_exists('Pathway analysis', gsea_logs_path, 
                  file.path(proj_dir, 'src', 'geomx_pathway_analysis.R'))


# conditionally run differential gene expression --------------------------
scrna_anno <<- 'mid_lvl_ct' # either 'cell_type' or 'mid_lvl_ct'

dge_inp_data_type <<- c('all', 'bp') # within c('all', 'bp')
# all - full geomx data (not-deconvoluted)
# bp - bayes prism deconvoluted data

ct_of_interest <<- c("tumor", "Tcells", "Fibroblasts", "Macrophages", "Endothelial cells")
# if running for 'bp' (bayes prism deconvolution results) 
# specifies for which cell types GSEA should be computed (as in scrna_anno column in scRNAseq reference ds)
# if ct_of_interest <<- NULL - GSEA will be computed for all cell types

# DGE parameters
comparison_type <<- 'within' 
# 'within' when you compare different ROI types within sample
# between - comparisons between slides
main_var_name <<- 'Annotation_cell' # main variable to make comparison between
main_var_is_bin <<- TRUE # should variable be compared with all others at once (TRUE) or with each other separately
# if FALSE all labels in main_var_name will be compared as they are

main_var_main_val <<- 'CD8' # if main_var_is_bin - TRUE - name of the main value (or regex - careful!)
dge_categories <<- c('Segment', 'NACT_status') # categories to divide to when making DGE separately


# don't change it - identifier of dge run
dge_name <<- paste0('dge_', comparison_type, '_slide_', main_var_name, 
                   '_bin_', main_var_is_bin, '_', gsub('\\*', '', main_var_main_val), '_',
                   paste0(dge_categories, collapse = '_'))

dge_logs_path <<- file.path(output_dir, 'dge', dge_name, 'dge_logs.txt')

run_unless_exists('Differential Gene Expression', dge_logs_path, 
                  file.path(proj_dir, 'src', 'geomx_dge.R'))



