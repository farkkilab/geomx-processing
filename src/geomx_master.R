# main packages for all scripts
library(plyr, quietly =T)
library(dplyr, quietly =T)
library(data.table, quietly =T)
library(tibble, quietly =T)
library(tools, quietly = T)

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
# TODO move to renv

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

# TODO make 1 parameter for important metadata column names reused in many scripts
# TODO optimise all output paths and logs to contqain all important infor about the run

# library(ggpubr)
# library(topGO)
# library(fgsea)
# library(fpc)
# library(dbscan)

# define variables and paths ----------------------------------------------

# TODO all the batches should be merged and qc-ed + processed together and bigbatch + smallbatch variable as batch effects
batch <<- 'batch1' # just for running slightly different batches separately

proj_dir <<- '~/Documents/phd/st'

#data_dir <<- '~/Documents/phd/st/data/geomx/geomx_batch2_1124/' # batch2 
data_dir <<- '~/Documents/phd/st/data/geomx/geomx_batch1_nact' # batch1

#output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch2-1802') # batch2
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch1-1903') # batch1

# input data
dcc_path <<- dir(file.path(data_dir, "dcc"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
pkc_path <<- file.path(data_dir, 'metadata', 'Hs_R_NGS_WTA_v1.0.pkc')

# anno file have to contain sheet named 'Sheet1' and following column names:
# 'Sample_ID', 'Slide_Name',  'Aoi', 'Roi' and 'Panel' 'dcc_filename' (main id of AOI)
# and dcc_name of proper NTC in 'NTC_ID' column if theres no 1NTC/batch
# anno_path <<- file.path(data_dir, 'metadata', 'dcc_metadata_all_batch2_1124.xlsx') #batch2
 anno_path <<- file.path(data_dir, 'metadata', 'dcc_metadata_all_cleaned.xlsx') #batch1

# path to reference scRNAseq dataset for deconvolution
# have to contain 'cell_type' column name
scrna_ref_path <<- file.path(proj_dir, 'data/scrna/vaharautio_scrnaseq_dataset_downsampled_for_iga_processed.RDS')

# path to csv file with custom gene signatures
custom_sign_path <<- file.path(proj_dir, 'geomx-processing', 'data', 'signatures',
                              'stromal_cell_subtype_signatures_symbols_ensembl_ids_revised.csv')


# set up metadata variables names -----------------------------------------

aoi_id <<- 'Sample_ID'
roi_id <<- 'Roi'
slide_id <- 'Slide_Name'

aoi_segment_var <<- "Segment"
main_roi_label <<- "Annotation_cell" 
main_experimental_condition <<- 'NACT_status'
sample_name <<- 'Sample'

batch_var <<- 'batch_nr'

# if analysing 1 batch separately
main_batch_var <<- batch_var
secondary_batch_var <<- NULL

# if analysisng many big batches together
# main_batch_var <- 'main_batch_nr'
# secondary_batch_var <- batch_var

other_vars_bio <<- c("Segment_geomx", "Patient", "Site", 'PFS', 'PFS_months')
other_vars_tech <<- c(slide_id, "batch_nr_sample_collection")

# load util functions and create dirs -------------------------------------

source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))

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
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_qc.R'))

# conditionally run normalisation -----------------------------------------

run_unless_exists('Normalisation', geomx_norm_path, 
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_normalisation.R'))


# conditonally run batch effect removal -----------------------------------

# TODO compute voom() weights for dge?
# The primary purpose of the voom() function is to compute precision 
# weights for the downstream differential expression analysis.

run_unless_exists('Batch effect removal', geomx_norm_batch_eff_rm_path, 
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_batch_effect_rmv.R'))

# conditionally run deconvolution -----------------------------------------

# TODO remove artifact - peak of low counts genes after vst
# column name of cell type label in scRNAseq metadata
scrna_anno <<- 'mid_lvl_ct' # either 'cell_type' or 'mid_lvl_ct'

deconv_bp_harm_path <<- file.path(output_dir,'deconvolution', 'bayes_prism', 
                                paste0('bp_res_', scrna_anno, '_expr_mtx_cleaned_vst_harmony_batch_corr.RDS'))

deconv_bp_limma_path <<- file.path(output_dir,'deconvolution', 'bayes_prism', 
                                  paste0('bp_res_', scrna_anno, '_expr_mtx_cleaned_vst_limma_batch_corr.RDS'))

deconv_sd_path <<- file.path(output_dir,'deconvolution', 'spatial_decon', 
                            paste0('sd_res_bg', scrna_anno, '_geomxfilt_ct_fraction.RDS'))

run_unless_exists('Deconvolution', deconv_sd_path, 
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_deconvolution.R'))

# conditionally run pathway analysis --------------------------------------

# TODO add limma fry calculation - another algorithm for pathway analysis not super important

pathway_inp_data_type <<- c('all', 'bp') # within c('all', 'bp')
# all - full geomx data (not-deconvoluted)
# bp - bayes prism deconvoluted data

signature_type <<- 'msigdb' # c('msigdb', 'custom')
# msigdb - on all pathways from msigdb (Hallmark + CP)
# custom - on custom signatures list

gsea_type <<- 'ssgsea' # 'gsva' or 'ssgsea'

gsea_logs_path <<- file.path(output_dir,'pathway_analysis', 'gsea', 
                             paste0(gsea_type, '_', signature_type,'.txt'))


run_unless_exists('Pathway analysis', gsea_logs_path, 
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_pathway_analysis.R'))


# conditionally run differential gene expression --------------------------

# TODO anova(full model, reduced model) - check if significantly improves the effect for interesting genes
# TODO add limma voom - not so important
# https://davislaboratory.github.io/GeoMXAnalysisWorkflow/articles/GeoMXAnalysisWorkflow.html#batch-correction

dge_inp_data_type <<- c('all', 'bp') # within c('all', 'bp')
# all - full geomx data (not-deconvoluted)
# bp - bayes prism deconvoluted data

# DGE parameters
comparison_type <<- 'within' 
# 'within' when you compare different ROI types within sample
# between - comparisons between slides
main_var_name <<- 'Annotation_cell' # main variable to make comparison between
main_var_is_bin <<- TRUE # should variable be compared with all others at once (TRUE) or with each other separately
# if FALSE all labels in main_var_name will be compared as they are
main_var_main_val <<- 'CD8_.*Iba1' # if main_var_is_bin - TRUE - name of the main value (or regex - careful!)
dge_categories <<- c('Segment', 'NACT_status') # categories to divide to when making DGE separately

# don't change it - identifier of dge run
dge_name <<- paste0('dge_', comparison_type, '_slide_', main_var_name, 
                   '_bin_', main_var_is_bin, '_', gsub('\\*', '', main_var_main_val), '_',
                   paste0(dge_categories, collapse = '_'))

dge_logs_path <<- file.path(output_dir, 'dge', dge_name, 'dge_logs.txt')

run_unless_exists('Differential Gene Expression', dge_logs_path, 
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_dge.R'))



