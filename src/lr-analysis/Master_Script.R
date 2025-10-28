# Master Script

#devtools::install_github('immunogenomics/presto')

# load libraries

library(BayesPrism)
library(DESeq2)
library(ggplot2)
library(tools)
library(parallel)

library(BulkSignalR, quietly =T)
library(igraph, quietly =T)
library(dplyr, quietly =T)
library(scales, quietly =T)
library(circlize, quietly =T)
library(pheatmap, quietly =T)
library(ComplexHeatmap, quietly =T)

library(NMF)
library(DESeq2, quietly =T)
library(CellChat, quietly =T)
library(patchwork, quietly =T)
library(ggh4x)
library(rlang)
library(presto)

#library(tidyverse) # check
library(stringr)

# To load Geomx object
library(GeomxTools, quietly =T)
library(Biobase, quietly =T)
library(BiocGenerics, quietly =T)
library(NanoStringNCTools, quietly =T)
library(S4Vectors, quietly =T)
library(stats4, quietly =T)
library(ggplot2, quietly =T)


# For MultiNicheNetr
library(SingleCellExperiment, quietly =T)
library(nichenetr, quietly =T)
library(multinichenetr, quietly =T)
library(parallel)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
library(data.table)
library(dplyr)
library(tibble)

#TODO make umap of pseudosc gene set

# define intermediate output folders and paths  ----------------------------------------

# TODO make plotdirs here/inside downstream scripts
dir.create(file.path(output_dir, 'lr_interactions', 'bulk_signalr') , recursive = T, showWarnings = F)
dir.create(file.path(output_dir, 'lr_interactions', 'cell_chat') , recursive = T, showWarnings = F)
dir.create(file.path(output_dir, 'lr_interactions', 'multi_niche_netr') , recursive = T, showWarnings = F)

geomx_BulkSignalR_path <<- file.path(output_dir, 'lr_interactions', 'bulk_signalr' ,'BulkSignalR_output.RDS')
geomx_CellChat_path <<- file.path(output_dir, 'lr_interactions', 'cell_chat','CellChat_output.RDS')
geomx_MultiNicheNet_path <<- file.path(output_dir, 'lr_interactions', 'multi_niche_netr', 'multinichenet_output.rds')

# load util functions  -------------------------------------

# TODO move all utils to 1 script
source(file.path(proj_dir, 'geomx-processing', 'src','lr-analysis', 'BulkSignalR_LR_utils.R'))
source(file.path(proj_dir, 'geomx-processing', 'src','lr-analysis', 'cellChat_util.R'))


# define params -----------------------------------------------------------

#TODO move to the main script
# column name of cell type label in scRNAseq metadata
scrna_anno <<- 'mid_lvl_ct_updated' #either 'cell_type' / 'mid_lvl_ct' / 'mid_lvl_ct_updated' / 'low_lvl_ct'

# common parameters
grouping_var_col_ids <- c("Segment") # define the meta data column names of the groups that needed to be compared separately eg: c("Segment","NACT_status")

# TODO BSR - handled with code - check for NN and CC
# define the groups from  "grouping_var_col_ids" that needed to be compared eg: c("pre_stroma","pre_tumor") order matters. 
# Can compare only two groups at a time
comparison <- c("stroma","tumor") 

# parameters for CellChat and MultiNicheNet : Single cell approaches
# names of cells to fin
cell_types = c("Tcells_CD8","Macrophages_Monocytes") # set to NULL to get all the cell types : ct_of_interest


# parameters for MultiNicheNet only

# covariates for EdgeR DEGs calculated by MultiNicheNet. 
# How to define the covariate: If the defined covaraite id not present in both groups of interest EdgeR will not run.
# therefore in such case leave the covariate to "NA" 
covariates =  "Sample"


# load data ----------------------------------------------

# TODO move to low-lvl scripts
geomx_obj <<- readRDS(geomx_norm_batch_eff_rm_path) # batch effect corrected Geomx Object

bprism_res <<- readRDS(file.path(output_dir, 'deconvolution', 'bayes_prism', 
                               paste0('bp_res_', scrna_anno, '.RDS'))) # raw BayesPrism result Object

cell_fractions_df <<- read.csv(file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                       paste0('bp_res_', scrna_anno, '_ct_fraction.csv'))) # cell fractions from BayesPrism

ligand_target_matrix <<- readRDS(file.path(nichenet_data_dir,"ligand_target_matrix_nsga2r_final.rds")) # NicheNet model
lr_network_all <<- readRDS(file.path(nichenet_data_dir,"lr_network_human_allInfo_30112033.rds")) # NicheNetR LR network

# set up metadata variables names -----------------------------------------
# already defined in a main script

#proj_dir <<- 'C:/Users/Sahas/Downloads/Masters_Thesis/Project_LR_prediction'
#output_dir <<- file.path(proj_dir, 'results', 'Batch01','LR_prediction') 

# aoi_id <- 'dcc_filename'
# sample_name <- 'Sample'
# aoi_segment_var <- "Segment"
# main_experimental_condition <- 'NACT_status' # eg 'NACT_status', but NULL for demo data
# other_vars_bio <- c("Patient", "Site")
# other_vars_tech <- c('Slide_Name')



# conditionally run BulkSignaR Analysis  -----------------------------------------

# TODO this have to be re-written
#paired_id <- "Patient" # If you want to predict for paired samples only in BUlkSignalR: I did not check this 

qval_threshold = 0.01 # thr for filtering significant LR pairs from LR output dfs
grouping_var_col_ids <- c("Segment") # define the meta data column names of the groups that needed to be compared separately eg: c("Segment","NACT_status")


lr_output_path <- file.path(output_dir, 'lr_interactions', 'bulk_signalr', 
                            paste0('LR_list_all_qval_', qval_threshold, '.RDS'))

run_unless_exists('BulkSignaR LR Analysis', geomx_BulkSignalR_path,
                  file.path(proj_dir, 'geomx-processing', 'src','lr-analysis', 'BulkSignaR_LR_Analysis.R'))


# BulkSignalR Visualization -----------------------------------------

# names of interesting pathways from REACTOME+GO:BP for BulkSignalR plotting: 
# A list of reactome pathways in a .csv file. This needed to be provided to plot the heatmap
pathway <<- read.csv(file.path(proj_dir, 'geomx-processing', 'data', 'signatures', 'immune_signatures_selected_forpaper_names_reactome_gobp.csv')) 

# for plots
# can be hardcoded inside
plot_dir = file.path(output_dir, 'lr_interactions', 'bulk_signalr', 'plots_and_csv_files')
dir.create(plot_dir , recursive = T, showWarnings = F)

# Params
# TODO make it different from qval from Analysis
qval_threshold = 0.001 # filter significant LR pairs
n = 50 # number of top LR pairs needed to visualize in the signature scoes heatmap
heatmap_col_ann = "Segment" # based on what you want to annotate the heatmap
manually_filtered_BulkSignalr_df = NULL# readRDS(file.path(output_dir,"df_combined_BulkSignalr_for_plot.RDS")) # to plot the bubble plot: If you want to visulaze your own filtered dataframe provide the dataframe  


source(file.path(proj_dir,'geomx_processing', 'src',  'lr-analysis', 'BulkSignaR_LR_Visualization.R')) # ??
source("/home/iganiemi/Documents/phd/st/geomx-processing/src/lr-analysis/BulkSignaR_LR_Visualization.R")

# conditionally run CellChat LR Analysis -----------------------------------------

scrna_anno <<- 'mid_lvl_ct_updated' #either 'cell_type' / 'mid_lvl_ct' / 'mid_lvl_ct_updated' / 'low_lvl_ct'
sample_name <- 'Sample'
# common parameters
grouping_var_col_ids <- c("Segment") # define the meta data column names of the groups that needed to be compared separately eg: c("Segment","NACT_status")

# parameters for CellChat and MultiNicheNet : Single cell approaches
# names of cells to fin
cell_types_selected = c("Tcells_CD8","Macrophages_Monocytes") # set to NULL to get all the cell types : ct_of_interest

# parameters for Cellchat
cell_frac_cutoff = 0.005 # 0.01 or 0.005
min_cells = 10 # Number of minimum cells in each cell group

outct <- ifelse(!is.null(cell_types_selected), paste(cell_types_selected, collapse = '_'), 'all')

# path to cellchat output file
geomx_CellChat_path <- file.path(output_dir, 'lr_interactions', 'cell_chat', 
                                 paste0('CellChat_output_',paste(grouping_var_col_ids, collapse = '_'),
                                        '_', outct, '.RDS'))

# path to output lr dataframe
cc_lr_df_path <- file.path(output_dir, 'lr_interactions', 'cell_chat', 
                           paste0('CellChat_df_',paste(grouping_var_col_ids, collapse = '_'),
                                  '_', outct, '_lr.csv'))

cc_path_df_path <- file.path(output_dir, 'lr_interactions', 'cell_chat', 
                             paste0('CellChat_df_', paste(grouping_var_col_ids, collapse = '_'),
                                    '_', outct, '_pathway.csv'))

run_unless_exists('CellChat LR Analysis', geomx_CellChat_path, 
                  file.path(proj_dir, 'LR_Analysis', 'CellChat_LR_Analysis.R'))


# CellChat Visualization -----------------------------------------

# for plots
# can be hardcoded inside
plot_dir = file.path(output_dir, 'lr_interactions', 'cell_chat', 'plots_and_csv_files')
dir.create(plot_dir , recursive = T, showWarnings = F)

# Params

manually_filtered_cellchat_df = NULL # to plot the bubble plot: If you want to visulaze your own filtered dataframe provide the dataframe  eg: readRDS(file.path(output_dir,"df_combined_cellchat_for_plot.RDS"))
pval_threshold = 0.01
prob_threshold = 0.1

source(file.path(proj_dir, 'LR_Analysis', 'CellChat_LR_Visualization.R'))


# conditionally run MultiNicheNet LR Analysis -----------------------------------------

scrna_anno <<- 'mid_lvl_ct_updated' #either 'cell_type' / 'mid_lvl_ct' / 'mid_lvl_ct_updated' / 'low_lvl_ct'
sample_name <- 'Sample'
# common parameters
grouping_var_col_ids <- c("Segment") # define the meta data column names of the groups that needed to be compared separately eg: c("Segment","NACT_status")

# parameters for CellChat and MultiNicheNet : Single cell approaches
# names of cells to fin
cell_types_selected = c("Tcells_CD8","Macrophages_Monocytes") # set to NULL to get all the cell types : ct_of_interest

# parameters for Cellchat
cell_frac_cutoff = 0.005 # 0.01 or 0.005
min_cells = 10 # Number of minimum cells in each cell group

# path where to download nichenet data
nichenet_data_dir <<- file.path(proj_dir, 'geomx-processing', 'data', 'nichenet')
dir.create(nichenet_data_dir, recursive = T, showWarnings = F)

plot_dir = file.path(output_dir, 'lr_interactions', 'multi_niche_netr','plots_and_csv_files')
dir.create(plot_dir , recursive = T, showWarnings = F)


run_unless_exists('MultiNicheNet LR Analysis', geomx_MultiNicheNet_path,
                  file.path(proj_dir, 'LR_Analysis', 'MultiNicheNet_LR_Analysis.R'))


# MultiNicheNet Visualization -----------------------------------------

# for plots
manually_filtered_LR_pairs_dfplot_median_bulk_expr = NULL # provide the dataframe  
manually_filtered_LR_pairs_dfplot_ligand_activity = NULL # provide the dataframe  

source(file.path(proj_dir, 'LR_Analysis', 'MultiNicheNet_LR_Visualization.R'))
