# Master Script

# load libraries

library(BulkSignalR, quietly =T)
library(igraph, quietly =T)
library(dplyr, quietly =T)
library(scales, quietly =T)
library(circlize, quietly =T)
library(pheatmap, quietly =T)
library(ComplexHeatmap, quietly =T)

# TODO check the libraries needed
library(DESeq2, quietly =T)
library(CellChat, quietly =T)
library(patchwork, quietly =T)
library(ggh4x)
library(rlang)

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


# define variables and paths ----------------------------------------------

proj_dir <<- 'C:/Users/Sahas/Downloads/Masters_Thesis/Project_LR_prediction'


data_dir <- file.path(proj_dir, 'Batch01_Data') 
output_dir <<- file.path(proj_dir, 'results', 'Batch01','LR_prediction') 


# load data ----------------------------------------------

# path to Geomx object 
# path to BayesPrism Object 
# path to cell fractions from BayesPrism
# path to NicheNet modal
# path to NicheNet lr network

geomx_obj = readRDS(file.path(data_dir,"geomx_qc_norm_batch_eff_rm.RDS")) # Geomx Object
bprism_res = readRDS(file.path(data_dir,"bp_res_mid_lvl_ct.RDS")) # BayesPrism Object
cell_fractions_df = read.csv(file.path(data_dir,"bp_res_mid_lvl_ct_ct_fraction.csv")) # cell fractions from BayesPrism
ligand_target_matrix = readRDS(file.path(data_dir,"Nichenet_Model","ligand_target_matrix_nsga2r_final.rds")) # NicheNet modal. Can be downloaded from MultiNicheNet repo : "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
lr_network_all = readRDS(file.path(data_dir,"Nichenet_Model","lr_network_human_allInfo_30112033.rds")) # NicheNetR LR network. Can be downloaded from MultiNicheNet repo : "https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds"
pathway <- read.csv(file.path(data_dir,"pathway_names.csv")) # pathways for plotting: A list of reactome pathways in a .csv file. This needed to be provided to plot the heatmap


# set up metadata variables names -----------------------------------------
aoi_id <- 'dcc_filename'
sample_name <- 'Sample'
aoi_segment_var <- "Segment"
main_experimental_condition <- 'NACT_status' # eg 'NACT_status', but NULL for demo data
other_vars_bio <- c("Patient", "Site")
other_vars_tech <- c('Slide_Name')



# params  -----------------------------------------

# common parameters


grouping_var_col_ids <- c("Segment") # define the meta data column names of the groups that needed to be compared eg: c("Segment","NACT_status")

# define the groups from  "grouping_var_col_ids" that needed to be compared eg: c("pre_stroma","pre_tumor") order matters. 
# Can compare only two groups at a time
comparison <- c("stroma","tumor") 


# TODO define above as a tibble

# parameters for CellChat and MultiNicheNet : Single cell approaches
cell_types = c("Tcells","Macrophages") # set to NULL to get all the cell types : ct_of_interest


# parameters for MultiNicheNet only
 
# covariates for EdgeR DEGs calculated by MultiNicheNet. 
# How to define the covariate: If the defined covaraite id not present in both groups of interest EdgeR will not run.
# therefore in such case leave the covariate to "NA" 
covariates =  "Sample"


# define intermediate output folders and paths  ----------------------------------------

BulkSignalR_folder_name = 'BulkSignalR_outputs'
CellChat_folder_name = 'CellChat_outputs'
MultiNicheNet_folder_name = 'MultiNicheNet_outputs'


geomx_BulkSignalR_path <<- file.path(output_dir,BulkSignalR_folder_name ,'BulkSignalR_output.RDS')
geomx_CellChat_path <<- file.path(output_dir, CellChat_folder_name,'CellChat_output.RDS')
geomx_MultiNicheNet_path <<- file.path(output_dir, MultiNicheNet_folder_name,'multinichenet_output.rds')


# create dirs -------------------------------------

dir.create(file.path(output_dir, BulkSignalR_folder_name) , recursive = T, showWarnings = F)
dir.create(file.path(output_dir, CellChat_folder_name) , recursive = T, showWarnings = F)
dir.create(file.path(output_dir, MultiNicheNet_folder_name) , recursive = T, showWarnings = F)



# load util functions  -------------------------------------

source(file.path(proj_dir, 'LR_Analysis', 'BulkSignalR_LR_utils.R'))
source(file.path(proj_dir, 'LR_Analysis', 'cellChat_util.R'))


# conditionally run BulkSignaR Analysis  -----------------------------------------

run_unless_exists('BulkSignaR LR Analysis', geomx_BulkSignalR_path,
                  file.path(proj_dir, 'LR_Analysis', 'BulkSignaR_LR_Analysis.R'))


# BulkSignalR Visualization -----------------------------------------

# for plots
plot_dir = file.path(output_dir, BulkSignalR_folder_name,'plots_and_csv_files_2')
dir.create(plot_dir , recursive = T, showWarnings = F)

# Params

qval_threshold = 0.001 # filter significant LR pairs
n = 50 # number of top LR pairs needed to visualize in the signature scoes heatmap
heatmap_col_ann = "Segment" # based on what you want to annotate the heatmap
LR_corr_threshold = 0.4 # For the bubble plot : correlation threshold
manually_filtered_BulkSignalr_df = NULL# readRDS(file.path(output_dir,"df_combined_BulkSignalr_for_plot.RDS")) # to plot the bubble plot: If you want to visulaze your own filtered dataframe provide the dataframe  

source(file.path(proj_dir, 'LR_Analysis', 'BulkSignaR_LR_Visualization.R'))

# conditionally run CellChat LR Analysis -----------------------------------------

run_unless_exists('CellChat LR Analysis', geomx_CellChat_path, 
                  file.path(proj_dir, 'LR_Analysis', 'CellChat_LR_Analysis.R'))


# CellChat Visualization -----------------------------------------

# for plots
plot_dir = file.path(output_dir, CellChat_folder_name,'plots_and_csv_files')
dir.create(plot_dir , recursive = T, showWarnings = F)

# Params

manually_filtered_cellchat_df = NULL # to plot the bubble plot: If you want to visulaze your own filtered dataframe provide the dataframe  eg: readRDS(file.path(output_dir,"df_combined_cellchat_for_plot.RDS"))
pval_threshold = 0.01
prob_threshold = 0.05

source(file.path(proj_dir, 'LR_Analysis', 'CellChat_LR_Visualization.R'))


# conditionally run MultiNicheNet LR Analysis -----------------------------------------

# for multiNicheNetR if you are providing DEGS externally follow the MultiNicheNet_LR_util.R script to prepare the DEGs dataframe. 
# Else MUltiNicheNet will not work
external_DE_info = FALSE # if TRUE,  provide the prepared DEGs dataframe to  'celltype_de_external'. Eg: celltype_de_external = readRDS(file.path(output_dir,MultiNicheNet_folder_name,"celltype_de_combined_calculated_externally.RDS"))
celltype_de_external = NULL 


run_unless_exists('MultiNicheNet LR Analysis', geomx_MultiNicheNet_path,
                  file.path(proj_dir, 'LR_Analysis', 'MultiNicheNet_LR_Analysis.R'))


# MultiNicheNet Visualization -----------------------------------------

# for plots

plot_dir = file.path(output_dir, MultiNicheNet_folder_name,'plots_and_csv_files_2')
dir.create(plot_dir , recursive = T, showWarnings = F)

manually_filtered_LR_pairs_dfplot_median_bulk_expr = NULL # provide the dataframe  
manually_filtered_LR_pairs_dfplot_ligand_activity = NULL # provide the dataframe  

source(file.path(proj_dir, 'LR_Analysis', 'MultiNicheNet_LR_Visualization.R'))
