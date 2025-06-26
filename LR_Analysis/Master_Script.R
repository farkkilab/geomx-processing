# Master Script

# load libraries

library(BulkSignalR, quietly =T)
library(igraph, quietly =T)
library(dplyr, quietly =T)
library(DESeq2, quietly =T)
library(CellChat, quietly =T)
library(patchwork, quietly =T)
library(rlang)
library(purrr)
#library(tidyverse) # check
library(parallel)

library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(scales)


# TODO check the libraries needed

library(Biobase, quietly =T)
library(NanoStringNCTools, quietly =T)
library(BiocGenerics, quietly =T)
library(S4Vectors, quietly =T)
library(stats4, quietly =T)
library(GeomxTools, quietly =T)
library(ggplot2, quietly =T)


# For MultiNicheNetr
library(SingleCellExperiment, quietly =T)
library(nichenetr, quietly =T)
library(multinichenetr, quietly =T)



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
ligand_target_matrix = readRDS(file.path(data_dir,"Nichenet_Model","ligand_target_matrix_nsga2r_final.rds")) # NicheNet modal
lr_network_all = readRDS(file.path(data_dir,"Nichenet_Model","lr_network_human_allInfo_30112033.rds")) # NicheNetR LR network
pathway <- read.csv(file.path(data_dir,"pathway_names.csv")) # pathways for plotting: A list of reactome pathways in a .csv file
# TODO If not provided plot for all pathways : Need to adjust the size of the pdf 


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
comparison <- c("tumor","stroma") # define the groups from  "grouping_var_col_ids" that needed to be compared eg: c("pre-stroma","pre-tumor") order does not matter



# TODO define above as a tibble

# parameters for CellChat and MultiNicheNet : Single cell approaches
cell_types = c("Tcells","Macrophages") # set to NULL to get all the cell types : ct_of_interest


# parameters for MultiNicheNet only
 
batches = NA # this did not work
covariates =  "Sample" #  "Patient" if paired



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

source(file.path(proj_dir, 'LR_Analysis', 'BulkSignaR_LR_Visualization.R'))

# conditionally run CellChat LR Analysis -----------------------------------------

run_unless_exists('CellChat LR Analysis', geomx_CellChat_path, 
                  file.path(proj_dir, 'LR_Analysis', 'CellChat_LR_Analysis.R'))


# CellChat Visualization -----------------------------------------

source(file.path(proj_dir, 'LR_Analysis', 'CellChat_LR_Visualization.R'))



# conditionally run MultiNicheNet LR Analysis -----------------------------------------

run_unless_exists('MultiNicheNet LR Analysis', geomx_MultiNicheNet_path,
                  file.path(proj_dir, 'LR_Analysis', 'MultiNicheNet_LR_Analysis.R'))


# MultiNicheNet Visualization -----------------------------------------


source(file.path(proj_dir, 'LR_Analysis', 'MultiNicheNet_LR_Visualization.R'))
