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
output_dir <<- file.path(proj_dir, 'results', 'Batch01') 


# load data ----------------------------------------------


geomx_obj = readRDS(paste0(data_dir,"/geomx_qc_norm_batch_eff_rm.RDS"))
bprism_res = readRDS(paste0(data_dir,"/bp_res_mid_lvl_ct.RDS"))
cell_fractions_df = read.csv(paste0(data_dir,"/bp_res_mid_lvl_ct_ct_fraction.csv"))


# set up metadata variables names -----------------------------------------
aoi_id <- 'dcc_filename'
sample_name <- 'Sample'
aoi_segment_var <<- "Segment"
main_experimental_condition <<- 'NACT_status' # eg 'NACT_status', but NULL for demo data
grouping_var_col_ids <- c("Segment") # define the meta data column names of the groups that needed to be compared eg: c("NACT_status","Segment")

# params  -----------------------------------------

# parameters for BulkSignalR

comparison <- c("stroma","tumor") # order of the group should matches with the order of the column names eg: c("pre-stroma","pre-tumor")

# TODO Run for combined data 
combined_Data = FALSE # Run for combined data


# parameters for CellChat

cell_types = c("Tcells","Macrophages") # set to NULL to get all the cell types : ct_of_interest


# parameters for MultiNicheNet




# load util functions and create dirs -------------------------------------

source(file.path(proj_dir, 'LR_Analysis', 'BulkSignalR_LR_utils.R'))
source(file.path(proj_dir, 'LR_Analysis', 'cellChat_util.R'))

dir.create(output_dir, recursive = T, showWarnings = F)

# define intermediate output paths ----------------------------------------



# if (combined_Data == TRUE) {
#   
#   geomx_BulkSignalR_path <<- file.path(output_dir,'BulkSignalR_objects' ,'BulkSignalR_combined_output.RDS')
#   
# } else {
#   
#   file_name = paste(parts_non_NA, collapse = "_")
#   geomx_BulkSignalR_path <<- file.path(output_dir, 'BulkSignalR_objects',paste0('BulkSignalR_',file_name,'_output.RDS'))
#   
# }

geomx_CellChat_path <<- file.path(output_dir, 'CellChat_objects','CellChat_output.RDS')
geomx_MultiNicheNet_path <<- file.path(output_dir, 'MultiNicheNet_objects','multinichenet_output.rds')

# conditionally run BulkSignaR Analysis  -----------------------------------------

run_unless_exists('BulkSignaR LR Analysis', geomx_BulkSignalR_path,
                  file.path(proj_dir, 'LR_Analysis', 'BulkSignaR_LR_Analysis.R'))


# BulkSignalR Visualization -----------------------------------------



# conditionally run CellChat LR Analysis -----------------------------------------

run_unless_exists('CellChat LR Analysis', geomx_CellChat_path, 
                  file.path(proj_dir, 'LR_Analysis', 'CellChat_LR_Analysis.R'))


# CellChat Visualization -----------------------------------------




# conditionally run MultiNicheNet LR Analysis -----------------------------------------

run_unless_exists('MultiNicheNet LR Analysis', geomx_MultiNicheNet_path,
                  file.path(proj_dir, 'LR_Analysis', 'MultiNicheNet_LR_Analysis.R'))




