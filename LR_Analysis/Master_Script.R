# Master Script

# load libraries

library(BulkSignalR, quietly =T)
library(igraph, quietly =T)
library(dplyr, quietly =T)
# library(Seurat, quietly =T)
# library(SeuratObject, quietly =T)
library(DESeq2, quietly =T)
library(CellChat, quietly =T)
library(patchwork, quietly =T)


# TODO check the libraries needed

library(Biobase, quietly =T)
library(NanoStringNCTools, quietly =T)
library(BiocGenerics, quietly =T)
library(S4Vectors, quietly =T)
library(stats4, quietly =T)
library(GeomxTools, quietly =T)
library(ggplot2, quietly =T)






# define variables and paths ----------------------------------------------

proj_dir <<- 'C:/Users/Sahas/Downloads/Masters_Thesis/Project_LR_prediction'


data_dir <- file.path(proj_dir, 'Batch01_Data') 
output_dir <<- file.path(proj_dir, 'results', 'Batch01') 


# params  -----------------------------------------

# parameters for BulkSignalR
# keep all as NA to get the LR prediction for the whole dataset

nact_status = NA # define as C("pre","post")
segment = c("stroma")
annotation = NA

# TODO Run for combined data 
combined_Data = FALSE # Run for combined data


# parameters for CellChat

cell_types = c("Tcells","Macrophages") # set to NULL to get all the cell types 


# parameters for MultiNicheNet




# load util functions and create dirs -------------------------------------

source(file.path(proj_dir, 'src', 'geomx_LR_utils.R'))
source(file.path(proj_dir, 'src', 'cellChat_util.R'))

dir.create(output_dir, recursive = T, showWarnings = F)

# define intermediate output paths ----------------------------------------

parts <- c(nact_status, segment, annotation)
parts_non_NA <- parts[!sapply(parts, is.na)]

if (length(parts_non_NA) == 0) {
  
  geomx_BulkSignalR_path <<- file.path(output_dir,'BulkSignalR_objects' ,'BulkSignalR_combined_output.RDS')
  
} else {
  
  file_name = paste(parts_non_NA, collapse = "_")
  geomx_BulkSignalR_path <<- file.path(output_dir, 'BulkSignalR_objects',paste0('BulkSignalR_',file_name,'_output.RDS'))
  geomx_CellChat_path <<- file.path(output_dir, 'CellChat_objects',paste0('CellChat_',file_name,'_output.RDS'))
}






# conditionally run BulkSignaR Analysis  -----------------------------------------

run_unless_exists('BulkSignaR LR Analysis', geomx_BulkSignalR_path, 
                  file.path(proj_dir, 'src', 'BulkSignaR_LR_Analysis.R'))


# BulkSignalR Visualization -----------------------------------------



# conditionally run CellChat LR Analysis -----------------------------------------

run_unless_exists('CellChat LR Analysis', geomx_CellChat_path, 
                  file.path(proj_dir, 'src', 'CellChat_LR_Analysis.R'))


# CellChat Visualization -----------------------------------------




# conditionally run MultiNicheNet LR Analysis -----------------------------------------






