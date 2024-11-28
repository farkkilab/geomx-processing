library(Biobase, quietly =T)
library(NanoStringNCTools, quietly =T)
library(GeomxTools, quietly =T)
library(GeoDiff, quietly =T)
library(Biobase, quietly =T)
library(DESeq2, quietly =T)
#install preprocessCore manually from source
# BiocManager::install("preprocessCore", configure.args = c(preprocessCore = "--disable-threading"), 
# force= TRUE, update=TRUE, type = "source")
library(preprocessCore, quietly =T)
library(umap, quietly =T)
library(Rtsne, quietly =T)

library(ggforce, quietly =T)
library(plyr, quietly =T)
library(dplyr, quietly =T)
library(ggplot2, quietly =T)
library(cowplot, quietly =T)
library(reshape2, quietly =T)

# define variables and paths ----------------------------------------------
# TODO move it finally to the 'geomx_runner'
proj_dir <- '~/Documents/phd/st'
data_dir <- '~/Documents/phd/st/data/geomx/geomx_batch2_1124/'

dcc_path <- dir(file.path(data_dir, "dcc"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
pkc_path <- file.path(data_dir, 'metadata', 'Hs_R_NGS_WTA_v1.0.pkc')
anno_path <- file.path(data_dir, 'metadata', 'dcc_metadata_all.xlsx')

output_dir <- file.path(proj_dir, 'geomx-processing', 'results', 'batch2')

# load utils functions ----------------------------------------------------

source(file.path(proj_dir, 'st-processing', 'src', 'visium_utils.R')) #TODO add needed functions to geomx_utils
source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))

# define intermediate output paths ----------------------------------------

dir.create(output_dir, recursive = T, showWarnings = F)

geomx_qc_path <<- file.path(output_dir, 'geomx_qc_neggeo_ntc.RDS')
geomx_norm_path <<- file.path(output_dir, 'geomx_qc_norm.RDS')

# start the pipeline ------------------------------------------------------

print('#############')
print('GeoMx pipeline starting :O')
print('#############')

# conditionally run preprocessing -----------------------------------------

run_unless_exists('Preprocessing', geomx_qc_path, 
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_qc.R'))

# conditionally run normalisation -----------------------------------------

run_unless_exists('Preprocessing', geomx_norm_path, 
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_normalisation.R'))


