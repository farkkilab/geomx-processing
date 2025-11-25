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
library(biomaRt)

# define variables --------------------------------------------------------

data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'
output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2'

input_rds_path <- file.path(output_dir, 'geomx_qc.RDS')

geomx <- readRDS(input_rds_path)

meta <- geomx@phenoData@data

# rm annotation
# sample, patient scodes into c codes
x <- 'Patient'

meta[[x]] <- gsub('S015', 'C063', meta[[x]])
meta[[x]] <- gsub('S027', 'C033', meta[[x]])
meta[[x]] <- gsub('S032', 'C686', meta[[x]])
meta[[x]] <- gsub('S053', 'C879', meta[[x]])
meta[[x]] <- gsub('S057', 'C538', meta[[x]])
meta[[x]] <- gsub('S065', 'C799', meta[[x]])
meta[[x]] <- gsub('S072', 'C917', meta[[x]])
meta[[x]] <- gsub('S073', 'C790', meta[[x]])
meta[[x]] <- gsub('S076', 'C423', meta[[x]])
meta[[x]] <- gsub('S084', 'C129', meta[[x]])
meta[[x]] <- gsub('S139', 'C535', meta[[x]])

unique(meta$Patient)

geomx@phenoData@data <- meta

saveRDS(geomx, 'geomx_qc_for_synapse.RDS')
