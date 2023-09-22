library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(plyr)
library(dplyr)

# define variables --------------------------------------------------------

data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'

dcc_path <- dir(file.path(data_dir, "dcc"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)

pkc_path <- file.path(data_dir, 'metadata/Hs_R_NGS_WTA_v1.0.pkc')
anno_path <- file.path(data_dir, 'metadata/dcc_metadata_all.xlsx')
# anno_path <-
#   dir(file.path(data_dir, "annotation"), pattern = ".xlsx$",
#       full.names = TRUE, recursive = TRUE)

output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact'
dir.create(output_dir, showWarnings = T, recursive = T)

# load geomx dataset ------------------------------------------------------

geomx_obj <- readNanoStringGeoMxSet(dccFiles = dcc_path, 
                                    pkcFiles = pkc_path, # this goes into fData() - features (probes) annotation
                                    phenoDataFile = anno_path, # this goes into pData() - protocol (samples) annotation
                                    phenoDataSheet = "Sheet1",
                                    phenoDataDccColName = "Sample_ID",
                                    protocolDataColNames = c("Aoi", "Roi"), #TODO adjust
                                    experimentDataColNames = c("Panel")) #TODO adjust


View(assayData(geomx_obj)$exprs)
dim(assayData(geomx_obj)$exprs)

View(pData(geomx_obj))
View(pData(protocolData(geomx_obj)))
View(fData(geomx_obj))
featureType(geomx_obj)
annotation(geomx_obj)

View(summary(geomx_obj, MARGIN = 1)) # for probes
View(summary(geomx_obj, MARGIN = 2)) # for segments (rois)

