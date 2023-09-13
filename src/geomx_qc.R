library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(plyr)
library(dplyr)
a new change

# define variables --------------------------------------------------------

data_dir <- system.file("extdata", "WTA_NGS_Example",
                       package="GeoMxWorkflows")

dcc_path <- dir(file.path(data_dir, "dccs"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)

pkc_path <- unzip(zipfile = dir(file.path(data_dir, "pkcs"), pattern = ".zip$",
                                full.names = TRUE, recursive = TRUE))
anno_path <-
  dir(file.path(data_dir, "annotation"), pattern = ".xlsx$",
      full.names = TRUE, recursive = TRUE)

output_dir <- '/media/iganiemi/T7-iga/st/st-processing/results/geomx/demo'
dir.create(output_dir, showWarnings = T, recursive = T)

# load geomx dataset ------------------------------------------------------

geomx_obj <- readNanoStringGeoMxSet(dccFiles = dcc_path, 
                                    pkcFiles = pkc_path, # this goes into fData() - features (probes) annotation
                                    phenoDataFile = anno_path, # this goes into pData() - protocol (samples) annotation
                                    phenoDataSheet = "Template", #TODO what about this param?
                                    phenoDataDccColName = "Sample_ID",
                                    protocolDataColNames = c("aoi", "roi"), #TODO adjust
                                    experimentDataColNames = c("panel")) #TODO adjust

View(assayData(geomx_obj)$exprs)
dim(assayData(geomx_obj)$exprs)

View(pData(geomx_obj))
View(pData(protocolData(geomx_obj)))
View(fData(geomx_obj))
featureType(geomx_obj)
annotation(geomx_obj)

View(summary(geomx_obj, MARGIN = 1)) # for probes
View(summary(geomx_obj, MARGIN = 2)) # for segments (rois)

