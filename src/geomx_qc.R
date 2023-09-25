library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(plyr)
library(dplyr)
library(ggforce)

# define variables --------------------------------------------------------

data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'

dcc_path <- dir(file.path(data_dir, "dcc"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
pkc_path <- file.path(data_dir, 'metadata/Hs_R_NGS_WTA_v1.0.pkc')
anno_path <- file.path(data_dir, 'metadata/dcc_metadata_all.xlsx')


output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact'
dir.create(output_dir, showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'qc'), showWarnings = T, recursive = T)

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')

# load geomx dataset ------------------------------------------------------

geomx_obj <- readNanoStringGeoMxSet(dccFiles = dcc_path, 
                                    pkcFiles = pkc_path, # this goes into fData() - features (probes) annotation
                                    phenoDataFile = anno_path, # this goes into pData() - protocol (samples) annotation
                                    phenoDataSheet = "Sheet1",
                                    phenoDataDccColName = "Sample_ID",
                                    protocolDataColNames = c("Aoi", "Roi"),
                                    experimentDataColNames = c("Panel")) # TODO dunno if this is needed

# explore
View(assayData(geomx_obj)$exprs)
dim(assayData(geomx_obj)$exprs)

View(pData(geomx_obj))
View(pData(protocolData(geomx_obj)))
View(fData(geomx_obj))
featureType(geomx_obj)
annotation(geomx_obj)

View(summary(geomx_obj, MARGIN = 1)) # for probes
View(summary(geomx_obj, MARGIN = 2)) # for segments (rois)


# make overall sankey plot ------------------------------------------------
#TODO to remove
# pkcs <- annotation(geomx_obj)
# modules <- gsub(".pkc", "", pkcs)

count_segments <- filter(geomx_obj@phenoData@data, Annotation_cell != 'NA' & !is.na(Annotation_cell))
variables_to_plot <- c("Slide Name", "NACT status", "Segment", "Annotation_cell")

plot_sankey(count_segments, variables_to_plot, "NACT status", 
            file.path(output_dir, 'qc/sankey_slides.png'))


# set and plot basic qc parameters ----------------------------------------

qc_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10, 1 in log scale)
       maxNTCCount = 9000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 20,         # Minimum # of nuclei estimated (100)
       minArea = 1000)         # Minimum segment area (5000)

geomx_obj <- setSegmentQCFlags(geomx_obj, qcCutoffs = qc_params)
qc_results <- protocolData(geomx_obj)[["QCFlags"]]

qc_summary <- qc_summarize(qc_results)

# plot qc histograms
QC_histogram(sData(geomx_obj), "Trimmed (%)", "Segment", qc_params[["percentTrimmed"]], 
             scale_trans = NULL, file.path(output_dir, 'qc/qc_hist_trim.png'))

QC_histogram(sData(geomx_obj), "Stitched (%)", "Segment", qc_params[["percentStitched"]],
             scale_trans = NULL, file.path(output_dir, 'qc/qc_hist_stich.png'))

QC_histogram(sData(geomx_obj), "Aligned (%)", "Segment", qc_params[["percentAligned"]],
             scale_trans = NULL, file.path(output_dir, 'qc/qc_hist_align.png'))

QC_histogram(sData(geomx_obj), "Saturated (%)", "Segment", qc_params[["percentSaturation"]],
             scale_trans = NULL, file.path(output_dir, 'qc/qc_hist_satur.png'))

QC_histogram(sData(geomx_obj), "Area", "Segment", qc_params[["minArea"]], 
             scale_trans = "log10", file.path(output_dir, 'qc/qc_hist_area.png'))

QC_histogram(sData(geomx_obj), "Nuclei", "Segment", qc_params[["minNuclei"]],
             scale_trans = NULL, file.path(output_dir, 'qc/qc_hist_nuclei.png'))


# negative probes modelling -----------------------------------------------

# TODO to check with Zhihan


# remove flagged segments -------------------------------------------------
#TODO code repetition from function
qc_results$qc_status <- apply(qc_results, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})

table(qc_results$qc_status)

geomx_obj <- geomx_obj[, qc_results$qc_status == "PASS"]


