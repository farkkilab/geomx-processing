library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(plyr)
library(dplyr)
library(ggforce)
library(data.table)
library(umap)
library(cowplot)
library(preprocessCore)
library(Rtsne)

# define variables --------------------------------------------------------

data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'

dcc_path <- dir(file.path(data_dir, "dcc"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
pkc_path <- file.path(data_dir, 'metadata/Hs_R_NGS_WTA_v1.0.pkc')
anno_path <- file.path(data_dir, 'metadata/dcc_metadata_all.xlsx')


output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact'
dir.create(output_dir, showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'qc'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'dcc_post_qc'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'umap_tsne'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'umap_tsne', 'tumor'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'umap_tsne', 'stroma'), showWarnings = T, recursive = T)

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
# View(assayData(geomx_obj)$exprs)
# dim(assayData(geomx_obj)$exprs)
# 
# View(pData(geomx_obj))
# View(pData(protocolData(geomx_obj)))
# View(fData(geomx_obj))
# featureType(geomx_obj)
# annotation(geomx_obj)
# 
# View(summary(geomx_obj, MARGIN = 1)) # for probes
# View(summary(geomx_obj, MARGIN = 2)) # for segments (rois)

dim(geomx_obj)
# !! no reads for DSP-1001660016656-A-H01.dcc

# make overall sankey plot ------------------------------------------------

count_segments <- filter(geomx_obj@phenoData@data, Annotation_cell != 'NA' & !is.na(Annotation_cell))
variables_to_plot <- c("Slide Name", "NACT status", "Segment", "Annotation_cell")

plot_sankey(count_segments, variables_to_plot, "NACT status", 
            file.path(output_dir, 'qc/sankey_slides.png'))


# set and plot basic qc parameters ----------------------------------------
# Shift 0 counts to one (needed for downstream analysis - ?)
#TODO this is needed for Q3 norm (but not 100% sure why)
geomx_obj <- shiftCountsOne(geomx_obj, useDALogic = TRUE)

qc_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10, 1 in log scale)
       maxNTCCount = 9000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 20,         # Minimum # of nuclei estimated (100) #TODO maybe bigger
       minArea = 1000)         # Minimum segment area (5000)

geomx_obj <- setSegmentQCFlags(geomx_obj, qcCutoffs = qc_params)
qc_results_segment <- protocolData(geomx_obj)[["QCFlags"]]

qc_summary <- qc_summarize(qc_results_segment)

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
qc_results_segment$qc_status <- apply(qc_results_segment, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})

table(qc_results_segment$qc_status)
segments_to_rmv <- geomx_obj@phenoData@data[qc_results_segment$qc_status != "PASS", ]

geomx_obj <- geomx_obj[, qc_results_segment$qc_status == "PASS"]

dim(geomx_obj)
# removed segments:
# DSP-1001660016656-A-H01.dcc - no reads at all
# DSP-1001660016656-A-G01.dcc - low saturation, low negatives
# DSP-1001660016658-B-A12.dcc - low stiched, low aligned
# DSP-1001660016658-B-H01.dcc - low stiched

###########################################################################
###########################################################################

# rmv of probes based on geometric mean and grubbs test -------------------
# the geometric mean of that probe’s counts from all segments divided by the geometric mean 
# of all probe counts representing the target from all segments is less than 0.1
# the probe is an outlier according to the Grubb’s test in at least 20% of the segments

geomx_obj <- setBioProbeQCFlags(geomx_obj, 
                                qcCutoffs = list(minProbeRatio = 0.1,
                                                 percentFailGrubbs = 20), 
                                removeLocalOutliers = TRUE)

qc_results_probe <- fData(geomx_obj)[["QCFlags"]]

# summarise probe qc results
qc_probe_df <- data.frame(Passed = sum(rowSums(qc_results_probe[, -1]) == 0),
                          Global = sum(qc_results_probe$GlobalGrubbsOutlier),
                          Local = sum(rowSums(qc_results_probe[, -2:-1]) > 0
                                      & !qc_results_probe$GlobalGrubbsOutlier))
dim(geomx_obj)

# retain only probes that passed qc (globally)
geomx_obj <- 
  subset(geomx_obj, 
         fData(geomx_obj)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(geomx_obj)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)

dim(geomx_obj)
# 1 probe removed
# TODO what about local removal of probes per segment?

# aggregate probes to features --------------------------------------------

# TODO if we want to do background modelling 
# from geoDiff it should be done before aggregating counts

# nr of unique targets
length(unique(featureData(geomx_obj)[["TargetName"]]))

# collapse features to targets
geomx_obj <- aggregateCounts(geomx_obj)

dim(geomx_obj)


# filter based on LOQ per segment and per gene ----------------------------

loq_cutoff <- 2
loq_min <- 2

# choosen based on the data (nothing is filtered out rn)
gene_detect_thr <- 0.1 # segment is removed if <10% of genes > LOQ
# TODO adjustment may be needed
segment_detect_rate_thr <- 0.01 # genes are removed if > LOQ in less than 1% of segments

# Calculate LOQ for each segment
# TODO adjust if > 1 modules
LOQ <- data.frame(row.names = colnames(geomx_obj))
module <- gsub(".pkc", "", annotation(geomx_obj))

LOQ[, module] <-
  pmax(loq_min,
       pData(geomx_obj)[, paste0("NegGeoMean_", module)] * # coming from aggregate_counts
         pData(geomx_obj)[, paste0("NegGeoSD_", module)] ^ loq_cutoff)

pData(geomx_obj)$LOQ <- LOQ

# calculate if expr > LOQ per each gene per segment
LOQ_Mat <- t(esApply(geomx_obj, MARGIN = 1,
                   FUN = function(x) {
                     x > LOQ[, module]
                   }))

LOQ_Mat <- LOQ_Mat[fData(geomx_obj)$TargetName, ] # ensure ordering

# Save detection rate information to pheno data
# how many genes have been detected in each segment  
pData(geomx_obj)$GenesDetected <- colSums(LOQ_Mat, na.rm = TRUE)
pData(geomx_obj)$GeneDetectionRate <- pData(geomx_obj)$GenesDetected / nrow(geomx_obj)

#TODO calculate signal/noise ratio = Count/LOQ per segment (similar to genedetectionrate)

plot_detection_rate(pData(geomx_obj), "NACT status", 
                    file.path(output_dir, 'qc/gene_detect_rate_nact.png'))
plot_detection_rate(pData(geomx_obj), "Segment", 
                    file.path(output_dir, 'qc/gene_detect_rate_segment.png'))
plot_detection_rate(pData(geomx_obj), "Annotation_cell", 
                    file.path(output_dir, 'qc/gene_detect_rate_anno_cell.png'))


# filter out segments with too low gene detection rate
geomx_obj <- geomx_obj[, pData(geomx_obj)$GeneDetectionRate >= gene_detect_thr]

dim(geomx_obj)

# in how many segments the given gene was detected
# save to probe data
LOQ_Mat <- LOQ_Mat[, colnames(geomx_obj)]
fData(geomx_obj)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(geomx_obj)$DetectionRate <- fData(geomx_obj)$DetectedSegments / nrow(pData(geomx_obj))
LOQ_Mat <- LOQ_Mat[fData(geomx_obj)$TargetName, ]

# plot detection rate per gene
plot_gene_detection_rate(fData(geomx_obj), file.path(output_dir, 'qc/gene_detection_rate.png'))


# manually include the negative control probe, for downstream use
# TODO check how many in the sample data. why just 1?
negativeProbefData <- subset(fData(geomx_obj), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)

# filter out genes detected > LOQ in less then thr nr of segments (1% for now)
geomx_obj <- 
  geomx_obj[fData(geomx_obj)$DetectionRate >= segment_detect_rate_thr |
              fData(geomx_obj)$TargetName %in% neg_probes, ]
dim(geomx_obj)

# save geomx dcc files after QC
#TODO sth wrong here
#writeNanoStringGeoMxSet(geomx_obj, dir = file.path(output_dir, 'dcc_post_qc'))

# Q3 normalisation --------------------------------------------------------

plot_q3_stats(geomx_obj, "Annotation_cell", file.path(output_dir, 'qc/q3_stats.png'))


geomx_obj <- normalize(geomx_obj ,
                       norm_method = "quant", 
                       desiredQuantile = .75,
                       toElt = "q3_norm")

# quantile normalisation --------------------------------------------------

norm.quantile = normalize.quantiles(as.matrix(geomx_obj@assayData$exprs))
dimnames(norm.quantile) = dimnames(geomx_obj@assayData$exprs)

# hacking GeoMx class object 
# TODO this is experimental - newassay is not identical and it may cause problems
# if so, store this in another mtx and use when needed
newassay <- new.env(parent=geomx_obj@assayData)
newassay$exprs <- geomx_obj@assayData$exprs
newassay$q3_norm <- geomx_obj@assayData$q3_norm
newassay$quant_norm <- norm.quantile

geomx_obj@assayData <- newassay

# TODO some problems with plotting, dunno for a while
# plot_norm_effect <- function(expr_data, norm_name, output_name){
#   norm_box <- boxplot(expr_data,
#           col = "#9EDAE5", main = norm_name,
#           log='y', names = seq(1:ncol(expr_data)), xlab = "Segment",
#           ylab = norm_name)
#   
#   png(filename=output_name, width=2000, height=1500, units="px")
#   plot(norm_box)
#   dev.off()
# }
# 
# plot_norm_effect(exprs(geomx_obj)[,1:10], 'Raw Counts', file.path(output_dir, 'qc/norm_raw.png'))
# 
# plot_norm_effect(assayDataElement(geomx_obj[,1:10], elt = "q3_norm"),
#                  'Q3 normalised', file.path(output_dir, 'qc/norm_q3.png'))
# 
# plot_norm_effect(assayDataElement(geomx_obj[,1:10], elt = "quant_norm"),
#                  'Quantile normalised', file.path(output_dir, 'qc/norm_quant.png'))

# plot effects of normalisation
boxplot(exprs(geomx_obj)[,1:10],
        col = "#9EDAE5", main = "Raw Counts",
        log='y', names = 1:10, xlab = "Segment",
        ylab = "Counts, Raw")

boxplot(assayDataElement(geomx_obj[,1:10], elt = "q3_norm"),
        col = "#2CA02C", main = "Q3 Norm Counts",
        log='y', names = 1:10, xlab = "Segment",
        ylab = "Counts, Q3 Normalized")

# TODO I don't like sth with this plot, why all outliers are the same in each segment?
boxplot(assayDataElement(geomx_obj[,1:10], elt = "quant_norm"),
        col = "#2CA02C", main = "Quantile Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Quantile Normalized")


# make UMAP and t-SNE -----------------------------------------------------

# divide for tumor and stroma and do dimentionality reduction for all
geomx_obj_tumor <- geomx_obj[, geomx_obj@phenoData@data$Segment == "tumor"]
geomx_obj_stroma <- geomx_obj[, geomx_obj@phenoData@data$Segment == "stroma"]

geomx_list <- list(all = geomx_obj, tumor = geomx_obj_tumor, stroma = geomx_obj_stroma)

geomx_list_dim_red <- lapply(1:length(geomx_list), function(n){
  
  geomx <- geomx_list[[n]]

  # run UMAP and tSNE on Q3 and quantile norm
  for(norm in c('q3_norm', 'quant_norm')){
    # update defaults for umap to contain a stable random_state (seed)
    custom_umap <- umap::umap.defaults
    custom_umap$random_state <- 42
    
    umap_out <-
      umap(t(log2(assayDataElement(geomx , elt = norm))),  
           config = custom_umap)
    
    # save UMAP1 and 2 results to pData
    pData(geomx)[, c(paste0("UMAP1_", norm), paste0("UMAP2_", norm))] <- umap_out$layout[, c(1,2)]
    
    # set the seed for tSNE as well
    set.seed(42) 
    tsne_out <-
      Rtsne(t(log2(assayDataElement(geomx , elt = norm))),
            perplexity = ncol(geomx)*.15)
    
    # save tSNE1 and 2 results to pData
    pData(geomx)[, c(paste0("tSNE1_", norm), paste0("tSNE2_", norm))] <- tsne_out$Y[, c(1,2)]
  }

  # generate umap and tsne plots and color by variables
  for(method in c('UMAP', 'tSNE')){
    for(norm in c('q3', 'quant')){
      for(color_var in c('Annotation_cell', 'Patient', 'NACT status', 'PFS', 'Site', 'Sample')){
        plot_umap_tsne(pData(geomx), method_type = method, 
                       norm_type = norm, color_var = color_var,
                       output_name = file.path(output_dir, 'umap_tsne2', names(geomx_list)[n], 
                                               paste0(method, '_', norm, '_', color_var, '.png')))
      }
    }
  }
  
  return(geomx)
})

# update objects
geomx_obj <- geomx_list_dim_red[[1]]
geomx_obj_tumor <- geomx_list_dim_red[[2]]
geomx_obj_stroma <- geomx_list_dim_red[[3]]

rm(geomx_list)
rm(geomx_list_dim_red)


# make DGE between selected ROI groups ------------------------------------


