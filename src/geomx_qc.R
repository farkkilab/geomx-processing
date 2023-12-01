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
library(clusterProfiler)
library(msigdbr)
library(progeny)
library(reshape2)
library(biomaRt)
library(GSVA)
library(ggpubr)

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
dir.create(file.path(output_dir, 'dge'), showWarnings = T, recursive = T)

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')

# signature files
sig_t_exhaustion <- '/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/t_cell_exhaustion.csv'
sig_macro <- '/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/macrophages.csv'
sig_mhc <- '/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/MHC.csv'

sig_all <- list(sig_t_exhaustion, sig_macro, sig_mhc)

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
                                               paste0(method, '_', norm, '_', color_var, '.pdf')))
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
##########################################################################
##########################################################################
##########################################################################
# TODO PREPROCESSING ENDS HERE UPPER PART OF A SCRIPT SEPARATED FROM THE ANALYSIS BELOW
# TODO SAVE GEOMX_OBJ AS EXTERNAJ OBJ AND LOAD FOR THE NEXT SCRIPTS
# make DGE between selected ROI groups ------------------------------------

# within slide analysis - with random slope in LLM
# comparison between ++ (posCD8_posIBA1) and other groups

# convert test variables to factors
pData(geomx_obj)[["Annotation_cell_factor"]] <- factor(pData(geomx_obj)[["Annotation_cell"]])
pData(geomx_obj)[["Sample_factor"]] <- factor(pData(geomx_obj)[["Sample"]])
pData(geomx_obj)[["NACT_status_factor"]] <- factor(pData(geomx_obj)[["NACT status"]])
pData(geomx_obj)[["PFS_factor"]] <- factor(pData(geomx_obj)[["PFS"]])

# convert normalized counts to log scale
assayDataElement(object = geomx_obj, elt = "log_q3_norm") <-
  assayDataApply(geomx_obj, 2, FUN = log, base = 2, elt = "q3_norm")

assayDataElement(object = geomx_obj, elt = "log_quant_norm") <-
  assayDataApply(geomx_obj, 2, FUN = log, base = 2, elt = "quant_norm")

# TODO same for quant norm
# run LMM:
# formula follows conventions defined by the lme4 package
results <- c()
for(segment in c("tumor", "stroma")){
  # careful! this is from all table with umap made for all ROIs - should be change to avoid confusion
  geomx_segment <- geomx_obj[, geomx_obj@phenoData@data$Segment == segment]
  for(status in c("pre", "post")) {
    ind <- geomx_segment@phenoData@data$`NACT status` == status
    
    mixedOutmc <-
      mixedModelDE(geomx_segment[, ind],
                   elt = "log_q3_norm",
                   modelFormula = ~ Annotation_cell_factor + (1 + Annotation_cell_factor | Sample_factor), # random slope
                   groupVar = "Annotation_cell_factor",
                   nCores = (parallel::detectCores() - 1),
                   multiCore = FALSE)
    
    
    # format results as data.frame
    r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests
    
    # use lapply in case you have multiple levels of your test factor to
    # correctly associate gene name with it's row in the results table
    r_test$Gene <- 
      unlist(lapply(colnames(mixedOutmc),
                    rep, nrow(mixedOutmc["lsmeans", ][[1]])))
    r_test$Subset <- status
    r_test$Segment <- segment
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
    r_test <- r_test[, c("Gene", "Subset", "Segment",  "Contrast", "Estimate", 
                         "Pr(>|t|)", "FDR")]
    results <- rbind(results, r_test)
  }
}

fwrite(results, file.path(output_dir, 'dge/dge_annotation_cell_pre_post_separately.csv'))

# without differentiation to pre and post

results2 <- c()
for(segment in c("tumor", "stroma")){
  # careful! this is from all table with umap made for all ROIs - should be change to avoid confusion
  geomx_segment <- geomx_obj[, geomx_obj@phenoData@data$Segment == segment]
  
  mixedOutmc <-
    mixedModelDE(geomx_segment,
                 elt = "log_q3_norm",
                 modelFormula = ~ Annotation_cell_factor + (1 + Annotation_cell_factor | Sample_factor), # random slope
                 groupVar = "Annotation_cell_factor",
                 nCores = (parallel::detectCores() - 1),
                 multiCore = FALSE)
  
  
  # format results as data.frame
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  
  # use lapply in case you have multiple levels of your test factor to
  # correctly associate gene name with it's row in the results table
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Segment <- segment
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Segment",  "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results2 <- rbind(results2, r_test)
}

fwrite(results2, file.path(output_dir, 'dge/dge_annotation_cell_all.csv'))


#####################################
results_signif <- results[results$FDR <= 0.05, ]
results2_signif <- results2[results2$FDR <= 0.05, ] 

fwrite(results_signif, file.path(output_dir, 'dge/dge_annotation_cell_pre_post_separately_signif.csv'))
fwrite(results2_signif, file.path(output_dir, 'dge/dge_annotation_cell_all_signif.csv'))

# TODO redo for 1group vs 3groups all together

######################################
# BETWEEN SLIDES COMPARISON

# run LMM without random slope:
# formula follows conventions defined by the lme4 package
# results_pre_post_per_cell <- c()
# for(segment in c("tumor", "stroma")){
#   # careful! this is from all table with umap made for all ROIs - should be change to avoid confusion
#   geomx_segment <- geomx_obj[, geomx_obj@phenoData@data$Segment == segment]
#   for(anno_cell in unique(geomx_obj@phenoData@data$Annotation_cell)) {
#     ind <- geomx_segment@phenoData@data$Annotation_cell == anno_cell
#     
#     # TODO trycatch if too litle nr of ROIs - return empty frame
#     mixedOutmc <-
#       mixedModelDE(geomx_segment[, ind],
#                    elt = "log_q3_norm",
#                    modelFormula = ~ NACT_status_factor + (1 | Sample_factor), # random slope
#                    groupVar = "NACT_status_factor",
#                    nCores = (parallel::detectCores() - 1),
#                    multiCore = FALSE)
#     
#     
#     # format results as data.frame
#     r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
#     tests <- rownames(r_test)
#     r_test <- as.data.frame(r_test)
#     r_test$Contrast <- tests
#     
#     # use lapply in case you have multiple levels of your test factor to
#     # correctly associate gene name with it's row in the results table
#     r_test$Gene <- 
#       unlist(lapply(colnames(mixedOutmc),
#                     rep, nrow(mixedOutmc["lsmeans", ][[1]])))
#     r_test$Annotation_cell <- anno_cell
#     r_test$Segment <- segment
#     r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
#     r_test <- r_test[, c("Gene", "Annotation_cell", "Segment",  "Contrast", "Estimate", 
#                          "Pr(>|t|)", "FDR")]
#     results_pre_post_per_cell <- rbind(results_pre_post_per_cell, r_test)
#   }
# }
# 
# fwrite(results_pre_post_per_cell, file.path(output_dir, 'dge/dge_pre_post_per_cell_separately.csv'))

#########################
##########################
# all cell anno mixed together doesnt give any meaningful results!!!

results_pre_post_doublepos <- c()
for(segment in c("tumor", "stroma")){
  # careful! this is from all table with umap made for all ROIs - should be change to avoid confusion
  geomx_segment <- geomx_obj[, geomx_obj@phenoData@data$Segment == segment]
  geomx_segment_doublepos <- geomx_segment[, geomx_segment@phenoData@data$Annotation_cell == "posCD8_posIBA1"]
  
  mixedOutmc <-
    mixedModelDE(geomx_segment_doublepos,
                 elt = "log_q3_norm",
                 modelFormula = ~ NACT_status_factor + (1 | Sample_factor), # random slope
                 groupVar = "NACT_status_factor",
                 nCores = (parallel::detectCores() - 1),
                 multiCore = FALSE)
  
  
  # format results as data.frame
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  
  # use lapply in case you have multiple levels of your test factor to
  # correctly associate gene name with it's row in the results table
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Segment <- segment
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Segment",  "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results_pre_post_doublepos <- rbind(results_pre_post_doublepos, r_test)
}

fwrite(results_pre_post_doublepos, file.path(output_dir, 'dge/dge_pre_post_doublepos.csv'))
results_pre_post_doublepos_signif <- results_pre_post_doublepos[results_pre_post_doublepos$FDR <= 0.05, ]
fwrite(results_pre_post_doublepos_signif, file.path(output_dir, 'dge/dge_pre_post_doublepos_signif.csv'))

################################
# for alltogether post samples long vs short pfs
geomx_post <- geomx_obj[, geomx_obj@phenoData@data$`NACT status` == 'post']

results_pfs_doublepos <- c()
for(segment in c("tumor", "stroma")){
  # careful! this is from all table with umap made for all ROIs - should be change to avoid confusion
  geomx_segment <- geomx_post[, geomx_post@phenoData@data$Segment == segment]
  geomx_segment_doublepos <- geomx_segment[, geomx_segment@phenoData@data$Annotation_cell == "posCD8_posIBA1"]
  
  mixedOutmc <-
    mixedModelDE(geomx_segment_doublepos,
                 elt = "log_q3_norm",
                 modelFormula = ~ PFS_factor + (1 | Sample_factor), # random slope
                 groupVar = "PFS_factor",
                 nCores = (parallel::detectCores() - 1),
                 multiCore = FALSE)
  
  
  # format results as data.frame
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  
  # use lapply in case you have multiple levels of your test factor to
  # correctly associate gene name with it's row in the results table
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Segment <- segment
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Segment",  "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results_pfs_doublepos <- rbind(results_pre_post_doublepos, r_test)
}

fwrite(results_pfs_doublepos, file.path(output_dir, 'dge/dge_pfs_doublepos.csv'))
results_pfs_doublepos_signif <- results_pfs_doublepos[results_pfs_doublepos$FDR <= 0.05, ]
fwrite(results_pfs_doublepos_signif, file.path(output_dir, 'dge/dge_pfs_doublepos_signif.csv'))


#############################################################
#############################################################
# separately per each patient pre and post, all 
paired_patients <- c("S015", "S027", "S032", "S084", "S139")

results_patient_pairs <- c()
for(segment in c("tumor", "stroma")){
  # careful! this is from all table with umap made for all ROIs - should be change to avoid confusion
  geomx_segment <- geomx_obj[, geomx_obj@phenoData@data$Segment == segment]
  
  for(patient in paired_patients){
    geomx_patient <- geomx_segment[, geomx_segment@phenoData@data$Patient == patient]
    
    mixedOutmc <-
      mixedModelDE(geomx_patient,
                   elt = "log_q3_norm",
                   modelFormula = ~ NACT_status_factor + (1 | Sample_factor), # random slope
                   groupVar = "NACT_status_factor",
                   nCores = (parallel::detectCores() - 1),
                   multiCore = FALSE)
    
    
    # format results as data.frame
    r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests
    
    # use lapply in case you have multiple levels of your test factor to
    # correctly associate gene name with it's row in the results table
    r_test$Gene <- 
      unlist(lapply(colnames(mixedOutmc),
                    rep, nrow(mixedOutmc["lsmeans", ][[1]])))
    r_test$Segment <- segment
    r_test$Patient <- patient
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
    r_test <- r_test[, c("Gene", "Segment", 'Patient',  "Contrast", "Estimate", 
                         "Pr(>|t|)", "FDR")]
    results_patient_pairs <- rbind(results_patient_pairs, r_test)
  }
}

fwrite(results_patient_pairs, file.path(output_dir, 'dge/dge_patient_pairs.csv'))
results_patient_pairs_signif <- results_patient_pairs[results_patient_pairs$FDR <= 0.05, ]
fwrite(results_patient_pairs_signif, file.path(output_dir, 'dge/dge_patient_pairs_signif.csv'))

############
# separately for patients, only doublepos
results_patient_pairs_doublepos <- c()
for(segment in c("tumor", "stroma")){
  # careful! this is from all table with umap made for all ROIs - should be change to avoid confusion
  geomx_segment <- geomx_obj[, geomx_obj@phenoData@data$Segment == segment]
  geomx_segment_doublepos <- geomx_segment[, geomx_segment@phenoData@data$Annotation_cell == "posCD8_posIBA1"]
  
  for(patient in paired_patients){
    geomx_patient <- geomx_segment_doublepos[, geomx_segment_doublepos@phenoData@data$Patient == patient]
    
    mixedOutmc <-
      mixedModelDE(geomx_patient,
                   elt = "log_q3_norm",
                   modelFormula = ~ NACT_status_factor + (1 | Sample_factor), # random slope
                   groupVar = "NACT_status_factor",
                   nCores = (parallel::detectCores() - 1),
                   multiCore = FALSE)
    
    
    # format results as data.frame
    r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests
    
    # use lapply in case you have multiple levels of your test factor to
    # correctly associate gene name with it's row in the results table
    r_test$Gene <- 
      unlist(lapply(colnames(mixedOutmc),
                    rep, nrow(mixedOutmc["lsmeans", ][[1]])))
    r_test$Segment <- segment
    r_test$Patient <- patient
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
    r_test <- r_test[, c("Gene", "Segment", 'Patient',  "Contrast", "Estimate", 
                         "Pr(>|t|)", "FDR")]
    results_patient_pairs_doublepos <- rbind(results_patient_pairs_doublepos, r_test)
  }
}

fwrite(results_patient_pairs_doublepos, file.path(output_dir, 'dge/dge_patient_pairs_doublepos.csv'))
results_patient_pairs_doublepos_signif <- results_patient_pairs_doublepos[results_patient_pairs_doublepos$FDR <= 0.05, ]
fwrite(results_patient_pairs_doublepos_signif, file.path(output_dir, 'dge/dge_patient_pairs_doublepos_signif.csv'))


###########################################
# check DEG results
list.files(file.path(output_dir, 'dge/old'))

dge_anno <- fread(file.path(output_dir, 'dge/old/dge_annotation_cell_all_signif.csv'))
dge_anno_prepost <- fread(file.path(output_dir, 'dge/old/dge_annotation_cell_pre_post_separately_signif.csv'))

# length = 0
#dge_pfs <- fread(file.path(output_dir, 'dge/old/dge_pfs_all_signif.csv'))
#dge_pfs_doublepos <- fread(file.path(output_dir, 'dge/old/dge_pfs_doublepos_signif.csv'))

#dge_prepost <- fread(file.path(output_dir, 'dge/old/dge_pre_post_all_signif.csv'))
#dge_prepost_doublepos <- fread(file.path(output_dir, 'dge/old/dge_pre_post_doublepos_signif.csv'))

dge_patient_pairs <- fread(file.path(output_dir, 'dge/old/dge_patient_pairs_signif.csv'))
dge_patient_pairs_doublepos <- fread(file.path(output_dir, 'dge/old/dge_patient_pairs_doublepos_signif.csv'))

dge_list <- list(dge_anno, dge_anno_prepost, dge_patient_pairs, dge_patient_pairs_doublepos)

##################################################
###################################################
# prepare volcano plots

library(ggrepel) 

####
# anno
res_anno <- fread(file.path(output_dir, 'dge/old/dge_annotation_cell_all.csv'))

for(contrast in unique(res_anno$Contrast)){
  
  cont_pos <- unlist(strsplit(contrast, ' - '))[1]
  cont_neg <- unlist(strsplit(contrast, ' - '))[2]
  
  plot_volcano_deg(filter(res_anno, Contrast == contrast),
                   'anno', 10, cont_pos, cont_neg)
}

#####
# anno pre-post
res_anno_prepost <- fread(file.path(output_dir, 'dge/old/dge_annotation_cell_pre_post_separately.csv'))

for(subset in unique(res_anno_prepost$Subset)){
  for(contrast in unique(res_anno_prepost$Contrast)){
    
    cont_pos <- unlist(strsplit(contrast, ' - '))[1]
    cont_neg <- unlist(strsplit(contrast, ' - '))[2]
    
    plot_volcano_deg(filter(res_anno_prepost, Contrast == contrast & Subset == subset),
                     paste0('anno_', subset), 10, cont_pos, cont_neg)
  }
}

######
# patient_pairs all
res_patient_pairs <- fread(file.path(output_dir, 'dge/old/dge_patient_pairs.csv'))

for(patient in unique(res_patient_pairs$Patient)){
  for(contrast in unique(res_patient_pairs$Contrast)){
    
    cont_pos <- unlist(strsplit(contrast, ' - '))[1]
    cont_neg <- unlist(strsplit(contrast, ' - '))[2]
    
    plot_volcano_deg(filter(res_patient_pairs, Contrast == contrast & Patient == patient),
                     paste0('patient_pairs_all_', patient), 10, cont_pos, cont_neg)
  }
}

######
# patient_pairs doublepos
res_patient_pairs_doublepos <- fread(file.path(output_dir, 'dge/old/dge_patient_pairs_doublepos.csv'))

for(patient in unique(res_patient_pairs_doublepos$Patient)){
  for(contrast in unique(res_patient_pairs_doublepos$Contrast)){
    
    cont_pos <- unlist(strsplit(contrast, ' - '))[1]
    cont_neg <- unlist(strsplit(contrast, ' - '))[2]
    
    plot_volcano_deg(filter(res_patient_pairs_doublepos, Contrast == contrast & Patient == patient),
                     paste0('patient_pairs_doublepos_', patient), 10, cont_pos, cont_neg)
  }
}

#######################################
#######################################

# ORA on all Hallmarks + CP -----------------------------------------------



# get bcg genes - all genes in dataset
bcg_genes <- rownames(geomx_obj)

# prepare mdigdb
msigdb_df <- msigdbr(species = "Homo sapiens")
msigdb_df <- filter(msigdb_df, gs_cat %in% c("H", "C2") & gs_subcat != "CGP")

################################
# ora for dge_anno
ora_anno <- data.frame()
for(segment in unique(dge_anno$Segment)){
  for(contrast in unique(dge_anno$Contrast)){
    dge_pos <- filter(dge_anno, Segment == segment & Contrast == contrast & Estimate > 0)
    dge_neg <- filter(dge_anno, Segment == segment & Contrast == contrast & Estimate < 0)
    
    ora_pos <- calculate_ora(dge_pos$Gene, bcg_genes, msigdb_df, padj = 0.05)
    ora_neg <- calculate_ora(dge_neg$Gene, bcg_genes, msigdb_df, padj = 0.05)
    
    if(nrow(ora_pos > 0)){
      ora_pos$direction <- 'up'
    }
    
    if(nrow(ora_neg > 0)){
      ora_neg$direction <- 'down'
    }
    
    ora <- rbind(ora_pos, ora_neg)
    
    if(nrow(ora) > 0){
      ora$contrast <- contrast
      ora$segment <- segment
      ora$GeneRatio_perc <- as.numeric(gsub("\\/[0-9]*", "", ora$GeneRatio))/
        as.numeric(gsub("[0-9]*\\/", "", ora$GeneRatio))
      
      ora_anno <- rbind(ora_anno, ora)
    }

  }
}

fwrite(ora_anno, file.path(output_dir, 'dge/ora/ora_anno.csv'))
#######################################################
# ora for dge_anno_prepost

ora_anno_prepost <- data.frame()
for(segment in unique(dge_anno_prepost$Segment)){
  for(subset in unique(dge_anno_prepost$Subset)){
    for(contrast in unique(dge_anno_prepost$Contrast)){
      dge_pos <- filter(dge_anno_prepost, Segment == segment & Subset == subset & 
                          Contrast == contrast & Estimate > 0)
      dge_neg <- filter(dge_anno_prepost, Segment == segment & Subset == subset &
                          Contrast == contrast & Estimate < 0)
      
      ora_pos <- calculate_ora(dge_pos$Gene, bcg_genes, msigdb_df, padj = 0.05)
      ora_neg <- calculate_ora(dge_neg$Gene, bcg_genes, msigdb_df, padj = 0.05)
      
      if(nrow(ora_pos > 0)){
        ora_pos$direction <- 'up'
      }
      
      if(nrow(ora_neg > 0)){
        ora_neg$direction <- 'down'
      }
      
      ora <- rbind(ora_pos, ora_neg)
      
      if(nrow(ora) > 0){
        ora$contrast <- contrast
        ora$segment <- segment
        ora$subset <- subset
        ora$GeneRatio_perc <- as.numeric(gsub("\\/[0-9]*", "", ora$GeneRatio))/
          as.numeric(gsub("[0-9]*\\/", "", ora$GeneRatio))
        
        ora_anno_prepost <- rbind(ora_anno_prepost, ora)
      }
    }
  }
}
fwrite(ora_anno_prepost, file.path(output_dir, 'dge/ora/ora_anno_prepost.csv'))
#######################################################
# ora for dge_patient_pairs

ora_patient_pairs <- data.frame()
for(segment in unique(dge_patient_pairs$Segment)){
  for(patient in unique(dge_patient_pairs$Patient)){
    for(contrast in unique(dge_patient_pairs$Contrast)){
      dge_pos <- filter(dge_patient_pairs, Segment == segment & Patient == patient & 
                          Contrast == contrast & Estimate > 0)
      dge_neg <- filter(dge_patient_pairs, Segment == segment & Patient == patient &
                          Contrast == contrast & Estimate < 0)
      
      ora_pos <- calculate_ora(dge_pos$Gene, bcg_genes, msigdb_df, padj = 0.05)
      ora_neg <- calculate_ora(dge_neg$Gene, bcg_genes, msigdb_df, padj = 0.05)
      
      if(nrow(ora_pos > 0)){
        ora_pos$direction <- 'up'
      }
      
      if(nrow(ora_neg > 0)){
        ora_neg$direction <- 'down'
      }
      
      ora <- rbind(ora_pos, ora_neg)
      
      if(nrow(ora) > 0){
        ora$contrast <- contrast
        ora$segment <- segment
        ora$patient <- patient
        ora$GeneRatio_perc <- as.numeric(gsub("\\/[0-9]*", "", ora$GeneRatio))/
          as.numeric(gsub("[0-9]*\\/", "", ora$GeneRatio))
        
        ora_patient_pairs <- rbind(ora_patient_pairs, ora)
      }
    }
  }
}

fwrite(ora_patient_pairs, file.path(output_dir, 'dge/ora/ora_patient_pairs.csv'))
#######################################################
# ora for dge_patient_pairs_doublepos

ora_patient_pairs_doublepos <- data.frame()
for(segment in unique(dge_patient_pairs_doublepos$Segment)){
  for(patient in unique(dge_patient_pairs_doublepos$Patient)){
    for(contrast in unique(dge_patient_pairs_doublepos$Contrast)){
      dge_pos <- filter(dge_patient_pairs_doublepos, Segment == segment & Patient == patient & 
                          Contrast == contrast & Estimate > 0)
      dge_neg <- filter(dge_patient_pairs_doublepos, Segment == segment & Patient == patient &
                          Contrast == contrast & Estimate < 0)
      
      ora_pos <- calculate_ora(dge_pos$Gene, bcg_genes, msigdb_df, padj = 0.05)
      ora_neg <- calculate_ora(dge_neg$Gene, bcg_genes, msigdb_df, padj = 0.05)
      
      if(nrow(ora_pos > 0)){
        ora_pos$direction <- 'up'
      }
      
      if(nrow(ora_neg > 0)){
        ora_neg$direction <- 'down'
      }
      
      ora <- rbind(ora_pos, ora_neg)
      
      if(nrow(ora) > 0){
        ora$contrast <- contrast
        ora$segment <- segment
        ora$patient <- patient
        ora$GeneRatio_perc <- as.numeric(gsub("\\/[0-9]*", "", ora$GeneRatio))/
          as.numeric(gsub("[0-9]*\\/", "", ora$GeneRatio))
        
        ora_patient_pairs_doublepos <- rbind(ora_patient_pairs_doublepos, ora)
      }
    }
  }
}

fwrite(ora_patient_pairs_doublepos, file.path(output_dir, 'dge/ora/ora_patient_pairs_doublepos.csv'))


# PROGENy scores ----------------------------------------------------------


dim(geomx_obj)
colnames(geomx_obj)[1:10]
rownames(geomx_obj)[1:10]

prog_noperm <- progeny(
  geomx_obj@assayData$q3_norm,
  scale = TRUE,
  organism = "Human",
  top = 100,
  perm = 1
  )


prog_perm <- progeny(
  geomx_obj@assayData$q3_norm,
  organism = "Human",
  top = 100,
  perm = 10,
  z_scores = FALSE,
  get_nulldist = FALSE
)

rownames(prog_perm) <- gsub('\\.', '\\-', rownames(prog_perm))
rownames(prog_perm) <- gsub('\\-dcc', '\\.dcc', rownames(prog_perm))


#####
#TODO do it in loop
prog_df <- prog_perm
prog_name <- 'perm'

# adjust df
prog_long <- melt(prog_df)
colnames(prog_long) <- c('dcc_filename', 'progeny_path', 'progeny_score')
prog_long <- left_join(prog_long, pData(geomx_obj)[c('dcc_filename', 'Segment', 'Annotation_cell', 'NACT status', 'PFS', 'Patient')])

fwrite(prog_long, file.path(output_dir, 'progeny', paste0('progeny_', prog_name, '.csv')))

##############
# make boxplots

# per Anno cell type
prog_boxpl <- ggplot(data = prog_long, aes(x = progeny_path, y = progeny_score, color = Annotation_cell)) +
  geom_boxplot() +
  facet_wrap(~Segment, scales = "fixed", dir="v") +
  geom_pwc(method = "t_test", label = "p.signif", hide.ns = TRUE) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ggtitle(paste('progeny', prog_name, 'scores'))

plot(prog_boxpl)
ggsave(file.path(output_dir, 'progeny', paste0('box_progeny_anno_', prog_name, '.png')), 
       height = 2000, width = 3000, unit = 'px')

################################
# per NACT status

prog_boxpl_nact_all <- ggplot(data = prog_long, aes(x = progeny_path, y = progeny_score, color = `NACT status`)) +
  geom_boxplot() +
  facet_wrap(~Segment, scales = "fixed", dir="v") +
  geom_pwc(method = "t_test", label = "p.signif", hide.ns = TRUE) +
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  ggtitle(paste('progeny', prog_name, 'scores'))

plot(prog_boxpl_nact_all)
ggsave(file.path(output_dir, 'progeny', paste0('box_progeny_prepost_all_', prog_name, '.png')), 
       height = 2000, width = 3000, unit = 'px')

prog_boxpl_nact_peranno <- ggplot(data = prog_long, aes(x = progeny_path, y = progeny_score, color = `NACT status`)) +
  geom_boxplot() +
  facet_wrap(Segment~Annotation_cell, scales = "fixed", ncol=4, nrow=2) +
  geom_pwc(method = "t_test", label = "p.signif", hide.ns = TRUE) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ggtitle(paste('progeny', prog_name, 'scores'))

plot(prog_boxpl_nact_peranno)
ggsave(file.path(output_dir, 'progeny', paste0('box_progeny_prepost_anno_', prog_name, '.png')), 
       height = 2000, width = 4000, unit = 'px')


################################
# per PFS in post samples
prog_long_post <- filter(prog_long, `NACT status` == 'post')

prog_boxpl_pfs_all <- ggplot(data = prog_long_post, aes(x = progeny_path, y = progeny_score, color = PFS)) +
  geom_boxplot() +
  facet_wrap(~Segment, scales = "fixed", dir="v") +
  geom_pwc(method = "t_test", label = "p.signif", hide.ns = TRUE) +
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  ggtitle(paste('progeny', prog_name, 'scores'))

plot(prog_boxpl_pfs_all)
ggsave(file.path(output_dir, 'progeny', paste0('box_progeny_pfs_all_', prog_name, '.png')), height = 2000, width = 3000, unit = 'px')

prog_boxpl_pfs_peranno <- ggplot(data = prog_long_post, aes(x = progeny_path, y = progeny_score, color = PFS)) +
  geom_boxplot() +
  facet_wrap(Segment~Annotation_cell, scales = "fixed", ncol=4, nrow=2) +
  geom_pwc(method = "t_test", label = "p.signif", hide.ns = TRUE) +
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  ggtitle(paste('progeny', prog_name, 'scores'))

plot(prog_boxpl_pfs_peranno)
ggsave(file.path(output_dir, 'progeny', paste0('box_progeny_pfs_anno_', prog_name, '.png')), 
       height = 2000, width = 4000, unit = 'px')




# UCell scores for selected pathways --------------------------------------


# changing signatures to HGNC symbols

# sig_list <- lapply(sig_all, function(x){
#   sig <- fread(x)
#   sig_list <- as.list(sig)
#   return(sig_list)
# })
# 
# sig_list <- unlist(sig_list, recursive = F)
# 
# sig_list_names <- lapply(sig_list, function(x){
#   print('XXXXX')
#   print(x)
#   x <- x[!is.na(x)]
#   print(x)
#   x <- unlist(gene_2names(x, conv='entrez'))
#   print(x)
#   })
# 
# library(qpcR)
# sig_df <- do.call(qpcR:::cbind.na, sig_list_names)
# fwrite(sig_df, '/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/texh_macro_mhc.csv')


# GSVA and ssGSEA ---------------------------------------------------------
library(GSVA)

# selected pathways
#sig_texh_macro_mhc_list <- as.list(fread('/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/texh_macro_mhc.csv'))
sig_texh_macro_mhc_list <- as.list(fread('/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/texh_macro_mhc_ifng_forpaper.csv'))
sig_caf_revised_list <- as.list(fread('/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/stromal_cell_subtype_signatures_symbols_ensembl_ids_revised.csv'))

# pathways from Glasgow
sig_pycr1 <- fread('~/Downloads/Gs_markers_02_04_10_23.csv')
sig_pycr1_05 <- filter(sig_pycr1, avg_log2FC >= 0.5)
pycr_list <- list('pycr1' = sig_pycr1$gene, 'pycr1_05' = sig_pycr1_05$gene)

pycr_sig <- data.frame('pycr1' = sig_pycr1$gene, 'pycr1_05' = sig_pycr1_05$gene)

# all Hallmark + CP from msigDB
msigdb_df <- msigdbr(species = "Homo sapiens")
msigdb_df <- filter(msigdb_df, gs_cat == 'H' | gs_subcat %in% c('CP:BIOCARTA', 'CP:KEGG', 'CP:REACTOME'))

hal_cp_list <- lapply(unique(msigdb_df$gs_name), function(x){
  gs <- filter(msigdb_df, gs_name == x)
  gs_genes <- unique(gs$human_gene_symbol)
})

names(hal_cp_list) <- unique(msigdb_df$gs_name)

#####

expr_mtx <- assayDataElement(geomx_obj, elt = "q3_norm")
View(expr_mtx[1:10, 1:10])

# do gsva on pathway lists
gsva_hal_cp <- gsva(expr_mtx, hal_cp_list, method = 'gsva', kcdf="Poisson", min.sz = 5)
gsva_texh_macro_mhc <- gsva(expr_mtx, sig_texh_macro_mhc_list, method = 'gsva', kcdf="Poisson", min.sz = 5)
gsva_caf <- gsva(expr_mtx, sig_caf_revised_list, method = 'gsva', kcdf="Poisson", min.sz = 5)
gsva_pycr1 <- gsva(expr_mtx, pycr_list, method = 'gsva', kcdf="Poisson", min.sz = 5)

# do ssgsea on selected lists
# ssgsea_hal_cp <- gsva(expr_mtx, hal_cp_list, method = 'ssgsea', kcdf="Poisson", min.sz = 5)
# ssgsea_texh_macro_mhc <- gsva(expr_mtx, sig_texh_macro_mhc_list, method = 'ssgsea', kcdf="Poisson", min.sz = 5)
# ssgsea_caf <- gsva(expr_mtx, sig_caf_revised_list, method = 'ssgsea', kcdf="Poisson", min.sz = 5)

######

# adjust df

# TODO do the same for hallmark+cp

gsva_texh_macro_mhc_long <- melt(gsva_texh_macro_mhc)
colnames(gsva_texh_macro_mhc_long) <- c('pathway','dcc_filename', 'gsva_score')
gsva_texh_macro_mhc_long <- left_join(gsva_texh_macro_mhc_long, pData(geomx_obj)[c('dcc_filename', 'Segment', 'Annotation_cell', 'NACT status', 'PFS', 'Patient')])

fwrite(gsva_texh_macro_mhc_long, file.path(output_dir, 'gsva', paste0('gsva_texh_macro_mhc_forpaper.csv')))

gsva_caf_long <- melt(gsva_caf)
colnames(gsva_caf_long) <- c('pathway','dcc_filename', 'gsva_score')
gsva_caf_long <- left_join(gsva_caf_long, pData(geomx_obj)[c('dcc_filename', 'Segment', 'Annotation_cell', 'NACT status', 'PFS')])

fwrite(gsva_caf_long, file.path(output_dir, 'gsva', paste0('gsva_caf.csv')))

gsva_list <- list(gsva_texh_macro_mhc_long, gsva_caf_long)
names(gsva_list) <- c('texh_macro_mhc', 'caf')

gsva_pycr1_long <- melt(gsva_pycr1)
colnames(gsva_pycr1_long) <- c('pathway','dcc_filename', 'gsva_score')
gsva_pycr1_long <- left_join(gsva_pycr1_long, pData(geomx_obj)[c('dcc_filename', 'Segment', 'Annotation_cell', 'NACT status', 'PFS')])

fwrite(gsva_pycr1_long, file.path(output_dir, 'gsva', paste0('gsva_pycr1.csv')))

# VISUALISATION IN  GSVA_VISUALISATION.R  


# individual genes expression ---------------------------------------------

goi <- c('PDCD1', 'HAVCR2', 'TIGIT', 'CD96', 'NECTIN2', 'CXCR3', 'CXCL9', 'IL2RG', 'IL2RB') 

expr_mtx <- assayDataElement(geomx_obj, elt = "q3_norm")
expr_mtx_goi <- expr_mtx[goi, ]

expr_mtx_goi_long <- melt(expr_mtx_goi)
colnames(expr_mtx_goi_long) <- c('gene','dcc_filename', 'expr')
expr_mtx_goi_long <- left_join(expr_mtx_goi_long, pData(geomx_obj)[c('dcc_filename', 'Segment', 'Annotation_cell', 'NACT status', 'PFS', 'Patient')])
fwrite(expr_mtx_goi_long, file.path(output_dir, 'ind_genes', paste0('ind_genes_exh_lr.csv')))
