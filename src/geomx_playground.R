library(Matrix) # HAVE TO BE 1.5.3 (1.6 doesn't work)
# install.packages('https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.5-3.tar.gz', repos = NULL, type ='source')
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(plyr)
library(dplyr)
library(ggforce)
library(ggplot2)
library(scales)
library(reshape2)
library(cowplot)
library(preprocessCore) # have to be installed with additional param 
# BiocManager::install("preprocessCore", configure.args="--disable-threading", force = TRUE)
library(umap)
library(Rtsne)
library(pheatmap)
library(parallel)
library(ggrepel)
library(GeoDiff)

# running vignette
# https://bioconductor.org/packages/devel/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html

# other vignettes:
# https://www.bioconductor.org/packages/release/bioc/vignettes/GeomxTools/inst/doc/Developer_Introduction_to_the_NanoStringGeoMxSet.html

# other important packages - check vignettes!
# SpatialDecon
# standR
# GeoDiff - another QC and normalisation!!


# ?
# what exactly is 'module' and 'segment' (segment - ROI/AOI, including geometric segment areas(ROI) and subarea with different markers(AOI))
# module - comes from PKC file (pkc's name?)

# pData - metadate for segments (pheno)
# fData - metadata for genes (feature)

# get path to DCC, PKC and anno files -------------------------------------

datadir <- system.file("extdata", "WTA_NGS_Example",
                       package="GeoMxWorkflows")

DCCFiles <- dir(file.path(datadir, "dccs"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)

PKCFiles <- unzip(zipfile = dir(file.path(datadir, "pkcs"), pattern = ".zip$",
                                full.names = TRUE, recursive = TRUE))
SampleAnnotationFile <-
  dir(file.path(datadir, "annotation"), pattern = ".xlsx$",
      full.names = TRUE, recursive = TRUE)


# load as GeomX dataset ---------------------------------------------------

geomx_obj <-
  readNanoStringGeoMxSet(dccFiles = DCCFiles, 
                         pkcFiles = PKCFiles, # this goes into fData() - features (probes) annotation
                         phenoDataFile = SampleAnnotationFile, # this goes into pData() - protocol (samples) annotation
                         phenoDataSheet = "Template", #TODO what about this param? - An optional character string representing the excel sheet name containing the phenotypic data.
                         phenoDataDccColName = "Sample_ID",
                         protocolDataColNames = c("aoi", "roi"), #TODO adjust
                         experimentDataColNames = c("panel")) #TODO adjust

# inspect -----------------------------------------------------------------

pkcs <- annotation(geomx_obj)
modules <- gsub(".pkc", "", pkcs)

anno <- geomx_obj@phenoData@data # same as pData(geomx_obj)

count_mat <- dplyr::count(pData(geomx_obj), `slide name`, class, region, segment)


count_mat$`slide name` <- gsub("disease", "d",
                               gsub("normal", "n", count_mat$`slide name`))

test_gr <- gather_set_data(count_mat, 1:4)


test_gr$x <- mapvalues(test_gr$x, 
                               from=c(1, 2, 3, 4), 
                               to=c("slide name", "class", "region", "segment"))
test_gr$x <- factor(test_gr$x,
                    levels = c("class", "slide name", "region", "segment"))

#####
# count_mat <- dplyr::count(anno, across(all_of(anno_vect)))
# 
# if('slide name' %in% anno_vect){
# count_mat$`slide name` <- gsub("disease", "d",
#                                gsub("normal", "n", count_mat$`slide name`))
# }
# 
# test_gr <- gather_set_data(count_mat, 1:length(anno_vect))
# 
# 
# test_gr$x <- mapvalues(test_gr$x, 
#                        from=seq(1, length(anno_vect)), 
#                        to=anno_vect)
# test_gr$x <- factor(test_gr$x,
#                     levels = anno_vect)
# 
# 

# plot Sankey
ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = region), alpha = 0.5, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.2) +
  geom_parallel_sets_labels(color = "white", size = 5) +
  theme_classic(base_size = 17) + 
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank()) +
  scale_y_continuous(expand = expansion(0)) + 
  scale_x_discrete(expand = expansion(0)) +
  labs(x = "", y = "") +
  annotate(geom = "segment", x = 4.25, xend = 4.25,
           y = 20, yend = 120, lwd = 2) +
  annotate(geom = "text", x = 4.19, y = 70, angle = 90, size = 5,
           hjust = 0.5, label = "100 segments")


# QC for segments ---------------------------------------------------------

# Every ROI/AOI segment will be tested for:
#   
# Raw sequencing reads: segments with >1000 raw reads are removed.
# % Aligned,% Trimmed, or % Stitched sequencing reads: segments below ~80% for one or more of these QC parameters are removed.
# % Sequencing saturation ([1-deduplicated reads/aligned reads]%): segments below ~50% require additional sequencing to capture full sample diversity and are not typically analyzed until improved.
# Negative Count: this is the geometric mean of the several unique negative probes in the GeoMx panel that do not target mRNA and establish the background count level per segment; segments with low negative counts (1-10) are not necessarily removed but may be studied closer for low endogenous gene signal and/or insufficient tissue sampling.
# No Template Control (NTC) count: values >1,000 could indicate contamination for the segments associated with this NTC; however, in cases where the NTC count is between 1,000- 10,000, the segments may be used if the NTC data is uniformly low (e.g. 0-2 counts for all probes).
# Nuclei: >100 nuclei per segment is generally recommended; however, this cutoff is highly study/tissue dependent and may need to be reduced; what is most important is consistency in the nuclei distribution for segments within the study.
# Area: generally correlates with nuclei; a strict cutoff is not generally applied based on area.

# Shift 0 counts to one (needed for downstream analysis - ?)
geomx_obj <- shiftCountsOne(geomx_obj, useDALogic = TRUE)

# TODO have to be adjusted per each new dataset
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

#####
# Collate QC Results
QCResults <- protocolData(geomx_obj)[["QCFlags"]]

QC_Summary <- data.frame(Pass = colSums(!QCResults[, colnames(QCResults)]),
                         Warning = colSums(QCResults[, colnames(QCResults)]))

QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})

QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))
#####

col_by <- "segment"

# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

# visualise segment qc metrics
QC_histogram(sData(geomx_obj), "Trimmed (%)", col_by, qc_params[["percentTrimmed"]])

QC_histogram(sData(geomx_obj), "Stitched (%)", col_by, qc_params[["percentStitched"]])

QC_histogram(sData(geomx_obj), "Aligned (%)", col_by, qc_params[["percentAligned"]])

QC_histogram(sData(geomx_obj), "Saturated (%)", col_by, qc_params[["percentSaturation"]])

QC_histogram(sData(geomx_obj), "area", col_by, qc_params[["minArea"]], scale_trans = "log10")

QC_histogram(sData(geomx_obj), "nuclei", col_by, qc_params[["minNuclei"]])


############
# calculate negative count
# calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(geomx_obj), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 

protocolData(geomx_obj)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData  
# this is only to make plot - later on is detached from pData
negCols <- paste0("NegGeoMean_", modules)
pData(geomx_obj)[, negCols] <- sData(geomx_obj)[["NegGeoMean"]]

# visualise negative counts
for(ann in negCols) {
  plt <- QC_histogram(pData(geomx_obj), ann, col_by, qc_params[['minNegativeCount']], scale_trans = "log10") # why 2 and not 1 as in qc_metrics?
  print(plt)
}

# detatch neg_geomean columns ahead of aggregateCounts call
# just for plot - see above
pData(geomx_obj) <- pData(geomx_obj)[, !colnames(pData(geomx_obj)) %in% negCols]

# count segments with neg counts
# NTC - no template control
table(NTC_Count = sData(geomx_obj)$NTC)

#TODO decide what to do with information about NegGeoMean and NCT - so far anything happens here
#########
# background modelling

paste("## of Negative Probes:", sum(fData(geomx_obj)$Negative))
# This model estimates a feature factor for each negative probe and a background size factor for each ROI.
geomx_obj <- fitPoisBG(geomx_obj)
summary(pData(geomx_obj)$sizefact)
summary(fData(geomx_obj)$featfact[fData(geomx_obj)$Negative])

set.seed(123)
geomx_diag <- diagPoisBG(geomx_obj)
notes(geomx_diag)$disper 
# dispersion - should be <2, if it's higher there may be problem with the modules
# if the dispersion is higher, some ROIs may be removed and Pois can be run again
which(assayDataElement(geomx_diag, "low_outlier") == 1, arr.ind = TRUE)
which(assayDataElement(geomx_diag, "up_outlier") == 1, arr.ind = TRUE)

# or if we assume batch effect we may want to group by eg slide, or other group we have and check 
# if the distribution is better
# geomx_obj <- fitPoisBG(geomx_obj, groupvar = "slide name")
# set.seed(123)
# geomx_diag <- diagPoisBG(geomx_obj, split = TRUE)
# notes(geomx_diag)$disper_sp # why not disper? - check in documentation
# TODO decide if any segments should be removed based on this score


##########
# remove flagged segments

table(QCResults$QCStatus)

geomx_obj <- geomx_obj[, QCResults$QCStatus == "PASS"]


# QC for probes -----------------------------------------------------------
#######
# rmv probes with 0 counts
# may not be needed

all0probeidx <- which(rowSums(exprs(geomx_obj))==0)

if (length(all0probeidx) > 0) {
  geomx_obj <- geomx_obj[-all0probeidx, ]
}
#######

# A probe is removed globally from the dataset if either of the following is true:
# the geometric mean of that probe’s counts from all segments divided by the geometric mean 
# of all probe counts representing the target from all segments is less than 0.1
# the probe is an outlier according to the Grubb’s test in at least 20% of the segments
# A probe is removed locally (from a given segment) if the probe is an outlier according to 
# the Grubb’s test in that segment.
# We do not typically adjust these QC parameters.

geomx_obj <- setBioProbeQCFlags(geomx_obj, 
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = TRUE)

ProbeQCResults <- fData(geomx_obj)[["QCFlags"]]

# summarise probe qc results
qc_probe_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))
dim(geomx_obj)
# retain only probes that passed qc
geomx_obj <- 
  subset(geomx_obj, 
         fData(geomx_obj)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(geomx_obj)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)

dim(geomx_obj)

# gene lvl count data -----------------------------------------------------

# nr of unique targets
length(unique(featureData(geomx_obj)[["TargetName"]]))

# collapse features to targets
geomx_obj <- aggregateCounts(geomx_obj)

exprs(geomx_obj)[1:5, 1:2]



# filter based on limit of quantification per segment ---------------------

# The LOQ is calculated based on the distribution of negative control probes and is intended to approximate 
# the quantifiable limit of gene expression per segment. Please note that this process is more stable in larger 
# segments. Likewise, the LOQ may not be as accurately reflective of true signal detection rates in 
# segments with low negative probe counts (ex: <2). The formula for calculating the LOQ in the ith
# segment is:
#   
#   LOQi=geomean(NegProbei)∗geoSD(NegProbei)n
# 
# We typically use 2 geometric standard deviations (n=2
# ) above the geometric mean as the LOQ, which is reasonable for most studies. 
# We also recommend that a minimum LOQ of 2 be used if the LOQ calculated in a segment is below this threshold.

loq_cutoff <- 2
loq_min <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(geomx_obj))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(geomx_obj)))) {
    LOQ[, module] <-
      pmax(loq_min,
           pData(geomx_obj)[, vars[1]] * 
             pData(geomx_obj)[, vars[2]] ^ loq_cutoff)
  }
}

pData(geomx_obj)$LOQ <- LOQ

###################
# filtering

# After determining the limit of quantification (LOQ) per segment, 
# filtering out either segments and/or genes with abnormally low signal

LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(geomx_obj)$Module == module
  Mat_i <- t(esApply(geomx_obj[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}

# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(geomx_obj)$TargetName, ]

# TODO check this LOQ values for a bias regarding sample type
# like in the paper

# Save detection rate information to pheno data
# how many genes have been detected in each segment  
pData(geomx_obj)$GenesDetected <- colSums(LOQ_Mat, na.rm = TRUE)
pData(geomx_obj)$GeneDetectionRate <- pData(geomx_obj)$GenesDetected / nrow(geomx_obj)

#TODO #calculate signal/noise ratio = Count/LOQ per segment (similar to genedetectionrate)

# filtering out segments with exceptionally low signal:
# small fraction of panel genes detected above the LOQ relative to the other segments in the study.

# visualization of the distribution of segments with respect to their % genes detected

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(geomx_obj)$DetectionThreshold <- 
  cut(pData(geomx_obj)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(pData(geomx_obj),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = region)) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")

# check for tissue type
table(pData(geomx_obj)$DetectionThreshold,
      pData(geomx_obj)$class)

# TODO appropriate thr should be choosen here based on these plots + stats!!!

# filter based on thr choosen
gene_detect_thr <- 0.1

geomx_obj <- geomx_obj[, pData(geomx_obj)$GeneDetectionRate >= gene_detect_thr]

dim(geomx_obj)

# TODO plot again Sankey plot here

# Calculate gene detection rate:
# in how many segments the given gene was detected
LOQ_Mat <- LOQ_Mat[, colnames(geomx_obj)]
fData(geomx_obj)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(geomx_obj)$DetectionRate <- fData(geomx_obj)$DetectedSegments / nrow(pData(geomx_obj))

LOQ_Mat <- LOQ_Mat[fData(geomx_obj)$TargetName, ]


# Gene of interest detection table
# TODO would be good to change for genes in pathways of interest/marker genes for immune cells
goi <- c("PDCD1", "CD274", "IFNG", "CD8A", "CD68", "EPCAM",
         "KRT18", "NPHS1", "NPHS2", "CALB1", "CLDN8")

goi_df <- data.frame(
  Gene = goi,
  Number = fData(geomx_obj)[goi, "DetectedSegments"],
  DetectionRate = percent(fData(geomx_obj)[goi, "DetectionRate"]))


######
# plot detection rate for all genes
# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(geomx_obj)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(geomx_obj))
rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")

# TODO think if this should be done. there may be lowcount immune cells populations
# only in some samples - do it after estimating cell type/stated diversity
# or ignore this step completely 

# gene filtering based on detection rate
# TODO segment cutoff have to be determined based on this AND biological diversity

# We typically set a % Segment cutoff ranging from 5-20% based on the biological 
# diversity of our dataset. For this study, we will select 10% as our cutoff. 
# In other words, we will focus on the genes detected in at least 10% of our segments; 
# we filter out the remainder of the targets.

#Note: if we know that a key gene is represented in only a small number of segments (<10%) 
# due to biological diversity, we may select a different cutoff or keep the target gene by 
# manually selecting it for inclusion in the data object.

segment_detect_rate_thr <- 0.1

# Subset to target genes detected in at least 10% of the samples.
#   Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(geomx_obj), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)

# filter out genes detected in less then thr nr of segments 
geomx_obj <- 
  geomx_obj[fData(geomx_obj)$DetectionRate >= segment_detect_rate_thr |
                    fData(geomx_obj)$TargetName %in% neg_probes, ]
dim(geomx_obj)

# retain only detected genes of interest
goi <- goi[goi %in% rownames(geomx_obj)]



# filter based on background modelling ------------------------------------
# alternative to filtering based on LOQ from the previous section

geomx_obj <- fitPoisBG(geomx_obj, groupvar = "slide name")
# probe aggregation once again and storage in the other object for GeoDiff lib usage
# !!!!! geomx_obj <- aggregateCounts(geomx_obj) shouldnt be run before!!!!
geomx_obj <- aggreprobe(geomx_obj, use = "cor")

# Using the background score test, we can determine which targets are expressed 
# above the background of the negative probes across this dataset. We can then filter 
# the data to only targets above background, using a suggested pvalue threshold of 1e-3.

geomx_obj <- BGScoreTest(geomx_obj)
sum(fData(geomx_obj)[["pvalues"]] < 1e-3, na.rm = TRUE)
# removeoutlier = TRUE ??


# To estimate the signal size factor, we use the fit negative binomial threshold function. 
# This size factor represents technical variation between ROIs like sequencing depth
# The feature_high_fitNBth labeled genes are ones well above background that will be used in later steps.

set.seed(123)
geomx_obj <- fitNBth(geomx_obj, split = TRUE)

features_high <- rownames(fData(geomx_obj))[fData(geomx_obj)$feature_high_fitNBth == 1]
length(features_high)

# We can compare this threshold to the mean of the background as a sanity check.
# TODO why it is so very different from vignette? - something made different to geomx_obj from other vignette
# TODO check the exact meaning
bgMean <- mean(fData(geomx_obj)$featfact, na.rm = TRUE)
notes(geomx_obj)[["threshold"]]
bgMean

#This is a sanity check to see that the signal size factor and background size factor are correlated but not redundant.

cor(geomx_obj$sizefact, geomx_obj$sizefact_fitNBth)
plot(geomx_obj$sizefact, geomx_obj$sizefact_fitNBth, xlab = "Background Size Factor",
     ylab = "Signal Size Factor")
abline(a = 0, b = 1)

# !!!!!!!
# In this dataset, this size factor correlate well with different quantiles, including 75%
# quantile which is used in Q3 normalization.

# get only biological probes
posdat <- geomx_obj[-which(fData(geomx_obj)$CodeClass == "Negative"), ]
posdat <- exprs(posdat)

quan <- sapply(c(0.75, 0.8, 0.9, 0.95), function(y)
  apply(posdat, 2, function(x) quantile(x, probs = y)))

corrs <- apply(quan, 2, function(x) cor(x, geomx_obj$sizefact_fitNBth))
names(corrs) <- c(0.75, 0.8, 0.9, 0.95)

corrs

quan75 <- apply(posdat, 2, function(x) quantile(x, probs = 0.75))

#Quantile range (quantile - background size factor scaled by the mean 
#feature factor of negative probes) has better correlation with the signal size factor.

geomx_obj <- QuanRange(geomx_obj, split = FALSE, probs = c(0.75, 0.8, 0.9, 0.95))

corrs <- apply(pData(geomx_obj)[, as.character(c(0.75, 0.8, 0.9, 0.95))], 2, function(x)
  cor(x, geomx_obj$sizefact_fitNBth))

names(corrs) <- c(0.75, 0.8, 0.9, 0.95)

corrs

# normalisation -----------------------------------------------------------

# 1) standard normalisation - Q3
# relies a normalization factor that aligns the third quartile gene count value for all samples
# does not address potential global variations in data distributions or expression ranges

# 2) # quantile normalisation (rank-based method)
# https://www.sciencedirect.com/science/article/pii/S2589004222020338
# https://github.com/LevivanHijfte/NanoString_normalization_methods
# better correction if signal-to-noise ratio is different between samples 
# this may be assessed from LOQ
# quantile normalisation assumes that differences in global variation
# between the samples do not represent biological data and correct for this
# aligns count values according to the rank of the genes and forces the data 
# of all samples into the same distribution
# difference in count value range or signal-to-noise ratio is no longer
# dependent on raw data distributions.
# quantile normalization, due to its course method of normalization, may limit
# the detection of more subtle differences in gene expression



# Q3 normalisation --------------------------------------------------------

# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "region"
Stat_data <- 
  data.frame(row.names = colnames(exprs(geomx_obj)),
             Segment = colnames(exprs(geomx_obj)),
             Annotation = pData(geomx_obj)[, ann_of_interest],
             Q3 = unlist(apply(exprs(geomx_obj), 2,
                               quantile, 0.75, na.rm = TRUE)),
             NegProbe = exprs(geomx_obj)[neg_probes, ])

Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")

plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
  geom_histogram(bins = 40) + theme_bw() +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) + 
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "Segments, #")

plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                     rel_widths = c(0.43,0.57))
plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))


# Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
geomx_obj <- normalize(geomx_obj ,
                             norm_method = "quant", 
                             desiredQuantile = .75,
                             toElt = "q3_norm")

# plot effects of normalisation
boxplot(exprs(geomx_obj)[,1:10],
        col = "#9EDAE5", main = "Raw Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Raw")

boxplot(assayDataElement(geomx_obj[,1:10], elt = "q3_norm"),
        col = "#2CA02C", main = "Q3 Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Q3 Normalized")


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

boxplot(assayDataElement(geomx_obj[,1:10], elt = "quant_norm"),
        col = "#2CA02C", main = "Quantile Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Quantile Normalized")


# UMAP + tSNE -------------------------------------------------------------

# update defaults for umap to contain a stable random_state (seed)
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42

# run UMAP on Q3 and quantile norm
umap_out <-
  umap(t(log2(assayDataElement(geomx_obj , elt = "q3_norm"))),  
       config = custom_umap)
umap_out_quant <-
  umap(t(log2(assayDataElement(geomx_obj , elt = "quant_norm"))),  
       config = custom_umap)

pData(geomx_obj)[, c("UMAP1_q3_norm", "UMAP2_q3_norm")] <- umap_out$layout[, c(1,2)]
pData(geomx_obj)[, c("UMAP1_quant_norm", "UMAP2_quant_norm")] <- umap_out_quant$layout[, c(1,2)]

ggplot(pData(geomx_obj),
       aes(x = UMAP1_q3_norm, y = UMAP2_q3_norm, color = region, shape = class)) +
  geom_point(size = 3) +
  theme_bw()

ggplot(pData(geomx_obj),
       aes(x = UMAP1_quant_norm, y = UMAP2_quant_norm, color = region, shape = class)) +
  geom_point(size = 3) +
  theme_bw()


# run tSNE
set.seed(42) # set the seed for tSNE as well
tsne_out <-
  Rtsne(t(log2(assayDataElement(geomx_obj , elt = "q3_norm"))),
        perplexity = ncol(geomx_obj)*.15)

tsne_out_quant <-
  Rtsne(t(log2(assayDataElement(geomx_obj , elt = "quant_norm"))),
        perplexity = ncol(geomx_obj)*.15)

pData(geomx_obj)[, c("tSNE1_q3_norm", "tSNE2_q3_norm")] <- tsne_out$Y[, c(1,2)]
pData(geomx_obj)[, c("tSNE1_quant_norm", "tSNE2_quant_norm")] <- tsne_out_quant$Y[, c(1,2)]


ggplot(pData(geomx_obj),
       aes(x = tSNE1_q3_norm, y = tSNE2_q3_norm, color = region, shape = class)) +
  geom_point(size = 3) +
  theme_bw()

ggplot(pData(geomx_obj),
       aes(x = tSNE1_quant_norm, y = tSNE2_quant_norm, color = region, shape = class)) +
  geom_point(size = 3) +
  theme_bw()

# TODO add clustering based on UMAP
# find genes with high coefficient variation ------------------------------

# create a log2 transform of the data for analysis
assayDataElement(object = geomx_obj, elt = "log_q3") <-
  assayDataApply(geomx_obj, 2, FUN = log, base = 2, elt = "q3_norm")

assayDataElement(object = geomx_obj, elt = "log_quant") <-
  assayDataApply(geomx_obj, 2, FUN = log, base = 2, elt = "quant_norm")

# create CV function
calc_CV <- function(x) {sd(x) / mean(x)}

CV_dat <- assayDataApply(geomx_obj,
                         elt = "log_q3", MARGIN = 1, calc_CV)
CV_dat_quant <- assayDataApply(geomx_obj,
                         elt = "log_quant", MARGIN = 1, calc_CV)

# show the highest CD genes and their CV values
sort(CV_dat, decreasing = TRUE)[1:10]
sort(CV_dat_quant, decreasing = TRUE)[1:10]

# Identify genes in the top 3rd of the CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8)]
GOI_quant <- names(CV_dat_quant)[CV_dat_quant > quantile(CV_dat_quant, 0.8)]

pheatmap(assayDataElement(geomx_obj[GOI, ], elt = "log_q3"),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = 
           pData(geomx_obj)[, c("class", "segment", "region")])

pheatmap(assayDataElement(geomx_obj[GOI_quant, ], elt = "log_quant"),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = 
           pData(geomx_obj)[, c("class", "segment", "region")])


# differential gene expression --------------------------------------------

# The LMM allows the user to account for the subsampling per tissue,
# we adjust for the fact that the multiple regions of interest placed per tissue section 
# are not independent observations, as is the assumption with other traditional statistical tests.

# two flavors of the LMM model when used with GeoMx data: i) with and ii) without random slope.
# 
# When comparing features that co-exist in a given tissue section (e.g. glomeruli vs tubules in DKD kidneys),
# a random slope is included in the LMM model. 
# When comparing features that are mutually exclusive in a given tissue section 
# (healthy glomeruli versus DKD glomeruli) the LMM model does not require a random slope. 

# TODO determine if we have to use LMM w or wo random slope for:
# pre vs post - wo random slope
# different macrophage/t-cell rich regions - w random slope


# within slide DEG --------------------------------------------------------
# with random slope

# TODO figure out about test variables + tissue (if more tissue/slide etc)
# test region and status reversed in the latter analysis (?)
# convert test variables to factors
pData(geomx_obj)$testRegion <- 
  factor(pData(geomx_obj)$region, c("glomerulus", "tubule"))
pData(geomx_obj)[["slide"]] <- 
  factor(pData(geomx_obj)[["slide name"]])

#TODO same for quant norm
# run LMM:
# formula follows conventions defined by the lme4 package
results <- c()
for(status in c("DKD", "normal")) {
  ind <- pData(geomx_obj)$class == status
  mixedOutmc <-
    mixedModelDE(geomx_obj[, ind],
                 elt = "log_q3",
                 modelFormula = ~ testRegion + (1 + testRegion | slide), # random slope
                 groupVar = "testRegion",
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
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results <- rbind(results, r_test)
}

# the log2 fold change value (Estimate), P-value (Pr(>|t|))
# false-discovery adjusted P-values (FDR).
# The contrast column is used to interpret the log2 fold change value as it specifies 
# which levels are compared (e.g. positive fold change values when comparing glomerulus - tubule 
# indicates an enrichment in the glomerulus; negative indicates enrichment in tubules).


# between slide DEG -------------------------------------------------------
# wo random slope

# convert test variables to factors
pData(geomx_obj)$testClass <-
  factor(pData(geomx_obj)$class, c("normal", "DKD"))

# TODO run also for log_quant
# run LMM:
# formula follows conventions defined by the lme4 package
results2 <- c()
for(region in c("glomerulus", "tubule")) {
  ind <- pData(geomx_obj)$region == region
  mixedOutmc <-
    mixedModelDE(geomx_obj[, ind],
                 elt = "log_q3",
                 modelFormula = ~ testClass + (1 | slide), # no random slope, comp: (~ testRegion + (1 + testRegion | slide))
                 groupVar = "testClass",
                 nCores = parallel::detectCores(),
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
  r_test$Subset <- region
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results2 <- rbind(results2, r_test)
}


# DGE visualisation -------------------------------------------------------

# Categorize Results based on P-value & FDR for plotting
# TODO please change it into sth normal
results$Color <- "NS or FC < 0.5"
results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results$Color[results$FDR < 0.05] <- "FDR < 0.05"
results$Color[results$FDR < 0.001] <- "FDR < 0.001"
results$Color[abs(results$Estimate) < 0.5] <- "NS or FC < 0.5"
results$Color <- factor(results$Color,
                        levels = c("NS or FC < 0.5", "P < 0.05",
                                   "FDR < 0.05", "FDR < 0.001"))

# pick top genes for either side of volcano to label
# order genes for convenience:
results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)
top_g <- c()
for(cond in c("DKD", "normal")) {
  ind <- results$Subset == cond
  top_g <- c(top_g,
             results[ind, 'Gene'][
               order(results[ind, 'invert_P'], decreasing = TRUE)[1:15]],
             results[ind, 'Gene'][
               order(results[ind, 'invert_P'], decreasing = FALSE)[1:15]])
}
top_g <- unique(top_g)
results <- results[, -1*ncol(results)] # remove invert_P from matrix WTF why not [, -c('invert_P')]

# Graph results
ggplot(results,
       aes(x = Estimate, y = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "Enriched in Tubules <- log2(FC) -> Enriched in Glomeruli",
       y = "Significance, -log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS or FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(results, Gene %in% top_g & FDR < 0.001),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom") +
  facet_wrap(~Subset, scales = "free_y")


# show expression for a single target
gene_name <- "PDHA1"

ggplot(pData(geomx_obj),
       aes(x = region, fill = region,
           y = assayDataElement(geomx_obj[gene_name, ],
                                elt = "q3_norm"))) +
  geom_violin() +
  geom_jitter(width = .2) +
  labs(y = "PDHA1 Expression") +
  scale_y_continuous(trans = "log2") +
  facet_wrap(~class) +
  theme_bw()

# show expression of 2 targets
gene_names <- c("PDHA1", "ITGB1")

glom <- pData(geomx_obj)$region == "glomerulus"

# show expression of PDHA1 vs ITGB1
# TODO fix this plot
ggplot(assayDataElement(geomx_obj, elt = "q3_norm")) +
  geom_vline(xintercept =
               max(assayDataElement(geomx_obj[gene_names[1], glom],
                                    elt = "q3_norm")),
             lty = "dashed", col = "darkgray") +
  geom_hline(yintercept =
               max(assayDataElement(geomx_obj[gene_names[2], !glom],
                                    elt = "q3_norm")),
             lty = "dashed", col = "darkgray") +
  geom_point(aes(x = assayDataElement(geomx_obj[gene_names[1], ],
                                      elt = "q3_norm"),
                 y = assayDataElement(geomx_obj[gene_names[2], ],
                                      elt = "q3_norm"),
                 color = region), size = 3) +
  theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  labs(x = "PDHA1 Expression", y = "ITGB1 Expression") +
  facet_wrap(~class)


# heatmap with significant genes
# select top significant genes based on significance, plot with pheatmap
GOI <- unique(subset(results, `FDR` < 0.001)$Gene)
pheatmap(log2(assayDataElement(geomx_obj[GOI, ], elt = "q3_norm")),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         cutree_cols = 2, cutree_rows = 2,
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = pData(geomx_obj)[, c("region", "class")])


# MA plot

results$MeanExp <-
  rowMeans(assayDataElement(geomx_obj,
                            elt = "q3_norm"))

top_g2 <- results$Gene[results$Gene %in% top_g &
                         results$FDR < 0.001 &
                         abs(results$Estimate) > .5 &
                         results$MeanExp > quantile(results$MeanExp, 0.9)]

ggplot(subset(results, !Gene %in% neg_probes),
       aes(x = MeanExp, y = Estimate,
           size = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
  geom_hline(yintercept = c(0.5, -0.5), lty = "dashed") +
  scale_x_continuous(trans = "log2") +
  geom_point(alpha = 0.5) + 
  labs(y = "Enriched in Glomeruli <- log2(FC) -> Enriched in Tubules",
       x = "Mean Expression",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS or FC < 0.5` = "gray")) +
  geom_text_repel(data = subset(results, Gene %in% top_g2),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2) +
  theme_bw(base_size = 16) +
  facet_wrap(~Subset, nrow = 2, ncol = 1)
