# Ligand Receptor Analysis  by BulkSignaR : Bulk Data,Geomx Full Transcriptomic signal



# load libraries

library(BulkSignalR)
library(igraph)
library(dplyr)

# file names

input_dir <- 'C:/Users/Sahas/Downloads/Masters_Thesis/Data/Batch01/'
output_dir <- 'C:/Users/Sahas/Downloads/Masters_Thesis/Ligand-receptor/BulkSignalR/TestRun_New_scripts/'


# parameters


null.model = NULL # c("automatic", "mixedNormal", "normal", "kernelEmpirical","empirical", "stable")

# To do

normalize_needed = TRUE # Default "UQ", FALSE if the provided data are normalized 

normalize_method = "UQ" # c("UQ",TC") or user defined if the provided data are normalized. ('UQ' for upper quartile or 'TC' for total count. UQ.pc = 0.75.
UQ_pc = 0.75
paired_only = FALSE # for paired samples only TRUE
nact_status = NULL
segment = NULL # tumor, stroma or NULL (to get all the data)
annotation = NULL
qval_threshold = 0.01 # filter significant LR pairs


# load data

geomx_obj <- readRDS(paste0(input_dir,'geomx_qc_norm_batch_eff_rm.RDS'))
count_geomx = data.frame(geomx_obj@assayData$exprs) # count data
#count_geomx  = data.frame(geomx_obj@assayData$q3_norm) # q3 normalized data


# create directory for output files

output_folder_name  = "/BulkSignalR_objects"

if (!dir.exists(paste0(output_dir,output_folder_name))) {
  new_dir = paste0(output_dir,output_folder_name)
  dir.create(new_dir,recursive = TRUE)
  message("Directory created")
  output_dir = new_dir
} else {
  output_dir = paste0(output_dir,output_folder_name)
}

# Filtering data based on nact_status, segment and annotation

meta_data = sData(geomx_obj )
meta_data = data.frame(meta_data[,c(2,5,6,7,24,25,28)])
meta_data$dcc_filename <- gsub('-', '.', meta_data$dcc_filename)

if (!is.null(nact_status)) {
  meta_data <- meta_data %>% filter(NACT_status == nact_status)
}

if (!is.null(segment)) {
  meta_data <- meta_data %>% filter(Segment == segment)
}

if (!is.null(annotation)) {
  meta_data <- meta_data %>% filter(Annotation == annotation)
}


if (paired_only == TRUE) { # check !!
  
  paired_samples <- meta_data %>%
    group_by(Patient) %>%
    filter(all(c("pre", "post") %in% NACT_status)) %>%
    summarise() %>%
    pull(Patient)

  meta_data = meta_data %>% filter(Patient %in% paired_samples)
  
}

col_ids = colnames(count_geomx)  %in% meta_data$dcc_filename
count_geomx = count_geomx[,col_ids]



### Analysis : written according to the BulkSignalR Vignette
# browseVignettes("BulkSignalR")

# step 01 : Prepare Dataset

bsrdm <- prepareDataset(counts = count_geomx, normalize = normalize_needed , method = normalize_method, UQ.pc = UQ_pc, log.transformed = FALSE, min.count = 10, prop = 0.1) 

# step 02 : learnParameters

set.seed(123)
bsrdm <- learnParameters(bsrdm, 
                         plot.folder = file.path(output_dir), filename = paste0("geomxUQ_",nact_status,'_',segment,'_',annotation), verbose = TRUE)


# step 03 : Building a BSRInference object


bsrinf <- initialInference(bsrdm)
LRinter.dataframe <- LRinter(bsrinf)
LRinter.dataframe <- LRinter.dataframe[order(LRinter.dataframe$qval <= qval_threshold),]



# reducing to best pathways before calculating signature scores

bsrinf.redBP    <- reduceToBestPathway(bsrinf) 
LRinter_pairs_best_pws = LRinter(bsrinf.redBP)


# Save Objects

BulkSignalR_output = list(
  bsrdm = bsrdm,
  bsrinf = bsrinf,
  bsrinf_redBP = bsrinf.redBP,
  LRinter_dataframe = LRinter.dataframe,
  LRinter_pairs_best_pws = LRinter_pairs_best_pws,
  meta_data = meta_data
)


parts <- c(nact_status, segment, annotation)
parts_non_null <- parts[!sapply(parts, is.null)]
if (length(parts_non_null) == 0) {
  saveRDS(BulkSignalR_output, file = paste0(output_dir,'/BulkSignalR_combined_output.RDS'))
} else {
  file_name = paste(parts_non_null, collapse = "_")
  saveRDS(BulkSignalR_output, file = paste0(output_dir,'/BulkSignalR_',file_name,'_output.RDS'))
}




