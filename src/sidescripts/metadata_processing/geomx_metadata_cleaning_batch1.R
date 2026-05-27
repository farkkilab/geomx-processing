library(data.table)
library(dplyr)
library(readxl)
library(tibble)
library(xlsx)

proj_dir <- '/media/iganiemi/T7-iga/st/'
dcc_dir <- file.path(proj_dir, 'data/geomx/nact_experiment/dcc')
fastq_dir <- file.path(proj_dir, 'data/geomx/nact_experiment/fastq')
stats_dir <- file.path(proj_dir, 'data/geomx/nact_experiment/stats')

clin_path <- file.path(proj_dir, 'data/geomx/nact_experiment/metadata/GeoMX sample set_full_info_12_07_2023.xlsx')
batch_path <- file.path(proj_dir, 'data/geomx/nact_experiment/metadata/batch_annotations.xlsx')
data_anno_path <- file.path(proj_dir, 'data/geomx/nact_experiment/metadata/download_data.xlsx')

########################################################################

dcc_files <- list.files(dcc_dir)
dcc_files <- data.frame(dcc_filename = dcc_files, Sample_ID = gsub('\\.dcc', '', dcc_files))

############################################################
# 1 less fastq and stats files than dcc!
fastq_files <- list.files(fastq_dir, pattern = 'c4gh')
stats_files <- list.files(stats_dir, pattern = 'c4gh')

st <- gsub('\\.stats\\.c4gh', '', stats_files)
ff <- unique(sapply(fastq_files, function(x){unlist(strsplit(x, split='_'))[1]}))

nn <- setdiff(dcc_files$Sample_ID, st)
nn2 <- setdiff(dcc_files$Sample_ID, ff)
############################################################
# load batch data and merge with dcc files

batch_all <- lapply(seq(1:8), function(x){
  batch <- as.data.frame(read_xlsx(batch_path, sheet = x))
  batch$batch_nr <- x
  return(batch)
})

batch_all <- do.call(rbind, batch_all)

length(unique(dcc_files$Sample_ID))
length(unique(batch_all$Sample_ID))
length(intersect(dcc_files$Sample_ID, batch_all$Sample_ID))
nrow(batch_all) - nrow(dcc_files)

# get duplicated ids
n_occur <- data.frame(table(batch_all$Sample_ID))
dup_ids <- n_occur[n_occur$Freq > 1, ]
batch_duplicated <- batch_all[batch_all$Sample_ID %in% as.character(dup_ids$Var1), ] %>%
  arrange(Sample_ID)

# 1 NTCc id is duplicated - contain the same info, but coming from batch 4 and 6

# remove duplicated rows
batch_all <- batch_all %>%
  distinct(across(-batch_nr), .keep_all = T)

length(intersect(dcc_files$Sample_ID, batch_all$Sample_ID))

dcc_data <- left_join(dcc_files, batch_all, by = 'Sample_ID')

# change Sample names without NACT status into 'post'
dcc_data <- mutate(dcc_data, Sample = ifelse(grepl('_', Sample), Sample, paste0(Sample, '_post'))) %>%
  mutate(Sample = ifelse(Sample == 'NA_post', 'NA', Sample))

###################################################\
# load clinical data and merge with dcc table

clin <- read_xlsx(clin_path, sheet = 'Sheet1')
clin <- clin %>%
  mutate(Notes = ifelse(PFS_months == 'no pregression', 'no progression', Notes)) %>%
  mutate(PFS_months = ifelse(PFS_months == 'no pregression', NA, PFS_months)) %>%
  mutate(PFS_months = as.numeric(PFS_months)) %>%
  mutate(`NACT status` = tolower(`NACT status`)) %>%
  mutate(Sample = paste0(Patient, '_', `NACT status`))

sort(unique(dcc_data$Sample))
sort(unique(clin$Sample))
length(unique(dcc_data$Sample))
length(unique(clin$Sample))
length(intersect(dcc_data$Sample, clin$Sample))

dcc_data <- left_join(dcc_data, clin, by = 'Sample')

# clean some variables
#TODO move it to annotation making script
# make cell annotation and fix some values
dcc_data$Annotation_cell <- gsub("S[0-9]*_", "", dcc_data$Annotation)
dcc_data$Annotation_cell <- gsub("pre_|post_", "", dcc_data$Annotation_cell)
dcc_data$Annotation_cell <- ifelse(dcc_data$`Slide Name` == 'No Template Control', NA, dcc_data$Annotation_cell)
dcc_data$Annotation_cell <- gsub("posBA1", "posIBA1", dcc_data$Annotation_cell)
dcc_data$Annotation_cell <- gsub("negBA1", "negIBA1", dcc_data$Annotation_cell)

dcc_data$Segment <- tolower(dcc_data$Segment)


fwrite(dcc_data, file.path(proj_dir, 'data/geomx/nact_experiment/metadata/dcc_metadata_all.csv'))
write.xlsx(dcc_data, file.path(proj_dir, 'data/geomx/nact_experiment/metadata/dcc_metadata_all.xlsx'))

# the one dcc file with no fastq
nofastq <- filter(dcc_data, Sample_ID == nn)
fwrite(nofastq, file.path(proj_dir, 'data/geomx/nact_experiment/metadata/missing_fastq.csv'))
