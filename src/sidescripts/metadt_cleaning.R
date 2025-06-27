library(data.table)
library(dplyr)
#library(readxls)
library(plyr)
library(openxlsx)
library(tibble)

#TODO S188 check if it was in batch2 or excluded

df <- fread('~/Documents/phd/st/data/geomx/clinical_data/9_eyemt_patient_clinical_data_SENSITIVE2.csv')

df_per_pt <- df[, c('Patient', 'PFS_months', 'OS_months')]
df_per_pt <- df_per_pt[!duplicated(df_per_pt), ]

batches_list <- list(b123 = c(1,2,3), b12 = c(1,2), b1 = c(1), b2 = c(2), b3 = c(3), b23 = c(2,3)) 

##########
for(batchname in names(batches_list)){
  batch_nr <- batches_list[[batchname]]
  
  df_per_pt_batch <- df_per_pt[df_per_pt$Patient %in% unique(df$Patient[df$main_batch_nr %in% batch_nr]),]
  
  pfs_quant <- as.numeric(quantile(df_per_pt_batch$PFS_months, probs = c(0.25, 0.5, 0.75)))
  os_quant <- as.numeric(quantile(df_per_pt_batch$OS_months, probs = c(0.25, 0.5, 0.75)))
  
  # hist(df_per_pt_batch$PFS_months, breaks = 40)
  # hist(df_per_pt_batch$OS_months, breaks = 40)
  
  df_per_pt_batch$PFS_quartile <- as.integer(cut(df_per_pt_batch$PFS_months, quantile(df_per_pt_batch$PFS_months, probs=0:4/4), include.lowest=TRUE))
  df_per_pt_batch$OS_quartile <- as.integer(cut(df_per_pt_batch$OS_months, quantile(df_per_pt_batch$OS_months, probs=0:4/4), include.lowest=TRUE))
  
  df_per_pt_batch <- df_per_pt_batch[, c('Patient', 'PFS_quartile', 'OS_quartile')]
  colnames(df_per_pt_batch) <- c('Patient', paste0(c('PFS_quartile', 'OS_quartile'), '_', batchname))
  
  df <- left_join(df, df_per_pt_batch, by = 'Patient')
}

# clean from sensitive data
df$HRP_status <- ifelse(df$HRP_status == 'HRP', 1, ifelse(df$HRP_status == 'HRD', 0, NA))
df$BRCA_status <- ifelse(df$BRCA_status == 'wt', 0, ifelse(grepl('BRCA', df$BRCA_status), 1, NA))

df <- df[, -c("PFS_months", "OS_months", "R0", "ovaHRD_score", "GIS_HRD_score", "selected block", "dead", "progression")]

fwrite(df, '~/Documents/phd/st/data/geomx/clinical_data/9_eyemt_patient_clinical_data.csv')


##############################

kk <- readRDS('~/Documents/phd/st/geomx-processing/results/batch2-1903/deconvolution/bayes_prism/bp_res_mid_lvl_ct_expr_mtx_cleaned_vst_harmony_batch_corr.RDS')
#
metab1 <- read_xlsx('~/Documents/phd/st/data/geomx/geomx_batch1_0823/metadata/dcc_metadata_batch1_0823.xlsx')

pt_data <- fread('~/Documents/phd/st/data/geomx/clinical data/9_eyemt_patient_clinical_data.csv')

metab1 <- left_join(metab1, pt_data, by = c('Sample'))
metab1 <- subset(metab1, select=-c(PFS_quartile, OS_quartile))


############################
# combine batches metadata and pt clinical data
# also combine with clinical for b1,2,3 separately
pt_data <- fread('~/Documents/phd/st/data/geomx/clinical_data/9_eyemt_patient_clinical_data.csv')

b1 <- read.xlsx('~/Documents/phd/st/data/geomx/geomx_batch1_0823/metadata/dcc_metadata_batch1_0823.xlsx')
b2 <- read.xlsx('~/Documents/phd/st/data/geomx/geomx_batch2_1124/metadata/dcc_metadata_batch2_1124.xlsx')
b3 <- read.xlsx('~/Documents/phd/st/data/geomx/geomx_batch3_0525/metadata/dcc_metadata_all_batch3_0525.xlsx')


b12 <- rbind.fill(b1, b2)
length(unique(b12$dcc_filename))
b12 <- left_join(b12, pt_data[, c('Sample', 'HRP_status', 'BRCA_status', 'PFS_quartile_b12', 'OS_quartile_b12')], by = c('Sample'))
b12$batch_nr <- paste0(b12$main_batch_nr, '_', b12$batch_nr)
b12$batch_nr_sample_collection <- paste0(b12$main_batch_nr, '_', b12$batch_nr_sample_collection)
write.xlsx(b12, '~/Documents/phd/st/data/geomx/batch12/metadata/dcc_metadata_batch12.xlsx')

b23 <- rbind.fill(b2, b3)
length(unique(b23$dcc_filename))
b23 <- left_join(b23, pt_data[, c('Sample', 'HRP_status', 'BRCA_status', 'PFS_quartile_b23', 'OS_quartile_b23')], by = c('Sample'))
b23$batch_nr <- paste0(b23$main_batch_nr, '_', b23$batch_nr)
b23$batch_nr_sample_collection <- ifelse(!is.na(b23$batch_nr_sample_collection),
                                         paste0(b23$main_batch_nr, '_', b23$batch_nr_sample_collection), b23$batch_nr)
write.xlsx(b23, '~/Documents/phd/st/data/geomx/batch23/metadata/dcc_metadata_batch23.xlsx')


b123 <- rbind.fill(b1, b2, b3)
length(unique(b123$dcc_filename))
b123 <- left_join(b123, pt_data[, c('Sample', 'HRP_status', 'BRCA_status', 'PFS_quartile_b123', 'OS_quartile_b123')], by = c('Sample'))
b123$batch_nr <- paste0(b123$main_batch_nr, '_', b123$batch_nr)
b123$batch_nr_sample_collection <- ifelse(!is.na(b123$batch_nr_sample_collection),
                                         paste0(b123$main_batch_nr, '_', b123$batch_nr_sample_collection), b123$batch_nr)

write.xlsx(b123, '~/Documents/phd/st/data/geomx/batch123/metadata/dcc_metadata_batch123.xlsx')


#############################
# extract metadata and expr
b12 <- readRDS('~/Documents/phd/st/geomx-processing/results/batch12-1004/geomx_qc_norm_batch_eff_rm.RDS')

expr <- as.data.frame(b12@assayData$harmony_batch_corr)
expr <- rownames_to_column(expr, 'gene_name')
mt <- as.data.frame(sData(b12))
mt <- mt[, -49]
mt <- mt[, 1:91]



fwrite(expr, '~/Documents/phd/st/geomx-processing/results/batch12-1004/geomx_qc_norm_batch_eff_rm_harmony_expr.csv')
fwrite(mt, '~/Documents/phd/st/geomx-processing/results/batch12-1004/geomx_qc_norm_batch_eff_rm_harmony_metadata.csv')
