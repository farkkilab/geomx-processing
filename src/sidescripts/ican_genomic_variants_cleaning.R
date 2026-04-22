library(data.table)
library(readxl)
#library(rvest)
library(plyr)
library(dplyr)
library(tidyr)
library(tibble)

# for S131 there was no proper SNV export:
# info from iCAN: total nr of coding mutations: 79
# 74 - tier4, 5 - tier3
# TP53 missense_variant p.Ile195Thr
# ERBB4 missense_variant p.Ser303Phe
# FAT4 p.Glu1673Ala
# PLXNB1 p.Pro850Thr
# CBL p.Cys384Tyr

eyemt_clin_path <- '/home/iganiemi/Documents/phd/st/data/geomx/clinical_data/9_eyemt_patient_clinical_data_SENSITIVE.csv' 
id_path <- '/media/ldrive/ltdk_farkkila/Projects/12-SPACE/HBP2019007_iCan_raportit_03_2025_EXT/studyId_Ican_WES_key_06_2025.csv'
snv_path <- '/media/ldrive/ltdk_farkkila/Projects/12-SPACE/HBP2019007_iCan_raportit_03_2025_EXT/iMTB_report/merged_global_SNV_variants.xlsx'
germ_clinvar_path <- '/media/ldrive/ltdk_farkkila/Projects/12-SPACE/HBP2019007_iCan_raportit_03_2025_EXT/iMTB_report/merged_CPSR_germline_ClinVar_variants.xlsx'
germ_biomarker_path <- '/media/ldrive/ltdk_farkkila/Projects/12-SPACE/HBP2019007_iCan_raportit_03_2025_EXT/iMTB_report/merged_CPSR_germline_biomarkers.xlsx'
imtb_report_dir <- '/media/ldrive/ltdk_farkkila/Projects/12-SPACE/HBP2019007_iCan_raportit_03_2025_EXT/iMTB_report'


eyemt_clin <- fread(eyemt_clin_path)
id_dt <- fread(id_path) %>%
  filter(id_dt, study_id %in% eyemt_clin$Patient)

# TODO check these patients? are they even on iCAN?
setdiff(unique(eyemt_clin$Patient), id_dt$study_id)

snv_eyemt <- read_excel(snv_path) %>%
  filter(folder_source %in% id_dt$return_pseudo) %>%
  filter(TIER != 'NONCODING')

# probably not needed
germ_clinvar_eyemt <- read_excel(germ_clinvar_path) %>%
  filter(folder_source %in% id_dt$return_pseudo) %>%
  filter(CLINVAR_CLASSIFICATION != 'VUS')

germ_biomarker_eyemt <- read_excel(germ_biomarker_path) %>%
  filter(folder_source %in% id_dt$return_pseudo)

# only 1 germline BRCA1 mutation - already mentioned in out clinical data

fwrite(snv_eyemt, '/media/ldrive/ltdk_farkkila/Projects/9_EyeMT/9_eyemt_somatic_snv_coding_imtb.csv')
fwrite(germ_clinvar_eyemt, '/media/ldrive/ltdk_farkkila/Projects/9_EyeMT/9_eyemt_germline_clinvar_significant_imtb.csv')
fwrite(id_dt, '/media/ldrive/ltdk_farkkila/Projects/9_EyeMT/9_eyemt_ican_keys.csv')

################################################
sort(table(snv_eyemt$folder_source))
table(snv_eyemt$TIER)

sort(table(germ_biomarker_eyemt$folder_source))

################################################
snv <- fread('~/Documents/phd/st/data/geomx/clinical_data/9_eyemt_somatic_snv_coding_imtb_SENSITIVE.csv') # same as snv_eyemt, just copied to the disc
# count all somatic variants
# calculate nr of mutations per tier + sum of all
snv_per_tier <- group_by(snv, study_id, TIER) %>%
  summarise(n = n()) %>%
  spread(key = 'TIER', value = 'n') %>%
  replace(is.na(.), 0) %>%
  mutate(all = rowSums(across(where(is.numeric)))) %>%
  select(-`TIER 4`)

colnames(snv_per_tier) <- c('Patient', 'snv_coding_nr_tier2', 'snv_coding_nr_tier3', 'snv_coding_number_total')

#manually ad patient which was not exported
s131_dt <- data.frame('S131', 0, 5, 79)
colnames(s131_dt) <- colnames(snv_per_tier)
snv_per_tier <- rbind(snv_per_tier, s131_dt) 

tp53_snv <- snv %>%
  filter(SYMBOL == 'TP53') %>%
  select(study_id, TIER) %>%
  mutate(TIER = gsub(' ', '_', TIER)) 

tp53_snv <- rbind(tp53_snv, data.frame(study_id = 'S131', TIER = 'TIER_3')) #manual add S131
colnames(tp53_snv) <- c('Patient', 'TP53_mut')

snv_summary <- left_join(snv_per_tier, tp53_snv) %>%
  replace(is.na(.), 'wt')

fwrite(snv_summary, '~/Documents/phd/st/data/geomx/clinical_data/9_eyemt_imtb_snv_summary.csv')

################################################
# fetch TMB information form htmpl reports

reports <- list.files(imtb_report_dir, pattern = 'html', full.names = T)
reports <- reports[grepl(paste(id_dt$return_pseudo, collapse = '|'), reports)]

tmb_all <- sapply(reports, function(rep_path){
  print('$$$$$$$$$')
  print(rep_path)
  report <- read_html(rep_path)
  paragraphs <- report %>% html_nodes("p")
  tmb <- unlist(strsplit(as.character(paragraphs[9]), split = '>'))
  print(tmb)
  tmb <- tmb[grepl('[0-9]*mutations', tmb)]
  print(tmb)
  tmb <- as.numeric(gsub('\\s[a-z]*.*', '', tmb))
  print(tmb)
  return(tmb)
})

ids <- sapply(reports, function(rep_path){
  id <- gsub('_pcgr_acmg.grch38.html', '', basename(rep_path))
})

tmb_df <- data.frame(return_pseudo = ids, TMB = tmb_all)

tmb_df <- left_join(tmb_df, id_dt)

# merge with clinical data
eyemt_clin <- left_join(eyemt_clin, tmb_df[, c('TMB', 'study_id')], by = c('Patient' = 'study_id'))

fwrite(eyemt_clin, eyemt_clin_path)
