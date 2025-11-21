library(data.table)
library(dplyr)
library(readxl)
library(rvest)

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
# clean somatic variant

snv_eyemt_important <- filter(snv_eyemt, TIER %in% c('TIER 1', 'TIER 2', 'TIER 3'))
sort(table(snv_eyemt_important$folder_source))

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
