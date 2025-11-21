library(data.table)
library(plyr)
library(dplyr)
library(readxl)
library(phsmethods)

template_meta <- read_xlsx('/home/iganiemi/Documents/phd/st/data/geomx/geomx_batch2_1124/metadata/dcc_metadata_batch2_1124.xlsx')

meta_dir <- '/home/iganiemi/Documents/phd/st/data/geomx/geomx_batch3_0525/metadata'

meta_b123_path <- '/home/iganiemi/Documents/phd/st/data/geomx/batch123/metadata/dcc_metadata_batch123.xlsx'
meta_b123_out_path <- '/home/iganiemi/Documents/phd/st/data/geomx/batch123/metadata/dcc_metadata_batch123_cleaned.csv'
meta_b123_out_no_tls_path <- '/home/iganiemi/Documents/phd/st/data/geomx/batch123/metadata/dcc_metadata_batch123_no_tls_cleaned.csv'
clin_path <- '/home/iganiemi/Documents/phd/st/data/geomx/clinical_data/9_eyemt_patient_clinical_data.csv'

# combine worksheets
worksheets_list <- list.files(file.path(meta_dir, 'lab_worksheets'), full.names = T, pattern = "csv")

worksheet_all <- lapply(worksheets_list, function(x){
  batchname <- unlist(strsplit(basename(x), split = "_"))[4]
  worksheet <- fread(x)
  worksheet$batch_nr <- batchname 
  
  return(worksheet)
})

worksheet_all <- do.call(rbind, worksheet_all)

worksheet_all$Roi <- gsub('=', '', worksheet_all$Roi)
worksheet_all$Roi <- gsub('"', '', worksheet_all$Roi, fixed = T)

# combine roi selection files
roi_list <- list.files(file.path(meta_dir, 'roi_selection'), full.names = T, pattern = "csv")

roi_all <- lapply(roi_list, function(x){
  samplename <- paste0(unlist(strsplit(basename(x), split = "_"))[1], '_', unlist(strsplit(basename(x), split = "_"))[2])
  roi <- fread(x)
  roi$Sample <- samplename 
  colnames(roi) <- tolower(colnames(roi))
  
  if('segment' %in% colnames(roi)){
    roi <- dplyr::rename(roi, geomx_segment = segment)
  }
  
  if('geomx_' %in% colnames(roi)){
    roi <- dplyr::rename(roi, geomx_roi = geomx_)
  }
  
  return(roi)
})


roi_all <- do.call(rbind.fill, roi_all, )


# roi annotation table cleaning -------------------------------------------

# clean

# manual changes to the tables:
# S088 - empty column name manually changed into 'geomx_segment'
# S311 S112_iOme empty column with '!' and 'reload' manually changed to 'geomx_note'
# S118 - updated roi changed from 567 into 456
# S247_iOme + S081_iOme 'geomx_roi_run2' column with different roi numbers for some roi added _1 for 2nd run ROI
# S118_iOme + S195_iOme 'updated_commenr' column with different roi numbers for some roi added _1 for 2nd run ROI
# S195_iOme1 - file ROIs wo 'geomx_segment' label - Agg - changed to 'stroma' 
# S355 - 5 ROIs without labels into Agg in tls_status (taken from note)
# S080 + 4 other samples - additional Agg marked as 'Agg_additional'


# assuming that new roi numbers are correct ones
roi_all$geomx_roi_run2 <- ifelse(roi_all$geomx_roi_run2 == '', NA, roi_all$geomx_roi_run2)
roi_all$updated_comment <- ifelse(roi_all$updated_comment == '', NA, roi_all$updated_comment)
  
roi_all$geomx_roi_fixed <- ifelse(!is.na(roi_all$geomx_roi_run2), roi_all$geomx_roi_run2, roi_all$geomx_roi)
roi_all$geomx_roi_fixed <- ifelse(!is.na(roi_all$updated_comment), roi_all$updated_comment, roi_all$geomx_roi_fixed)


roi_geomx <- roi_all[!is.na(roi_all$geomx_roi_fixed), c('sample', 'geomx_roi_fixed', 
                                                        'geomx_segment', 'tls_id', 'tls_status',
                                                        'roi index', 'final_label', 'nk_status', 'note',
                                                        'start x', 'end x', 'start y', 'end y', 'width',
                                                        'height', 'label', 'manual_selection')]

# additional cleaning
roi_geomx$tls_status <- ifelse(!(roi_geomx$tls_status %in% c('Agg','Agg_additional', 'GC', 'S', 'TB')), NA,
                                roi_geomx$tls_status)

roi_geomx$sample_roi <- paste0(roi_geomx$sample, '_', roi_geomx$geomx_roi_fixed)
length(unique(roi_geomx$sample_roi)) # no duplicated values now

# rm GS with no tls_id
roi_geomx$tls_status <- ifelse(roi_geomx$tls_status == 'GC' & is.na(roi_geomx$tls_id), NA, roi_geomx$tls_status)

# rmv tls id for Agg
roi_geomx$tls_id <- ifelse(roi_geomx$tls_status == 'Agg', NA, roi_geomx$tls_id)

roi_geomx <- dplyr::rename(roi_geomx, tcycif_roi = `roi index`, geomx_roi = geomx_roi_fixed)

# clean final label to our nomenclature
clean_labs <- function(x, nk=F){
  y <- ifelse(grepl('CD4', x), 'CD4', '')
  y <- ifelse(grepl('CD8', x), paste(y, 'CD8', sep = '_'), y)
  y <- ifelse(grepl('CD11', x), paste(y, 'CD11', sep = '_'), y)
  y <- ifelse(grepl('Iba1', x), paste(y, 'Iba1', sep = '_'), y)
  if(nk){y <- ifelse(grepl('NK|hub', x), paste(y, 'NK', sep = '_'), y)}
  y <- gsub('^_', '', y)
  return(y)
}

roi_geomx$nk_status <- ifelse(is.na(roi_geomx$nk_status) | roi_geomx$nk_status == '', FALSE, roi_geomx$nk_status)

roi_geomx$final_label <- clean_labs(roi_geomx$final_label)
roi_geomx$final_label_nk <- ifelse(roi_geomx$nk_status == TRUE, paste0(roi_geomx$final_label, '_', 'NK'), 
                                   roi_geomx$final_label)


#change final label - without tls - eyemt_nolabel with tls - TLS
roi_geomx$final_label <- ifelse(roi_geomx$final_label == '' & is.na(roi_geomx$tls_status), 
                                'label_unknown', roi_geomx$final_label)

roi_geomx$final_label <- ifelse(roi_geomx$final_label == '' & !is.na(roi_geomx$tls_status), 
                                'TLS', roi_geomx$final_label)

roi_geomx$final_label_nk <- ifelse(roi_geomx$final_label == 'label_unknown', 
                                'label_unknown', roi_geomx$final_label_nk)

roi_geomx$final_label_nk <- ifelse(roi_geomx$final_label == 'TLS', 
                                'TLS', roi_geomx$final_label_nk)

fwrite(roi_geomx, file.path(meta_dir, 'geomx_batch3_0525_roi_selection_all_cleaned.csv'))

# clean worksheet file and select sample from slide -----------------------

# fix roi names from 2nd run for 2 slides (add _1)
worksheet_all$roi_clean <- ifelse(worksheet_all$`Scan Name` %in% c('S081_S247_130325_1', 'S118_S195_200325_1'), 
                             paste0(worksheet_all$Roi, '_1'), worksheet_all$Roi)

worksheet_all$roi_clean <- gsub('^0|^00', '', worksheet_all$roi_clean)
  
# divide slidename and fix samplenames
worksheet_all$Sample1 <- apply(worksheet_all, 1, function(x){
  
  s <- ifelse(x[['Slide Name']] == 'No Template Control', NA, 
              unlist(strsplit(x[['Slide Name']], split = "_"))[1])
  
  if(!is.na(s)){
    s <- paste0(s, '_', ifelse(grepl('PDS', x[['Slide Name']]), 'p', 'i'))
    s <- grep(s, unique(roi_geomx$sample), value = T)
  } else{ s <- NA}
  return(s)
})

worksheet_all$Sample2 <- apply(worksheet_all, 1, function(x){
  
  s <- ifelse(x[['Slide Name']] == 'No Template Control', NA, 
              unlist(strsplit(x[['Slide Name']], split = "_"))[2])
  s <- ifelse(s %in% c('PDS', 'IDS'), NA, s)
  
  if(!is.na(s)){
    s <- paste0(s, '_', ifelse(grepl('PDS', x[['Slide Name']]), 'p', 'i'))
    s <- grep(s, unique(roi_geomx$sample), value = T)
  } else{ s <- NA}
  return(s)
})


# make sample_roi column with roi matching
worksheet_all$Sample <- apply(worksheet_all, 1, function(x){
  # na handling

  if(paste0(x[['Sample1']], '_', x[['roi_clean']]) %in% roi_geomx$sample_roi){
  s <- x[['Sample1']]
  } else if (paste0(x[['Sample2']], '_', x[['roi_clean']]) %in% roi_geomx$sample_roi){
    s <- x[['Sample2']]
  } else{
    s <- NA
  }
  return(s)
})

worksheet_all$sample_roi <- ifelse(!is.na(worksheet_all$Sample), 
                                   paste0(worksheet_all$Sample, '_', worksheet_all$roi_clean), NA)

worksheet_all$sample_roi_aoi <- ifelse(!is.na(worksheet_all$Sample), 
                                       paste0(worksheet_all$sample_roi, '_', worksheet_all$Segment), NA)

which(duplicated(worksheet_all$sample_roi_aoi)) # only NTC and ROI25 in S081_2247 

# remove roi 25 from S081_2247 which was collected with an error
# DSP-1001660039810-H-D03 - actually it is empty as well
worksheet_all <- worksheet_all[!(worksheet_all$`Scan Name` == 'S081_S247_130325' & worksheet_all$Roi == '025'), ]

ll <- worksheet_all[which(duplicated(worksheet_all$sample_roi_aoi)), ]

setdiff(roi_geomx$sample_roi, worksheet_all$sample_roi)
setdiff(worksheet_all$sample_roi, roi_geomx$sample_roi)

#  clean segment in worksheet - S118_S195 'Full ROI' = stroma
worksheet_all$Segment <- ifelse(worksheet_all$Segment == 'Full ROI', 'stroma', worksheet_all$Segment)
worksheet_all$Segment <- ifelse(worksheet_all$Segment == 'Stroma', 'stroma', worksheet_all$Segment)
worksheet_all$Segment <- ifelse(worksheet_all$Segment == 'Tumor', 'tumor', worksheet_all$Segment)

worksheet_all <- worksheet_all[, -c('Sample1', 'Sample2')]

fwrite(worksheet_all, file.path(meta_dir, 'geomx_batch3_0525_worksheet_all_cleaned.csv'))

# combine worksheet and roi table -----------------------------------------

# match DCC per sample+roi
worksheet_roi_all <- left_join(worksheet_all[, -c('sample_roi_aoi')], roi_geomx, by = 'sample_roi')

# adjust colnames from roi to template_meta
worksheet_roi_all <- worksheet_roi_all[, -c('sample', 'roi_clean', 'sample_roi', 'nk_status')]
worksheet_roi_all$dcc_filename <- paste0(worksheet_roi_all$Sample_ID, '.dcc')
worksheet_roi_all$main_batch_nr <- 3
worksheet_roi_all$NACT_status <- ifelse(grepl('p', worksheet_roi_all$Sample), 'pre', 
                                        ifelse(grepl('i', worksheet_roi_all$Sample), 'post', NA))
worksheet_roi_all$Patient <- gsub('_.*', '', worksheet_all$Sample)
worksheet_roi_all$Site <- ifelse(grepl('Ome', worksheet_roi_all$Sample), 'Omentum', 
                                 ifelse(grepl('Per', worksheet_roi_all$Sample), 'Peritroneum', 
                                        ifelse(grepl('Adn', worksheet_roi_all$Sample), 'Adnex',
                                               ifelse(grepl('Ova', worksheet_roi_all$Sample), 'Ovary', NA))))
worksheet_roi_all$Clinical_Notes <- NA
worksheet_roi_all$Correct_label <- NA
worksheet_roi_all$batch_nr_sample_collection <- NA
worksheet_roi_all$Annotation_cell_first_labels <- NA
worksheet_roi_all$PanCK_positive <- NA 



worksheet_roi_all <- dplyr::rename(worksheet_roi_all, Roi_original = Roi, Roi = geomx_roi, 
                                   ROI_Index_tCycIF = tcycif_roi, Annotation_cell = final_label, 
                                   Annotation_cell_nk = final_label_nk, Segment_geomx = geomx_segment,
                                   Label = label, Start_X = `start x`, End_X = `end x`, 
                                   Start_Y = `start y`, End_Y = `end y`, Note_tCycIF = note, Width = width,
                                   Height = height, Manual_selection = manual_selection)

colnames(worksheet_roi_all) <- gsub(' ', '_', colnames(worksheet_roi_all))

template_cols <- colnames(template_meta)[-c(1,44)]

worksheet_roi_all <- as.data.frame(worksheet_roi_all)[, c(template_cols, "Roi_original", 
                             "tls_id", "tls_status", "Annotation_cell_nk")]

fwrite(worksheet_roi_all, file.path(meta_dir, 'dcc_metadata_all_batch3_0525.csv'))


# TODO add negative probe name



# examine an unify for all batches df -------------------------------------

b1 <- read_excel('/home/iganiemi/Documents/phd/st/data/geomx/geomx_batch1_0823/metadata/dcc_metadata_batch1_0823.xlsx')
b2 <- read_excel('/home/iganiemi/Documents/phd/st/data/geomx/geomx_batch2_1124/metadata/dcc_metadata_batch2_1124.xlsx')
b3 <- read_excel('/home/iganiemi/Documents/phd/st/data/geomx/geomx_batch3_0525/metadata/dcc_metadata_all_batch3_0525.xlsx')
# load all batches df
meta <- as.data.frame(read_xlsx(meta_b123_path))

meta$Roi_geomx <- meta$Roi # just to have it saved while loading geomx_obj
meta$ROI_Coordinate_Y <- ifelse(is.na(meta$ROI_Coordinate_Y), meta$ROI.Coordinate.Y, meta$ROI_Coordinate_Y) #b1 has wrong name

# fix nk annotations for b2
meta$Annotation_cell_nk <- ifelse(meta$main_batch_nr != 3, meta$Annotation_cell, meta$Annotation_cell_nk)
meta$Annotation_cell_nk <- ifelse(meta$main_batch_nr == 2 & grepl('NK', meta$Correct_label), paste0(meta$Annotation_cell_nk, '_NK'), meta$Annotation_cell_nk)

meta <- dplyr::rename(meta, tCycIF_preselection_ROI_Start_X = Start_X, tCycIF_preselection_ROI_End_X = End_X, 
                      tCycIF_preselection_ROI_Start_Y = Start_Y, tCycIF_preselection_ROI_End_Y = End_Y, 
                      tCycIF_preselection_ROI_Height = Height, tCycIF_preselection_ROI_Width = Width,
                      tCycIF_preselection_ROI_Index = ROI_Index_tCycIF,
                      ROI_Coordinate_X_geomx = ROI_Coordinate_X, ROI_Coordinate_Y_geomx = ROI_Coordinate_Y,
                      tCycIF_preselection_Note = Note_tCycIF, tCycIF_preselection_Manual_selection = Manual_selection,
                      tCycIF_preselection_initial_label = Label, Roi_geomx_original = Roi_original) %>%
  dplyr::select(!c(CD8, IBA1, Annotation, Tags, ROI.Coordinate.Y, Notes,
                   Clinical_Notes, PanCK_positive, Correct_label, Annotation_cell_first_labels, HRP_status,
                   BRCA_status, PFS_quartile_b123, OS_quartile_b123))

# add clinical data
clin <- fread(clin_path) %>%
  filter(Patient %in% meta$Patient) %>%
  select(Patient, Sample, NACT_status, HRP_status, BRCA_status, PFS_quartile_b123, OS_quartile_b123)

# add paired status
paired <- distinct(clin, Patient, NACT_status) %>%
  group_by(Patient) %>%
  summarise(n = n()) %>%
  filter(n == 2)

meta$paired_status <- ifelse(meta$Patient %in% paired$Patient, 'paired', 'unpaired')

# merge once again with clinical data (some NAs previously..)
clin_per_pt <- clin %>% 
  select(Patient, HRP_status, BRCA_status, PFS_quartile_b123, OS_quartile_b123) %>%
  distinct()

meta <- left_join(meta, clin)

# calculate additional PFS/OS quartiles

meta$PFS_median_b123 <- ifelse(meta$PFS_quartile_b123 %in% c(1, 2), 1, 2)
meta$OS_median_b123 <- ifelse(meta$OS_quartile_b123 %in% c(1, 2), 1, 2)

meta$PFS_quartile_paired <- mapvalues(meta$Patient, from=c("S015", "S027", "S032", "S069", "S084", "S139",
                                                            "S229", "S333"), 
                                                     to=c(3, 2, 4, 2, 1, 1, 4, 3))

meta$PFS_quartile_paired <- ifelse(meta$paired_status == 'paired', meta$PFS_quartile_paired, NA)

meta$OS_quartile_paired <- mapvalues(meta$Patient, from=c("S015", "S027", "S032", "S069", "S084", "S139",
                                                          "S229", "S333"), 
                                     to=c(4, 3, 4, 2, 1, 1, 3, 2))

meta$OS_quartile_paired <- ifelse(meta$paired_status == 'paired', meta$OS_quartile_paired, NA)

meta$PFS_median_paired <- ifelse(meta$PFS_quartile_paired %in% c(1, 2), 1,
                                                ifelse(meta$PFS_quartile_paired %in% c(3, 4), 2, NA))

meta$OS_median_paired <- ifelse(meta$OS_quartile_paired %in% c(1, 2), 1,
                                               ifelse(meta$OS_quartile_paired %in% c(3, 4), 2, NA))

fwrite(meta, meta_b123_out_path)

# remove TLS
meta_no_tls <- meta[!(meta$tls_status %in% c('GC', 'S', 'TB')), ]
fwrite(meta_no_tls, meta_b123_out_no_tls_path)

# check empty DCC files ---------------------------------------------------

b31 <- file_size(filepath = '/home/iganiemi/Documents/phd/st/data/geomx/geomx_batch3_0525/dcc/Batch3_1')
b32 <- file_size(filepath = '/home/iganiemi/Documents/phd/st/data/geomx/geomx_batch3_0525/dcc/Batch3_2')
b33 <- file_size(filepath = '/home/iganiemi/Documents/phd/st/data/geomx/geomx_batch3_0525/dcc/Batch3_3')

b3_size <- rbind(b31, b32, b33)
b3_size <- b3_size[!(grepl('zip', b3_size$name)), ]
b3_size$empty <- ifelse(b3_size$size == '354 B', TRUE, FALSE)
duplicates <- b3_size$name[which(duplicated(b3_size$name))]

# check if all empty has a duplicate
# 1 file empty wo duplicate DSP-1001660039810-H-D03.dcc
# 1 file almost empty wo duplicate DSP-1001660039812-B-A01.dcc
table(b3_size$empty)
length(which(b3_size$name[b3_size$empty == T] %in% duplicates))

######
meta_all <- fread(meta_b123_out_path)
meta_all <- meta_all[, c('dcc_filename','Slide_Name', 'Scan_Name', 'Roi', 'Segment', 'Sample', 'main_batch_nr', 'batch_nr')]
meta_all3 <- meta_all[meta_all$main_batch_nr == 3, ]
meta_all3 <- meta_all3[meta_all3$Sample != '', ]

for(sample in unique(meta_all3$Sample)){
  print(sample)
  print(unique(meta_all3$Roi[meta_all3$Sample == sample]))
  print('$$$$$$$$$$')
}

