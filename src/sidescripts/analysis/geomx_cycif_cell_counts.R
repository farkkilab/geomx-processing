library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)



# load variables ----------------------------------------------------------

data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'
output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2'

cell_counts_dir <- file.path(data_dir, 'image_cell_counts') 
cell_counts_dir <- '/home/ad/P-drive/h345/afarkkilab/Data/9-tCycIF-GeoMx-PreandPost/data/geomx/quantification/'

geomx_path <- file.path(output_dir, 'geomx_qc_norm.RDS')
prism_deconv_path <- file.path(output_dir, 'prism','results', 'geomx_weights_low_lvl_ct.tsv')


#################

geomx_obj <- readRDS(geomx_path)
geomx_meta <- sData(geomx_obj)

cell_counts <- lapply(list.files(cell_counts_dir, full.names = T), fread)
cell_counts <- do.call(rbind, cell_counts)

prism_deconv <- fread(prism_deconv_path)

prism_deconv <- left_join(prism_deconv, geomx_meta[, c('dcc_filename', 'Segment', 'Patient', 'NACT status', 'Annotation_cell')],
                          by = c('V1' = 'dcc_filename'))

########################################


geomx_meta$Sample_ROI <- paste0(geomx_meta$Patient, '-', geomx_meta$Roi)

# remove wrong ROIs
# TODO it's freaking madness. 
set1 <- c('001', '002', '003', '004', '005', '006', '007', '008', '009', '010')
set2 <- c('011', '012', '013', '014', '015', '016', '017', '018', '019', '020')

cell_counts <- filter(cell_counts, !(`ROI-Sample` %in% paste0("S032_Post-", set1)))
cell_counts$`ROI-Sample` <- gsub('_Post|_Pre', '', cell_counts$`ROI-Sample`)

cell_counts$`ROI-Sample` <- ifelse(cell_counts$`ROI-Sample` %in% paste0("S053_S076-", c(set1, '021')),
                                   gsub('_S076', '', cell_counts$`ROI-Sample`), cell_counts$`ROI-Sample`)
cell_counts$`ROI-Sample` <- ifelse(cell_counts$`ROI-Sample` %in% paste0("S053_S076-", set2),
                                   gsub('S053_', '', cell_counts$`ROI-Sample`), cell_counts$`ROI-Sample`)

cell_counts$`ROI-Sample` <- ifelse(cell_counts$`ROI-Sample` %in% paste0("S057_S072-", set1),
                                   gsub('_S072', '', cell_counts$`ROI-Sample`), cell_counts$`ROI-Sample`)
cell_counts$`ROI-Sample` <- ifelse(cell_counts$`ROI-Sample` %in% paste0("S057_S072-", set2),
                                   gsub('S057_', '', cell_counts$`ROI-Sample`), cell_counts$`ROI-Sample`)

cell_counts$`ROI-Sample` <- ifelse(cell_counts$`ROI-Sample` %in% paste0("S065_S073-", set1),
                                   gsub('S065_', '', cell_counts$`ROI-Sample`), cell_counts$`ROI-Sample`)
cell_counts$`ROI-Sample` <- ifelse(cell_counts$`ROI-Sample` %in% paste0("S065_S073-", set2),
                                   gsub('_S073', '', cell_counts$`ROI-Sample`), cell_counts$`ROI-Sample`)


sort(unique(geomx_meta$Sample_ROI))
sort(unique(cell_counts$`ROI-Sample`))

setdiff(geomx_meta$Sample_ROI, cell_counts$`ROI-Sample`)
setdiff(cell_counts$`ROI-Sample`, geomx_meta$Sample_ROI) #S084-006 got lost with geomx QC



sapply(unique(cell_counts$`ROI-Sample`), function(x){
  a <- ggplot(cell_counts[cell_counts$`ROI-Sample` == x, ]) +
    geom_histogram(aes(x = CD45), bins = 200)
  
  b <- ggplot(cell_counts[cell_counts$`ROI-Sample` == x, ]) +
    geom_histogram(aes(x = CD45), bins = 200) +
    xlim(0, 1000)
  
  ggsave(
    filename = file.path(output_dir, 'cell_count', paste0(x, '.pdf')), 
    plot = marrangeGrob(list(a, b), nrow=1, ncol=2), 
    width = 15, height = 9
  )
})

x <- 'S015-001'

ggplot(cell_counts[cell_counts$`ROI-Sample` == x, ]) +
  geom_histogram(aes(x = CD45), bins = 200) +
  xlim(0, 10000)




kk <- data.frame('Sample-ROI' = unique(cell_counts$`ROI-Sample`), 'CD45_thr' = NA)
fwrite(kk, file.path(output_dir, 'cell_count', paste0('cd45_thr.csv')))
########################################

ggplot(prism_deconv, aes(x = Segment, y = immune)) +
  geom_violin() +
  geom_point( size= 0.2, alpha = 0.6, position = 'jitter') +

ggplot(prism_deconv, aes(x = Segment, y = tumor)) +
  geom_violin() +
  geom_point(size= 0.2, alpha = 0.6, position = 'jitter') +

ggplot(prism_deconv, aes(x = Segment, y = stroma)) +
  geom_violin() +
  geom_point(size= 0.2, alpha = 0.6, position = 'jitter') 


############################################




##########################