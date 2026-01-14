# integrates geomx and cycif roi coordinates, updates cleaned metadata

library(data.table)
library(plyr)
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(ggplot2)
library(tibble)
library(GeomxTools)
library(readxl)
library(tools)

# TODO make sure X px and Y px for b1 are for c1 (top-left corner) - NO!!!!
# TODO calculate 2 missing ROIs: 45,46 of S106 iOme (not in input frames)
# c1: top-left, c2: top-right, c3: down-right, c4: down-left


# S080_iOme2, S091_iOme1, S195_iOme1 - corrected in the b3 filenames
# S225 + S309 roi naming was changed to _1 in input files to match metadata

# define paths ------------------------------------------------------------
# from master script
batch <- 'batch123'
proj_dir <<- '~/Documents/phd/st'
data_dir <<- '~/Documents/phd/st/data/geomx/batch123/' # batch1 2 and 3
anno_path <<- file.path(data_dir, 'metadata', 'dcc_metadata_batch123_cleaned.csv') #batch1 and 2 and 3
anno_path_notls <<- file.path(data_dir, 'metadata', 'dcc_metadata_batch123_no_tls_cleaned.csv')
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808') # batch123

# coordinates for cycif images pathways
eyemt_geomx_dir <- "/home/ad/P-drive/h30492/farkkilab2/9_EyeMT/Data/geomx"
coords_dir <- "roi_coordinates_cycif/cycif_roi_arrays" # coors dir within batch data dir
coords_geomx_dir <- "roi_coordinates_cycif/roi_coordinates_geomx" # coors dir within batch data dir

dir.create(file.path(output_dir, "cycif_integration"))

output_coords_b123_path <- file.path(output_dir, 'cycif_integration', 'batch123_geomx_cycif_coordinates.csv')

# load full cleaned metadata ----------------------------------------------

meta_cleaned <- as.data.frame(fread(anno_path))
meta_cleaned$Roi_geomx_original <- as.character(meta_cleaned$Roi_geomx_original)
meta_cleaned$main_batch_nr <- as.character(meta_cleaned$main_batch_nr)
meta_cleaned$sample_roi <- paste0(meta_cleaned$Sample, '_', meta_cleaned$Roi_geomx)

# clean and calculate cycif coordinates for b2 and 3 ----------------------

col_order <- c("main_batch_nr", "Sample", "Roi_geomx", "sample_roi", "roi_name",
               "roi_c1_X_geomx", "roi_c1_Y_geomx", "roi_width_geomx", "roi_height_geomx",
               "roi_c1_X_cycif", "roi_c1_Y_cycif", "roi_c2_X_cycif", "roi_c2_Y_cycif",
               "roi_c3_X_cycif", "roi_c3_Y_cycif", "roi_c4_X_cycif", "roi_c4_Y_cycif",
               "roi_width_cycif", "roi_height_cycif", "roi_center_X_cycif", "roi_center_Y_cycif")

# b2 and b3 ROIs are not rectangles but polygons.
# width+height are calculated as longer one from 2 possible height/width edges

for(batch_name in c('batch2', 'batch3')){
  roi_coords_dir <- file.path(eyemt_geomx_dir, batch_name, coords_dir)
  roi_coords_geomx_dir <- file.path(eyemt_geomx_dir, batch_name, coords_geomx_dir)
  sample_names <- file_path_sans_ext(basename(list.files(roi_coords_dir, pattern = ".csv")))
  
  roi_coords_all <- lapply(sample_names, function(sample_name){
    # clean roi coords df
    roi_coords <- as.data.frame(fread(file.path(roi_coords_dir, paste0(sample_name, '.csv')), drop = c(6, 11)))
    colnames(roi_coords) <- c('roi_name', 'c1_X', 'c2_X', 'c3_X', 'c4_X', 'c1_Y', 'c2_Y', 'c3_Y', 'c4_Y')
    roi_coords$c1_X <- gsub('\\[|\\]', '', roi_coords$c1_X)
    roi_coords$c1_Y <- gsub('\\[|\\]', '', roi_coords$c1_Y)
    roi_coords <- mutate_at(roi_coords, vars(matches("c[0-9]")), function(x){gsub(" ", "", x)})
    roi_coords <- mutate_at(roi_coords, vars(matches("c[0-9]")), as.numeric)
    roi_coords$roi_name <- as.character(roi_coords$roi_name)
    
    #calculate additional dimentions
    roi_coords <-  apply(roi_coords, 1, function(row){
      row$roi_width <- round(as.numeric(max(as.numeric(row[['c2_X']]), as.numeric(row[['c3_X']])) - min(as.numeric(row[['c1_X']]), as.numeric(row[['c4_X']]))),3)
      row$roi_height <- round(as.numeric(max(as.numeric(row[['c1_Y']]), as.numeric(row[['c2_Y']])) - min(as.numeric(row[['c3_Y']]), as.numeric(row[['c4_Y']]))),3)
      row$roi_center_X <- round(as.numeric(min(as.numeric(row[['c1_X']]), as.numeric(row[['c4_X']])) + (row$roi_width * 0.5)),3)
      row$roi_center_Y <- round(as.numeric(min(as.numeric(row[['c4_Y']]), as.numeric(row[['c3_Y']])) + (row$roi_height * 0.5)),3)
      return(as.data.frame(row))
    })
    
    roi_coords <- as.data.frame(do.call(rbind, roi_coords))
    roi_coords <- mutate_at(roi_coords, vars(matches("c[0-9]")), as.numeric)
    roi_coords$Sample <- sample_name
    
    # add original geomx roi coords
    if(batch_name == 'batch2'){
      coords_path <- file.path(roi_coords_geomx_dir, paste0(sample_name, '_GeoMx_ROIs.csv'))
    } else if(batch_name == 'batch3'){
      coords_path <- file.path(roi_coords_geomx_dir, paste0(sample_name, '.csv'))
    }

    roi_coords_geomx <- as.data.frame(fread(coords_path, select = c('Roi', 'ROI_Coordinate_X', 'ROI_Coordinate_Y', 'Width', 'Height'), 
                                            col.names = c('roi_name', "roi_c1_X_geomx", "roi_c1_Y_geomx", 
                                                          "roi_width_geomx", "roi_height_geomx")))

    roi_coords_geomx$roi_name <- as.character(roi_coords_geomx$roi_name)
    roi_coords <- left_join(roi_coords, roi_coords_geomx)


    return(roi_coords)
  })
  
  print('merged')
  
  roi_coords_all <- do.call(rbind, roi_coords_all)
  roi_coords_all$roi_name <- gsub(" ", "", roi_coords_all$roi_name)
  roi_coords_all$main_batch_nr <- gsub("batch", "", batch_name)
  
  # merge with other metadata to ensure roi_geomx naming
  if(batch_name == 'batch2'){
    roi_coords_all <- left_join(roi_coords_all, meta_cleaned[, c('Roi_geomx_original', 'Sample', 'main_batch_nr', 'Roi_geomx')], 
                                by = c('roi_name' = 'Roi_geomx_original', 'Sample', 'main_batch_nr')) %>%
      distinct()
  } else if(batch_name == 'batch3'){
    roi_coords_all$Roi_geomx <- paste0('roi-', roi_coords_all$roi_name)
  }
  
  roi_coords_all$sample_roi <- paste0(roi_coords_all$Sample, '_', roi_coords_all$Roi_geomx) 
  colnames(roi_coords_all)[2:9] <- paste0('roi_', colnames(roi_coords_all)[2:9], '_cycif')
  colnames(roi_coords_all)[10:13] <- paste0(colnames(roi_coords_all)[10:13], '_cycif')
  roi_coords_all <- roi_coords_all[, col_order] # reorder

  # double check if names are identical with the one in metadata
  stopifnot(all(roi_coords_all$sample_roi %in% 
                  meta_cleaned$sample_roi[meta_cleaned$main_batch_nr == gsub("batch", "", batch_name)]))

  fwrite(roi_coords_all, file.path(output_dir, "cycif_integration", paste0(batch_name, "_geomx_cycif_coordinates.csv")))
}


# combine cycif coordinates for b1 ----------------------------------------

# b1 are just rectangles calculated from c1 coords + width/height

batch_name <- 'batch1'

roi_coords_dir <- file.path(eyemt_geomx_dir, batch_name, coords_dir)
sample_names <- file_path_sans_ext(basename(list.files(roi_coords_dir, pattern = ".csv")))

roi_coords_all <- lapply(sample_names, function(sample_name){
  
  # read correct columns and rename
  roi_coords <- as.data.frame(fread(file.path(roi_coords_dir, paste0(sample_name, '.csv')), select = c(1, 7:10, 2:5),
                                    col.names = c("roi_name", "roi_c1_X_cycif", "roi_c1_Y_cycif",
                                                  "roi_width_cycif", "roi_height_cycif", 
                                                  "roi_c1_X_geomx", "roi_c1_Y_geomx",
                                                  "roi_width_geomx", "roi_height_geomx")))
  
  # calculate other corners
  roi_coords$roi_c2_X_cycif <- roi_coords$roi_c1_X_cycif + roi_coords$roi_width_cycif
  roi_coords$roi_c2_Y_cycif <- roi_coords$roi_c1_Y_cycif
  roi_coords$roi_c3_X_cycif <- roi_coords$roi_c2_X_cycif
  roi_coords$roi_c3_Y_cycif <- roi_coords$roi_c1_Y_cycif - roi_coords$roi_height_cycif
  roi_coords$roi_c4_X_cycif <- roi_coords$roi_c3_X_cycif - roi_coords$roi_width_cycif
  roi_coords$roi_c4_Y_cycif <- roi_coords$roi_c3_Y_cycif
  
  roi_coords$roi_center_X_cycif <- roi_coords$roi_c1_X_cycif + (roi_coords$roi_width_cycif * 0.5)
  roi_coords$roi_center_Y_cycif <- roi_coords$roi_c1_Y_cycif - (roi_coords$roi_height_cycif * 0.5)

  roi_coords$roi_name <- gsub("ROI", "", roi_coords$roi_name)
  roi_coords$Patient <- sample_name
  
  return(roi_coords)
})

roi_coords_all <- do.call(rbind, roi_coords_all)
roi_coords_all$roi_name <- gsub(" ", "", roi_coords_all$roi_name)
roi_coords_all$main_batch_nr <- gsub("batch", "", batch_name)

# merge patient+roi with metadata to find sample name
roi_coords_all <- left_join(roi_coords_all, meta_cleaned[, c('Roi_geomx_original', 'Patient', 'Sample', 'main_batch_nr', 'Roi_geomx')], 
                            by = c('roi_name' = 'Roi_geomx_original', 'Patient', 'main_batch_nr')) %>%
  distinct()

# unique identifier of roi
roi_coords_all$sample_roi <- paste0(roi_coords_all$Sample, '_', roi_coords_all$Roi_geomx) 

#reorder and save
roi_coords_all <- roi_coords_all[, col_order]
fwrite(roi_coords_all, file.path(output_dir, "cycif_integration", paste0(batch_name, "_geomx_cycif_coordinates.csv")))

# merge all files together and merge with cleaned metadata ----------------

roi_coords_b123 <- lapply(c("batch1", "batch2", "batch3"), function(batch_name){
  roi_coords_b <- fread(file.path(output_dir, "cycif_integration", paste0(batch_name, "_geomx_cycif_coordinates.csv")))
})

roi_coords_b123 <- do.call(rbind, roi_coords_b123)
roi_coords_b123$main_batch_nr <- as.integer(roi_coords_b123$main_batch_nr)
fwrite(roi_coords_b123, output_coords_b123_path)

# join with metadata
meta <- fread(anno_path)
meta <- left_join(meta, roi_coords_b123[, c(1, 2, 3, 6:21)], by = c('main_batch_nr','Sample', 'Roi_geomx'))
meta <- meta[, c(1:24, 55:70, 25:54)]

# 2 ROIs have not been computed 45,46 of S106 iOme (not in input frames)
meta_empty <- meta[is.na(meta$roi_c1_X_geomx), ]

# examine non-matching geomx coords
# batch1 coords are not at all in line..
non_matching <- which(round(as.numeric(meta$ROI_Coordinate_X_geomx), 3) != round(as.numeric(meta$roi_c1_X_geomx), 3))
meta_non_matching <- meta[non_matching, ]
