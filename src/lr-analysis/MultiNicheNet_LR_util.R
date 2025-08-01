# MultiNicheNet util

library(readr)
library(stringr)
library(dplyr, quietly =T)
# preparing externally provided DEGS


# Load all DEG CSVs into a named list
#deg_files <- list.files(file.path(data_dir, "DEGs"), pattern = "\\.csv$", full.names = TRUE)

#deg_list <- lapply(deg_files, read.csv)
#names(deg_list) <- c("Macrophages", "Tcells")  # check the file name and manually assign cell type names



df_Tcells = read.csv(file.path(data_dir,'DEGs','dge_deconv_Tcells_dge_within_slide_Segment_bin_FALSE__.csv'))
df_Macrophages = read.csv(file.path(data_dir,'DEGs','dge_deconv_Macrophages_dge_within_slide_Segment_bin_FALSE__.csv'))


deg_list = list(Tcells = df_Tcells ,
                Macrophages = df_Macrophages
                )


# Processing function for one DEG dataframe
process_deg_file <- function(df, celltype) {
  df1 <- data.frame(
    gene = df$Gene,
    cluster_id = celltype,
    logFC = df$Estimate,
    p_val = df$`Pr...t..`,
    p_adj = df$FDR,
    contrast = paste(trimws(strsplit(unique(df$Contrast), split = "-")[[1]]), collapse = "-")
  )

  contrast2 = paste(rev(trimws(strsplit(unique(df$Contrast), split = "-")[[1]])), collapse = "-")



  df2 <- df1
  df2$logFC <- -df2$logFC
  df2$contrast <- contrast2

  return(rbind(df1, df2))
}

# Apply the function to all items in the list
celltype_de_combined <- do.call(rbind, Map(
  process_deg_file,
  df = deg_list,
  celltype = names(deg_list)
))



rownames(celltype_de_combined) <- NULL

saveRDS(celltype_de_combined, file.path(output_dir,MultiNicheNet_folder_name,"celltype_de_combined_calculated_externally.RDS"))
        
        
        
        