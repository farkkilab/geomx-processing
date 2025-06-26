# MultiNicheNet util


# preparing externally provided DEGS

# Provide the DEGS seperetely for each cell type

library(dplyr)
library(readr)
library(stringr)


comparison <- c("tumor","stroma") 

# Create a contrast table
contrast_tbl <- tibble(contrast = c(paste(comparison[1], comparison[2], sep = "-"),paste(comparison[2], comparison[1], sep = "-")), 
                       group = c(comparison[1], comparison[2]))


# Load all DEG CSVs into a named list
deg_files <- list.files(file.path(data_dir, "DEGs"), pattern = "\\.csv$", full.names = TRUE)

deg_list <- lapply(deg_files, read.csv)
names(deg_list) <- c("Macrophages", "Tcells")  # check the file name and manually assign cell type names


# Processing function for one DEG dataframe
process_deg_file <- function(df, celltype, contrast1, contrast2) {
  df1 <- data.frame(
    gene = df$Gene,
    cluster_id = celltype,
    logFC = df$Estimate,
    p_val = df$`Pr...t..`,
    p_adj = df$FDR,
    contrast = contrast1
  )
  
  df2 <- df1
  df2$logFC <- -df2$logFC
  df2$contrast <- contrast2
  
  rbind(df1, df2)
}

# Apply the function to all items in the list
celltype_de_combined <- do.call(rbind, Map(
  process_deg_file,
  df = deg_list,
  celltype = names(deg_list),
  MoreArgs = list(
    contrast1 = contrast_tbl$contrast[1],
    contrast2 = contrast_tbl$contrast[2]
  )
))

rownames(celltype_de_combined) <- NULL

saveRDS(celltype_de_combined, file.path(output_dir,MultiNicheNet_folder_name,"celltype_de_combined_calculated_externally.RDS"))
        
        
        
        