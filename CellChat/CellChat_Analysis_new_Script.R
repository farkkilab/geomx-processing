# Ligand Receptor Analysis  by CellChat : Geomx Deconvoluted data


# load libraries

library(Seurat)
library(SeuratObject)
library(DESeq2)
library(CellChat)
library(patchwork)
#library(gridExtra)


# file names

input_dir <- 'C:/Users/Sahas/Downloads/Masters_Thesis/Data/Batch01/'
output_dir <- 'C:/Users/Sahas/Downloads/Masters_Thesis/Ligand-receptor/TestRun_New_scripts/'



# parameters for Cellchat

cell_types = c("Tcells","Macrophages") # set to NULL to get all the cell types  "Bcells","DCs"
cell_frac_cutoff = 0.01 # 0.01 or 0.005
min_cells = 10
# row_variance_cutoff = NULL I did not filter based on rwo variance

nact_status = NA
segment = c("stroma","tumor") # stroma , tumor
annotation = NA




# load data


geomx_obj = readRDS(paste0(input_dir,"geomx_qc_norm_batch_eff_rm.RDS"))
bprism_res = readRDS(paste0(input_dir,"bp_res_mid_lvl_ct.RDS"))
cell_fractions_df = read.csv(paste0(input_dir,"bp_res_mid_lvl_ct_ct_fraction.csv"))


# create directory for output files

output_folder_name  = "Cellchat_objects"

if (!dir.exists(paste0(output_dir,output_folder_name))) {
  new_dir = paste0(output_dir,output_folder_name)
  dir.create(new_dir,recursive = TRUE)
  message("Directory created")
  output_dir = new_dir
} else {
  output_dir = paste0(output_dir,output_folder_name)
}




# preprocessing and generating log normalized expr matrix


###################  Extracting expression data of the desired cell types  ################ 


if (is.null(cell_types)) { # If the cell types are not defined take all the cell types in the prism object
  ct_names <- colnames(bprism_res@posterior.theta_f@theta.cv)
  
}else{
  ct_names <- cell_types
}



deconv_ct_list <- setNames(lapply(ct_names, function(ct_name) {
  cell_frac_cv <- as.data.frame(bprism_res@posterior.theta_f@theta.cv)
  
  # Mask cell fractions if CV > 0.2
  cell_to_rm <- rownames(cell_frac_cv)[cell_frac_cv[[ct_name]] > 0.2]
  
  deconv_ct <- BayesPrism::get.exp(bp = bprism_res,
                                   state.or.type = "type",
                                   cell.name = ct_name)
  
  # Filter out cells with high CV
  deconv_ct <- deconv_ct[!(rownames(deconv_ct) %in% cell_to_rm), ]
  
  return(deconv_ct)
}), ct_names) 


####################     combine expression data        ####################

deconv_ct_list_int <- lapply(names(deconv_ct_list), function(ct_name){
  mat <- deconv_ct_list[[ct_name]]
  rownames(mat) <- gsub('\\.dcc', paste0('_', ct_name), rownames(mat))
  
  mat <- apply(mat, c(1, 2), function(x) {(as.integer(x))})
  
  return(mat)
})

names(deconv_ct_list_int) <- names(deconv_ct_list)
deconv_ct_list_int_all <- do.call(rbind, deconv_ct_list_int)



####################  combine meta data ########################


deconv_metadt_list <- lapply(names(deconv_ct_list_int), function(ct_name){
  meta_data = sData(geomx_obj)
  
  new_meta <- data.frame(meta_data[,c(2,5,6,7,24,25,28)])
  
  new_meta$dcc_filename = gsub('\\.dcc', paste0('_', ct_name), new_meta$dcc_filename)
  mat = deconv_ct_list_int[[ct_name]]
  
  matched_entries = new_meta[new_meta$dcc_filename %in% rownames(mat), ]
  rownames(matched_entries) = matched_entries$dcc_filename  
  
  return(matched_entries)
})

metadt_all <- do.call(rbind, deconv_metadt_list)



###############################     normalization           ################



# Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = t(deconv_ct_list_int_all),
                              colData = metadt_all,
                              design = formula(~ Segment + NACT_status))

# TODO examine eg if add Annotation_cell or NACT status?

dds <- estimateSizeFactors(dds, type = "poscounts")

#dds <- estimateSizeFactors(dds, type = "iterate")

#dds <- estimateSizeFactors(dds) # did not work

deseq2_norm_counts <- counts(dds, normalized=TRUE)
expr_deseq2_norm_log <- log2(deseq2_norm_counts + 1) # log transformation



################   filtering based on cell fraction #####################



formatted_ct_names <- gsub(" ", ".", ct_names) #because the cell fraction dataframe colnames are different: dot instead of space

filtered_sample_list <- setNames(lapply(formatted_ct_names, function(ct_name) {
  
  # Filter out cells less than the cell_frac_cutoff
  if(!is.null(cell_frac_cutoff)){
    
    cells_greater_than_cutoff <- cell_fractions_df[cell_fractions_df[, ct_name] >= cell_frac_cutoff,]$dcc_filename
    
    cells_greater_than_cutoff <- gsub('\\.dcc', paste0('_', ct_name), cells_greater_than_cutoff)
    
    return(cells_greater_than_cutoff)
    
  }
  
  return(NULL)
  
}), ct_names) 

# Get a unique list of filtered samples
filtered_samples <- unlist(filtered_sample_list)

# Subset the expression matrix
expr_deseq2_norm_log <- expr_deseq2_norm_log[, colnames(expr_deseq2_norm_log) %in% filtered_samples]





# cell type abundance plot

plot = cell_type_abundance_plot(metadt_all,min_cells)

pdf(file.path(output_dir, "cell_chat_cell_type_abundance.pdf"), width = 10, height = 6)
plot
dev.off()

# TODO

# if (length(low_count_labels) > 0) {
#   warning("The following cell types have fewer than 10 cells: ", paste(low_count_labels, collapse = ", "))
#   stop("Stopping execution due to low cell count.check the cell type abundance plots and set the min_cells counts")
# }

# calculating the cellchat probabiities

cellchat_results <- list()

for (seg in segment) {
  
  if (is.na(seg)){seg = NULL}
  
  for (nact in nact_status) {
    
    if (is.na(nact)){nact = NULL}
    
    for (ann in annotation) {
      
      if (is.na(ann)){ann = NULL}
      
      key <- paste0(seg, "_", nact, "_", ann)
      cellchat_results[[key]] <- cellchat_predict_prob(
        metadt_all,
        expr_deseq2_norm_log,
        nact_status = nact,
        segment = seg,
        annotation = ann
      )

    }
  }
}


# save final LR predictions


cellchat_output = list(
  expr_deseq2_norm_log = expr_deseq2_norm_log,
  metadt_all = metadt_all,
  cellchat_results = cellchat_results
  )

saveRDS(cellchat_output, file = paste0(output_dir,'/cellchat_output.RDS'))



