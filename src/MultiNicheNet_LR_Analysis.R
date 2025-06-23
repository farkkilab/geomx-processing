# Ligand Receptor Analysis  by CellChat : Geomx Deconvoluted data


# param


organism = "human"
cell_types = c("Tcells","Macrophages") 
celltype_id = "labels" 
cell_idents = c("Tcells","Macrophages")
sample_id = "new_sample_ID"
group_id =  "Segment" 
batches = NA
covariates =  "samples" #  "Patient" #  Change the s to S later

min_cells = 4
cell_frac_cutoff = 0.01 # cutoff was null in run 02 should set to 0.01
min_sample_prop = 0.50 
fraction_cutoff = 0.05
logFC_threshold = 0.5 # 0.5
p_val_threshold = 0.05 # 0.05
p_val_adj = TRUE # TRUE 
empirical_pval = FALSE
ligand_activity_down = FALSE


top_n_target = 250 #250
verbose = TRUE
cores_system = 2


# Set contrasts
contrasts_oi <- c("'stroma-tumor','tumor-stroma'")
# Create a contrast table
contrast_tbl <- tibble(contrast = c("stroma-tumor","tumor-stroma"), 
                       group = c("stroma", "tumor"))




# create directory for output files

output_folder_name  = "/MultiNicheNet_objects"

if (!dir.exists(paste0(output_dir,output_folder_name))) {
  new_dir = paste0(output_dir,output_folder_name)
  dir.create(new_dir,recursive = TRUE)
  message("Directory created")
  output_dir = new_dir
} else {
  output_dir = paste0(output_dir,output_folder_name)
}


#  Extracting expression data of the desired cell types and combining


if (is.null(cell_types)) { # If the cell types are not defined take all the cell types in the prism object
  ct_names <- colnames(bprism_res@posterior.theta_f@theta.cv)
  
}else{
  ct_names <- cell_types
}


# TODO check if temp exsist 

temp <- extract_and_combine_expression_data(cell_types, bprism_res, ct_names)

deconv_ct_list_int_all = temp$deconv_ct_list_int_all
metadt_all = temp$metadt_all



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
deconv_ct_list_int_all <- deconv_ct_list_int_all[, colnames(deconv_ct_list_int_all) %in% filtered_samples]
deconv_ct_list_combined = t(deconv_ct_list_int_all)



# TODO check Filter based on meta data


meta = data.frame(samples = metadt_all$Sample, NACT_status = metadt_all$NACT_status, Patient = metadt_all$Patient, Segment = metadt_all$Segment , Annotation = metadt_all$Annotation_cell, row.names = rownames(metadt_all))
meta = data.frame(meta , new_sample_ID = paste0(meta$samples,'_',meta$Segment))
#meta = meta %>% filter(NACT_status == NACT_status_filtered_group) 
#meta = meta %>% filter(Segment  == Segment_filtered_group)

# for paired samples : should edit this code

# paired_samples <- meta %>%
#      group_by(Patient) %>%
#      filter(all(c("pre", "post") %in% NACT_status)) %>%
#      summarise() %>%
#      pull(Patient)
# #
# meta = meta %>% filter(Patient %in% paired_samples)

#################

cell_types <- sapply(strsplit(rownames(meta), "_"), function(x) x[2])
meta = data.frame(labels = cell_types, meta)
meta$samples <- as.factor(meta$samples)
meta$labels <- as.factor(meta$labels)
meta_filter = rownames(meta) %in% colnames(deconv_ct_list_combined)
meta = meta[meta_filter,]
#meta = data.frame(meta, new_sample_id = paste0(meta$samples,"_",meta$Segment)) # have to change this depending on the group

deconv_ct_list_combined_filter = colnames(deconv_ct_list_combined) %in% rownames(meta)
deconv_ct_list_combined = deconv_ct_list_combined[,deconv_ct_list_combined_filter]


# creating single cell experiment object


sce <- SingleCellExperiment(
  assays = list(counts = deconv_ct_list_combined),
  colData = meta
)


# make sure that gene symbols used in the expression data are updated
sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()

# Define sender and receiver cell types
senders_oi <- SummarizedExperiment::colData(sce)[, celltype_id] %>% unique() %>% .[.%in% cell_idents]
receivers_oi <- SummarizedExperiment::colData(sce)[, celltype_id] %>% unique() %>% .[.%in% cell_idents]
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
]




# TODO moving to utils and changing the model dir 


dir_name = "C:/Users/Sahas/Downloads/Masters_Thesis/Ligand-receptor/NichNet/Nichenet_Model/" # replace with the url at the end

options(timeout = 120)

if(organism == "human"){
  
  lr_network_all = 
    readRDS(url(
      "https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds"
    )) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(paste0(dir_name,"ligand_target_matrix_nsga2r_final.rds"))
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
} 


# diagnostic plots

# To check whether each cell type have enough number of cells in each sample

abundance_expression_info <- get_abundance_info(sce = sce, 
                                                sample_id = sample_id, 
                                                group_id = group_id, 
                                                celltype_id = celltype_id, 
                                                min_cells = min_cells, 
                                                senders_oi = senders_oi, 
                                                receivers_oi = receivers_oi, 
                                                #lr_network = lr_network, 
                                                batches = batches)



pdf(file = paste0(output_dir,"abund_plot.pdf"), width = 17, height = 10)
abundance_expression_info$abund_plot_sample
abundance_expression_info$abund_plot_group
abundance_expression_info$abund_barplot
dev.off()


################################ Analysis ###############################

# Step 01. Cell-type filtering

abundance_df_summarized = abundance_expression_info$abundance_data %>% 
  mutate(keep = as.logical(keep)) %>% 
  group_by(group_id, celltype_id) %>% 
  summarise(samples_present = sum((keep)))

celltypes_absent_one_condition = abundance_df_summarized %>% 
  filter(samples_present == 0) %>% pull(celltype_id) %>% unique() 
# find truly condition-specific cell types by searching for cell types 
# truely absent in at least one condition

celltypes_present_one_condition = abundance_df_summarized %>% 
  filter(samples_present >= 2) %>% pull(celltype_id) %>% unique() 
# require presence in at least 2 samples of one group so 
# it is really present in at least one condition

condition_specific_celltypes = intersect(
  celltypes_absent_one_condition, 
  celltypes_present_one_condition)

total_nr_conditions = SummarizedExperiment::colData(sce)[,group_id] %>% 
  unique() %>% length() 

absent_celltypes = abundance_df_summarized %>% 
  filter(samples_present < 2) %>% 
  group_by(celltype_id) %>% 
  dplyr::count() %>% 
  filter(n == total_nr_conditions) %>% 
  pull(celltype_id)



# print("condition-specific celltypes:")
# print(condition_specific_celltypes)
# print("absent celltypes:")
# print(absent_celltypes)



############## Step 02. Gene Filtering ############

frq_list = get_frac_exprs(
  sce = sce, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, 
  min_sample_prop = min_sample_prop)

# Now only keep genes that are expressed by at least one cell type:

genes_oi = frq_list$expressed_df %>% 
  filter(expressed == TRUE) %>% pull(gene) %>% unique() 

sce = sce[genes_oi, ]



################ Step 03. Pseudobulk expression calculation ####################

abundance_expression_info = process_abundance_expression_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_expression_info)


################ Step 04. Differential expression (DE) analysis: ####################

DE_info = get_DE_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  batches = batches, covariates = covariates, 
  contrasts_oi = contrasts_oi, 
  min_cells = min_cells, 
  expressed_df = frq_list$expressed_df)


if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
  celltype_de = DE_info_emp$de_output_tidy_emp %>% select(-p_val, -p_adj) %>% 
    rename(p_val = p_emp, p_adj = p_adj_emp)

} else {
  celltype_de = DE_info$celltype_de$de_output_tidy

} 


# Combine DE information for ligand-senders and receptors-receivers

sender_receiver_de = multinichenetr::combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)


# plot to check the p value distribution

pdf(file = paste0(output_dir,"hist_pvals.pdf"), width = 17, height = 10)
DE_info$hist_pvals
dev.off()


# TODO save geneset_assessment

geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 

#saveRDS(geneset_assessment, file = paste0(output_dir,"geneset_assessment.RDS"))


###### Step 05. Ligand activity prediction #####

# Perform the ligand activity analysis and ligand-target inference
# increase the number of scores to run the code faster

n.cores = min(cores_system, celltype_de$cluster_id %>% unique() %>% length()) 


ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = verbose, 
    n.cores = n.cores
    #verbose = TRUE
  )
))



###### Step 06. Prioritization: rank cell-cell communication patterns through multi-criteria prioritization


sender_receiver_tbl = sender_receiver_de %>% distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group")
}


# This table gives the final prioritization score of each interaction, and the values of the individual prioritization criteria.

prioritization_tables = suppressMessages(multinichenetr::generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular", # all prioritization criteria will be weighted equally
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down
))



# Calculate correlation

# TODO chcek thr output 

lr_target_prior_cor = lr_target_prior_cor_inference(
  receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(), 
  abundance_expression_info = abundance_expression_info, 
  celltype_de = celltype_de, 
  grouping_tbl = grouping_tbl, 
  prioritization_tables = prioritization_tables, 
  ligand_target_matrix = ligand_target_matrix, 
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj
)




# Save the Output

multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
) 
multinichenet_output = make_lite_output(multinichenet_output)


saveRDS(multinichenet_output, paste0(output_dir, "multinichenet_output.rds"))
  
