# Ligand Receptor Analysis  by CellChat : Geomx Deconvoluted data
# TODO NicheNet code is absolutely dreadful
# TODO the whole concept is based on the ideas of comparing 2 conditions
# TODO for DGE theres an assumption that 1 sample belongs to 1 condition 
# TODO in general it might not be useful for the future - abandoned for now
# TODO maybe with adding externally calculated DGE..

# define params -----------------------------------------------------------

# for multiNicheNetR if you are providing DEGS externally follow the MultiNicheNet_LR_util.R script to prepare the DEGs dataframe. 
# Else MUltiNicheNet will not work
celltype_de_external_path = NULL # either NULL or path to the prepared DEGs dataframe to  'celltype_de_external'. Eg: celltype_de_external = readRDS(file.path(output_dir,MultiNicheNet_folder_name,"celltype_de_combined_calculated_externally.RDS"))

# which groups from combination of grouping_var_col_ids should be compared (if NULL: everything with everything)
# note that even for 4 groups, there will be already 12 combinations so choose wisely!
# see documentation for get_DE_info()
# format eg: c('stroma_pre-stroma_post', 'stroma_post-stroma_pre') each comparison in both directions
groups_to_compare <- NULL

sample_name <- 'Sample'
grouping_var_col_ids <- c('Segment')

# which from grouping_var_col_ids are connected to ROI types within sample (such as Segment)
# and not directly to sample (such as NACT_status)
# if there are no such vars - NULL
grouping_var_col_ids_within_sample <- c('Segment')

cell_frac_cutoff = 0.005 # 0.01 or 0.005 ct specific expr from dcc with ct fraction lower than cutoff will be removed
min_cells = 3 # minimum number of rois containing given cell > cell_frac_cutoff per cell type per sample.Samples that have less than min_cells cells will be excluded from the analysis for that specific cell type

min_sample_prop = 0.1 # genes expressed if they are expressed in at least a min_sample_prop fraction of samples in the condition with the lowest number of samples
fraction_cutoff = 0.05 # genes as expressed if they have non-zero expression values in a fraction_cutoff fraction of cells of that cell type in that sample


# variables to merge the final csv with
meta_names <- c(aoi_id, roi_id, aoi_segment_var, sample_name, main_experimental_condition, 
                grouping_var_col_ids, main_batch_var, secondary_batch_var)

# often not enough cells in all batches and throws errors
main_batch_var <- "main_batch_nr"
main_batch_var <- NA

# path to the prepared normalized pseudo scRNaseq dataset from all bpres
bp_pseudosc_path <- file.path(output_dir, 'lr_interactions', 
                              paste0('bp_res_pseudosc_', scrna_anno, 'ct_frac_', cell_frac_cutoff, '_norm.csv'))

outct <- ifelse(!is.null(cell_types_selected), paste(cell_types_selected, collapse = '_'), 'all')

# download necessary files ------------------------------------------------

# !! this may change in the future! 
# keep up with https://github.com/saeyslab/nichenetr  and https://github.com/saeyslab/multinichenetr for updates
if(!file.exists(file.path(nichenet_data_dir, "lr_network_human_allInfo_30112033.rds"))){
  download.file('https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds', 
                destfile = file.path(nichenet_data_dir, "lr_network_human_allInfo_30112033.rds"), method = "wget", extra = "-r -p --random-wait")
}

if(!file.exists(file.path(nichenet_data_dir, "ligand_target_matrix_nsga2r_final.rds"))){
  download.file('https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds', 
                destfile = file.path(nichenet_data_dir, "ligand_target_matrix_nsga2r_final.rds"), method = "wget", extra = "-r -p --random-wait")
  
}

if(!file.exists(file.path(nichenet_data_dir, "weighted_networks_nsga2r_final.rds"))){
  download.file('https://zenodo.org/record/10229222/files/weighted_networks_nsga2r_final.rds', 
                destfile = file.path(nichenet_data_dir, "weighted_networks_nsga2r_final.rds"), method = "wget", extra = "-r -p --random-wait")
}

# lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
# ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
# weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

lr_network_all <- readRDS(file.path(nichenet_data_dir, "lr_network_human_allInfo_30112033.rds"))
ligand_target_matrix <- readRDS(file.path(nichenet_data_dir, "ligand_target_matrix_nsga2r_final.rds"))

# TODO cannot find this file
weighted_networks <- readRDS(file.path(nichenet_data_dir, "weighted_networks_nsga2r_final.rds"))

# load geomx metadata -----------------------------------------------------

geomx_obj <<- readRDS(geomx_norm_batch_eff_rm_path) # batch effect corrected Geomx Object
meta_data_all <- sData(geomx_obj)[, unique(meta_names)]

# TODO code repetition from BSR Analysis
# make variable with all categories from grouping_var_col_ids
meta_data_all$comparison_group <- apply(meta_data_all, 1, function(row){
  group <- sapply(grouping_var_col_ids, function(var){
    paste(row[var])
  })
  group <- paste(group, collapse = '_')
  return(group)
})

rm(geomx_obj)
gc()

# create normalised pseudo scRNAseq dataset -------------------------------

# create a combined 'artificial pseudo-bulk scRNAseq' dataset with all ct specific counts 
if(file.exists(bp_pseudosc_path)){
  bprism_res_norm <- as.matrix(fread(bp_pseudosc_path), rownames=1)
} else{
  bprism_res_norm <- create_norm_pseudosc_from_deconv(bp_res_path, bp_ct_frac_path, scrna_anno, cell_frac_cutoff, bp_pseudosc_path)
}

# make metadata -----------------------------------------------------------

# make metadata with dcc_cell type
dcc_ct <- data.frame('dcc_filename' = gsub('_.*', '', colnames(bprism_res_norm)), 
                     'dcc_ct' = colnames(bprism_res_norm))
meta_data_ct <- left_join(dcc_ct, meta_data_all)
meta_data_ct$ct_label <- gsub('^[^_]*', '', meta_data_ct$dcc_ct)
meta_data_ct$ct_label <- gsub('^_', '', meta_data_ct$ct_label)
rownames(meta_data_ct) <- meta_data_ct$dcc_ct
meta_data_ct$samples <- meta_data_ct[[sample_name]]
meta_data_ct$celltype_id <- meta_data_ct$ct_label

# hacking NichenetR - 1 sample can only be in 1 group
if(is.null(grouping_var_col_ids_within_sample)){
  meta_data_ct$sample_id <- meta_data_ct[[sample_name]]
} else{
  # add additional ROI type labels to sample name to treat it as separate sample
  meta_data_ct$sample_id <- apply(meta_data_ct, 1, function(row){
    sample_add <- paste(row[grouping_var_col_ids_within_sample], collapse = '_')
    sample_id <- paste(c(row[sample_name], sample_add), collapse = '_')
  })
  
}

bprism_res_norm <- bprism_res_norm[, !grepl('Mast_cells', colnames(bprism_res_norm))]
meta_data_ct <- meta_data_ct[meta_data_ct$ct_label != 'Mast_cells', ]

if(!is.null(cell_types_selected)){
  cell_idents <- cell_types_selected
} else{
  cell_idents <- unique(meta_data_ct$ct_label)
}

# which groups should be compared (by default: everything)
if(!is.null(groups_to_compare)){
  contrasts_oi <- groups_to_compare
} else{
  # find all combinations of comparison_groups
  comp_pairs <- combn(unique(meta_data_ct$comparison_group),2)
  
  # get vector with combinations, including reverse
  # see documentation for get_DE_info()
  contrasts_oi <- unique(as.vector(apply(comp_pairs, 2, function(x){
    comb <- paste(x, collapse = '-')
    comb_rev <- paste(rev(x), collapse = '-')
    return(c(comb, comb_rev))
  })))
}

contrasts_oi <- sapply(contrasts_oi, function(x){paste0("'", x, "'" )})
contrasts_oi <- paste(contrasts_oi, collapse = ',')

# prepare input expr mtx --------------------------------------------------

# creating single cell experiment object
sce <- SingleCellExperiment(
  assays = list(counts = bprism_res_norm),
  colData = meta_data_ct
)

# make sure that gene symbols used in the expression data are updated
sce = makenames_SCE(alias_to_symbol_SCE(sce, "human")) 

# Define sender and receiver cell types
# filter to interesting cell types
senders_oi <- SummarizedExperiment::colData(sce)[, 'ct_label'] %>% unique() %>% .[.%in% cell_idents]
receivers_oi <- SummarizedExperiment::colData(sce)[, 'ct_label'] %>% unique() %>% .[.%in% cell_idents]
sce = sce[, SummarizedExperiment::colData(sce)[,'ct_label'] %in% 
            c(senders_oi, receivers_oi)
]

# filter lr network files -------------------------------------------------

# TODO move it to outside function
colnames(ligand_target_matrix) = make.names(convert_alias_to_symbols(colnames(ligand_target_matrix), organism = 'human'))
rownames(ligand_target_matrix) = make.names(convert_alias_to_symbols(rownames(ligand_target_matrix), organism = 'human'))

lr_network_all = lr_network_all %>% 
  mutate(ligand = make.names(convert_alias_to_symbols(ligand, organism = 'human')), 
         receptor = make.names(convert_alias_to_symbols(receptor, organism = 'human'))) 

lr_network = lr_network_all %>% 
  distinct(ligand, receptor) %>%
  filter(ligand %in% colnames(ligand_target_matrix))

ligand_target_matrix = ligand_target_matrix[, colnames(ligand_target_matrix) %in% lr_network$ligand]


# check cell type abundance -----------------------------------------------
# diagnostic plots
# To check whether each cell type have enough number of cells in each sample

abundance_expression_info <- get_abundance_info(sce = sce, 
                                                sample_id = 'sample_id', 
                                                group_id = 'comparison_group', 
                                                celltype_id = 'ct_label', 
                                                min_cells = min_cells, 
                                                senders_oi = senders_oi, 
                                                receivers_oi = receivers_oi,
                                                batches = main_batch_var)


# plotted by comparison + batch
pdf(file.path(plot_dir,'abund_plot.pdf'), width = 17, height = 10)
print(abundance_expression_info$abund_plot_sample)
print(abundance_expression_info$abund_plot_group)
print(abundance_expression_info$abund_barplot)
dev.off()

# NicheNet Analysis 01 cell type filtering --------------------------------

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

total_nr_conditions = length(unique(SummarizedExperiment::colData(sce)[,'comparison_group']))

absent_celltypes = abundance_df_summarized %>%
  filter(samples_present < 2) %>%
  group_by(celltype_id) %>%
  dplyr::count() %>%
  filter(n == total_nr_conditions) %>%
  pull(celltype_id)

print("condition-specific celltypes:")
print(condition_specific_celltypes)
print("absent celltypes:")
print(absent_celltypes)

#TODO filter out absent cell types


# NicheNet analysis 02 - gene filtering ------------------------------------

frq_list = get_frac_exprs(
  sce = sce, 
  sample_id = 'sample_id', 
  group_id = 'comparison_group', 
  celltype_id = 'ct_label', 
  batches = main_batch_var,
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, 
  min_sample_prop = min_sample_prop)


# Now only keep genes that are expressed by at least one cell type:
genes_oi = frq_list$expressed_df %>% 
  filter(expressed == TRUE) %>% pull(gene) %>% unique() 

sce = sce[genes_oi, ]

# NicheNet analysis 03 pseudobulk expression calculation ------------------

abundance_expression_info = process_abundance_expression_info(
  sce = sce, 
  sample_id = 'sample_id',
  group_id = 'comparison_group',
  celltype_id = 'ct_label', 
  min_cells = min_cells, 
  senders_oi = senders_oi, 
  receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = main_batch_var, 
  frq_list = frq_list, 
  abundance_info = abundance_expression_info)

#TODO Warning message:
# In get_pseudobulk_logCPM_exprs(sce, sample_id = sample_id, celltype_id = celltype_id,  :
# Not all possible group-batch/batch combinations are present in your data. 
# This will result in errors during the batch effect correction process of Combat and/or Muscat DE analysis. 
# Please reconsider the groups and batches you defined.

# NicheNet Analysis 04 Differential Expression ----------------------------


if(!is.null(celltype_de_external_path)){
  celltype_de = readRDS(celltype_de_external_path)
} else{
  
  DE_info = get_DE_info(
    sce = sce, 
    sample_id = 'sample_id', 
    group_id = 'comparison_group',
    celltype_id = 'ct_label', 
    batches = main_batch_var, 
    covariates = 'sample_id', 
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

pdf(file.path(plot_dir,'hist_pvals.pdf'), width = 17, height = 10)
print(DE_info$hist_pvals)
dev.off()


# TODO save geneset_assessment

geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 

write.csv(geneset_assessment, file.path(output_dir,MultiNicheNet_folder_name,"geneset_assessment.csv"))


if (all(geneset_assessment$n_geneset_up == 0)) {
  stop("Execution halted: there are no up regulated genes between groups. check the geneset_assessment table saved in MultiNicheNet output folder")
}


############################################################################################
##############################################################################################
##############################################################################################
# param
comparison <- c('stroma', 'tumor')


cell_frac_cutoff = 0.01 # cutoff to filter out AOIs based on the bayesprsim cell fraction values

organism = "human"
batches = NA # TODO: gives an error when specified
min_cells = 4 # minimum number of cells per cell type per sample.Samples that have less than min_cells cells will be excluded from the analysis for that specific cell type
min_sample_prop = 0.50 # genes expressed if they are expressed in at least a min_sample_prop fraction of samples in the condition with the lowest number of samples
fraction_cutoff = 0.05 # genes as expressed if they have non-zero expression values in a fraction_cutoff fraction of cells of that cell type in that sample
logFC_threshold = 0.5 
p_val_threshold = 0.05 
p_val_adj = TRUE  
empirical_pval = FALSE # TODO : this i did not check. In case p-value distributions look irregular in p-value histogram, you can estimate empirical p-values
ligand_activity_down = FALSE # to focus specifically on upregulating ligands keep it FALSE
n = 50 # Number of top n LRpairs for visualization from each group 


# for ligand-target inference procedure, need to select which top n of the predicted target genes will be considered (here: top 250 targets per ligand). 
# This parameter will not affect the ligand activity predictions. 
# It will only affect ligand-target visualizations and construction of the intercellular regulatory network during the downstream analysis.  
top_n_target = 250

verbose = TRUE
cores_system = detectCores()-4

##################
############################
#############################

# TODO if null, take all types from metadf
cell_idents = cell_types_selected


# TODO generate contrasts and contrast table from grouping_var_col_ids

# TODO else create a different group_id column

# Set contrasts

contrasts_oi <- paste0(
  "'", comparison[2], "-", comparison[1], "'", ",",
  "'", comparison[1], "-", comparison[2], "'"
)
# Create a contrast table
contrast_tbl <- tibble(contrast = c(paste(comparison[1], comparison[2], sep = "-"),paste(comparison[2], comparison[1], sep = "-")), 
                       group = c(comparison[1], comparison[2]))



#  Extracting expression data of the desired cell types and combining

if (is.null(cell_types)) { # If the cell types are not defined take all the cell types in the prism object
  ct_names <- colnames(bprism_res@posterior.theta_f@theta.cv)
  
}else{
  ct_names <- cell_types
}


file_name = paste(ct_names, collapse = "_")
combined_expression_data_path = file.path(output_dir,paste0(file_name,'_expr_and_meta_list.RDS'))


  if(!file.exists(combined_expression_data_path)){
  
    # expr_and_meta_list
    
    combined_expression_data <- extract_and_combine_expression_data(bprism_res, ct_names, geomx_obj,aoi_id, sample_name, aoi_segment_var, main_experimental_condition, grouping_var_col_ids)
    deconv_ct_list_int_all <- combined_expression_data$deconv_ct_list_int_all
    metadt_all <- combined_expression_data$metadt_all
    
    print(paste0(' combined expression matrix and meta data generation succeeded!'))
    saveRDS(combined_expression_data, file = combined_expression_data_path)
    
  } else{
    
    
    combined_expression_data = readRDS(combined_expression_data_path)
    deconv_ct_list_int_all <- combined_expression_data$deconv_ct_list_int_all
    metadt_all <- combined_expression_data$metadt_all
    
    # print(paste0(' already produce the combined expression matrix and meta data'))
  }



combined_group_names = unique(c(sample_name,grouping_var_col_ids))
metadt_all$new_sample_ID <- apply(metadt_all[, combined_group_names], 1, function(x) paste0(x, collapse = "_"))


if(length(grouping_var_col_ids) == 1){
  
  group_id =  grouping_var_col_ids 
  
} else {
  
  metadt_all$new_group_ID <- apply(metadt_all[, grouping_var_col_ids], 1, function(x) paste0(x, collapse = "_"))
  group_id = "new_group_ID"
}


metadt_all = data.frame(labels = sapply(strsplit(rownames(metadt_all), "_"), function(x) x[2]), metadt_all)
metadt_all$labels <- as.factor(metadt_all$labels)
metadt_all <- metadt_all %>%
  mutate(
    across(all_of(sample_name), as.factor)
  )

celltype_id = "labels"
sample_id = "new_sample_ID" 

# filtering based on cell fraction 

deconv_ct_list_combined = t(deconv_ct_list_int_all)
expr_deseq2_norm_log_cf_filtered = filter_based_on_cell_fraction(ct_names, cell_frac_cutoff, cell_fractions_df, deconv_ct_list_combined)
metadt_all <- metadt_all[rownames(metadt_all) %in% colnames(expr_deseq2_norm_log_cf_filtered),]
expr_deseq2_norm_log_cf_filtered = expr_deseq2_norm_log_cf_filtered[,colnames(expr_deseq2_norm_log_cf_filtered) %in% rownames(metadt_all)]


#################################################33
######################################################

# creating single cell experiment object

sce <- SingleCellExperiment(
  assays = list(counts = expr_deseq2_norm_log_cf_filtered),
  colData = metadt_all
)


# make sure that gene symbols used in the expression data are updated
sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()

# Define sender and receiver cell types
senders_oi <- SummarizedExperiment::colData(sce)[, celltype_id] %>% unique() %>% .[.%in% cell_idents]
receivers_oi <- SummarizedExperiment::colData(sce)[, celltype_id] %>% unique() %>% .[.%in% cell_idents]
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
]


# TODO check options(timeout = 120)

if(organism == "human"){
  
  lr_network_all = lr_network_all %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  
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
                                                batches = batches)



pdf(file.path(plot_dir,'abund_plot.pdf'), width = 17, height = 10)
print(abundance_expression_info$abund_plot_sample)
print(abundance_expression_info$abund_plot_group)
print(abundance_expression_info$abund_barplot)
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

# IF there are externally provided DEGS

if(external_DE_info == TRUE){
  
  celltype_de = celltype_de_external
  
  
} else{
  
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

pdf(file.path(plot_dir,'hist_pvals.pdf'), width = 17, height = 10)
print(DE_info$hist_pvals)
dev.off()


# TODO save geneset_assessment

geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 

write.csv(geneset_assessment, file.path(output_dir,MultiNicheNet_folder_name,"geneset_assessment.csv"))


if (all(geneset_assessment$n_geneset_up == 0)) {
  stop("Execution halted: there are no up regulated genes between groups. check the geneset_assessment table saved in MultiNicheNet output folder")
}


###### Step 05. Ligand activity prediction #####

# Perform the ligand activity analysis and ligand-target inference
# increase the number of scores to run the code faster

# TODO celltype_de$cluster_id

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


# get top n LR pairs from each cell type communication combinations and save it as multiple csv files
# theses df can be directly used for plotting or can do manually filtering and provide seperetly


top_n_LR_pairs = list()

for (group in comparison) {
  for (receiver in cell_types) {
    for (sender in cell_types) {
      
      print(paste0(group,"-",receiver,"-",sender))
      prioritized_tbl_oi = get_top_n_lr_pairs(
        prioritization_tables, 
        n, 
        groups_oi = group, 
        receivers_oi = receiver,
        senders_oi = sender)
      
      
      # ligand-receptor pseudobulk product expression panel
      sample_data = prioritization_tables$sample_prioritization_tbl %>% 
        dplyr::filter(id %in% prioritized_tbl_oi$id) %>% 
        dplyr::mutate(
          sender_receiver = paste(sender, receiver, sep = " --> "), 
          lr_interaction = paste(ligand, receptor, sep = " - ")) %>%
        dplyr::arrange(receiver) %>% 
        dplyr::group_by(receiver) %>%  
        dplyr::arrange(sender, .by_group = TRUE)
      
      sample_data = sample_data %>% 
        dplyr::mutate(sender_receiver = factor(
          sender_receiver, 
          levels = sample_data$sender_receiver %>% unique()
        ))
      
      ################################################################
      
      keep_sender_receiver_values = c(0.25, 0.9, 1.75, 4) # TODO check
      names(keep_sender_receiver_values) = levels(sample_data$keep_sender_receiver)
      
      ######## calculate the median bulk expression for each group
      
      # calculate the median
      
      group_medians <- sample_data %>%
        group_by(group,lr_interaction) %>%
        summarize(median_scaled_LR = median(scaled_LR_pb_prod, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = group, values_from = median_scaled_LR)
      
      
      # Compute median difference (stroma - tumor)
      group_medians <- group_medians %>%
        mutate(
          diff_median = .[[comparison[1]]] - .[[comparison[2]]]
        )
      
      # Wilcoxon test per interaction
      wilcox_results <- sample_data %>%
        group_by(lr_interaction) %>%
        filter(group %in% comparison) %>%
        summarize(
          test = list(wilcox.test(scaled_LR_pb_prod ~ group)),
          .groups = "drop"
        ) %>%
        mutate(
          p_value = map_dbl(test, "p.value"),
          neg_log10_p = -log10(p_value)
        ) %>%
        select(lr_interaction, p_value, neg_log10_p)
      
      
      adj_pvals <- p.adjust(wilcox_results$p_value, method = "BH")
      
      # Merge with fold change data
      final_data <- group_medians %>%
        left_join(wilcox_results, by = "lr_interaction")
      
      
      final_data$adj_p_value = adj_pvals
      final_data$neg_log10_p_adj = -log10(adj_pvals)
      
      
      sender_receiver <- paste(sender, receiver, sep = " --> ")
      final_data$sender_receiver <- rep(sender_receiver, nrow(final_data))
      final_data$group <- rep(paste(comparison, collapse = "-"), nrow(final_data))
      df_plot1 = final_data
      
      
      
      #########################################################################
      
      group_data = multinichenet_output$prioritization_tables$group_prioritization_table_source  %>% 
        dplyr::mutate(
          sender_receiver = paste(sender, receiver, sep = " --> "), 
          lr_interaction = paste(ligand, receptor, sep = " - "))  %>% 
        dplyr::distinct(id, sender, receiver, sender_receiver, ligand, receptor, lr_interaction, group, activity_scaled, direction_regulation, prioritization_score) %>% 
        dplyr::filter(id %in% sample_data$id) %>% 
        dplyr::arrange(receiver) %>% 
        dplyr::group_by(receiver) %>% 
        dplyr::arrange(sender, .by_group = TRUE)
      
      df_plot2 = group_data %>% dplyr::mutate(
        sender_receiver = factor(
          sender_receiver, 
          levels = group_data$sender_receiver %>% unique()
        ))
      
      
      
      #################################################################
      
      
      table_name <- paste(group, receiver, sender, sep = "_")
      file_name_p1 <- file.path(plot_dir, paste0(table_name, "_LR_pairs_dfplot_median_bulk_expr.csv"))
      write.csv(df_plot1, file_name_p1, row.names = FALSE)
      file_name_p2 <- file.path(plot_dir, paste0(table_name, "_LR_pairs_dfplot_ligand_activity.csv"))
      write.csv(df_plot2, file_name_p2, row.names = FALSE)
      top_n_LR_pairs[[table_name]] <- list(df_plot1 = df_plot1, df_plot2 = df_plot2)
      
    }
    
  }
  
}





# Save the Output

multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor,
  top_n_LR_pairs = top_n_LR_pairs
) 
multinichenet_output = make_lite_output(multinichenet_output)


saveRDS(multinichenet_output, geomx_MultiNicheNet_path)


# check the warning

# Warning messages:
#   1: In get_frac_exprs(sce = sce, sample_id = sample_id, celltype_id = celltype_id,  :
#                          There are some genes with NA/NaN fraction of expression. This is the result of the muscat function `calcExprFreqs` which will give NA/NaN when there are no cells of a particular cell type in a particular group or no cells of a cell type in one sample. As a temporary fix, we give all these genes an expression fraction of 0 in that group for that cell type
#                        2: In DGEList.default(pb@assays@data[[celltype_oi]]) :
#                          At least one library size is zero


# TODO to write logs

# # write logs --------------------------------------------------------------
# 
# # save logs
# writeLines(c('deconvolution logs:',
#              '; spatial decon normalisation type : ', norm_type,
#              '; minimum cell type number : ', ct_nr_thr,
#              '; cell type annotation  : ', scrna_anno,
#              '; limma primary batch effect variable : ', primary_batch_var,
#              '; limma secondary batch effect variable : ', secondary_batch_var,
#              '; limma experimental design : ', as.character(exp_design)[2],
#              '; limma covariate : ', as.character(cov_design)[2]), deconv_logs_path)