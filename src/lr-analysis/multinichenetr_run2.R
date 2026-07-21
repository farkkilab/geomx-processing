library(GeomxTools, quietly =T)
library(ggplot2, quietly =T)

# For MultiNicheNetr
library(SingleCellExperiment, quietly =T)
library(nichenetr, quietly =T)
library(multinichenetr, quietly =T)
library(parallel)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
library(data.table)
library(dplyr)
library(tibble)


# from master file
proj_dir <<- '~/Documents/phd/st'
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808')
geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS') 

aoi_id <<- 'dcc_filename'
roi_id <<- 'Roi_geomx'
aoi_segment_var <<- "Segment"
sample_name <- "Sample"
main_experimental_condition <<- 'NACT_status'

custom_metadt_path <<- "~/Documents/phd/st/data/geomx/metadata_full_SENSITIVE.csv"
dge_dir <- "/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2808/dge/dge_within_slide_roi_cluster_label_gmm_bin_TRUE_Macro_domin_Segment_NACT_status"


plot_dir = file.path(output_dir, 'lr_interactions', 'multi_niche_netr','plots_and_csv_files')
dir.create(plot_dir , recursive = T, showWarnings = F)

# path where to download nichenet data
nichenet_data_dir <<- file.path(proj_dir, 'geomx-processing', 'data', 'nichenet')
dir.create(nichenet_data_dir, recursive = T, showWarnings = F)

source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))
source(file.path(proj_dir, 'geomx-processing', 'src', 'lr-analysis', 'MultiNicheNet_LR_Visualization.R'))

# define params -----------------------------------------------------------

scrna_anno <<- 'mid_lvl_ct_updated' #either 'cell_type' / 'mid_lvl_ct' / 'mid_lvl_ct_updated' / 'low_lvl_ct'

# names of cells to fin
# c("tumor", "Macrophages_Monocytes", "Tcells_CD8", "Tcells_CD4", "DCs", "Fibroblasts_Mesothelial", "Bcells")

cell_types_selected <- c("Macrophages_Monocytes", "Tcells_CD8")

# for multiNicheNetR if you are providing DEGS externally follow the MultiNicheNet_LR_util.R script to prepare the DEGs dataframe. 
# Else MUltiNicheNet will not work
celltype_de_external_path = NULL # either NULL or path to the prepared DEGs dataframe to  'celltype_de_external'. Eg: celltype_de_external = readRDS(file.path(output_dir,MultiNicheNet_folder_name,"celltype_de_combined_calculated_externally.RDS"))

# which groups from combination of grouping_var_col_ids should be compared (if NULL: everything with everything)
# note that even for 4 groups, there will be already 12 combinations so choose wisely!
# see documentation for get_DE_info()
# format eg: c('stroma_pre-stroma_post', 'stroma_post-stroma_pre') each comparison in both directions
groups_to_compare <- NULL
# groups_to_exclude <- c('stroma_pre_Bcell_domin', 'stroma_pre_mixed_w_CD4', 'stroma_pre_mixed_w_others', 'tumor_pre_mixed_w_CD4')
groups_to_exclude <- NULL

sample_id <- 'dcc_filename'
covariates <- NA
batches <- NA

#grouping_var_col_ids <- 'Segment' 
grouping_var_col_ids <- c('Segment', 'NACT_status', 'roi_cluster_label_gmm')

# which from grouping_var_col_ids are connected to ROI types within sample (such as Segment)
# and not directly to sample (such as NACT_status)
# if there are no such vars - NULL
grouping_var_col_ids_within_sample <- NULL

cell_frac_cutoff = 0.01 # 0.01 or 0.005 ct specific expr from dcc with ct fraction lower than cutoff will be removed
min_cells = 1 # minimum number of rois containing given cell > cell_frac_cutoff per cell type per sample.Samples that have less than min_cells cells will be excluded from the analysis for that specific cell type

min_sample_prop = 0.25 # genes expressed if they are expressed in at least a min_sample_prop fraction of samples in the condition with the lowest number of samples
fraction_cutoff = 0.05 # genes as expressed if they have non-zero expression values in a fraction_cutoff fraction of cells of that cell type in that sample


# variables to merge the final csv with
meta_names <- unique(c(aoi_id, roi_id, aoi_segment_var, sample_name, main_experimental_condition, 
                       grouping_var_col_ids, batches, covariates, 'roi_cluster_label_gmm'))
meta_names <- meta_names[!is.na(meta_names)]

# path to the prepared normalized pseudo scRNaseq dataset from all bpres
bp_pseudosc_path <- file.path(output_dir, 'lr_interactions', 
                              paste0('bp_res_pseudosc_', scrna_anno, '_ct_frac_', cell_frac_cutoff, '_raw.csv'))


outct <- ifelse(!is.null(cell_types_selected), paste(cell_types_selected, collapse = '_'), 'all')

geomx_MultiNicheNet_path <<- file.path(output_dir, 'lr_interactions', 'multi_niche_netr', 
                                       paste0('multinichenet_output_', outct, '.rds'))

# load files --------------------------------------------------------------

organism = "human"

if(!file.exists(file.path(nichenet_data_dir, "lr_network_human_allInfo_30112033.rds"))){
  download.file('https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds', 
                destfile = file.path(nichenet_data_dir, "lr_network_human_allInfo_30112033.rds"), method = "wget", extra = "-r -p --random-wait")
}

if(!file.exists(file.path(nichenet_data_dir, "ligand_target_matrix_nsga2r_final.rds"))){
  download.file('https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds', 
                destfile = file.path(nichenet_data_dir, "ligand_target_matrix_nsga2r_final.rds"), method = "wget", extra = "-r -p --random-wait")
  
}

lr_network_all <- readRDS(file.path(nichenet_data_dir, "lr_network_human_allInfo_30112033.rds"))
ligand_target_matrix <- readRDS(file.path(nichenet_data_dir, "ligand_target_matrix_nsga2r_final.rds"))

lr_network_all = lr_network_all %>%
  mutate(
    ligand = convert_alias_to_symbols(ligand, organism = organism), 
    receptor = convert_alias_to_symbols(receptor, organism = organism)) %>%
  mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 

lr_network = lr_network_all %>% 
  distinct(ligand, receptor)

colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
  convert_alias_to_symbols(organism = organism) %>% make.names()
rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
  convert_alias_to_symbols(organism = organism) %>% make.names()

lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]


# load geomx metadata and raw pseudosc dataset ----------------------------

geomx_obj <<- readRDS(geomx_norm_batch_eff_rm_path) # batch effect corrected Geomx Object
custom_metadt <- fread(custom_metadt_path)
meta_data_all <- left_join(pData(geomx_obj), custom_metadt, by = aoi_id, suffix = c("_orig", "")) %>%
  dplyr::select(c('dcc_filename', unique(meta_names)))

#TODO make it automatic later
meta_data_all$roi_cluster_label_gmm <- ifelse(meta_data_all$roi_cluster_label_gmm == 'Macro_domin', 'Macro_domin', 'other_roi_type')

# TODO code repetition from BSR Analysis
# make variable with all categories from grouping_var_col_ids
meta_data_all$groups_to_compare <- apply(meta_data_all, 1, function(row){
  group <- sapply(grouping_var_col_ids, function(var){
    paste(row[var])
  })
  group <- paste(group, collapse = '_')
  return(group)
})

rm(geomx_obj)
gc()

bprism_res_raw <- fread(bp_pseudosc_path) %>%
  column_to_rownames('V1')

# make metadata -----------------------------------------------------------
group_id <- 'groups_to_compare'
celltype_id <- 'ct_label'

dcc_ct <- data.frame('dcc_filename' = gsub('_.*', '', colnames(bprism_res_raw)), 
                     'dcc_ct' = colnames(bprism_res_raw))
meta_data_ct <- left_join(dcc_ct, meta_data_all)
meta_data_ct$ct_label <- gsub('^[^_]*', '', meta_data_ct$dcc_ct)
meta_data_ct$ct_label <- gsub('^_', '', meta_data_ct$ct_label)
rownames(meta_data_ct) <- meta_data_ct$dcc_ct

# dashes cannot be in the names
meta_data_ct[[sample_id]] <- gsub('-', '_', meta_data_ct[[sample_id]])

# add 'batch' to batchnr
#meta_data_ct[[batches]] <- paste0('batch', as.character(meta_data_ct[[batches]]))


# load and format externally calculated DGE  ------------------------------

# load DGE which were already calculated and change format to mimic the NicheNetR results
dge_all <- lapply(cell_types_selected, function(ct_name){
  dge_filepath <- list.files(dge_dir, pattern = paste0("dge_deconv_", ct_name, ".*csv"), full.names = T)
  dge_ct <- fread(dge_filepath)
  dge_ct$cluster_id <- ct_name
  return(dge_ct)
})

dge_all <- do.call(rbind, dge_all)

# combine data group and contrast
dge_all$cont1 <- paste0(dge_all$data_group, '_', gsub(' -.*', '', dge_all$Contrast))
dge_all$cont2 <- paste0(dge_all$data_group, '_', gsub('.*- ', '', dge_all$Contrast))

if(is.null(groups_to_compare)){
  groups_to_compare <- unique(c(unique(dge_all$cont1), unique(dge_all$cont2)))
}

# in original nn theres tum-str and str-tum contrast with negatives of each value
celltype_de_pos <- data.frame(gene = dge_all$Gene, cluster_id = dge_all$cluster_id, logFC = dge_all$Estimate,
                          logCPM = 0, Fval = 0, p_val = dge_all$`Pr(>|t|)`, p.adj.loc = 0,
                          p_adj = dge_all$FDR, contrast = paste0(dge_all$cont1, '-', dge_all$cont2))
celltype_de_neg <- data.frame(gene = dge_all$Gene, cluster_id = dge_all$cluster_id, logFC = -(dge_all$Estimate),
                          logCPM = 0, Fval = 0, p_val = dge_all$`Pr(>|t|)`, p.adj.loc = 0,
                          p_adj = dge_all$FDR, contrast = paste0(dge_all$cont2, '-', dge_all$cont1))

celltype_de <- rbind(celltype_de_neg, celltype_de_pos)

# filter to groups to compare
celltype_de <- celltype_de[grepl(paste0(groups_to_compare, collapse = '|'), celltype_de$contrast), ]


rm(celltype_de_neg)
rm(celltype_de_pos)
rm(dge_all)

# make contrasts ----------------------------------------------------------

# which groups should be compared (by default: everything)
# if(!is.null(groups_to_compare)){
#   contrasts_oi <- groups_to_compare
#   comp_groups <- groups_to_compare
# } else{
#   # find all combinations of comparison_groups
#   if(!is.null(groups_to_exclude)){
#     comp_groups <- setdiff(unique(meta_data_ct$groupes_to_compare), groups_to_exclude)
#   } else{
#     comp_groups <- unique(meta_data_ct$groupes_to_compare)
#   }
#   comp_pairs <- combn(unique(meta_data_ct$groupes_to_compare),2)
#   
#   # get vector with combinations, including reverse
#   # see documentation for get_DE_info()
#   contrasts_oi <- unique(as.vector(apply(comp_pairs, 2, function(x){
#     comb <- paste(x, collapse = '-')
#     comb_rev <- paste(rev(x), collapse = '-')
#     return(c(comb, comb_rev))
#   })))
# }


# contrasts_oi_list <- sapply(contrasts_oi, function(x){paste0("'", x, "'" )})

contrasts_oi_list <- unname(sapply(unique(celltype_de$contrast), function(x){paste0("'", x, "'" )}))
contrasts_oi <- paste(contrasts_oi_list, collapse = ',')

# Create a contrast table
contrast_tbl <- tibble(contrast = contrasts_oi_list, 
                       group = unique(groups_to_compare))

# fix contrast table
contrast_tbl$contrast <- gsub("'", "", contrast_tbl$contrast) # in NichenetR 2.1 it goes wo quites

# filter metadata and expression to group, cells and samples ------------

# keep only cell types and conditions of interest
meta_data_ct <- meta_data_ct[meta_data_ct$ct_label %in% cell_types_selected & 
                               meta_data_ct$groups_to_compare %in% groups_to_compare, ]

# # keep only AOIs where all cells are present
# samples_w_all_cells <- meta_data_ct %>%
#   group_by(get(sample_id)) %>%
#   summarise(ncells = n_distinct(ct_label)) %>%
#   filter(ncells == length(cell_types_selected))
# 
# meta_data_ct <- meta_data_ct[meta_data_ct[[sample_id]] %in% samples_w_all_cells$`get(sample_id)`, ]

bprism_res_raw <- bprism_res_raw[, meta_data_ct$dcc_ct]

# prepare single cell experiment object -----------------------------------

# creating single cell experiment object
sce <- SingleCellExperiment(
  assays = list(counts = bprism_res_raw),
  colData = meta_data_ct
)

# make sure that gene symbols used in the expression data are updated
sce = makenames_SCE(alias_to_symbol_SCE(sce, "human")) 

# make sure names are valid
SummarizedExperiment::colData(sce)[[sample_id]] = SummarizedExperiment::colData(sce)[[sample_id]] %>% make.names()
SummarizedExperiment::colData(sce)[[celltype_id]] = SummarizedExperiment::colData(sce)[[celltype_id]] %>% make.names()
SummarizedExperiment::colData(sce)[[group_id]] = SummarizedExperiment::colData(sce)[[group_id]] %>% make.names()
#SummarizedExperiment::colData(sce)[[batches]] = SummarizedExperiment::colData(sce)[[batches]] %>% make.names()
#SummarizedExperiment::colData(sce)[[covariates]] = SummarizedExperiment::colData(sce)[[covariates]] %>% make.names()

senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
]

rm(bprism_res_raw)
gc()

# calculate abundance info ------------------------------------------------

abundance_info = get_abundance_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  batches = batches
)

# plotted by comparison + batch
pdf(file.path(plot_dir,'abund_plot.pdf'), width = 17, height = 10)
print(abundance_info$abund_plot_sample)
print(abundance_info$abund_plot_group)
print(abundance_info$abund_barplot)
dev.off()

# ct filtering based on abundances ----------------------------------------
analyse_condition_specific_celltypes = TRUE

abundance_df_summarized = abundance_info$abundance_data %>% 
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
  count() %>% 
  filter(n == total_nr_conditions) %>% 
  pull(celltype_id)

print("condition-specific celltypes:")
print(condition_specific_celltypes)

print("absent celltypes:")
print(absent_celltypes)

if(analyse_condition_specific_celltypes == TRUE){
  senders_oi = senders_oi %>% setdiff(absent_celltypes)
  receivers_oi = receivers_oi %>% setdiff(absent_celltypes)
} else {
  senders_oi = senders_oi %>% 
    setdiff(union(absent_celltypes, condition_specific_celltypes))
  receivers_oi = receivers_oi %>% 
    setdiff(union(absent_celltypes, condition_specific_celltypes))
}

sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
]

# filter sufficiently expressed genes -------------------------------------

frq_list = get_frac_exprs(
  sce = sce, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)

# keep genes that are expressed by at least one cell type
genes_oi = frq_list$expressed_df %>% 
  filter(expressed == TRUE) %>% pull(gene) %>% unique() 
sce = sce[genes_oi, ]


# pseudobulk dge expression calculation -----------------------------------

abundance_expression_info = process_abundance_expression_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)


# differential expression -------------------------------------------------

# DE_info = get_DE_info(
#   sce = sce, 
#   sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
#   batches = NA, covariates = covariates, 
#   contrasts_oi = contrasts_oi, 
#   min_cells = min_cells, 
#   expressed_df = frq_list$expressed_df,
#   contrast_tbl = contrast_tbl)
# 
# celltype_de_orig = DE_info$celltype_de$de_output_tidy
# 
# fwrite(celltype_de_orig, file.path(nichenet_data_dir, 'DE_info_tumor_stroma_Macro_CD8.RDS'))


# combine DGE for ligand-senders and receptor-receivers -------------------

sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)

# assess gene to bcg ratio ------------------------------------------------
logFC_threshold = 0.50
p_val_threshold = 0.05
p_val_adj = TRUE

geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 


# ligand activity analysis and ligand-target inference --------------------
top_n_target = 250
verbose = TRUE
cores_system = 20
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
  )
))


# prioritise cell-cell communication patterns through multi-criteria --------

ligand_activity_down = TRUE

sender_receiver_tbl = sender_receiver_de %>% distinct(sender, receiver)
metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% 
  tibble::as_tibble() %>% distinct()
colnames(grouping_tbl) = c("sample","group")

# These should have perfect matches (same strings, same order)
all(contrast_tbl$contrast %in% unique(sender_receiver_de$contrast))  # Should return TRUE

prioritization_tables = generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular", # all prioritization criteria will be weighted equally
  fraction_cutoff = 0.05, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down
)

# across-samples expression correlation between LR and target genes --------

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

# save all the output -----------------------------------------------------

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

saveRDS(multinichenet_output, geomx_MultiNicheNet_path)


# get top LR pairs for each communication ---------------------------------
# get top n LR pairs from each cell type communication combinations and save it as multiple csv files
# theses df can be directly used for plotting or can do manually filtering and provide seperetly

top_n_LR_pairs = list()
top_n <- 100


for (group in groups_to_compare) {
  for (receiver in cell_types_selected) {
    for (sender in cell_types_selected) {
      
      print(paste0(group,"-",receiver,"-",sender))
      prioritized_tbl_oi = get_top_n_lr_pairs(
        prioritization_tables, 
        top_n, 
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


# make plots --------------------------------------------------------------


# TODO check if works for more groups
for (group in groups_to_compare){ 
  for (receiver in cell_types_selected) {
    for (sender in cell_types_selected) {
      
      table_name <- paste(group, sender, receiver, sep = "_")
      sample_data = top_n_LR_pairs[[table_name]]
      
      
      df_plot1 = sample_data$df_plot1
      df_plot2 = sample_data$df_plot2
      
      p1 = plot_bulk_expression(df_plot1)
      p2 = plot_igand_activity(df_plot2)
      
      p = patchwork::wrap_plots(
        p1,p2,
        nrow = 1,guides = "collect",
        widths = c(6,6)
      )
      
      pdf(file = file.path(plot_dir,paste0("MultiNicheNet_plot_",table_name,".pdf")), width = 17, height = 10)
      print(p)
      dev.off()

    }
  }
}

########################
# multinichenet_output <- readRDS(geomx_MultiNicheNet_path)
# celltype_info <- multinichenet_output$celltype_info
# celltype_de <- multinichenet_output$celltype_de
# sender_receiver_info <- multinichenet_output$sender_receiver_info
# sender_receiver_de <- multinichenet_output$sender_receiver_de
# ligand_activities_targets_DEgenes <- multinichenet_output$ligand_activities_targets_DEgenes
# prioritization_tables2 <- multinichenet_output$prioritization_tables
# grouping_tbl <- multinichenet_output$grouping_tbl
# lr_target_prior_cor <- multinichenet_output$lr_target_prior_cor
# 
