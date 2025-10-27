# cellChat Functions

# CellChat Preprocessing

# extract and combine expression data and meta data (Used by both CellChat and MultiNicheNet)




extract_and_combine_expression_data <- function(bprism_res, ct_names, geomx_obj, aoi_id, sample_name, aoi_segment_var, main_experimental_condition, grouping_var_col_ids){

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

  # combine expression data

  deconv_ct_list_int <- lapply(names(deconv_ct_list), function(ct_name){
    mat <- deconv_ct_list[[ct_name]]
    rownames(mat) <- gsub('\\.dcc', paste0('_', ct_name), rownames(mat))

    mat <- apply(mat, c(1, 2), function(x) {(as.integer(x))})

    return(mat)
  })

  names(deconv_ct_list_int) <- names(deconv_ct_list)
  deconv_ct_list_int_all <- do.call(rbind, deconv_ct_list_int)

  # Meta Data

  deconv_metadt_list <- lapply(names(deconv_ct_list_int), function(ct_name){

    meta_data = sData(geomx_obj)

    # TODO here I am selecting only few of the selected columns


    meta_data = meta_data %>% select(!!sym(aoi_id), !!sym(sample_name), !!sym(aoi_segment_var), !!sym(main_experimental_condition), all_of(grouping_var_col_ids))
    # meta_data = data.frame(meta_data[,c(2,5,6,7,24,25,28)])

    meta_data$dcc_filename = gsub('\\.dcc', paste0('_', ct_name), meta_data$dcc_filename)
    mat = deconv_ct_list_int[[ct_name]]

    matched_entries = meta_data[meta_data$dcc_filename %in% rownames(mat), ]
    rownames(matched_entries) = matched_entries$dcc_filename

    return(matched_entries)
  })

  metadt_all <- do.call(rbind, deconv_metadt_list)

  return(list(
    deconv_ct_list_int_all = deconv_ct_list_int_all,
    metadt_all = metadt_all
    ))

}



# function to filter based on cell fraction (Used by both CellChat and MultiNicheNet)

filter_based_on_cell_fraction <- function(ct_names, cell_frac_cutoff, cell_fractions_df, expr){

  # TODO check this formatted_ct_names
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
  expr <- expr[, colnames(expr) %in% filtered_samples]


  return(expr)

}



# filter expression data based on colnames eg: NACT_status, segment and Annotations

filter_expr_deseq2_norm_log = function(metadt_all, expr_deseq2_norm_log, grouping_var_col_ids, group){
  

  groups <- strsplit(group, "_")[[1]]
  print(groups)
  
  
  
  for (column_name in grouping_var_col_ids){
    
    meta_data = meta_data %>% filter(!!sym(column_name) %in% groups)
    
  }
  
  
  cell_types <- sapply(strsplit(rownames(meta), "_"), function(x) x[2])
  meta = data.frame(labels = cell_types, meta)
  meta$labels <- as.factor(meta$labels)
  meta_filter = rownames(meta) %in% colnames(expr_deseq2_norm_log)
  meta = meta[meta_filter,]
  temp_filter = colnames(expr_deseq2_norm_log) %in% rownames(meta)
  expr_deseq2_norm_log_filtered = expr_deseq2_norm_log[,temp_filter]
  
  
  return(list(
    expr_deseq2_norm_log_filtered = expr_deseq2_norm_log_filtered,
    meta_filtered = meta
  ))

}



# cell type abundance plot 


cell_type_abundance_plot = function(metadt_all, min_cells, aoi_segment_var, main_experimental_condition, grouping_var_col_ids, output_dir){
  
 
  cell_types <- sapply(strsplit(rownames(metadt_all), "_"), function(x) x[2])
  meta = data.frame(labels = cell_types, metadt_all)
  meta$labels <- as.factor(meta$labels)
  
  group_vars <- c("labels",grouping_var_col_ids)  # always group by labels
  df_grouped_samples <- meta %>%
    group_by(!!!syms(group_vars)) %>%
    summarise(count = n(), .groups = "drop")
  
  
  make_plot <- function(x_var, facet_var) {
    ggplot(df_grouped_samples, aes(x = .data[[x_var]], y = count, fill = labels)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_grid(as.formula(paste(". ~", facet_var))) +
      theme_minimal() +
      geom_hline(yintercept = min_cells, linetype = "solid", color = "black") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(angle = 0)
      )
      # + ggtitle(paste("X:", x_var, "| Facet:", facet_var))
  }
  
  col_names = unique(c(aoi_segment_var, main_experimental_condition, aoi_segment_var))
  
  combinations <- combn(col_names, 2, simplify = FALSE)
  list_of_plots <- map(combinations, ~ make_plot(.x[1], .x[2]))
  plot_combined = wrap_plots(list_of_plots, guides = "collect")
  
  pdf(file.path(output_dir, "cell_type_abundance.pdf"), width = 10, height = 6)
  print(plot_combined)
  dev.off()
  

}



# cellChat probability prediction function

cellchat_predict_prob <- function(metadt_all, expr_deseq2_norm_log, DE_genes = NULL, DEG_info = NULL, sample_name, grouping_var_col_ids, group, thresh_fc = 0, thresh_p = 0.05, min_cells = 10, output_dir){
  

  
  expr_deseq2_norm_log_filtered_list <- filter_expr_deseq2_norm_log(
    metadt_all,
    expr_deseq2_norm_log,
    grouping_var_col_ids,
    group
  )
  
  expr_deseq2_norm_log_filtered = expr_deseq2_norm_log_filtered_list$expr_deseq2_norm_log_filtered
  meta = expr_deseq2_norm_log_filtered_list$meta_filtered
  meta$samples = meta[,sample_name] # for cellChat
  meta$samples <- as.factor(meta$samples)
  
  # check min number of cells else stop the execution
  # # Build grouping variable names
  group_vars <- c("labels",grouping_var_col_ids)  # always group by labels
  
  
  df_grouped_samples <- meta %>%
    group_by(!!!syms(group_vars)) %>%
    summarise(count = n(), .groups = "drop")
  


  low_count_labels <- unique(df_grouped_samples[df_grouped_samples$count < min_cells,]$labels)
 
  # TODO cell type abundance plot: error in the plotting function
  # cell_type_abundance_plot(metadt_all, min_cells, aoi_segment_var, main_experimental_condition, grouping_var_col_ids, output_dir)


  if (length(low_count_labels) > 0) {

    warning("The following cell types have fewer than 10 cells: ", paste(low_count_labels, collapse = ", "))
    stop("Stopping execution due to low cell count.check the cell type abundance plots and set the min_cells counts")
  }

  
  ##########    CellChat Analysis
  
  
  # create cellChat object
  
  cellchat_obj = createCellChat(object = expr_deseq2_norm_log_filtered, meta = meta, group.by = "labels")
  
  # Import cellChat object
  
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  #showDatabaseCategory(CellChatDB)
  CellChatDB.use <- subsetDB(CellChatDB) # use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
  cellchat_obj@DB <- CellChatDB.use
  cellchat_obj <- subsetData(cellchat_obj) # This step is necessary even if using the whole database
  
  
  
  # Preprocessing the expression data for cell-cell communication analysis
  # Identify over-expressed signaling genes associated with each cell group
  
  # TODO check future::plan("multisession", workers = 4) 
  
 
  
  if (is.null(DE_genes)) {
    
    # TODO check both DE_genes and DEG_info does not exsist else error !
    
    cellchat_obj <- identifyOverExpressedGenes(cellchat_obj, thresh.fc = thresh_fc, thresh.p = thresh_p) # default thresholds thresh.fc = 0, thresh.p = 0.05
  } else {
    
    # TODO check both DE_genes and DEG_info exsist else error !
    
    # to provide DEGs externally
    cellchat_obj@var.features[["features"]] = DE_genes
    cellchat_obj@var.features[[paste0("features", ".info")]] = DEG_info
  }
  
  
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  
  
  ##  Part II: Inference of cell-cell communication network
  
  cellchat_obj <- computeCommunProb(cellchat_obj, type = "triMean")
  
  
  # compute communication probabilities
  
  
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  cellchat_obj = filterCommunication(cellchat_obj, min.cells = min_cells) 
  df.net <- subsetCommunication(cellchat_obj) # check !!!
  
  
  return(list(
    cell_type_abundance_plot = plot,
    cellchat_obj = cellchat_obj,
    df.net = df.net
  ))
  
  
}








      
