# cellChat Functions

# CellChat Preprocessing

# extract and combine expression data and meta data (Used by both CellChat and MultiNicheNet)

extract_and_combine_expression_data <- function(cell_types, bprism_res, ct_names, geomx_obj){
  
  
  
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
    
    new_meta <- data.frame(meta_data[,c(2,5,6,7,24,25,28)])
    
    new_meta$dcc_filename = gsub('\\.dcc', paste0('_', ct_name), new_meta$dcc_filename)
    mat = deconv_ct_list_int[[ct_name]]
    
    matched_entries = new_meta[new_meta$dcc_filename %in% rownames(mat), ]
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



# filter expression data based on NACT_status, segment and Annotations

filter_expr_deseq2_norm_log = function(metadt_all,expr_deseq2_norm_log,nact_status,segment,annotation){
  
  meta = data.frame(samples = metadt_all$Sample, 
                    NACT_status = metadt_all$NACT_status, 
                    Patient = metadt_all$Patient, 
                    Segment = metadt_all$Segment , 
                    Annotation = metadt_all$Annotation_cell, 
                    row.names = rownames(metadt_all)
                    )
  
  
  if (!is.null(nact_status)) {
    # TODO cehck the legth of the list 
    
    meta <- meta %>% filter(NACT_status == nact_status)
  }

  if (!is.null(segment)) {
    meta <- meta %>% filter(Segment == segment)
  }

  if (!is.null(annotation)) {
    meta <- meta %>% filter(Annotation == annotation)
  }


  
  cell_types <- sapply(strsplit(rownames(meta), "_"), function(x) x[2])
  meta = data.frame(labels = cell_types, meta)
  meta$samples <- as.factor(meta$samples)
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


cell_type_abundance_plot = function(metadt_all, min_cells){
  

  meta = data.frame(samples = metadt_all$Sample,
                    NACT_status = metadt_all$NACT_status,
                    Patient = metadt_all$Patient,
                    Segment = metadt_all$Segment ,
                    Annotation = metadt_all$Annotation_cell,
                    row.names = rownames(metadt_all))

  cell_types <- sapply(strsplit(rownames(meta), "_"), function(x) x[2])
  meta = data.frame(labels = cell_types, meta)
  
  df_grouped_samples <- meta %>%
    group_by(Segment, NACT_status, Annotation, labels) %>%
    summarise(count = n(), .groups = "drop")
  
  # Define a custom ggplot function to avoid repetition
  create_bar_plot <- function(data, x, facet_formula) {
    ggplot(data, aes_string(x = x, y = "count", fill = "labels")) +
    # TODO check ggplot(data, aes(x = .data[[var1]], y = .data[[var2]])) + geom_point()  
      geom_bar(stat = "identity", position = "dodge") +
      facet_grid(facet_formula) +
      theme_minimal() +
      geom_hline(yintercept = min_cells, linetype = "solid", color = "black") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(angle = 0)
      )
  }
  
  # Create individual plots
  p1 <- create_bar_plot(df_grouped_samples, "Annotation", ". ~ Segment")
  p2 <- create_bar_plot(df_grouped_samples, "Annotation", ". ~ NACT_status")
  p3 <- create_bar_plot(df_grouped_samples, "NACT_status", ". ~ Segment")
  
  plot_combined = patchwork::wrap_plots(
    p1,p2,p3,
    nrow = 2, guides = "collect",
    widths = c(6,6)
  )
  
  
  
  return(list(
    plot_combined = plot_combined,
    df_grouped_samples = df_grouped_samples

  ))
  

  
}


# cellChat probability prediction function

cellchat_predict_prob <- function(metadt_all, expr_deseq2_norm_log, DE_genes = NULL, DEG_info = NULL, nact_status, segment, annotation, thresh_fc = 0, thresh_p = 0.05, min_cells = 10, output_dir){
  
  
  parts <- c(nact_status, segment, annotation)
  parts_non_NA <- parts[!sapply(parts, is.na)]
  file_name = paste(parts_non_NA, collapse = "_")
  
  expr_deseq2_norm_log_filtered_list <- filter_expr_deseq2_norm_log(
    metadt_all,
    expr_deseq2_norm_log,
    nact_status,
    segment,
    annotation 
  )
  
  expr_deseq2_norm_log_filtered = expr_deseq2_norm_log_filtered_list$expr_deseq2_norm_log_filtered
  meta = expr_deseq2_norm_log_filtered_list$meta_filtered
  
  # check min number of cells else stop the execution
  
  
  # Build grouping variable names
  group_vars <- c("labels")  # always group by labels
  
  if (!is.null(segment)) group_vars <- c("Segment", group_vars)
  if (!is.null(nact_status)) group_vars <- c("NACT_status", group_vars)
  if (!is.null(annotation)) group_vars <- c("Annotation", group_vars)
  
  # Now use across + all_of
  df_grouped_samples <- meta %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(count = n(), .groups = "drop")
  
  
  low_count_labels <- unique(df_grouped_samples[df_grouped_samples$count < min_cells,]$labels)
  
  # cell type abundance plot
  
  plot = cell_type_abundance_plot(metadt_all,min_cells)
  
  pdf(file.path(output_dir, "cell_type_abundance.pdf"), width = 10, height = 6)
  print(plot$plot_combined)
  dev.off()
  
  
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
  } 
  
  else {
    
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
