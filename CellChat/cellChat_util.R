# cellChat Functions

# CellChat Preprocessing


filter_expr_deseq2_norm_log = function(metadt_all,expr_deseq2_norm_log,nact_status,segment,annotation){
  
  meta = data.frame(samples = metadt_all$Sample, 
                    NACT_status = metadt_all$NACT_status, 
                    Patient = metadt_all$Patient, 
                    Segment = metadt_all$Segment , 
                    Annotation = metadt_all$Annotation_cell, 
                    row.names = rownames(metadt_all)
                    )
  
  
  if (!is.null(nact_status)) {
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


cell_type_abundance_plot = function(meta, min_cells){
  
  
  meta = data.frame(samples = metadt_all$Sample, 
                    NACT_status = metadt_all$NACT_status, 
                    Patient = metadt_all$Patient, 
                    Segment = metadt_all$Segment , 
                    Annotation = metadt_all$Annotation_cell, 
                    row.names = rownames(metadt_all))
  
  meta = data.frame(labels = cell_types, meta)
  
  df_grouped_samples <- meta %>%
    group_by(Segment, NACT_status, Annotation, labels) %>%
    summarise(count = n(), .groups = "drop")
  
  # Define a custom ggplot function to avoid repetition
  create_bar_plot <- function(data, x, facet_formula) {
    ggplot(data, aes_string(x = x, y = "count", fill = "labels")) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_grid(facet_formula) +
      theme_minimal() +
      geom_hline(yintercept = 10, linetype = "solid", color = "black") +
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
  
  
  #low_count_labels <- unique(df_grouped_samples[df_grouped_samples$count < min_cells,]$labels)

  return(plot_combined)
  
}


# cellChat probability prediction function

cellchat_predict_prob <- function(metadt_all,expr_deseq2_norm_log,nact_status,segment,annotation){
  
  expr_deseq2_norm_log_filtered_list <- filter_expr_deseq2_norm_log(
    metadt_all,
    expr_deseq2_norm_log,
    nact_status,
    segment,
    annotation 
  )
  
  expr_deseq2_norm_log_filtered = expr_deseq2_norm_log_filtered_list$expr_deseq2_norm_log_filtered
  meta = expr_deseq2_norm_log_filtered_list$meta_filtered
  ### CellChat Analysis
  
  
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
  
  
  future::plan("multisession", workers = 4) # do parallel
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj) # thresh.fc = 0, thresh.p = 0.05 fold chage and pvalue: Can use other thresholds
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  
  
  ##  Part II: Inference of cell-cell communication network
  
  cellchat_obj <- computeCommunProb(cellchat_obj, type = "triMean")
  
  
  # compute communication probabilities
  
  
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  cellchat_obj = filterCommunication(cellchat_obj, min.cells = min_cells) 
  df.net <- subsetCommunication(cellchat_obj) # check !!!
  
  
  return(list(
    cellchat_obj = cellchat_obj,
    df.net = df.net
  ))
  
  
}
