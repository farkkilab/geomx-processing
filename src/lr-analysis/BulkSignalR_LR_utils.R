

run_unless_exists <- function(step_name, expected_output, script){
  if(!file.exists(expected_output)){
    source(script, local = TRUE)
    gc()
    print('$$$$$$$$$$')
    print(paste0(step_name, ' succeeded!'))
  } else{
    print(paste0(step_name, ' have been already run'))
  }
}



# return a list of paired_samples

filter_paired_data = function(geomx_obj, main_experimental_condition, paired_id){
  meta_data = sData(geomx_obj)
  meta_data = meta_data %>% select(!!sym(aoi_id), !!sym(sample_name), !!sym(paired_id), !!sym(main_experimental_condition)) 
  
  
  main_experimental_condition_groups = unique(meta_data[,main_experimental_condition])
  
  paired_patients <- meta_data %>%
       group_by(!!sym(paired_id)) %>%
       filter(all(main_experimental_condition_groups %in% !!sym(main_experimental_condition))) %>%
       summarise() %>%
       pull(!!sym(paired_id))
  

  paired_samples = unique(meta_data %>% filter(!!sym(paired_id) %in% paired_patients) %>% pull(!!sym(sample_name)))
  return(paired_samples)
  
}


# fuction to filter data for BulkSignalR prediction
# TODO simplify by passing metadata colnames as 1 argument
Filter_for_BulkSignaR_LR_prediction <- function(geomx_obj, aoi_id, sample_name, aoi_segment_var, main_experimental_condition, grouping_var_col_ids, count_geomx,group,paired_only, paired_id = NULL){
  
    meta_data = sData(geomx_obj)
    meta_data = meta_data %>% select(!!sym(aoi_id), !!sym(sample_name), !!sym(aoi_segment_var), !!sym(main_experimental_condition), all_of(grouping_var_col_ids)) 
    
    
    # TODO check this code
    # if (paired_only == TRUE){
    #   
    #   paired_samples = filter_paired_data(geomx_obj, main_experimental_condition, paired_id)
    #   meta_data = meta_data %>% filter(!!sym(sample_name) %in% paired_samples)
    # 
    #   
    # }  
    
    
    groups <- strsplit(group, "_")[[1]]
    print(groups)
    
    for (column_name in grouping_var_col_ids){
      
      meta_data = meta_data %>% filter(!!sym(column_name) %in% groups)
      
    }
    
   
  meta_data[,aoi_id] = gsub('-', '.', meta_data[,aoi_id])
  col_ids = colnames(count_geomx)  %in% meta_data$dcc_filename
  count_geomx = count_geomx[,col_ids]
  
  return(list(count_geomx = count_geomx,
              meta_data = meta_data
              ))
  
}


# function to run BulkSignaR LR prediction 

BulkSignaR_LR_prediction <- function(count_geomx_list, normalize_needed, normalize_method, UQ_pc, output_dir, qval_threshold, group, null_model){
  

  
  if(is.null(group)){
    plot_name = "combined"
  } else{
    plot_name = group

  }
  
  count_geomx = count_geomx_list$count_geomx
  meta_data = count_geomx_list$meta_data
  
  ### Analysis : written according to the BulkSignalR Vignette
  # browseVignettes("BulkSignalR")
  
  # step 01 : Prepare Dataset
  # bsrdm <- BSRDataModel(speSubset,
  #                       min.count = 1,
  #                       prop = 0.01,
  #                       method = "TC",
  #                       symbol.col = 2,
  #                       x.col = 4,
  #                       y.col = 5, 
  #                       barcodeID.col = 1)
  
  # TODO this function cannot be found :0
  #bsrdm <- prepareDataset(counts = count_geomx, normalize = normalize_needed , method = normalize_method, UQ.pc = UQ_pc, log.transformed = FALSE, min.count = 10, prop = 0.1) 
  bsrdm <- BSRDataModel(counts = count_geomx, normalize = normalize_needed , method = normalize_method, UQ.pc = UQ_pc, log.transformed = FALSE, min.count = 10, prop = 0.1) 
  
  # step 02 : learnParameters

  set.seed(123)
  
  if(is.null(null_model)){
    bsrdm <- learnParameters(bsrdm, 
                             plot.folder = file.path(output_dir), 
                             filename = paste0("geomxUQ_",plot_name), 
                             verbose = TRUE
                             )
  } else{
    
    bsrdm <- learnParameters(bsrdm,
                             null.model = null_model,
                             plot.folder = file.path(output_dir), 
                             filename = paste0("geomxUQ_",plot_name), 
                             verbose = TRUE
                             )
    
    }
  
  # step 03 : Building a BSRInference object
  
  bsrinf <- BSRInference(bsrdm, min.cor = 0.3, reference="REACTOME-GOBP")
  LRinter.dataframe <- LRinter(bsrinf)
  LRinter.dataframe <- LRinter.dataframe[order(LRinter.dataframe$qval <= qval_threshold),]
  
  # reducing to best pathways before calculating signature scores
  
  bsrinf.redBP    <- reduceToBestPathway(bsrinf) 
  LRinter_pairs_best_pws = LRinter(bsrinf.redBP)
  
  
  return(list(
    bsrdm = bsrdm,
    bsrinf = bsrinf,
    bsrinf_redBP = bsrinf.redBP,
    LRinter_dataframe = LRinter.dataframe,
    LRinter_pairs_best_pws = LRinter_pairs_best_pws,
    count_geomx = count_geomx,
    meta_data = meta_data
    ))

}



# functions for plotting in BulkSignalr


# Creating a function to plot the heatmap

plot_heatmap <- function(bsrinf_redBP, bsrdm, meta_data, pathway_names = NULL, qval_threshold=0.01, n=50, heatmap_col_ann){
  
  pairs <- LRinter(bsrinf_redBP)
  pairs$index <- seq_len(nrow(pairs))
  selected_pairs <- filter(pairs, qval < qval_threshold) 
  
  if(!is.null(pathway_names)){
    selected_pairs <- filter(selected_pairs, pw.name  %in% pathway_names)
  }
  
  if(nrow(selected_pairs) == 0){
    print('no pairs selected. try changing qval thr or add more pathway names')
    return(NULL)
  }
  
  # TODO if selected_pairs == 0 then an error message
  
  top_n_pairs <- arrange(selected_pairs, desc(LR.corr))
  top_n_pairs = head(top_n_pairs, n)

  top_n_pairs_index <- top_n_pairs$index
  
  
  ligands   <- ligands(bsrinf_redBP)[top_n_pairs_index]
  receptors   <- receptors(bsrinf_redBP)[top_n_pairs_index]
  pathways  <- pairs$pw.name
  t.genes   <- tgGenes(bsrinf_redBP)[top_n_pairs_index]
  t.corrs   <- tgCorr(bsrinf_redBP)[top_n_pairs_index]

  for (i in seq_len(nrow(top_n_pairs))){
    tg <- t.genes[[i]]
    t.genes[[i]] <- tg[top_n_pairs$rank[i]:length(tg)]
    
    tc <- t.corrs[[i]]
    t.corrs[[i]] <- tc[top_n_pairs$rank[i]:length(tc)]
  }
  
  bsrinf_redBP@LRinter = top_n_pairs
  bsrinf_redBP@ligands = ligands
  bsrinf_redBP@receptors = receptors
  bsrinf_redBP@tg.genes = t.genes
  bsrinf_redBP@tg.corr = t.corrs
  
  
  # step 4 : Building a BSRSignature object
  
  bsrsig.redBP <- BSRSignature(bsrinf_redBP, qval.thres = qval_threshold)
  
  scoresLR <- scoreLRGeneSignatures(bsrdm, bsrsig.redBP,
                                    name.by.pathway=FALSE)
  
  
  LRinter = bsrinf_redBP@LRinter # check !!
  LRinter$lr_inter = paste0("{",LRinter$L,"} / {",LRinter$R,"}")
  
  # Filter the dataframe based on matches with the column names
  
  row_names <- rownames(scoresLR)
  matched_pw_names <- LRinter %>%
    filter(lr_inter %in% row_names) %>%
    arrange(match(lr_inter, row_names)) %>%  # Preserve the column order
    pull(pw.name)
  
  # row annotations based on pathways
  
  row_annot <- data.frame(Pathway = matched_pw_names)
  rownames(row_annot) <- rownames(scoresLR)  # ensure they align with row_annot 
  
  
  # column annotations
  
  col_annot <- data.frame(SampleGroup = factor(meta_data[,heatmap_col_ann]))

  # if both 
  # col_annot <- data.frame(SampleGroup = factor(paste0(meta_data$Segment,"_",meta_data$NACT_status)))
  
  rownames(col_annot) <- gsub('-', '.', meta_data$dcc_filename)
  
  # colors for the row and col annotations
  
  pathways_anno <- as.character(unique(matched_pw_names)) # Convert to character vector (if not already)
  Sample_Group <- as.character(unique(col_annot$SampleGroup))
  
  colors_pw <- hue_pal()(length(pathways_anno))
  colors_pw <- setNames(colors_pw, pathways_anno)
  colors_SGroup <- hue_pal()(length(Sample_Group))
  colors_SGroup <- setNames(colors_SGroup, Sample_Group)
  
  ann_colors <- list(
    SampleGroup = colors_SGroup,
    Pathway = colors_pw)
  
  # TO DO
  # viridisLite::viridis() or colorspace::qualitative_hcl() better colors
  
  # Finally the heatmap
  
  col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  
  heatmap <- pheatmap(scoresLR, 
                      annotation_row = row_annot,
                      annotation_col = col_annot,
                      annotation_colors = ann_colors,
                      fontsize_col = 11,
                      cellheight = 13,
                      fontsize_row = 11,
                      color = col_fun, 
                      main = paste0("Ligand receptor signatures per samples"),
                      name = "signature score",
                      labels_col = rep("", ncol(scoresLR)))
  
  
  return(heatmap)
  
  
}



