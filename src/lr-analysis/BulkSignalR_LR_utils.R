

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

# function to run BulkSignaR LR prediction 
BulkSignalR_LR_prediction <- function(count_geomx, meta_data, normalize_needed, norm_is_log,
                                     min_expr_counts, null_model,
                                     qval_threshold, group, output_dir){

  ### Analysis : written according to the BulkSignalR Vignette
  # browseVignettes("BulkSignalR")
  
  # step 01 : Prepare Dataset
  bsrdm <- BSRDataModel(counts = count_geomx, 
                        normalize = normalize_needed , 
                        method = "UQ", 
                        UQ.pc = 0.75,
                        log.transformed = norm_is_log, 
                        min.count = min_expr_counts, 
                        prop = 0.05) 
  
  # step 02 : learnParameters
  set.seed(123)
  
  bsrdm <- learnParameters(bsrdm,
                           null.model = null_model,
                           plot.folder = file.path(output_dir), 
                           filename = paste0("geomxUQ_",group), 
                           verbose = TRUE)
  

  # step 03 : Building a BSRInference object
  bsrinf <- BSRInference(bsrdm, reference="REACTOME-GOBP")
  
  # reducing to best pathways
  bsrinf.redBP <- reduceToBestPathway(bsrinf)

  # reducing to pathways and to  best pathways 
  bsrinf.redP <- reduceToPathway(bsrinf)
  bsrinf.redPBP <- reduceToBestPathway(bsrinf.redP)
  
  # reducing to ligands and receptors and best bathways
  bsrinf.redL <- reduceToLigand(bsrinf)
  bsrinf.redLBP <- reduceToBestPathway(bsrinf.redL)
  bsrinf.redR <- reduceToReceptor(bsrinf)
  bsrinf.redRBP <- reduceToBestPathway(bsrinf.redR)

  BSR_all <- list(
    bsrdm = bsrdm,
    bsrinf_all = bsrinf,
    bsrinf_redP = bsrinf.redP,
    bsrinf_redPBP = bsrinf.redPBP,
    bsrinf_redBP = bsrinf.redBP,
    bsrinf_redL = bsrinf.redL,
    bsrinf_redR = bsrinf.redR,
    bsrinf_redLBP = bsrinf.redLBP,
    bsrinf_redRBP = bsrinf.redRBP,
    count_geomx = count_geomx,
    meta_data = meta_data
  )
  
  saveRDS(BSR_all, file.path(output_dir, paste0('BSR_results_', group, '.RDS')))

  return(BSR_all)
}



# functions for plotting in BulkSignalr


# Creating a function to plot the heatmap

plot_heatmap <- function(bsrinf_red, bsrdm, reduction_name, meta_data, pathway_names = NULL, 
                         qval_threshold=0.001, top_n=50,
                         heatmap_col_ann, aoi_id = 'dcc_filename',
                         out_path){

  # step 4 : Building a BSRSignature object
  bsrsig_red <- BSRSignature(bsrinf_red, qval.thres = qval_threshold)
  
  # different name.by.pathway param dependig on reduction type
  scoresLR <- scoreLRGeneSignatures(bsrdm, bsrsig_red, 
                                    name.by.pathway=ifelse(reduction_name == 'redPBP', TRUE, FALSE))
  
  # filtering LR interactions df
  LRinter <- LRinter(bsrinf_red) %>%
    filter(qval <= qval_threshold) %>%
    mutate(lr_inter = paste0("{",L,"} / {",R,"}"))
  
  # fix names for different reduction types
  LRinter$lr_inter <- gsub('{{', '{', LRinter$lr_inter, fixed = T)
  LRinter$lr_inter <- gsub('}}', '}', LRinter$lr_inter, fixed = T)
  
  if(!is.null(pathway_names)){
    LRinter <- filter(LRinter, pw.name %in% pathway_names)
  }
  
  if(nrow(LRinter) < 2){
    if(nrow(LRinter) > 0){
      print('only 1 pathway found. try add more pathway names')
      print(LRinter$pw.name)
      return(NULL)
    }
    print('no pairs selected. try changing qval thr or add more pathway names')
    return(NULL)
  }
    
  LRinter_top <- LRinter %>%
    arrange(desc(LR.corr)) %>%
    head(top_n)
  
  # Filter LRscores to match top pairs or pathways
  if(reduction_name == 'redPBP'){
    scoresLR_top <- scoresLR[rownames(scoresLR) %in% LRinter_top$pw.name, ]
    scoresLR_top <- scoresLR_top[LRinter_top$pw.name, ] # ensure ordering
    row_annot <- NA
    colors_pw <- NA
  } else{
    scoresLR_top <- scoresLR[rownames(scoresLR) %in% LRinter_top$lr_inter, ]
    scoresLR_top <- scoresLR_top[LRinter_top$lr_inter, ] # ensure ordering
    
    # row annotations based on pathways
    row_annot <- data.frame(Pathway = LRinter_top$pw.name)
    rownames(row_annot) <- rownames(scoresLR_top)  
    
    # colors for col annotations
    pathways_anno <- as.character(unique(as.character(unique(row_annot$Pathway))))
    colors_pw <- hue_pal()(length(pathways_anno))
    colors_pw <- setNames(colors_pw, pathways_anno)
  }


  # column annotations
  col_annot <- data.frame(SampleGroup = factor(meta_data[,heatmap_col_ann]), row.names = meta_data[, aoi_id])

  # TODO col_annot <- data.frame(SampleGroup = factor(paste0(meta_data$Segment,"_",meta_data$NACT_status))) # if both
  
  # colors for col annotations

  Sample_Group <- as.character(unique(col_annot$SampleGroup))
  colors_SGroup <- hue_pal()(length(Sample_Group))
  colors_SGroup <- setNames(colors_SGroup, Sample_Group)
  
  ann_colors <- list(
    SampleGroup = colors_SGroup,
    Pathway = colors_pw)
  
  # TODO viridisLite::viridis() or colorspace::qualitative_hcl()  # better colors
  
  # Finally the heatmap
  col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  wid <- ifelse(reduction_name == 'redPBP', 12, 24)
  
  heatmap <- pheatmap(scoresLR_top, 
                      annotation_row = row_annot,
                      annotation_col = col_annot,
                      annotation_colors = ann_colors,
                      fontsize_col = 11,
                      cellheight = 13,
                      fontsize_row = 11,
                      color = col_fun, 
                      main = paste0("Ligand receptor signatures per samples"),
                      name = "signature score",
                      labels_col = rep("", ncol(scoresLR_top)))
  
  pdf(out_path, width = wid, height = 12)  # Width and height in inches
  print(heatmap)
  dev.off()
}



