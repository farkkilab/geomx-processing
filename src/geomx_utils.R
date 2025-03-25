####################################################
# Description: 
#   Executes script if a specified file (expected_output) does not already exist.
# Parameters:
#   step_name: (string) The name of the current step, used for printing a success message.
#   expected_output: (string) The path to the expected output file that determines whether the script should be run.
#   script: (string) Path to the script file to be sourced and executed if the expected output does not exist.
# Return value:
#   None. The function will print status messages to the console.
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

##############################################
# Description:
#   Generates a Sankey plot from provided data, grouping by specified variables, and saves the plot to a file.
# Parameters:
#   data: (data.frame) The dataset to be used for plotting.
#   variables_to_plot: (vector of strings) The names of the columns in data to group by in the Sankey plot.
#   fill_var: (string) The name of the column in data used to fill the colors in the plot.
#   output_name: (string) The path to the output file where the plot will be saved.
# Return value:
#   None. The function saves the generated plot to output_name.
plot_sankey <- function(data, variables_to_plot, fill_var, output_name){
  count_mat <-   data %>%
    group_by_at(variables_to_plot) %>% 
    summarise(n = n())
  
  test_gr <- gather_set_data(count_mat, 1:length(variables_to_plot))
  
  
  test_gr$x <- mapvalues(test_gr$x, 
                         from=seq(1:length(variables_to_plot)), 
                         to=variables_to_plot)
  test_gr$x <- factor(test_gr$x,
                      levels = variables_to_plot)
  
  # plot Sankey
  ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
    geom_parallel_sets(aes(fill = get(fill_var)), alpha = 0.5, axis.width = 0.1) +
    geom_parallel_sets_axes(axis.width = 0.2) +
    geom_parallel_sets_labels(color = "white", size = 3) +
    theme_classic(base_size = 17) + 
    theme(legend.position = "bottom",
          axis.ticks.y = element_blank(),
          axis.line = element_blank(),
          axis.text.y = element_blank()) +
    scale_y_continuous(expand = expansion(0)) + 
    scale_x_discrete(expand = expansion(0)) +
    labs(x = "", y = "") +
    annotate(geom = "segment", x = 4.25, xend = 4.25,
             y = 20, yend = 120, lwd = 2) +
    annotate(geom = "text", x = 4.19, y = 70, angle = 90, size = 5,
             hjust = 0.5, label = "100 segments")
  
  ggsave(output_name, width = 2000, height = 2000, unit = 'px')
}

###############################################
# Description:
#   Summarizes the quality control (QC) results, creating a summary table with pass/fail statistics.
# Parameters:
#   QCResults: (data.frame) A dataframe containing results from a QC process, with logical values indicating pass or fail statuses.
# Return value:
#   (data.frame) A summary dataframe with columns indicating the number of passes and warnings for the QC checks.
qc_summarize <- function(QCResults){
  QC_Summary <- data.frame(Pass = colSums(!QCResults[, colnames(QCResults)]),
                           Warning = colSums(QCResults[, colnames(QCResults)]))
  
  QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
    ifelse(sum(x) == 0L, "PASS", "WARNING")
  })
  
  QC_Summary["TOTAL FLAGS", ] <-
    c(sum(QCResults[, "QCStatus"] == "PASS"),
      sum(QCResults[, "QCStatus"] == "WARNING"))
  
  return(QC_Summary)
}

################################################
# Description:
#   Generates and saves a histogram of QC statistics with optional threshold line and scale transformation.
# Parameters:
#   assay_data: (data.frame) Dataset containing assay results to be plotted.
#   annotation: (string) Name of the plot and x-axis.
#   fill_by: (string) Column name used for filling colors in the histogram.
#   thr: (numeric) The threshold value to be marked as a vertical line in the plot.
#   scale_trans: (string) The type of scale transformation to apply to the x-axis (eg log)
#   output_name: (string) Path to the file where the plot will be saved.
# Return value:
#   None. The function saves the plot to output_name.
QC_histogram <- function(assay_data,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL,
                         output_name = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
  
  print(output_name)
  ggsave(output_name, width = 2000, height = 1000, unit='px')
}

#######################################
# Description:
#   Creates and saves a stacked bar plot of nr of AOI above defined gene detection rate thresholds.
# Parameters:
#   segment_data: (data.frame) Dataframe containing data for each AOI 
#                 (typically pData(geomx_obj)) with a column GeneDetectionRate.
#   fill_var: (string) Column in the segment_data used for filling the plot colors.
#   output_name: (string) Path where the output plot will be saved.
# Return value:
#   None. The function saves the plot as specified in output_name.
plot_detection_rate <- function(segment_data, fill_var, output_name){
  segment_data$DetectionThreshold <- 
    cut(segment_data$GeneDetectionRate,
        breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
        labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))
  
  # stacked bar plot of different cut points (1%, 5%, 10%, 15%)
  ggplot(segment_data,
         aes(x = DetectionThreshold)) +
    geom_bar(aes(fill = get(fill_var))) +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Gene Detection Rate",
         y = "Segments, #",
         fill = "Segment Type")
  
  ggsave(output_name, width = 2000, height = 2000, unit='px')
}

############################################
# Description:
#   Generates a bar plot showing how many genes were detected in given % of AOIs
# Parameters:
#   gene_data: (data.frame) Dataframe containing gene detection rate data (typically fData(geomx_obj))
#   output_name: (string) The name of the file where the plot will be saved.
# Return value:
#   None. Saves the bar plot to the file specified by output_name.
plot_gene_detection_rate <- function(gene_data, output_name){
  plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
  plot_detect$Number <-
    unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                  function(x) {sum(gene_data$DetectionRate >= x)}))
  plot_detect$Rate <- plot_detect$Number / nrow(gene_data)
  rownames(plot_detect) <- plot_detect$Freq
  
  ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
              vjust = 1.6, color = "black", size = 4) +
    scale_fill_gradient2(low = "orange2", mid = "lightblue",
                         high = "dodgerblue3", midpoint = 0.65,
                         limits = c(0,1),
                         labels = scales::percent) +
    theme_bw() +
    scale_y_continuous(labels = scales::percent, limits = c(0,1),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "% of Segments",
         y = "Genes Detected, % of Panel > LOQ")
  
  ggsave(output_name, width = 2000, height = 1500, unit='px')
}

#######################################################
# plot_q3_stats
# Description:
#   Creates multiple plots of Q3 statistics from the given geomx_obj object and saves them to a single file.
# Parameters:
#   geomx_obj: (S4 object) A GeoMx object containing expression data.
#   ann_of_interest: (string) The annotation column to focus on in the plot.
#   output_name: (string) The path to the file where the combined plot will be saved.
# Return value:
#   None. The function saves plots to the file specified by output_name.
plot_q3_stats <- function(geomx_obj, ann_of_interest, output_name){
  
  negativeProbefData <- subset(fData(geomx_obj), CodeClass == "Negative") # 1 bcs already collapsed to targets
  neg_probes <- unique(negativeProbefData$TargetName)
  
  Stat_data <- 
    data.frame(row.names = colnames(exprs(geomx_obj)),
               Segment = colnames(exprs(geomx_obj)),
               Annotation = pData(geomx_obj)[, ann_of_interest],
               Q3 = unlist(apply(exprs(geomx_obj), 2,
                                 quantile, 0.75, na.rm = TRUE)),
               NegProbe = exprs(geomx_obj)[neg_probes, ])
  
  Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                      variable.name = "Statistic", value.name = "Value")
  
  plt1 <- ggplot(Stat_data_m,
                 aes(x = Value, fill = Statistic)) +
    geom_histogram(bins = 40) + theme_bw() +
    scale_x_continuous(trans = "log2") +
    facet_wrap(~Annotation, nrow = 1) + 
    scale_fill_brewer(palette = 3, type = "qual") +
    labs(x = "Counts", y = "Segments, #")
  
  plt2 <- ggplot(Stat_data,
                 aes(x = NegProbe, y = Q3, color = Annotation)) +
    geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
    geom_point() + guides(color = "none") + theme_bw() +
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2") +
    theme(aspect.ratio = 1) +
    labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")
  
  plt3 <- ggplot(Stat_data,
                 aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
    geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
    geom_point() + theme_bw() +
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2") +
    theme(aspect.ratio = 1) +
    labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")
  
  btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                       rel_widths = c(0.43,0.57))
  plt_all <- plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))
  
  ggsave(output_name, width=2000, height=1500, unit='px')
}

###########################################################
# Description:
#   Plots and saves a boxplot visualizing counts values for first 10 AOI to visualise 
#   the effects of normalisation
# Parameters:
#   expr_data: (matrix) Expression matrix of raw/normalised data.
#   norm_name: (string) Name of the normalization method or dataset.
#   output_name: (string) Path where the output plot will be saved.
#   is_log: (logical) A flag indicating if the data are in log scale.
# Return value:
#   None. The function outputs the plot to the specified file.
plot_norm_effect <- function(expr_data, norm_name, output_name, is_log = F){
  png(filename=output_name, width=1000, height=750, units="px")
  
  if(!is_log){
    boxplot(expr_data,
            col = "#9EDAE5", main = norm_name,
            log='y', names = seq(1:ncol(expr_data)), xlab = "Segment",
            ylab = norm_name)
  } else{
    boxplot(expr_data,
            col = "#9EDAE5", main = norm_name,
            names = seq(1:ncol(expr_data)), xlab = "Segment",
            ylab = norm_name)
  }

  
  dev.off()
}

###########################################################
# Description:
#   Visualizes the distribution of expression data and optionally its log-transformed values.
# Parameters:
#   expr_data: (matrix) Expression data to be plotted.
#   norm_name: (string) Label for the dataset (eg normalisation type), used in the plot title.
#   output_name: (string) Path to the file where the plots will be saved.
#   is_log: (logical) Boolean indicating if data is already in log scale. 
#   if not - the log-transformed plot will also be produced
# Return value:
#   None. The plots are saved using ggsave().
plot_expr_distribution <- function(expr_data, norm_name, output_name, is_log = F){
  
  expr_df <- as.data.frame(as.vector(expr_data))
  colnames(expr_df) <- 'expr'
  
  ggplot(data = expr_df) +
    geom_histogram(aes(x = expr), bins = 100) +
    xlim(0, as.numeric(quantile(expr_df$expr, probs = 0.99))) +
    ggtitle(paste0(norm_name, ' [.99 percentile]'))
  
  ggsave(output_name)
  
  if(!is_log){
    expr_df$expr_log2 <- log2(expr_df$expr + 1)
    
    ggplot(data = expr_df) +
      geom_histogram(aes(x = expr_log2), bins = 100) +
      xlim(0, as.numeric(quantile(expr_df$expr_log2, probs = 0.99)))+
      ggtitle(paste0('log2 ', norm_name, ' [.99 percentile]'))
    
    ggsave(paste0(file_path_sans_ext(output_name), '_log2.', file_ext(output_name)))

  }
}

###########################################################
# Description:
#   Performs UMAP and t-SNE dimentionality reduction on given expression data and adds results to geomx objest metadata
# Parameters:
#   geomx: (S4 object) GeoMx object with assay data and pData (metadata)
#   assay_name: (string) Name of the assay or expression matrix to use (eg normalisation type) for reduction
#   assay_is_log: (logical) If the expr data is already in log scale.
# Return value:
#   (S4 object) The geomx object  with UMAP and t-SNE results added to metadata
make_umap_tsne <- function(geomx, assay_name, assay_is_log = F){
  
  # set the seed for UMAP
  custom_umap <- umap::umap.defaults
  custom_umap$random_state <- 42
  
  # log2 have to be used if the data are not in the log scale
  if(assay_is_log){
    inp_expr <- assayDataElement(geomx , elt = assay_name)
  } else{
    inp_expr <- log2(assayDataElement(geomx , elt = assay_name))
  }
  
  # make umap
  umap_out <- umap(t(inp_expr), config = custom_umap)
  
  # save UMAP1 and 2 results to pData
  pData(geomx)[, c(paste0("UMAP1_", assay_name), paste0("UMAP2_", assay_name))] <- umap_out$layout[, c(1,2)]
  
  # set the seed for tSNE 
  set.seed(42) 
  
  # make tsne
  tsne_out <- Rtsne(t(inp_expr), perplexity = ncol(geomx)*.15)
  
  # save tSNE1 and 2 results to pData
  pData(geomx)[, c(paste0("tSNE1_", assay_name), paste0("tSNE2_", assay_name))] <- tsne_out$Y[, c(1,2)]
  
  return(geomx)
}

############################################################
# Description:
#   Generates plots for UMAP or t-SNE results from data frame (eg geomx pData())
# Parameters:
#   pheno_data: (data.frame) Dataframe containing UMAP/t-SNE results.
#   method_type: (character vector) Either 'UMAP' or 'tSNE' indicating the dimensionality reduction method.
#   norm_type: (string) input expression matrix name (eg normalisation type) 
#   color_var: (string) Column name for color grouping in the plot.
#   shape_var: (string) Column used for shape grouping in the plot, default is 'Segment'.
#   output_name: (string) The path to the file where the plot will be saved.
# Return value:
#   None. The plot is saved to output_name.
plot_umap_tsne <- function(pheno_data, method_type = c('UMAP', 'tSNE'), 
                           norm_type, color_var, shape_var = 'Segment',
                           output_name){
  
  pheno_data[[color_var]] <- as.character(pheno_data[[color_var]])
  
  ggplot(pheno_data,
         aes(x = get(paste0(method_type, '1_', norm_type)), 
             y = get(paste0(method_type, '2_', norm_type)), 
             color = get(color_var), shape = get(shape_var))) +
    geom_point(size = 3) +
    xlab(paste0(method_type, '1_', norm_type)) +
    ylab(paste0(method_type, '2_', norm_type)) +
    scale_color_discrete(name = color_var) + 
    scale_shape_discrete(name = shape_var) + 
    theme_bw()
  
  ggsave(output_name, width = 2000, height = 1500, unit='px', device='pdf')
}

############################################################
# Description:
#   Plots the proportion of variance explained for each variable in PVCA analysis.
# Parameters:
#   pvca_obj: (list) PVCA object containing the PVCA analysis results returned by pvcaBatchAssess()
#   plot_name: (string) Name used in the file title of the saved plot.
#   output_dir: (string) Directory where the plot will be saved.
# Return value:
#   None. The function saves a bar plot to the specified directory.
plot_pvca <- function(pvca_obj, plot_name, output_dir){
  pvca_dt <- data.frame(effect_name = pvca_obj$label, var = t(pvca_obj$dat))
  pvca_dt$effect_name <- gsub('_factor', '', pvca_dt$effect_name)
  
  pvca_plot <- ggplot(pvca_dt, aes(x = effect_name, y = var)) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6)) +
    geom_text(aes(label = round(var, 3)), angle = 90, vjust = 0.5, hjust = 0, size = 2) +
    ylab('Weighted average proportion variance') +
    xlab('Effects') +
    ylim(0, max(pvca_dt$var)+0.1)
  
  ggsave(file.path(output_dir, paste0('pvca_', plot_name, '.png')))
}

############################################################
# Description:
#   Conducts an Over-Representation Analysis (ORA) for the given set of genes 
#   against a background gene set library.
# Parameters:
#   gene_vect: (vector) List of genes to be tested
#   bcg_gene_vect: (vector) Background gene vector used for comparison.
#   msigdb_df: (data.frame) The gene set database to be used for ORA.
#   padj: (numeric) Adjusted p-value cutoff for the analysis, default 0.1.
# Return value:
#   (data.frame) A dataframe with the results of the ORA.
calculate_ora <- function(gene_vect, bcg_gene_vect, msigdb_df, padj = 0.1){
  ora <- enricher(
    gene = gene_vect,
    pvalueCutoff = padj, # Can choose a FDR cutoff
    pAdjustMethod = "BH", 
    universe = bcg_gene_vect, 
    TERM2GENE = dplyr::select(msigdb_df, gs_name, human_gene_symbol)
  )
  
  if(!is.null(ora)){
    ora_df <- data.frame(ora@result) %>%
      filter(p.adjust <= padj)
  } else{
    ora_df <- data.frame()
  }
  return(ora_df)
}

##########################################################
# have to be used for each data group and contrast separately!!
# remeber to always use Segment as a data grouping variable in DEG
# Description:
#   Creates a volcano plot to visualize Differential Gene Expression (DGE) results.
# Parameters:
#   results: (data.frame) The results of the DGE analysis.
#   plot_name: (string) Title used for the plot.
#   top_n_lab: (integer) Number of top genes to label in the plot.
#   group_pos: (string) Name of the positive comparison group.
#   group_neg: (string) Name of the negative comparison group.
#   output_dir: (string) Directory where the plot file will be saved.
# Return value:
#   None. Saves a volcano plot as a PNG file.
plot_volcano_deg <- function(results, plot_name, top_n_lab, group_pos, group_neg, output_dir){
  # Categorize Results based on P-value & FDR for plotting
  results$Color <- "NS or FC < 0.5"
  results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
  results$Color[results$FDR < 0.05] <- "FDR < 0.05"
  results$Color[results$FDR < 0.001] <- "FDR < 0.001"
  results$Color[abs(results$Estimate) < 0.5] <- "NS or FC < 0.5"
  results$Color <- factor(results$Color,
                          levels = c("NS or FC < 0.5", "P < 0.05",
                                     "FDR < 0.05", "FDR < 0.001"))
  
  # pick top genes for either side of volcano to label
  # order genes for convenience:
  results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)
  # top_g <- c()
  # for(cond in c("tumor", "stroma")) {
  #   ind <- results$Segment == cond
  #   top_g <- c(top_g,
  #              results[ind, 'Gene'][
  #                order(results[ind, 'invert_P'], decreasing = TRUE)[1:top_n_lab]],
  #              results[ind, 'Gene'][
  #                order(results[ind, 'invert_P'], decreasing = FALSE)[1:top_n_lab]])
  # }
  # top_g <- unique(unlist(top_g))
  top_g <- unique(c(results$Gene[
               order(results$invert_P, decreasing = TRUE)[1:top_n_lab]],
             results$Gene[
               order(results$invert_P, decreasing = FALSE)[1:top_n_lab]]))
  
  # Graph results
  volc <- ggplot(results,
                 aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                     color = Color, label = Gene)) +
    geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
    geom_hline(yintercept = -log10(0.05), lty = "dashed") +
    geom_point() +
    labs(x = paste("Enriched in ",  group_neg, " <- log2(FC) -> Enriched in ", group_pos),
         y = "Significance, -log10(P)",
         color = "Significance") +
    scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                  `FDR < 0.05` = "lightblue",
                                  `P < 0.05` = "orange2",
                                  `NS or FC < 0.5` = "gray"),
                       guide = guide_legend(override.aes = list(size = 4))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
    geom_text_repel(data = subset(results, (Gene %in% top_g) & (FDR < 0.05) & (Estimate > 0.5 | Estimate < -0.5)),
                    size = 4, point.padding = 0.15, color = "black",
                    min.segment.length = .1, box.padding = .2, lwd = 2,
                    max.overlaps = 50) +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom") +
    #facet_wrap(~Segment, scales = "fixed") +
    ggtitle(paste(plot_name, group_pos, group_neg))
  
  ggsave(file.path(output_dir, paste0('volc_', plot_name, '_', group_pos, '_', group_neg,  '.png')), width = 4000, height = 2000, unit='px')
}

################################################################
# Description:
#   Converts Ensembl or Entrez gene IDs to gene names for list of vecors of gene sets
# Parameters:
#   gene_inp_list: (list) Nested list of gene IDs to be converted.
#   conv: (character) The type of conversion from 'ens' (Ensembl) or 'entrez'.
#   type: (string) The expected input data structure type, default is 'list'.
# Return value:
#   (list) list of vectors with gene sets with converted gene names.
gene_2names <- function(gene_inp_list, conv = c('ens', 'entrez'), type = 'list'){
  
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  
  #TODO make it for df if needed and also the other way around
  # gene_df_names <- getBM(attributes=c('external_gene_name', 'ensembl_gene_id'),
  #                           filters = 'ensembl_gene_id',
  #                           values = as.character(unlist(gene_list)),
  #                           mart = ensembl)
  # 
  # 
  # 
  # marker_ind <- left_join(marker_ind, marker_ind_names)
  # rm(marker_ind_names)
  # 
  if(conv == 'ens'){
    conv_name <- 'ensembl_gene_id'
  } else if(conv == 'entrez'){
    conv_name <- 'entrezgene_id'
  } else{stop()}
  
  if(type == 'list'){
    
    gene_list <- lapply(gene_inp_list, function(x){
      gene_names <- getBM(attributes=c('external_gene_name', conv_name),
                          filters = conv_name,
                          values = x,
                          mart = ensembl)
      
      gene_names <- gene_names$external_gene_name
    })
    return(gene_list)
  }
}

############################################################
# Description:
#   Creates  a violin plot for pathway scores across groups, with wilcox test statistical annotations.
# Parameters:
#   df: (data.frame) Dataframe containing pathway scores and group information.
#   pathway_colname: (string) Column name for pathways.
#   score_colname: (string) Column name indicating the scores.
#   color_colname: (string) Column used to color the plot.
#   facet_var: (string or list) Variables for facet wrapping.
#   plot_title: (string) Title of the plot.
#   output_path: (string) File path where the plot is saved.
#   ymin: (numeric) Minimum y-axis value for the plot.
#   ymax: (numeric) Maximum y-axis value.
#   manual_colours: (vector) Manual color values for the plot.
# Return value:
#   None. The plot is saved as a PDF.
pathway_boxplot <- function(df, pathway_colname, score_colname, color_colname, facet_var,
                            plot_title, output_path, ymin=-1, ymax=1.4,
                            manual_colours = c("#F8766D", "#00BA38", "#619CFF", "#C77CFF")){
  # per Anno cell type
  gsva_boxpl <- ggplot(data = df, aes(x = get(pathway_colname), y = get(score_colname), fill = get(color_colname))) +
    #geom_boxplot() +
    geom_violin() +
    # geom_point(position= position_jitterdodge(dodge.width = 1, jitter.width= .3, jitter.height = 0),
    #            size= 0.2, alpha = 0.6) +
    stat_summary(fun = "mean", geom = "point", colour = "red", position = position_dodge(0.9), size=0.3) +
    geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = TRUE, size = 0.2, label.size = 2.8) +
    theme(axis.text.x = element_text(angle=45, hjust=1, size = 5)) +
    ggtitle(plot_title)+
    xlab(pathway_colname) +
    ylab(paste0(score_colname)) +
    guides(fill=guide_legend(title=color_colname)) +
    scale_fill_manual(values=manual_colours) #+
    #scale_y_continuous(trans='log10')
    #ylim(ymin, ymax)
  
  if(length(facet_var) == 1){
    gsva_boxpl <- gsva_boxpl +
      facet_wrap(~get(facet_var), scales = "fixed", dir="v")
  } else if(length(facet_var) == 2){
    gsva_boxpl <- gsva_boxpl +
      facet_wrap(get(facet_var[1])~get(facet_var[2]), scales = "fixed", dir="v", nrow=2)
  } else if(length(facet_var) > 2){
    stop('only 1 or 2 variables for facet')
  }
  
  pdf(file= output_path, width=8, height=5)
  plot(gsva_boxpl)
  dev.off()
  #ggsave(output_path, height = 2000, width = 3000, unit = 'px', device='pdf')
}

############################################
#############################################
# Description:
#   Adjusts gene names in a vector (eg scRNAseq dataste) to match those in a given 
#   reference vector (eg in geomx_obj), using synonyms from Ensembl.
# Parameters:
#   geomx_gene_names: (vector) Reference vector of HGSC gene names. 
#     Names in the gene_vector will be adjusted to the ones in geomx_gene_names
#     it doesnt have to be geomx, any vector is fine, but keep it to avoid confusion
#   gene_vector: (vector) Vector of gene names to be adjusted.
# Return value:
#   (vector) Adjusted vector of gene names with synonyms converted to match the reference.
adjust_synonym_genes <- function(geomx_gene_names, gene_vector){
  gene_vector <- unlist(gene_vector)
  geo_non_ex <- setdiff(gene_vector, geomx_gene_names)
  
  if(length(geo_non_ex) > 0){
    ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    geo_non_ex_syn <- getBM(attributes = c('external_gene_name', 'external_synonym'),
                            filters = 'external_gene_name',
                            values = geo_non_ex,
                            mart = ensembl)
    
    
    geo_syn_in_gene_vector <- filter(geo_non_ex_syn, external_synonym %in% geomx_gene_names & 
                                       !(external_synonym %in% gene_vector)) %>%
      distinct(external_gene_name, .keep_all = T) %>% # it'll remove a handful of weird genes with multiple synonyms simultaneously present in scrna, may be ignored
      distinct(external_synonym, .keep_all = T)
    
    if(nrow(geo_syn_in_gene_vector) > 0){
      common_genes <- sapply(gene_vector, function(x){
        if(x %in% geo_syn_in_gene_vector$external_gene_name){
          gname <- geo_syn_in_gene_vector$external_synonym[geo_syn_in_gene_vector$external_gene_name == x]
        } else{
          gname <- x
        }
        return(gname)
      })
      
      stopifnot(length(common_genes) == length(gene_vector))
      return(common_genes)
    } else{
      return(gene_vector)
    }
  } else{
    return(gene_vector)
  }
}

###################################################
###################################################
# adjust_synonym - useful to rescue couple hundred synonym genes, but often ensembl does not work 

# Description:
#   Prepares a list of signatures from the MSigDB database, optionally adjusting synonyms to match a reference.
# Parameters:
#   adjust_synonym: (logical) Whether to adjust synonyms in the gene list, default is TRUE.
#   geomx_obj: (S4 object) GeoMx object with its rownames used for synonym adjustment if needed.
#   hal: (logical) Whether to include hallmark gene sets, default is TRUE.
#   db_subcat_list: (vector) List of database subcategories to filter
# Return value:
#   (list) list of vectors with pathway signatures with adjusted genes.
prepare_msigdb_sign_list <- function(adjust_synonym = T, geomx_obj = NULL, hal = T, 
                                      db_subcat_list = c('CP:BIOCARTA', 'CP:KEGG', 'CP:REACTOME', 'CP:PID', 'CP:WIKIPATHWAYS', 'GO:BP')){
  
  # prepare msigdb signatures list
  msigdb_df <- msigdbr(species = "Homo sapiens")
  if(hal){
    msigdb_df <- filter(msigdb_df, gs_cat == 'H' | gs_subcat %in% db_subcat_list)
  } else{
    msigdb_df <- filter(msigdb_df, gs_subcat %in% db_subcat_list)
  }
  
  if(adjust_synonym){
    msigdb_df$gene_symbol_adj <- adjust_synonym_genes(rownames(geomx_obj), msigdb_df$gene_symbol)
    gene_colname <- 'gene_symbol_adj'
  } else{
    gene_colname <- 'gene_symbol'
  }

  hal_cp_list <- lapply(unique(msigdb_df$gs_name), function(x){
    gs <- filter(msigdb_df, gs_name == x)
    gs_genes <- unique(gs[[gene_colname]])
  })
  
  names(hal_cp_list) <- unique(msigdb_df$gs_name)
  
  return(hal_cp_list)
}

###################################################
###################################################
# Description:
#   Prepares a signature list from a custom data frame, with optional synonym adjustment.
# Parameters:
#   custom_sign_df: (data.frame) DataFrame with custom gene sets. 
#      With gene set name as column name and genes ar rows in a given column
#   adjust_synonym: (logical) Whether to adjust synonym names, default is TRUE.
#   geomx_obj: (S4 object) GeoMx object for benchmarking synonyms if needed.
# Return value:
#   (list) list of vectors with custom signatures, adjusted for synonyms if specified.
prepare_custom_sign_list <- function(custom_sign_df, adjust_synonym = T, geomx_obj = NULL){
  sign_list <- as.list(custom_sign_df)
  sign_list <- lapply(sign_list, function(l){l[l !=""]})
  
  if(adjust_synonym){
    sign_list <- lapply(sign_list, function(x){
      adjust_synonym_genes(rownames(geomx_obj), x)})
  }
  
  return(sign_list)
}
  

###############################################
#################################################
# Description:
#   Processes an expression matrix, log-transforming if necessary and optionally removing low complexity genes. 
#   For usage in downstream analysis (eg GSEA)
# Parameters:
#   geomx_obj_path: (string) Path to the GeoMx object file.
#   norm_type: (string) Normalization type or assay name in the geomx_obj.
#   scrna_ref_cleaned_path: (string) Path to reference object for filtering low complexity genes.
#   norm_is_log: (logical) Whether the input data is already log-transformed.
# Return value:
#   (matrix) Log-tranformed expression matrix with low complexity genes removed if specified.
prepare_expr_mtx <- function(geomx_obj_path, norm_type, norm_is_log, scrna_ref_cleaned_path = NULL){
  
  geomx_obj <- readRDS(geomx_obj_path)
  
  # make log expression mtx if needed
  if(!norm_is_log){
    # make log2 transformed normalised counts if norm_type not in log scale
    expr_norm_log <- log2(geomx_obj@assayData[[norm_type]] + 1)
  } else{
    expr_norm_log <- geomx_obj@assayData[[norm_type]]
  }
  
  # filter out from low complexity genes
  if(file.exists(scrna_ref_cleaned_path)){
    scrna_ref_obj <- readRDS(scrna_ref_cleaned_path)
    scrna_mtx_name <- ifelse('RNA_common_genes' %in% colnames(scrna_ref_obj@meta.data), 'RNA_common_genes', 'RNA')
    
    geomx_stat <- plot.bulk.outlier(
      bulk.input=t(geomx_obj@assayData$exprs),#make sure the colnames are gene symbol or ENSMEBL ID
      sc.input=t(scrna_ref_obj@assays[[scrna_mtx_name]]@data), #make sure the colnames are gene symbol or ENSMEBL ID
      cell.type.labels=scrna_ref_obj@meta.data$cell_type,
      species="hs", 
      return.raw=TRUE,
      pdf.prefix= NULL
    )
    
    geomx_stat_to_rm <- geomx_stat[ rowSums(geomx_stat[, -c(1,2)]) >= 1, ]
    expr_norm_log <- expr_norm_log[!(rownames(expr_norm_log) %in% rownames(geomx_stat_to_rm)),  ]
    
  }
  return(expr_norm_log)
}


###########################################
# Description:
#   Processes metadata for differential gene expression analysis, 
#   preparing factors and groups for comparison (eg binary vs multiclass).
# Parameters:
#   metadt: (data.frame) Metadata dataframe, typically phenotype data.
#   main_var_name: (string) Main variable name to test within.
#   main_var_is_bin: (logical) Whether the main variable is binary.
#   main_var_main_val: (string) Main value of interest for the binary variable.
#   dge_categories: (vector) Categories to include for DGE analysis.
#   cofounder_name: (string) Name of the column used for cofounder effect adjustment.
# Return value:
#   (data.frame) Modified metadata with adjusted groupings and factors for DGE.
# 
prepare_dge_metadata <- function(metadt, main_var_name, main_var_is_bin, main_var_main_val,
                                 dge_categories, cofounder_name){
  if(main_var_is_bin){
    # make binary vector - either main variable has the desired value or not
    metadt$main_var <- ifelse(grepl(main_var_main_val, metadt[, main_var_name]),
                              main_var_main_val, 'other_roi_type')
  } else{
    metadt$main_var <- metadt[, main_var_name]
  }
  
  print('groups which will be compared:')
  print(table(metadt[, c(main_var_name, 'main_var')]))
  
  # convert test variables to factors
  for(col in c(dge_categories, 'main_var')){
    metadt[[paste0(col, "_factor")]] <- factor(metadt[[col]])
  }
  
  metadt$cofounder_factor <- factor(metadt[[cofounder_name]])
  
  # make variable with all dge categories
  metadt$dge_group <- apply(metadt, 1, function(row){
    group <- sapply(dge_categories, function(var){
      paste(row[var])
    })
    group <- paste(group, collapse = '_')
    return(group)
  })
  
  return(metadt)
}

###########################################################33
############################################################
# Description:
#   Removes samples and/or groups from a GeoMx object that do not meet a specified 
#   minimum number of AOI/group for within slide DGE
# Parameters:
#   geomx_obj_dge_group: (S4 object) GeoMx object with samples for given DGE comparison
#   min_aoi_nr: (integer) Minimum number of areas of interest required for a group to be kept.
#   main_var_is_bin: (logical) Indicates if the main variable is binary.
#   comparison_type: (string) Type of comparison, either 'within' or 'between' (for 'between' no group removal)
# Return value:
#   (S4 object) The cleaned GeoMx object with small groups removed.
rm_too_small_groups <- function(geomx_obj_dge_group, min_aoi_nr, main_var_is_bin, comparison_type){
  samples_freq <- data.frame(table(pData(geomx_obj_dge_group)$main_var_factor,
                                   pData(geomx_obj_dge_group)$cofounder_factor))
  
  msg1 <- 'frequency of AOI in given group per sample:'
  print(msg1)
  print(samples_freq)
  
  if(comparison_type == 'within'){
    groups_keep <- samples_freq[samples_freq$Freq >= min_aoi_nr, ]
    
    # rmv samples with only 1 group with enough nr of ROI
    groups_keep_per_sample <- data.frame(table(groups_keep$Var2))
    sample_to_rm <- as.character(groups_keep_per_sample$Var1[groups_keep_per_sample$Freq < min_aoi_nr])
    
    groups_keep2 <- groups_keep[!(groups_keep$Var2 %in% sample_to_rm), ]
    
    msg2 <- 'only this groups will be keeped for DGE:'
    print(msg2)
    print(groups_keep2)
    
    keep_ind <- inner_join(pData(geomx_obj_dge_group), groups_keep2, 
                           by = c('main_var_factor' = 'Var1', 'cofounder_factor' ='Var2'))
    
    keep_ind <-  pData(geomx_obj_dge_group)$dcc_filename %in% keep_ind$dcc_filename
    
    # test
    #kk <- pData(geomx_obj_dge_group)[keep_ind, c('main_var_factor', 'cofounder_factor')]
    
    geomx_obj_dge_group_cleaned <- geomx_obj_dge_group[, keep_ind]
    
    return(list(geomx_obj = geomx_obj_dge_group_cleaned, logs = c(msg1, as.character(samples_freq), msg2, as.character(groups_keep2))))
  } else{
    return(list(geomx_obj = geomx_obj_dge_group, logs = c(msg1, as.character(samples_freq))))
  }
}
