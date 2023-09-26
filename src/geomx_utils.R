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
    geom_parallel_sets_labels(color = "white", size = 5) +
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
# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
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

plot_q3_stats <- function(geomx_obj, ann_of_interest, output_name){
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

############################################################
plot_umap_tsne <- function(pheno_data, method_type = c('UMAP', 'tSNE'), 
                           norm_type = c('q3', 'quant'), color_var, shape_var = 'Segment',
                           output_name){
  ggplot(pData(geomx_obj),
         aes(x = get(paste0(method_type, '1_', norm_type, '_norm')), 
             y = get(paste0(method_type, '2_', norm_type, '_norm')), 
             color = get(color_var), shape = get(shape_var))) +
    geom_point(size = 3) +
    xlab(paste0('UMAP2_', norm_type, '_norm')) +
    ylab(paste0('UMAP2_', norm_type, '_norm')) +
    scale_color_discrete(name = color_var) + 
    scale_shape_discrete(name = shape_var) + 
    theme_bw()
  
  ggsave(output_name, width = 2000, height = 1500, unit='px')
}