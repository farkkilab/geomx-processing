# MultiNicheNet Visualization
# by_default take from MultiNichenet outputs
# TODO if prioritized_tbl_oi provided use it else use the unfiltered 

## Plot functions

plot_bulk_expression = function(df_plot1){
  
  p1 =  df_plot1 %>%
    ggplot(aes(group, lr_interaction, color = diff_median, size = neg_log10_p_adj)) +
    geom_point() +
    facet_grid(sender_receiver~group, scales = "free", space = "free", switch = "y")+
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      #axis.text.y = element_blank(),
      axis.text.y = element_text(face = "bold.italic", size = 7),
      axis.text.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.40, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x.top = element_text(size = 8, color = "black", face = "bold", angle = 0),
      strip.text.y.left = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "Median difference\nin scaled L-R\npseudobulk\nexpression\nproduct", size = "-log10(pval_adj)") 
  
  max_diff_median = abs(df_plot1$diff_median) %>% max()
  
  custom_scale_fill = scale_color_gradientn(
    colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),
    values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  
    limits = c(-1*max_diff_median, max_diff_median))
  
  p1 = p1+ custom_scale_fill + scale_size_binned_area(max_size = 4) 
  return(p1)
  
  
}

# If you want the bulkexpression data use sample_data directly and change the color accordignly



plot_igand_activity = function(df_plot2){
  
  p2 = df_plot2 %>%
    ggplot(aes(direction_regulation , lr_interaction, fill = activity_scaled)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(sender_receiver~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      #axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 8,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 8, color = "black", face = "bold"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(fill = "Scaled Ligand\nActivity RNA")
  
  max_activity = abs(df_plot2$activity_scaled) %>% max(na.rm = TRUE)
  
  custom_scale_fill = scale_fill_gradientn(
    colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]),
    values = c(0, 0.51, 0.575, 0.625, 0.675, 0.725, 1),  
    limits = c(-1*max_activity, max_activity))
  
  p2 = p2 + custom_scale_fill
  
  
  return(p2) 
  
}



######################################

# TODO

if (!is.null(manually_filtered_LR_pairs_dfplot_median_bulk_expr) && !is.null(manually_filtered_LR_pairs_dfplot_ligand_activity) &&
    is.data.frame(manually_filtered_LR_pairs_dfplot_median_bulk_expr) && is.data.frame(manually_filtered_LR_pairs_dfplot_ligand_activity)) {
  
  
  df_plot1 = manually_filtered_LR_pairs_dfplot_median_bulk_expr
  df_plot2 = manually_filtered_LR_pairs_dfplot_ligand_activity
  
  df_plot1 = df_plot1 %>% filter(lr_interaction %in% df_plot2$lr_interaction)
  df_plot2 = df_plot2 %>% filter(lr_interaction %in% df_plot1$lr_interaction)
  
  p1 = plot_bulk_expression(df_plot1)
  p2 = plot_igand_activity(df_plot2)
  
  
  p = patchwork::wrap_plots(
    p1,p2,
    nrow = 1,guides = "collect",
    widths = c(6,6)
  )
  
  pdf(file = file.path(plot_dir,paste0("MultiNicheNet_plot_manually_filtered.pdf")), width = 17, height = 10)
  print(p)
  dev.off()
  

  
} else {
  
  multinichenet_output = readRDS(geomx_MultiNicheNet_path)
  top_n_LR_pairs = multinichenet_output$top_n_LR_pairs
  
  for (group in comparison){
    for (receiver in cell_types) {
      for (sender in cell_types) {
        
        table_name <- paste(group, sender, receiver, sep = "_")
        sample_data = top_n_LR_pairs[[table_name]]
        
        
        df_plot1 = sample_data$df_plot1
        df_plot2 = sample_data$df_plot2
        
        p1 = plot_bulk_expression(df_plot1)
        p2 = plot_igand_activity(df_plot2)
        
        ###########################################################################
        
        p = patchwork::wrap_plots(
          p1,p2,
          nrow = 1,guides = "collect",
          widths = c(6,6)
        )
        
        pdf(file = file.path(plot_dir,paste0("MultiNicheNet_plot_",table_name,".pdf")), width = 17, height = 10)
        print(p)
        dev.off()
        
        
        
      }}}
  
  
  
  
}


