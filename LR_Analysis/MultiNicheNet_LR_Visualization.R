# MultiNicheNet Visualization


# for plots

plot_dir = file.path(output_dir, MultiNicheNet_folder_name,'plots_and_csv_files')
dir.create(plot_dir , recursive = T, showWarnings = F)


# by_default take from MultiNichenet outputs
# TODO if prioritized_tbl_oi provided use it else use the unfiltered 




for (group in comparison){
  for (receiver in cell_types) {
    for (sender in cell_types) {
      
      table_name <- paste(group, receiver, sender, sep = "_")
      
      multinichenet_output = readRDS(geomx_MultiNicheNet_path)
      top_n_LR_pairs = multinichenet_output$top_n_LR_pairs
      sample_data = top_n_LR_pairs[[table_name]]
      
      
      keep_sender_receiver_values = c(0.25, 0.9, 1.75, 4)
      names(keep_sender_receiver_values) = levels(sample_data$keep_sender_receiver)
      
      ######## calculate the median bulk expression for each group
      
      # calculate the median
      
      group_medians <- sample_data %>%
        group_by(group,lr_interaction) %>%
        summarize(median_scaled_LR = median(scaled_LR_pb_prod, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = group, values_from = median_scaled_LR)
      
      
      # Compute log2 fold change (stroma / tumor)
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
      
      ##########################################################################
      
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







