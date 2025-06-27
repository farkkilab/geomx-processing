# CellChat Visualization


if (!is.null(manually_filtered_cellchat_df) && is.data.frame(manually_filtered_cellchat_df)) {
  
  df_for_plotting = manually_filtered_cellchat_df
  
} else {
  
  cellchat_output = readRDS(file = geomx_CellChat_path)
  df_for_plotting = cellchat_output$unfiltered_LR_df_for_plotting
  
  
}



for (receiver in cell_types) {
  for (sender in cell_types) {
    
    print(paste0(sender,"-",receiver))
    
    df_plot = df_for_plotting %>% filter(source == sender, 
                                           target == receiver, 
                                           annotation != "ECM-Receptor",
                                           pval < pval_threshold,
                                           prob > prob_threshold)
    
    
    p1 =  df_plot %>%
      ggplot(aes(group, interaction_name_2, fill = prob)) +
      #geom_point() +
      geom_tile(color = "whitesmoke") +
      facet_nested(sender_receiver + pathway_name ~ group, scales = "free", space = "free", switch = "y")+
      #facet_grid(sender_receiver~group, scales = "free", space = "free", switch = "y")+
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
        panel.spacing.x = unit(0.15, "lines"),
        panel.spacing.y = unit(0.05, "lines"),
        strip.text.x.top = element_text(size = 8, color = "black", face = "bold", angle = 0),
        strip.text.y.left = element_text(size = 7, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=0.5, linetype="solid")
      ) + labs(fill = "Probability") 
    
    max_prob = df_plot$prob %>% max()
    
    custom_scale_fill = scale_fill_gradientn(
      colours = RColorBrewer::brewer.pal(9,"YlGnBu"),
      values = c(0, 0.2,0.5,0.8,1),  
      limits = c(0, max_prob))
    
    p1 = p1+ custom_scale_fill + scale_size_binned_area(max_size = 4) 
    
    
    file_name <- paste(sender,receiver, sep = "_")
    pdf(file.path(plot_dir, paste0(file_name,"_cellchat_plot.pdf")), width = 7.5, height = 9)
    print(p1)
    dev.off()
    
  }
}



