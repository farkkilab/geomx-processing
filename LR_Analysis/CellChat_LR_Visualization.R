# CellChat Visualization

# for plots

plot_dir = file.path(output_dir, CellChat_folder_name,'plots_and_csv_files')
dir.create(plot_dir , recursive = T, showWarnings = F)



pval_threshold = 0.01
prob_threshold = 0.05


cellchat_output = readRDS(file = geomx_CellChat_path)
cellchat_results = cellchat_output$cellchat_results

lr_df_list = list()

for (group in comparison){
  
  df = cellchat_results[[group]]$df.net
  df$group = group
  lr_df_list[[group]] = df
  
}

df_combined = do.call(rbind, lr_df_list)
df_combined$sender_receiver = paste(df_combined$source, df_combined$target, sep = " -> ")

# TODO need to group based on pathway

for (receiver in cell_types) {
  for (sender in cell_types) {
    
    df_plot = df_combined %>% filter(source == sender, 
                                           target == receiver, 
                                           annotation != "ECM-Receptor",
                                           pval < pval_threshold,
                                           prob > prob_threshold)
    
    
    p1 =  df_plot %>%
      ggplot(aes(group, interaction_name_2, fill = prob)) +
      #geom_point() +
      geom_tile(color = "whitesmoke") +
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
        panel.spacing.x = unit(0.25, "lines"),
        panel.spacing.y = unit(0.15, "lines"),
        strip.text.x.top = element_text(size = 8, color = "black", face = "bold", angle = 0),
        strip.text.y.left = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=0.8, linetype="solid")
      ) + labs(fill = "Probability") 
    
    max_prob = df_plot$prob %>% max()
    
    custom_scale_fill = scale_fill_gradientn(
      colours = RColorBrewer::brewer.pal(9,"YlGnBu"),
      values = c(0, 0.2,0.5,0.8,1),  
      limits = c(0, max_prob))
    
    p1 = p1+ custom_scale_fill + scale_size_binned_area(max_size = 4) 
    
    
    file_name <- paste(sender,receiver, sep = "_")
    pdf(file.path(plot_dir, paste0(file_name,"_cellchat_plot.pdf")), width = 6, height = 6)
    print(p1)
    dev.off()
    
  }
}



