# CellChat Visualization

# TODO super correct would be to first identify LR/path > thr in any group
# and then add value from another group to compare (now other group has NA)

# TODO for pathways max prob is very high and other vals are not visible - now scale for each

# load CellChat results ---------------------------------------------------

if (!is.null(manually_filtered_cellchat_df) && is.data.frame(manually_filtered_cellchat_df)) {
  df_for_plotting = manually_filtered_cellchat_df
} else {
  df_for_plotting = fread(cc_lr_df_path)
}

df_path_for_plotting <- fread(cc_path_df_path)

# max prob for consistent color scale
max_prob = df_for_plotting$prob %>% max()


# plot lr interactions between each ct ------------------------------------

for (receiver in unique(df_for_plotting$target)) {
  for (sender in unique(df_for_plotting$source)) {
    
    print(paste0(sender,"-",receiver))
    
    df_plot = df_for_plotting %>% 
      filter(source == sender, 
             target == receiver, 
             pval < pval_threshold,
             prob > prob_threshold)
    
    print(nrow(df_plot))
    
    p1 =  df_plot %>%
      ggplot(aes(group, interaction_name_2, fill = prob)) +
      #geom_point() +
      geom_tile(color = "whitesmoke") +
      facet_nested(pathway_name ~ group, scales = "free", space = "free", switch = "y")+
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
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=0.5, linetype="solid")) + 
      labs(fill = "Probability") +
      ggtitle(paste(sender,"->",receiver))
    
    custom_scale_fill = scale_fill_gradientn(
      colours = RColorBrewer::brewer.pal(9,"YlGnBu"),
      values = c(0, 0.2,0.5,0.8,1),  
      limits = c(0, max_prob))
    
    p1 = p1+ custom_scale_fill + scale_size_binned_area(max_size = 4) 
    
    file_name <- paste(sender,receiver, sep = "_")
    pdf(file.path(plot_dir, paste0(file_name,"_cellchat_plot_prob_", gsub('\\.', '', as.character(prob_threshold)), ".pdf")),
        width = 7.5, height = 9)
    print(p1)
    dev.off()
  }
}


# plot results on pathway lvl ---------------------------------------------

for (receiver in unique(df_path_for_plotting$target)) {
  for (sender in unique(df_path_for_plotting$source)) {
    print(paste0(sender,"-",receiver))
    
    df_path_plot = df_path_for_plotting %>% 
      filter(source == sender, 
             target == receiver, 
             pval < pval_threshold,
             prob > prob_threshold)
    
    
    p2 =  df_path_plot %>%
      ggplot(aes(group, pathway_name, fill = prob)) +
      #geom_point() +
      geom_tile(color = "whitesmoke") +
      facet_nested(~ group, scales = "free", space = "free", switch = "y")+
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
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=0.5, linetype="solid")) + 
      labs(fill = "Probability") +
      ggtitle(paste(sender,"->",receiver))
    
    max_prob_path = df_path_plot$prob %>% max()
    
    custom_scale_fill = scale_fill_gradientn(
      colours = RColorBrewer::brewer.pal(9,"YlGnBu"),
      values = c(0, 0.2,0.5,0.8,1),  
      limits = c(0, max_prob_path))
    
    p2 = p2+ custom_scale_fill + scale_size_binned_area(max_size = 4) 
    
    file_name <- paste(sender,receiver, sep = "_")
    pdf(file.path(plot_dir, paste0(file_name,"_cellchat_plot_pathway_prob_", gsub('\\.', '', as.character(prob_threshold)), ".pdf")),
        width = 7.5, height = 9)
    print(p2)
    dev.off()
  }
}

