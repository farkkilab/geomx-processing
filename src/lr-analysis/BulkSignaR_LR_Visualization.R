# Ligand Receptor Analysis by BulkSignaR : Bulk Data,Geomx Full Transcriptomic signal


# Reading BulkSignaR objects
# TODO loop over all objects in a list
BulkSignaR_Output_list = readRDS(geomx_BulkSignalR_path)
BulkSignaR_Output = BulkSignaR_Output_list$BulkSignaR_Output
#BulkSignaR_Output_combined = BulkSignaR_Output[["combined"]]

pathway_names <- pathway$pathway


for(outname in names(BulkSignaR_Output)){
  print(outname)
  
  bsrdm = BulkSignaR_Output[[outname]]$bsrdm
  bsrinf.redBP = BulkSignaR_Output[[outname]]$bsrinf_redBP 
  meta_data = BulkSignaR_Output[[outname]]$meta_data # collecting Meta data for the plot
  
  # heatmap with top-n pairs, without pathways selection
  heatmap_no_sel = plot_heatmap(bsrinf.redBP, bsrdm, meta_data, pathway_names = NULL, qval_threshold, n, heatmap_col_ann)
  
  if(!is.null(heatmap_no_sel)){
    pdf(file.path(plot_dir,paste0("heatmap_LR_",qval_threshold, '_', outname, ".pdf")), width = 12, height = 12)  # Width and height in inches
    print(heatmap_no_sel)
    dev.off()
  }

  # heatmap with top-n pairs with pathways selection
  heatmap_sel = plot_heatmap(bsrinf.redBP, bsrdm, meta_data, pathway_names = pathway$pathway, qval_threshold, n, heatmap_col_ann)
  
  if(!is.null(heatmap_sel)){
    pdf(file.path(plot_dir,paste0("heatmap_LR_",qval_threshold, '_', outname, "_selected_pathways.pdf")), width = 12, height = 12)  # Width and height in inches
    print(heatmap_sel)
    dev.off()
  }
  
}


############ bubble plot #####################

# Bubble plots if you are comparing between groups
# First need to generate LR pairs separately for each group
# provide pathways you want to visualize as a .csv file in the master script 


if (!is.null(manually_filtered_BulkSignalr_df) && is.data.frame(manually_filtered_BulkSignalr_df)) {
  
  df_for_plotting = manually_filtered_BulkSignalr_df
  
} else {
  
  
  df_for_plotting = BulkSignaR_Output_list$unfiltered_LR_df_for_plotting
  df_for_plotting = df_for_plotting %>% 
    filter(qval < qval_threshold, LR.corr > LR_corr_threshold) %>%
    arrange(desc(LR.corr))
  #df_for_plotting = filter(df_for_plotting, pw.name %in% pathway_names)
  df_for_plotting = head(df_for_plotting, n)
  
}





p1 =  df_for_plotting %>%
  ggplot(aes(group, lr_interaction, color = LR.corr, size = neg_log10_p_adj)) +
  geom_point() +
  facet_grid(pw.name~group, scales = "free", space = "free", switch = "y")+
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

max_corr = abs(df_for_plotting$LR.corr) %>% max()

custom_scale_fill = scale_color_gradientn(
  colours = RColorBrewer::brewer.pal(n = 7, name = "PuRd"),
  values = c(0, 0.5,0.6,0.7,0.8,0.9,1),  
  limits = c(0,max_corr))

p1 = p1+ custom_scale_fill + scale_size_binned_area(max_size = 4) 



pdf(file.path(plot_dir, "bubble_plot.pdf"), width = 9, height = 7)
print(p1)
dev.off()


