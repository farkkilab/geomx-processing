# Ligand Receptor Analysis by BulkSignaR : Bulk Data,Geomx Full Transcriptomic signal

# TODO fix legend placement, row names width etc
# TODO bubble plot for redPBP, red LBP and red RBP
# set up parameters -------------------------------------------------------
# params
qval_threshold = 0.001 # filter significant LR pairs
top_n = 50 # number of top LR pairs needed to visualize in the signature scoes heatmap
# TODO better setting
#heatmap_col_ann = "Segment" # based on what you want to annotate the heatmap

reduction_name <- 'redBP' # from c('redPBP', 'redBP', 'redLBP', 'redRBP')

pathway_names <- pathway$pathway

aoi_id <- 'dcc_filename' # from main
# read BSR outputs --------------------------------------------------------

# Reading BulkSignaR objects
BulkSignalR_Output = readRDS(geomx_BulkSignalR_path)

LR_output <- readRDS(lr_output_path)

# make heatmap ------------------------------------------------------------

outname <- 'combined'

for(outname in names(BulkSignalR_Output)){
  print(outname)
  
  heatmap_col_ann = ifelse(outname == 'combined', "Segment", "NACT_status")
  
  bsrdm = BulkSignalR_Output[[outname]]$bsrdm
  bsrinf_red = BulkSignalR_Output[[outname]][[paste0('bsrinf_', reduction_name)]] 
  meta_data = BulkSignalR_Output[[outname]]$meta_data # collecting Meta data for the plot
  
  # heatmap with top-n pairs, without pathways selection
  plot_heatmap(bsrinf_red, bsrdm, reduction_name, meta_data, pathway_names = NULL, 
               qval_threshold, top_n, heatmap_col_ann, aoi_id, 
               file.path(plot_dir,paste0("heatmap_LR_",qval_threshold, '_', outname, "_", reduction_name, ".pdf")))

  # heatmap with top-n pairs with pathways selection
  plot_heatmap(bsrinf_red, bsrdm, reduction_name, meta_data, pathway_names = pathway_names, 
               qval_threshold, top_n, heatmap_col_ann, aoi_id,
               file.path(plot_dir,paste0("heatmap_LR_",qval_threshold, '_', outname, "_", reduction_name, "_selpath.pdf")))
}


############ bubble plot #####################

# Bubble plots if you are comparing between groups
# First need to generate LR pairs separately for each group
# provide pathways you want to visualize as a .csv file in the master script 


# if (!is.null(manually_filtered_BulkSignalr_df) && is.data.frame(manually_filtered_BulkSignalr_df)) {
#   
#   df_for_plotting = manually_filtered_BulkSignalr_df
#   
# } else {
 
# }

if(reduction_name == 'redBP'){
  df_for_plotting <- LR_output[[paste0('bsrinf_', reduction_name)]]
  
  df_for_plotting = df_for_plotting %>% 
    filter(group != 'combined') %>%
    filter(qval < qval_threshold) %>%
    mutate(lr_inter = paste0("{",L,"} / {",R,"}")) 
  
  # fix names for different reduction types
  df_for_plotting$lr_inter <- gsub('{{', '{', df_for_plotting$lr_inter, fixed = T)
  df_for_plotting$lr_inter <- gsub('}}', '}', df_for_plotting$lr_inter, fixed = T)
  
  top_lr <- df_for_plotting %>%
    arrange(desc(LR.corr)) %>%
    head(n)
  
  # to retain values in all groups
  df_for_plotting_top <- df_for_plotting[df_for_plotting$lr_inter %in% top_lr$lr_inter, ]
  
  # make glued pah name if for the same pair, diff best pathway in diff group
  df_for_plotting_top <- df_for_plotting_top %>%
    group_by(lr_inter) %>%
    mutate(pw_both = paste(pw.name, collapse = ' / ')) %>%
    ungroup()
  
  
  p1 =  df_for_plotting_top %>%
    ggplot(aes(group, lr_inter, color = LR.corr, size = neg_log10_p_adj)) +
    geom_point() +
    facet_grid(pw_both~group, scales = "free", space = "free", switch = "y")+
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_text(face = "bold.italic", size = 7),
      axis.text.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.40, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x.top = element_text(size = 8, color = "black", face = "bold", angle = 0),
      strip.text.y.left = element_text(size = 6, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "LR correlation", size = "-log10(pval_adj)")
  
  max_corr = abs(df_for_plotting$LR.corr) %>% max()
  
  custom_scale_fill = scale_color_gradientn(
    colours = RColorBrewer::brewer.pal(n = 7, name = "PuRd"),
    values = c(0, 0.5,0.6,0.7,0.8,0.9,1),
    limits = c(0,max_corr))
  
  p1 = p1+ custom_scale_fill + scale_size_binned_area(max_size = 4)
  
  pdf(file.path(plot_dir, paste0("bubble_plot_", reduction_name, ".pdf")), width =12, height = 12)
  print(p1)
  dev.off()
  
  
  
}


