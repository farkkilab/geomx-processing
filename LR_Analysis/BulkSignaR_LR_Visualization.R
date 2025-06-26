# Ligand Receptor Analysis by BulkSignaR : Bulk Data,Geomx Full Transcriptomic signal


# Params

qval_threshold = 0.001 # filter significant LR pairs
n = 50 # number of top LR pairs needed to visualize in the signature scoes heatmap


# for plots

plot_dir = file.path(output_dir, BulkSignalR_folder_name,'plots_and_csv_files')
dir.create(plot_dir , recursive = T, showWarnings = F)


# Reading BulkSignaR objects

BulkSignaR_Output = readRDS(geomx_BulkSignalR_path)
BulkSignaR_Output = BulkSignaR_Output$BulkSignaR_Output
BulkSignaR_Output_combined = BulkSignaR_Output[["combined"]]


bsrdm = BulkSignaR_Output_combined$bsrdm
bsrinf.redBP = BulkSignaR_Output_combined$bsrinf_redBP 
meta_data = BulkSignaR_Output_combined$meta_data # collecting Meta data for the plot
pathway_names = pathway$Pathway.names



# Generate a heatmap
plot_heatmap = plot_heatmap(bsrinf.redBP, bsrdm, meta_data, pathway_names, qval_threshold, n, grouping_var_col_ids)


# save in a pdf
pdf(file.path(plot_dir,paste0("heatmap_LR_",qval_threshold,".pdf")), width = 12, height = 12)  # Width and height in inches
print(plot_heatmap)
dev.off()



############ bubble plot #####################

# Bubble plots if you are comparing between groups
# First need to generate LR pairs separately for each group
# Define what groups you need to compare : eg stroma vs tumor  












