library(data.table)
library(dplyr)
library(CellChat)

# following the vignette
# https://rdrr.io/github/sqjin/CellChat/f/tutorial/CellChat-vignette.Rmd

prob_thr <- 0.1

###########################

# tum vs stroma
# cc_res_dir <- "/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2808/lr_interactions/cell_chat/Segment_NACT_status_deseq2log_harmony/"
# lr_res_all_path <- file.path(cc_res_dir, "CellChat_df_Segment_NACT_status_tumor_Bcells_Tcells_CD4_Tcells_other_Tcells_CD8_Fibroblasts_Mesothelial_Macrophages_Monocytes_DCs_lr.csv")
# cc_obj_path <- file.path(cc_res_dir, "CellChat_output_Segment_NACT_status_tumor_Bcells_Tcells_CD4_Tcells_other_Tcells_CD8_Fibroblasts_Mesothelial_Macrophages_Monocytes_DCs.RDS")

# ct1 <- 'Macrophages_Monocytes'
# ct2 <- 'tumor'
# dt_group <- 'tumor_post'

##################################################33
# groups
# "stroma_post_mixed_w_others"  "stroma_post_mixed_w_CD4"     "stroma_post_Macro_domin"    
# "stroma_post_CD8_Macro_domin" "stroma_post_Bcell_domin"  

# gmm labels
cc_res_dir <- "/home/iganiemi/Documents/phd/st/geomx-processing/results/batch123-2808/lr_interactions/cell_chat/stroma_post_roi_cluster_gmm_deseq2log_harmony/"
lr_res_all_path <- file.path(cc_res_dir, "CellChat_df_Segment_NACT_status_roi_cluster_label_gmm_tumor_Bcells_Tcells_CD4_Tcells_other_Tcells_CD8_Fibroblasts_Mesothelial_Macrophages_Monocytes_DCs_lr.csv")
cc_obj_path <- file.path(cc_res_dir, "CellChat_output_Segment_NACT_status_roi_cluster_label_gmm_tumor_Bcells_Tcells_CD4_Tcells_other_Tcells_CD8_Fibroblasts_Mesothelial_Macrophages_Monocytes_DCs.RDS")

ct1 <- 'Macrophages_Monocytes' #
ct2 <- 'Macrophages_Monocytes'
dt_group <- 'stroma_post_mixed_w_others'


###################
cc_plot_dir <- file.path(cc_res_dir, paste0('cc_plots_', dt_group, '_', 'prob', as.character(prob_thr)))
dir.create(cc_plot_dir)

#######################

paths_oi <- c("MHC-II", "ApoE", "FN1", "APP", "SPP1", "CCL", "IL10",
              "CD86", "CXCL", "MHC-I", "TGFb")

# load lr results and cellchat obj ----------------------------------------

lr_res_all <- as.data.frame(fread(lr_res_all_path))
lr_res_all <- lr_res_all[lr_res_all$prob >= prob_thr, ]

cc_all <- readRDS(cc_obj_path)
cellchat <- cc_all[[dt_group]]

# calculate the aggregated cell-cell communication network 
# by counting the number of links or summarizing the communication probability
cellchat <- aggregateNet(cellchat)

# load LR results for given ct pair and select best  ----------------------

lr_res_sel1 <- lr_res_all %>%
  filter(source == ct1 & target == ct2 & group == dt_group) %>%
  arrange(ligand, -prob)
  
lr_res_sel2 <- lr_res_all %>%
  filter(source == ct2 & target == ct1 & group == dt_group) %>%
  arrange(ligand, -prob)  

lr_res_sel1_best <- lr_res_sel1 %>%
  group_by(ligand) %>%
  filter(prob == max(prob)) %>%
  arrange(-prob)

lr_res_sel2_best <- lr_res_sel2 %>%
  group_by(ligand) %>%
  filter(prob == max(prob)) %>%
  arrange(-prob)

sort(table(lr_res_sel1_best$pathway_name))
unique(lr_res_sel1_best$pathway_name[1:21])


# selected visualisations -------------------------------------------------
ct_names <- colnames(cellchat@net$weight)
groupSize <- as.numeric(table(cellchat@idents)) # in geomx it just means nr of ROIs with given ct nr > thr 


for(path_name in paths_oi){
  
  # aggregated to cells
  pdf(file.path(cc_plot_dir, paste0("path_aggregated_chord_", path_name, '_', dt_group, '_', as.character(prob_thr), '.pdf')), width = 10, height =10)
  netVisual_aggregate(cellchat, signaling = path_name, vertex.receiver = seq(1:length(ct_names)), # groups to show
                      vertex.size = groupSize, thresh = prob_thr, layout = "chord")
  dev.off()
  
  # hmap with path strenght across all cells
  pdf(file.path(cc_plot_dir, paste0("path_hmap_", path_name, '_', dt_group, '.pdf')), width = 6, height =6)
  print(netVisual_heatmap(cellchat, signaling = path_name, color.heatmap = "Reds", , measure = "weight"))
  dev.off()
  
  # Compute the contribution of each ligand-receptor pair to the overall signaling 
  pdf(file.path(cc_plot_dir, paste0("path_contrib_", path_name, '_', dt_group, '_', as.character(prob_thr), '.pdf')), width = 10, height =10)
  print(netAnalysis_contribution(cellchat, signaling = path_name, thresh = prob_thr))
  dev.off()
  
  # violinplot with expression of path genes across cells
  pdf(file.path(cc_plot_dir, paste0("path_violin_expr_", path_name, "_", dt_group, ".pdf")))
  print(plotGeneExpression(cellchat, signaling = path_name, enriched.only = TRUE))
  dev.off()
  
  # bubbleplot with all  significant interactions (L-R pairs) for a pathway
  pairLR.use <- extractEnrichedLR(cellchat, signaling = path_name)
  
  pdf(file.path(cc_plot_dir, paste0("path_bubble_all_", path_name, '_', dt_group, '_', as.character(prob_thr), '.pdf')), width = 10, height =10)
  print(netVisual_bubble(cellchat, sources.use = seq(1,7), targets.use = seq(1,7), 
                   pairLR.use = pairLR.use, remove.isolate = TRUE, thresh = prob_thr))
  dev.off()
  
  # chordplot with all the significant interactions (L-R pairs) associated with certain signaling pathways
  pdf(file.path(cc_plot_dir, paste0("path_chord_all_", path_name, '_', dt_group, '_', as.character(prob_thr), '.pdf')), width = 10, height =10)
  print(netVisual_chord_gene(cellchat, sources.use = seq(1,7), targets.use = seq(1,7), signaling = path_name, 
                       legend.pos.x = 4, thresh = prob_thr))
  dev.off()
  
  pdf(file.path(cc_plot_dir, paste0("path_chord_sender_", ct1, "_",  path_name, '_', dt_group, '_', as.character(prob_thr), '.pdf')), width = 10, height =10)
  print(netVisual_chord_gene(cellchat, sources.use = match(ct1, ct_names), targets.use = seq(1,7), signaling = path_name, 
                       legend.pos.x = 4, thresh = prob_thr))
  dev.off()
  
  pdf(file.path(cc_plot_dir, paste0("path_chord_receiver_", ct1, "_",  path_name, '_', dt_group, '_', as.character(prob_thr), '.pdf')), width = 10, height =10)
  print(netVisual_chord_gene(cellchat, sources.use = seq(1,7), targets.use = match(ct1, ct_names), signaling = path_name, 
                       legend.pos.x = 4, thresh = prob_thr))
  dev.off()
  

  
}


############################
# for selected pathways:
pairLR.use.sel <- extractEnrichedLR(cellchat, signaling = paths_oi, thresh = prob_thr)
pdf(file.path(cc_plot_dir, paste0("selpaths_bubble_all_", dt_group, '_', as.character(prob_thr), '.pdf')), width = 10, height =10)
netVisual_bubble(cellchat, sources.use = seq(1,7), targets.use = seq(1,7), 
                 pairLR.use = pairLR.use.sel, remove.isolate = TRUE, thresh = prob_thr, font.size = 8)
dev.off()

pdf(file.path(cc_plot_dir, paste0("selpaths_chord_sender_", ct1, "_",  dt_group, '_', as.character(prob_thr), '.pdf')), width = 10, height =10)
netVisual_chord_gene(cellchat, sources.use = match(ct1, ct_names), targets.use = seq(1,7), signaling = paths_oi, 
                     legend.pos.x = 4, thresh = prob_thr, lab.cex = 0.6)
dev.off()

pdf(file.path(cc_plot_dir, paste0("selpaths_chord_receiver_", ct1, "_",  dt_group, '_', as.character(prob_thr), '.pdf')), width = 10, height =10)
netVisual_chord_gene(cellchat, sources.use = seq(1,7), targets.use = match(ct1, ct_names), signaling = paths_oi, 
                     legend.pos.x = 4, thresh = prob_thr, lab.cex = 0.6)
dev.off()


############################
# for all pathways (only for > prob)
if(prob_thr >= 0.1){
  pathways.all <- cellchat@netP$pathways
  pairLR.use.all <- extractEnrichedLR(cellchat, signaling = pathways.all, thresh = prob_thr, enriched.only = TRUE)
  
  pdf(file.path(cc_plot_dir, paste0("allpaths_bubble_all_", dt_group, '_', as.character(prob_thr), '.pdf')), width = 10, height =10)
  netVisual_bubble(cellchat, sources.use = seq(1,7), targets.use = seq(1,7), 
                   signaling = pathways.all, remove.isolate = TRUE, thresh = prob_thr, font.size = 4)
  dev.off()
  
  pdf(file.path(cc_plot_dir, paste0("allpaths_chord_sender_", ct1, "_",  dt_group, '_', as.character(prob_thr), '.pdf')), width = 20, height =20)
  netVisual_chord_gene(cellchat, sources.use = match(ct1, ct_names), targets.use = seq(1,7), signaling = pathways.all, 
                       legend.pos.x = 4, thresh = prob_thr, lab.cex = 0.4, slot.name = "netP")
  dev.off()
  
  pdf(file.path(cc_plot_dir, paste0("allpaths_chord_receiver_", ct1, "_",  dt_group, '_', as.character(prob_thr), '.pdf')), width = 30, height =30)
  netVisual_chord_gene(cellchat, sources.use = seq(1,7), targets.use = match(ct1, ct_names), signaling = pathways.all, 
                       legend.pos.x = 4, thresh = prob_thr, lab.cex = 0.4, slot.name = "netP", small.gap = 0.5, big.gap = 5)
  dev.off()

}

# 
# # visualisations with in-build functions ----------------------------------
# 
# # number of interactions or the total interaction strength (weights) 
# # between every cell groups with circle plot
# # not very useful for us since we artificially combine signal from given cell type/ROI into 1
# 
# groupSize <- as.numeric(table(cellchat@idents)) # in geomx it just means nr of ROIs with given ct nr > thr 
# # png(file.path(cc_plot_dir, "circleplot_interactions_all_nr.png"))
# # netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
# # dev.off()
# # 
# # png(file.path(cc_plot_dir, "circleplot_interactions_all_weights.png"))
# # netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# # dev.off()
# # 
# # # each cell separately
# mat <- cellchat@net$weight
# # for (i in 1:nrow(mat)) {
# #   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
# #   mat2[i, ] <- mat[i, ]
# #   png(file.path(cc_plot_dir, paste0("circleplot_interactions_", rownames(mat)[i], "_weights.png")))
# #   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
# #   dev.off()
# # }
# 
# ##################################################
# ##################################################
# path_name <- "MHC-II"
# gene_name <- "CD4"
# 
# dt_group <- names(cc_all)[1]
# 
# cellchat <- cc_all[[dt_group]]
# cellchat <- aggregateNet(cellchat)
# 
# png(file.path(cc_plot_dir, paste0("violin_expr_path_", path_name, "_", dt_group, ".png")))
# plotGeneExpression(cellchat, signaling = path_name, enriched.only = FALSE)
# dev.off()
# 
# png(file.path(cc_plot_dir, paste0("violin_expr_", gene_name, "_", dt_group, ".png")))
# plotGeneExpression(cellchat, features = gene_name, enriched.only = FALSE)
# dev.off()
# 
# #####################################################3
# ###################################################3##
# # visualise pathways
# allpaths <- cellchat@netP$pathways
# paths_oi <- c("MHC-II", "ApoE", "FN1", "APP", "SPP1", "CCL", "IL10",
#               "CD86", "CXCL", "MHC-I", "TGFb")
# 
# pathways.show <- paths_oi[1]
# # Hierarchy plot
# # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
# # vertex.receiver = seq(4,4) # a numeric vector. 
# # netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# # Circle plot
# # par(mfrow=c(1,1))
# # netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# # Chord diagram
# # par(mfrow=c(1,1))
# # netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# # Heatmap
# # par(mfrow=c(1,1))
# pdf(file.path(cc_plot_dir, paste0("hmap_path_", pathways.show, '_', dt_group, '.pdf')), width = 10, height =10)
# netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
# dev.off()
# 
# # vertex.receiver = seq(4,4)
# # for (i in 1:length(paths_oi)) {
# #   # Visualize communication network associated with both signaling pathway and individual L-R pairs
# #   #netVisual(cellchat, signaling = paths_oi[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
# #   # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
# #   gg <- netAnalysis_contribution(cellchat, signaling = paths_oi[i])
# #   ggsave(filename= file.path(cc_plot_dir, paste0('hmap_', paths_oi[i], "_L-R_contribution_Macrophages.pdf")),
# #                              plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
# # }
# 
# # Compute the contribution of each ligand-receptor pair to the overall signaling 
# # pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
# 
# netAnalysis_contribution(cellchat, signaling = pathways.show)
# 
# pairLR.oi <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
# LR.show <- pairLR.oi[1,] # show one ligand-receptor pair
# # Hierarchy plot
# vertex.receiver = seq(1,7) # a numeric vector
# netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# # Circle plot
# netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# # Chord diagram
# netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
# 
# 
# #####################################################3
# ###################################################3##
# # visualise cell-cell communication
# 
# # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# # netVisual_bubble(cellchat, sources.use = seq(1,7), targets.use = seq(1,7), remove.isolate = FALSE)
# # show all the significant interactions (L-R pairs) associated with certain signaling pathways
# netVisual_bubble(cellchat, sources.use = seq(1,7), targets.use = seq(1,7), signaling = paths_oi, remove.isolate = FALSE)
# # show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
# pairLR.use <- extractEnrichedLR(cellchat, signaling = paths_oi)
# netVisual_bubble(cellchat, sources.use = seq(1,7), targets.use = seq(1,7), pairLR.use = pairLR.use, remove.isolate = TRUE)
# 
# 
# # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# # show all the interactions sending from Inflam.FIB
# pdf(file.path(cc_plot_dir, paste0("chord_communication_Macrophages_Monocytes.pdf")), width = 20, height =16)
# netVisual_chord_gene(cellchat, sources.use = 4, targets.use = 4, lab.cex = 0.5,legend.pos.y = 30)
# dev.off()
# 
# # show all the significant interactions (L-R pairs) associated with certain signaling pathways
# netVisual_chord_gene(cellchat, sources.use = seq(1,7), targets.use = seq(1,7), signaling = paths_oi[2], 
#                      legend.pos.x = 4, thresh = 0.1)
# netVisual_chord_gene(cellchat, sources.use = seq(1,7), targets.use = 4, signaling = paths_oi[2], 
#                      legend.pos.x = 4, thresh = 0.1)
# netVisual_chord_gene(cellchat, sources.use = 4, targets.use = seq(1,7), signaling = paths_oi[2], 
#                      legend.pos.x = 4, thresh = 0.1)
# # show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# # netVisual_chord_gene(cellchat, sources.use = 4, targets.use = 4, slot.name = "netP", legend.pos.x = 10)
# 
