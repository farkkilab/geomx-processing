# Ligand Receptor Analysis by BulkSignaR : Bulk Data,Geomx Full Transcriptomic signal


# load libraries

library(dplyr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(tidyr)   
library(ggplot2)
library(scales)


# Params

nact_status = NULL
segment = NULL
annotation = NULL
heatmap_col_ann = "segment"

qval_threshold = 0.001 # filter significant LR pairs
n = 50 # number of top LR pairs needed to visualize in the signature scoes heatmap


# file names

input_dir <- 'C:/Users/Sahas/Downloads/Masters_Thesis/Ligand-receptor/BulkSignalR/TestRun_New_scripts/BulkSignalR_objects/'
output_dir <- 'C:/Users/Sahas/Downloads/Masters_Thesis/Ligand-receptor/BulkSignalR/TestRun_New_scripts/'


# pathways for plotting

# If not provided plot for all pathways : Need to adjust the size of the pdf
pathway_names <- read.csv("C:/Users/Sahas/Downloads/Masters_Thesis/Ligand-receptor/BulkSignalR/TestRun_New_scripts/pathway_names.csv")

# creating the folder for the plots

output_folder_name  = "/BulkSignalR_plots"

if (!dir.exists(paste0(output_dir,output_folder_name))) {
  new_dir = paste0(output_dir,output_folder_name)
  dir.create(new_dir,recursive = TRUE)
  message("Directory created")
  output_dir = new_dir
} else {
  output_dir = paste0(output_dir,output_folder_name)
}


# Reading BulkSignaR objects 

# TO DO : check the length of the segment and get the correct object 


parts <- c(nact_status, segment, annotation)
parts_non_null <- parts[!sapply(parts, is.null)]
if (length(parts_non_null) == 0) {
  mid_name = "combined"
  BulkSignalR_output = readRDS(paste0(input_dir,'/BulkSignalR_',mid_name,'_output.RDS'))
} else {
  mid_name = paste(parts_non_null, collapse = "_")
  BulkSignalR_output = readRDS(paste0(input_dir,'/BulkSignalR_',mid_name,'_output.RDS'))
}

bsrdm = BulkSignalR_output$bsrdm
bsrinf.redBP = BulkSignalR_output$bsrinf_redBP 


# Extracting LR pairs in the provided pathways

obj = bsrinf.redBP
pairs = LRinter(obj)
pairs$index <- seq_len(nrow(pairs))
selected_pairs = pairs %>% filter(qval < qval_threshold) %>% filter(pw.name  %in% pathway_names$Pathway.names)
 

top_n_pairs <- selected_pairs %>%
  arrange(desc(LR.corr)) %>%
  slice(1:n)

top_n_pairs_index <- top_n_pairs$index


ligands   <- ligands(obj)[top_n_pairs_index]
receptors   <- receptors(obj)[top_n_pairs_index]
pathways  <- pairs$pw.name
t.genes   <- tGenes(obj)[top_n_pairs_index]
t.corrs   <- tgCorr(obj)[top_n_pairs_index]



for (i in seq_len(nrow(top_n_pairs))){
  tg <- t.genes[[i]]
  t.genes[[i]] <- tg[top_n_pairs$rank[i]:length(tg)]
  
  tc <- t.corrs[[i]]
  t.corrs[[i]] <- tc[top_n_pairs$rank[i]:length(tc)]
}


bsrinf.redBP@LRinter = top_n_pairs
bsrinf.redBP@ligands = ligands
bsrinf.redBP@receptors = receptors
bsrinf.redBP@t.genes = t.genes
bsrinf.redBP@tg.corr = t.corrs


# step 4 : Building a BSRSignature object

bsrsig.redBP <- getLRGeneSignatures(bsrinf.redBP, qval.thres = qval_threshold)
scoresLR_pw <- scoreLRGeneSignatures(bsrdm, bsrsig.redBP,
                                     name.by.pathway=TRUE)

scoresLR <- scoreLRGeneSignatures(bsrdm, bsrsig.redBP,
                                  name.by.pathway=FALSE)


LRinter = obj@LRinter # check for obj
LRinter$lr_inter = paste0("{",LRinter$L,"} / {",LRinter$R,"}")


# collecting Meta data for the plot

meta_data = BulkSignalR_output$meta_data


# Filter the dataframe based on matches with the column names

row_names <- rownames(scoresLR)
matched_pw_names <- LRinter %>%
  filter(lr_inter %in% row_names) %>%
  arrange(match(lr_inter, row_names)) %>%  # Preserve the column order
  pull(pw.name)

# row annotations based on pathways

row_annot <- data.frame(Pathway = matched_pw_names)
rownames(row_annot) <- rownames(scoresLR)  # ensure they align with mat1


# deciding the column annotations

if (heatmap_col_ann == "segment") {
  col_annot <- data.frame(SampleGroup = factor(meta_data$Segment))
}
if (heatmap_col_ann == "NACT_status") {
  col_annot <- data.frame(SampleGroup = factor(meta_data$NACT_status))
}

# if both 
# col_annot <- data.frame(SampleGroup = factor(paste0(meta_data$Segment,"_",meta_data$NACT_status)))

rownames(col_annot) <- gsub('-', '.', meta_data$dcc_filename)

# colors for the row and col annotations

pathways <- as.character(unique(matched_pw_names)) # Convert to character vector (if not already)
Sample_Group <- as.character(unique(col_annot$SampleGroup))

colors_pw <- hue_pal()(length(pathways))
colors_pw <- setNames(colors_pw, pathways)
colors_SGroup <- hue_pal()(length(Sample_Group))
colors_SGroup <- setNames(colors_SGroup, Sample_Group)

ann_colors <- list(
  SampleGroup = colors_SGroup,
  Pathway = colors_pw)


# viridisLite::viridis() or colorspace::qualitative_hcl() better colors

# Finally the heatmap

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

pdf(paste0(output_dir,"/heatmap_LR_",mid_name,".pdf"), width = 12, height = 12)  # Width and height in inches

pheatmap(scoresLR, 
         annotation_row = row_annot,
         annotation_col = col_annot,
         annotation_colors = ann_colors,
         fontsize_col = 11,
         cellheight = 13,
         fontsize_row = 11,
         color = col_fun, 
         main = paste0("Ligand receptor signatures per samples"),
         name = "signature score",
         labels_col = rep("", ncol(scoresLR)))


dev.off()

