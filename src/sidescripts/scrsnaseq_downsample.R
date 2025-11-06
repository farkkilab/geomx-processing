library(Seurat)
library(data.table)
library(plyr)
library(dplyr)
library(umap)
library(Rtsne)

#TODO copied from reference_scrnaseq_exploration - remove one of the scripts later

output_dir <- '/home/iganiemi/Documents/phd/st/data/scrna/'

mtx_path <- '/home/iganiemi/Documents/phd/st/data/scrna/GSE165897_UMIcounts_HGSOC.tsv'
meta_path <-  '/home/iganiemi/Documents/phd/st/data/scrna/GSE165897_cellInfo_HGSOC.tsv'


sc_ref_full_outpath <- '/home/iganiemi/Documents/phd/st/data/scrna/GSE165897_qc_full.RDS'
sc_ref_down_5k_outpath <- '/home/iganiemi/Documents/phd/st/data/scrna/GSE165897_qc_downsampled_5k.RDS'
sc_ref_down_10k_outpath <- '/home/iganiemi/Documents/phd/st/data/scrna/GSE165897_qc_downsampled_10k.RDS'


########################

mtx_path <- '/home/iganiemi/Documents/phd/st/data/scrna/GSE266577_counts_raw.mtx'
bar_path <- '/home/iganiemi/Documents/phd/st/data/scrna/GSE266577_barcodes.txt'
ft_path <- '/home/iganiemi/Documents/phd/st/data/scrna/GSE266577_seurat_features.txt'

meta_path <- '/home/iganiemi/Documents/phd/st/data/scrna/GSE266577_metadata.txt'

sc_ref_full_outpath <- '/home/iganiemi/Documents/phd/st/data/scrna/GSE266577_qc_full.RDS'
sc_ref_down_5k_outpath <- '/home/iganiemi/Documents/phd/st/data/scrna/GSE266577_qc_downsampled_5k.RDS'
sc_ref_down_keepfreq_outpath <- '/home/iganiemi/Documents/phd/st/data/scrna/GSE266577_qc_downsampled_keepfreq.RDS'

# from Seurat vignette
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html


# load scRNAseq as seurat object ------------------------------------------

# for GSE266577
sc_ref <- CreateSeuratObject(
  ReadMtx(mtx_path, bar_path, ft_path, feature.column = 1),
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = fread(meta_path),
  min.cells = 10,
  min.features = 200)

dim(sc_ref)

# for GSE165897
sc_ref <- CreateSeuratObject(
  counts = as.matrix(fread(mtx_path), rownames=1),
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = fread(meta_path),
  min.cells = 10,
  min.features = 200)

dim(sc_ref)

# basic QC ----------------------------------------------------------------

# Visualize QC metrics as a violin plot
VlnPlot(sc_ref, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

sc_ref <- subset(sc_ref, subset = percent.mt < 7.5) # 10 for GSE266577

dim(sc_ref)

# rename cell types -------------------------------------------------------

# for GSE266577
sc_ref@meta.data$mid_lvl_ct_updated <- mapvalues(sc_ref@meta.data$cell_type, 
                                                 from=c("Regulatory T cells", "Tem/Trm cytotoxic T cells", "Tcm/Naive helper T cells",
                                                        "Type 17 helper T cells", "ILC", "CD16- NK cells", "CD16+ NK cells", 
                                                        "NK cells", "Memory B cells", "Naive B cells", "Plasma cells", "Macrophages", 
                                                        "Classical monocytes", "Migratory DCs", "DC1", "DC2", "pDC", "Mast cells",
                                                        "Fibroblasts", "Endothelial cells", "Epithelial cells", "Mesothelial cells"), 
                                                 to=c("Tcells_other","Tcells_CD8","Tcells_CD4", "Tcells_other", "NKcells",
                                                      "NKcells", "NKcells", "NKcells", "Bcells", "Bcells", "Bcells", "Macrophages_Monocytes", "Macrophages_Monocytes",
                                                      "DCs", "DCs", "DCs", "DCs", "Mast_cells", "Fibroblasts_Mesothelial", "Endothelial_cells", "tumor", "Fibroblasts_Mesothelial"))

sc_ref@meta.data$low_lvl_ct <- mapvalues(sc_ref@meta.data$mid_lvl_ct_updated, 
                                         from=c("Tcells_other","Tcells_CD8","Tcells_CD4",
                                                "NKcells", "Bcells", "Macrophages_Monocytes",
                                                "DCs", "Mast_cells", "Fibroblasts_Mesothelial", "Endothelial_cells", "tumor"),
                                         to=c("Tcells_NK", "Tcells_NK", "Tcells_NK", "Tcells_NK", "Bcells",
                                              "Myeloids", "Myeloids","Mast_cells", "Fibroblasts_Mesothelial",
                                              "Endothelial_cells", "tumor"))

table(sc_ref@meta.data$cell_type, sc_ref@meta.data$mid_lvl_ct_updated)
table(sc_ref@meta.data$cell_type, sc_ref@meta.data$low_lvl_ct)

saveRDS(sc_ref, sc_ref_full_outpath)

##########################################
##########################################

# for GSE165897
sc_ref@meta.data$cell_type_main <- sc_ref@meta.data$cell_type
sc_ref@meta.data$cell_type <- sc_ref@meta.data$cell_subtype

sc_ref@meta.data$cell_type <- gsub('EOC', 'tumor', sc_ref@meta.data$cell_type)
sc_ref@meta.data$cell_type <- gsub('-', '_', sc_ref@meta.data$cell_type)
sc_ref@meta.data$mid_lvl_ct_updated <- ifelse(grepl('tumor', sc_ref@meta.data$cell_type), 'tumor', sc_ref@meta.data$cell_type)
sc_ref@meta.data$mid_lvl_ct_updated <- ifelse(grepl('CAF', sc_ref@meta.data$mid_lvl_ct_updated), 'Fibroblasts', sc_ref@meta.data$mid_lvl_ct_updated)
sc_ref@meta.data$mid_lvl_ct_updated <- ifelse(grepl('DC_', sc_ref@meta.data$mid_lvl_ct_updated), 'DC', sc_ref@meta.data$mid_lvl_ct_updated)

sc_ref@meta.data$low_lvl_ct <- mapvalues(sc_ref@meta.data$mid_lvl_ct_updated, 
                                                 from=c("tumor", "Fibroblasts", "Mesothelial", "Endothelial", "T_cells", "Plasma_cells",
                                                        "NK", "DC", "B_cells", "Macrophages", "pDC", "Mast_cells", "ILC"), 
                                                 to=c("tumor", "stroma", "stroma", "stroma",
                                                      "Tcells_NK", "B_cells", "Tcells_NK", "Myeloids", "B_cells", 
                                                      "Myeloids", "B_cells", "Mast_cells", "Tcells_NK"))

table(sc_ref@meta.data$cell_type, sc_ref@meta.data$mid_lvl_ct_updated)
table(sc_ref@meta.data$cell_type, sc_ref@meta.data$low_lvl_ct)
table(sc_ref@meta.data$low_lvl_ct, sc_ref@meta.data$mid_lvl_ct_updated)

saveRDS(sc_ref, sc_ref_full_outpath)

# downsample with min 5k per ct -------------------------------------------
max_ct_thr <- 5000 # downsample to this nr of cells per mid_lvl_ct
min_tumor_per_pt <- 10 # min nr of tumor cells in patient 50 for GSE577

mt <- sc_ref@meta.data 
table(mt$mid_lvl_ct_updated)

pt_name <- 'publication_patient_code_final' #GSE577
pt_name <- 'patient_id' #GSE 897

mt_downsampled <- lapply(unique(mt$mid_lvl_ct_updated), function(ct){
  ct_df <- mt[mt$mid_lvl_ct_updated == ct, ]
  freq_ratio <- max_ct_thr/nrow(ct_df)
  
  if(ct == 'tumor'){
    # remove patients with too little tumor cells
    ct_df <- ct_df %>%
      group_by(get(pt_name)) %>%
      mutate(n = n()) %>%
      filter(n > (min_tumor_per_pt/freq_ratio))
    
    freq_ratio <- max_ct_thr/nrow(ct_df) # updated freq ratio
  } 
  
  if(nrow(ct_df) > max_ct_thr){
    # downsample to max_ct_thr sample and fine grained cell type wise
    ct_df <- ct_df %>%
      group_by(get(pt_name), cell_type) %>%
      sample_n(round(n()*freq_ratio))
  }
  return(ct_df)
  })

mt_downsampled <- do.call(rbind, mt_downsampled)

table(mt_downsampled$mid_lvl_ct_updated)

#filter reference only to downsampled cells

sc_ref_down_5k <- subset(sc_ref, cells = mt_downsampled$cell) # cell_name

dim(sc_ref)
dim(sc_ref_down_5k)

saveRDS(sc_ref_down_5k, sc_ref_down_5k_outpath)

# downsample while keeping cell type proportions --------------------------
min_ct_nr <- 100 # downsample to min nr of cells in mid_lvl_ct
min_tumor_per_pt <- 10 # min nr of tumor cells in patient

mt <- sc_ref@meta.data %>%
  group_by(mid_lvl_ct_updated) %>%
  mutate(n = n())

freq_ratio_all <- min_ct_nr/min(mt$n)

downsample_df <- mt %>%
  group_by(publication_sample_code_final, cell_type) %>%
  sample_n(round(n()*freq_ratio_all))

# remove tumor cells not abundant per patient
tum_to_rm <- downsample_df %>%
  filter(mid_lvl_ct_updated == 'tumor') %>%
  group_by(publication_patient_code_final) %>%
  mutate(n = n()) %>%
  filter(n < min_tumor_per_pt)

downsample_df <- filter(downsample_df, !(cell_name %in% tum_to_rm$cell_name))

k <- downsample_df[downsample_df$mid_lvl_ct_updated == 'tumor', ]
table(k$publication_patient_code_final)

table(downsample_df$mid_lvl_ct_updated)
table(old_ref$mid_lvl_ct_updated)
table(sc_ref$mid_lvl_ct_updated)


sc_ref_down_keepfreq <- subset(sc_ref, cells = downsample_df$cell_name)

dim(sc_ref)
dim(sc_ref_down_5k)
dim(sc_ref_down_keepfreq)

saveRDS(sc_ref_down_keepfreq, sc_ref_down_keepfreq_outpath)


rm(sc_ref)
gc()

# make umap ---------------------------------------------------------------

sc_ref_path <- "/home/iganiemi/Documents/phd/st/data/scrna/GSE165897_qc_downsampled_5k.RDS"

sc_ref <- readRDS(sc_ref_path)
down_name <- '5k'

sc_ref <- NormalizeData(sc_ref)

sc_ref <- FindVariableFeatures(sc_ref, selection.method = "vst", nfeatures = 2000)

# scaling
sc_ref <- ScaleData(sc_ref, features = rownames(sc_ref))

# PCA
sc_ref <- RunPCA(sc_ref, features = VariableFeatures(object = sc_ref))

ElbowPlot(sc_ref)

# UMAP
sc_ref <- RunUMAP(sc_ref, dims = 1:10)

u1 <- DimPlot(sc_ref, reduction = "umap", group.by = "mid_lvl_ct_updated")
u2 <- DimPlot(sc_ref, reduction = "umap", group.by = "low_lvl_ct")
u3 <- DimPlot(sc_ref, reduction = "umap", group.by = "cell_type")


pdf(file= file.path(output_dir, paste0('umap_GSE165897_qc_downsampled_', down_name, '_mid_lvl_ct_updated.pdf')), width=8, height=5)
plot(u1)
dev.off()

pdf(file= file.path(output_dir, paste0('umap_GSE165897_qc_downsampled_', down_name, '_low_lvl_ct.pdf')), width=8, height=5)
plot(u2)
dev.off()

pdf(file= file.path(output_dir, paste0('umap_GSE165897_qc_downsampled_', down_name, '_cell_type.pdf')), width=8, height=5)
plot(u3)
dev.off()

