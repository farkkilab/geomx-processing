library(data.table)
library(plyr)
library(dplyr)
library(Seurat)
library(umap)
library(Rtsne)

output_dir <- "/home/iganiemi/Documents/phd/st/data/scrna/"
dir.create(output_dir)

sc_ref_orig_path <- "/home/iganiemi/Documents/phd/st/data/scrna/vaharautio_scrnaseq_dataset_downsampled_for_iga.RDS"
sc_ref_path <- "/home/iganiemi/Documents/phd/st/data/scrna/vaharautio_scrnaseq_dataset_downsampled_for_iga_processed.RDS"

sc_ref <- readRDS(sc_ref_path)
sc_ref_orig <- readRDS(sc_ref_orig_path)

###############################################
# re-cluster cell types

unique(sc_ref@meta.data$cell_type)
unique(sc_ref@meta.data$mid_lvl_ct)

sc_ref@meta.data$mid_lvl_ct_updated <- mapvalues(sc_ref@meta.data$cell_type, 
                               from=c("Regulatory T cells", "Tem/Trm cytotoxic T cells", "Tcm/Naive helper T cells",
                                      "Type 17 helper T cells", "ILC", "Late erythroid", "CD16- NK cells", "CD16+ NK cells", 
                                      "NK cells", "Memory B cells", "Naive B cells", "Plasma cells", "Macrophages", 
                                      "Classical monocytes", "Migratory DCs", "DC1", "DC2", "pDC", "Mast cells",
                                      "Fibroblasts", "Endothelial cells", "Epithelial cells"), 
                               to=c("Tcells_reg","Tcells_CD8","Tcells_CD4", "Tcells_other", "Tcells_other", "Tcells_other",
                                    "NKcells", "NKcells", "NKcells", "Bcells", "Bcells", "Bcells", "Macrophages_Monocytes", "Macrophages_Monocytes",
                                    "DCs", "DCs", "DCs", "DCs", "Mast_cells", "Fibroblasts", "Endothelial_cells", "tumor"))

sc_ref@meta.data$low_lvl_ct <- mapvalues(sc_ref@meta.data$mid_lvl_ct_updated, 
                                        from=c("Tcells_reg","Tcells_CD8","Tcells_CD4", "Tcells_other",
                                             "NKcells", "Bcells", "Macrophages_Monocytes",
                                             "DCs", "Mast_cells", "Fibroblasts", "Endothelial_cells", "tumor"),
                                        to=c("Tcells_NK", "Tcells_NK", "Tcells_NK", "Tcells_NK", "Tcells_NK", "Bcells",
                                             "Myeloids", "Myeloids","Mast_cells", "Fibroblasts_Endothelial",
                                             "Fibroblasts_Endothelial", "tumor"))


table(sc_ref@meta.data$cell_type, sc_ref@meta.data$mid_lvl_ct_updated)
table(sc_ref@meta.data$cell_type, sc_ref@meta.data$low_lvl_ct)

saveRDS(sc_ref, sc_ref_path)


####################################
# make umap
# from Seurat vignette
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

sc_ref <- FindVariableFeatures(sc_ref, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sc_ref), 10)

# scaling
sc_ref <- ScaleData(sc_ref, features = rownames(sc_ref))

# PCA
sc_ref <- RunPCA(sc_ref, features = VariableFeatures(object = sc_ref))

# UMAP
sc_ref <- RunUMAP(sc_ref, dims = 1:10)

u1 <- DimPlot(sc_ref, reduction = "umap", group.by = "mid_lvl_ct")
u2 <- DimPlot(sc_ref, reduction = "umap", group.by = "mid_lvl_ct_updated")
u3 <- DimPlot(sc_ref, reduction = "umap", group.by = "low_lvl_ct")
u4 <- DimPlot(sc_ref, reduction = "umap", group.by = "cell_type")

pdf(file= file.path(output_dir, 'umap_scrnaseq_ref_vaharautio_downsampled_mid_lvl_ct_orig.pdf'), width=8, height=5)
plot(u1)
dev.off()

pdf(file= file.path(output_dir, 'umap_scrnaseq_ref_vaharautio_downsampled_mid_lvl_ct_updated.pdf'), width=8, height=5)
plot(u2)
dev.off()

pdf(file= file.path(output_dir, 'umap_scrnaseq_ref_vaharautio_downsampled_low_lvl_ct.pdf'), width=8, height=5)
plot(u3)
dev.off()

pdf(file= file.path(output_dir, 'umap_scrnaseq_ref_vaharautio_downsampled_cell_type.pdf'), width=8, height=5)
plot(u4)
dev.off()

# Feature plot for Fibro and tumor markers expression

fibro_markers <- c('ACTA2', 'BGN', 'CAV1', 'COL1A1', 'COL1A2', 'COL6A1', 'COL6A2',
                   'DCN', 'DDR2', 'FAP', 'FBLN1', 'LUM', 'PDGFRA', 'PDGFRB', 'PDPN', 'POSTN')

tumor_markers <- c('BRCA2', 'MUC16', 'CD24', 'KRT7', 'CDH1', 'EPCAM', 'FAS', 'WFDC2', 'KLF6',
                   'KLF7', 'KLF8', 'KRT18', 'MUC16', 'PAX8', 'PIK3CA', 'SOX18', 'WFDC2', 'WT1')

fibro_plot <- FeaturePlot(sc_ref, features = fibro_markers)
tumor_plot <- FeaturePlot(sc_ref, features = tumor_markers)

pdf(file= file.path(output_dir, 'umap_scrnaseq_ref_vaharautio_downsampled_fibro_markers.pdf'), width=15, height=15)
plot(fibro_plot)
dev.off()

pdf(file= file.path(output_dir, 'umap_scrnaseq_ref_vaharautio_downsampled_tumor_markers.pdf'), width=20, height=12)
plot(tumor_plot)
dev.off()
