library(data.table)
library(tibble)

geomx_obj <- readRDS('/home/iganiemi/Documents/phd/st/geomx-processing/results/batch1-1903/geomx_qc_norm_batch_eff_rm.RDS')
deconv_ct_frac <- fread('/home/iganiemi/Documents/phd/st/geomx-processing/results/batch1-1903/deconvolution/bayes_prism/bp_res_mid_lvl_ct_ct_fraction.csv')

expr <- data.frame(geomx_obj@assayData$harmony_batch_corr)
expr <- rownames_to_column(expr, 'gene_name')


meta <- sData(geomx_obj)

cols <- c("dcc_filename", "Slide_Name", "Sample", "Segment", 'Roi', 'Aoi', "area", "nuclei",
          "NACT_status", "Patient", "Site", "PFS_months",  "Annotation_cell",
          "UMAP1_harmony_batch_corr", "UMAP2_harmony_batch_corr", "tSNE1_harmony_batch_corr",
          "tSNE2_harmony_batch_corr") # "R0", "OS",

meta <- meta[, cols]

fwrite(expr, '/home/iganiemi/Documents/phd/st/geomx-processing/results/batch1-1903/geomx_harmony_corr_expr.csv')
fwrite(meta, '/home/iganiemi/Documents/phd/st/geomx-processing/results/batch1-1903/geomx_metadata.csv')
