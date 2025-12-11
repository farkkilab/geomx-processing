# geomx_pcg_list <- c('devtools','BiocManager', 'plyr', 'dplyr', 'data.table', 'tibble', 'tools', 'parallel',
#                     'ggforce', 'ggplot2', 'cowplot', 'ggrepel', 'reshape2', 
#                     'Biobase', 'NanoStringNCTools', 'GeomxTools', 'GeoDiff', 'DESeq2', 'SpatialDecon',
#                     'preprocessCore', 'PCAtools', 'umap', 'Rtsne', 'limma', 'pvca', 'harmony', 
#                     'Seurat', 'BayesPrism', 'biomaRt', 'msigdbr', 'GSVA', 'clusterProfiler',
#                     'BulkSignalR', 'igraph', 'scales', 'pheatmap', 'ComplexHeatmap','NMF', 'CellChat',
#                     'BiocNeighbors', 'patchwork', 'ggh4x', 'rlang', 'stringr', 'BiocGenerics', 'S4Vectors', 'stats4', 'SingleCellExperiment',
#                     'nichenetr', 'multinichenetr', 'tidyr', 'purrr', 'readr')
# 
# to_install <- geomx_pcg_list[which(!(geomx_pcg_list %in% ins))]
# for(pcg in to_install){BiocManager::install(pcg)}

# devtools::install_github("jinworks/CellChat")
# devtools::install_github("saeyslab/nichenetr")
# devtools::install_github("saeyslab/multinichenetr")
# devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")


# TODO move loading pck to certain scripts
# main packages for all scripts
library(plyr, quietly =T)
library(dplyr, quietly =T)# geomx_pcg_list <- c('devtools','BiocManager', 'plyr', 'dplyr', 'data.table', 'tibble', 'tools', 'parallel',
#                     'ggforce', 'ggplot2', 'cowplot', 'ggrepel', 'reshape2', 
library(data.table, quietly =T)
library(tibble, quietly =T)
library(tools, quietly = T)
library(parallel, quietly = T)

library(ggforce, quietly =T)
library(ggplot2, quietly =T)
library(cowplot, quietly =T)
library(ggrepel, quietly =T)
library(reshape2, quietly =T)

# problems with Matrix package - i has to be lower that 1.7 to work with lmer
# devtools::install_version("Matrix","1.6.4")
# dependencies:  c("grDevices", "graphics", "grid", "lattice", "methods",  "stats", "utils")  
library(Matrix)

# script specific packages
library(Biobase, quietly =T)
library(NanoStringNCTools, quietly =T)
library(GeomxTools, quietly =T)
library(GeoDiff, quietly =T)
library(DESeq2, quietly =T)
library(SpatialDecon, quietly =T)
#install preprocessCore manually from source
# BiocManager::install("preprocessCore", configure.args = c(preprocessCore = "--disable-threading"), 
# force= TRUE, update=TRUE, type = "source")
library(preprocessCore, quietly =T)
library(PCAtools) # NEW
library(umap, quietly =T)
library(Rtsne, quietly =T)

library(limma)
library(pvca)
library(harmony)

library(Seurat, quietly =T)
library(BayesPrism, quietly =T)
library(biomaRt, quietly =T)
library(msigdbr, quietly =T)

library(GSVA, quietly =T)
library(clusterProfiler, quietly =T)
#library(progeny, quietly =T)

# library(ggpubr)
# library(topGO)
# library(fgsea)
# library(fpc)
# library(dbscan)


# TODO --------------------------------------------------------------------


# TODO simplify logs by putting all console info from source() to logs
# TODO move loading libraries to each script separately

#all the batches should be merged and qc-ed + processed together and bigbatch + smallbatch variable as batch effects
batch <<- 'batch123' # for correct paths and batch eff vars 
merge_per_roi <- TRUE # whether or not signal from all AOIs within ROI should be added

# define variables and paths ----------------------------------------------

proj_dir <<- '~/Documents/phd/st'

# anno file have to contain sheet named 'Sheet1' and following column names:
# 'Sample_ID', 'Slide_Name',  'Aoi', 'Roi' and 'Panel'
# and dcc_name of proper NTC in 'NTC_ID' column if theres no 1NTC/batch

if(batch == 'batch1'){
  data_dir <<- '~/Documents/phd/st/data/geomx/geomx_batch1_0823' # batch1
  output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch1-1903') # batch1
  anno_path <<- file.path(data_dir, 'metadata', 'dcc_metadata_batch1_0823.xlsx') #batch1
} else if(batch == 'batch2'){
  data_dir <<- '~/Documents/phd/st/data/geomx/geomx_batch2_1124/' # batch2 
  output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch2-1903') # batch2
  anno_path <<- file.path(data_dir, 'metadata', 'dcc_metadata_batch2_1124.xlsx') #batch2
} else if(batch %in% c('batch3')){
  data_dir <<- '~/Documents/phd/st/data/geomx/geomx_batch3_0525/'
  output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch3-2606') # batch3
  anno_path <<- file.path(data_dir, 'metadata', 'dcc_metadata_all_batch3_0525_no_tls.xlsx') #batch1 and 2
} else if(batch %in% c('batch3-tls')){
  data_dir <<- '~/Documents/phd/st/data/geomx/geomx_batch3_0525/'
  output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch3-tls-2808') # batch3
  anno_path <<- file.path(data_dir, 'metadata', 'dcc_metadata_all_batch3_0525_tls.xlsx') #batch1 and 2
} else if(batch == 'batch12'){
  data_dir <<- '~/Documents/phd/st/data/geomx/batch12/' # batch1 and 2
  # output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch12-1004') # batch12
  output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch12-1205-no-counts-shift2') # batch12
  anno_path <<- file.path(data_dir, 'metadata', 'dcc_metadata_batch12.xlsx') #batch1 and 2
} else if(batch == 'batch23'){
  # metadata havent been changed and contains all roi with tls
  data_dir <<- '~/Documents/phd/st/data/geomx/batch23/' # batch2 and 3
  output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch23-2706') # batch23
  anno_path <<- file.path(data_dir, 'metadata', 'dcc_metadata_batch23.xlsx') #batch1 and 2
} else if(batch == 'batch123'){
  data_dir <<- '~/Documents/phd/st/data/geomx/batch123/' # batch1 2 and 3
  output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2711-roibased') # batch123
  anno_path <<- file.path(data_dir, 'metadata', 'dcc_metadata_batch123_no_tls_cleaned.xlsx') #batch1 and 2 and 3
}else{
  stop('wrong batch nr')
}

# input data
dcc_path <<- dir(file.path(data_dir, "dcc"), pattern = ".dcc$",
                 full.names = TRUE, recursive = TRUE)
pkc_path <<- file.path(data_dir, 'metadata', 'Hs_R_NGS_WTA_v1.0.pkc')

# path to reference scRNAseq dataset for deconvolution
# have to contain 'cell_type' column name in metadata
#scrna_ref_path <<- file.path(proj_dir, 'data/scrna/GSE165897_qc_downsampled_5k.RDS') # hautaniemi 
scrna_ref_path <<- file.path(proj_dir, 'data/scrna/GSE266577_qc_downsampled_keepfreq.RDS') # vaharautio

# path to csv file with custom gene signatures
custom_sign_path <<- file.path(proj_dir, 'geomx-processing', 'data', 'signatures',
                               'ct_markers.csv')


# set up metadata variables names -----------------------------------------

aoi_id <<- 'dcc_filename'
roi_id <<- 'Roi_geomx'

main_batch_var <- 'main_batch_nr'
batch_var <<- 'batch_nr'


aoi_segment_var <<- "Segment"
main_roi_label <<- "Annotation_cell" 
main_experimental_condition <<- 'NACT_status'
sample_name <<- 'Sample'

other_vars_bio <<- c("Segment_geomx", "Patient", "Site") #
other_vars_tech <<- c('Slide_Name')

# variables for batch effect removal
# if analysing each batch separately, only batch_var is considered
# if analysisng many big batches together, both main_batch_var and batch_var are considered

# !!! check throughfully the 1st PVCA plots from batch effect removal step
# if another variables are responsible for variance 
# primary_batch_var and secondary_batch_var values should be changed
# secondary batch variable has to be INDEPENDENT from the primary_batch_var

if(batch %in% c('batch1', 'batch2', 'batch3', 'batch3-tls')){
  primary_batch_var <<- batch_var
  secondary_batch_var <<- NULL
} else{
  primary_batch_var <<- main_batch_var
  secondary_batch_var <<- batch_var
}

# load util functions and create dirs -------------------------------------

source(file.path(proj_dir, 'geomx-processing', 'src', 'geomx_utils.R'))

dir.create(output_dir, recursive = T, showWarnings = F)

# define intermediate output paths ----------------------------------------

geomx_qc_path <<- file.path(output_dir, 'geomx_qc.RDS')
geomx_qc_roibased_path <<- file.path(output_dir, 'geomx_qc_roibased.RDS')
geomx_norm_path <<- file.path(output_dir, 'geomx_qc_norm.RDS')
geomx_norm_batch_eff_rm_path <<- file.path(output_dir, 'geomx_qc_norm_batch_eff_rm.RDS') 

# start the pipeline ------------------------------------------------------

print('#############')
print('GeoMx pipeline starting :O')
print('#############')

# conditionally run preprocessing -----------------------------------------

run_unless_exists('Preprocessing', geomx_qc_path, 
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_qc.R'))


# conditionally run merging signal per ROI --------------------------------

if(merge_per_roi){
  run_unless_exists('Merging signal per roi', geomx_qc_roibased_path, 
                    file.path(proj_dir, 'geomx-processing', 'src', 'geomx_merge_per_roi.R'))
  
  # change paths and main sample parameter
  geomx_qc_path <- geomx_qc_roibased_path
  aoi_id <<- 'sample_roi'
}

# conditionally run normalisation -----------------------------------------

run_unless_exists('Normalisation', geomx_norm_path, 
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_normalisation.R'))


# conditonally run batch effect removal -----------------------------------

run_unless_exists('Batch effect removal', geomx_norm_batch_eff_rm_path, 
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_batch_effect_rmv.R'))

# conditionally run deconvolution -----------------------------------------

# column name of cell type label in scRNAseq metadata
scrna_anno <<- 'mid_lvl_ct_updated' # either 'cell_type' / 'mid_lvl_ct' / 'mid_lvl_ct_updated' / 'low_lvl_ct'

deconv_logs_path <<- file.path(output_dir,'deconvolution', 
                               paste0('deconv_', scrna_anno, '_logs.txt'))

run_unless_exists('Deconvolution', deconv_logs_path, 
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_deconvolution.R'))

# conditionally run pathway analysis --------------------------------------

scrna_anno <<- 'mid_lvl_ct_updated' # either 'cell_type' / 'mid_lvl_ct' / 'mid_lvl_ct_updated' / 'low_lvl_ct'

pathway_inp_data_type <<- c('all', 'bp') # within c('all', 'bp')
# all - full geomx data (not-deconvoluted)
# bp - bayes prism deconvoluted data

# low_lvl_ct
# ct_of_interest <- c("Tcells_NK", "Bcells", "Myeloids","Mast_cells",
#                     "Fibroblasts_Endothelial", "tumor")

# mid_lvl_ct_updated
ct_of_interest <- c("Tcells_other","Tcells_CD8","Tcells_CD4", "Bcells", 'NKcells', 'Mast_cells',
                    "Macrophages_Monocytes", "DCs", "Fibroblasts_Mesothelial", "Endothelial_cells", "tumor")
ct_of_interest <- NULL
# if running for 'bp' (bayes prism deconvolution results) 
# specifies for which cell types GSEA should be computed (as in scrna_anno column in scRNAseq reference ds)
# if ct_of_interest <<- NULL - GSEA will be computed for all cell types

signature_type <<- 'custom' # c('msigdb', 'custom')
# msigdb - on all pathways from msigdb (Hallmark + CP)
# custom - on custom signatures list specified in custom_sign_path

signature_name <<- ifelse(signature_type == 'custom', gsub('.csv', '', basename(custom_sign_path)), '')

gsea_type <<- 'ssgsea' # 'gsva' or 'ssgsea'

gsea_logs_path <<- file.path(output_dir,'pathway_analysis', 'gsea', 
                             paste0(gsea_type, '_', signature_type, '_', signature_name, '_logs.txt'))


run_unless_exists('Pathway analysis', gsea_logs_path, 
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_pathway_analysis.R'))


# conditionally run differential gene expression --------------------------
# TODO compute voom() weights for dge - not so important
# (voom computes precision weights for the downstream dge)
# https://davislaboratory.github.io/GeoMXAnalysisWorkflow/articles/GeoMXAnalysisWorkflow.html#batch-correction
# TODO anova(full model, reduced model) - check if significantly improves the effect for interesting genes

# column name of cell type label in scRNAseq metadata
scrna_anno <<- 'mid_lvl_ct_updated' #either 'cell_type' / 'mid_lvl_ct' / 'mid_lvl_ct_updated' / 'low_lvl_ct'

dge_inp_data_type <<- c('all', 'bp') # within c('all', 'bp')

#dge_inp_data_type <<- c('bp') # within c('all', 'bp')
# all - full geomx data (not-deconvoluted)
# bp - bayes prism deconvoluted data

#ct_of_interest <<- c("tumor", "Tcells", "Fibroblasts", "Macrophages", "Endothelial cells", "DCs")
ct_of_interest <<- c("tumor", "Macrophages_Monocytes", "Tcells_CD8", "Tcells_CD4", "DCs", "Bcells", "Fibroblasts_Mesothelial")
# if running for 'bp' (bayes prism deconvolution results) 
# specifies for which cell types GSEA should be computed (as in scrna_anno column in scRNAseq reference ds)
# if ct_of_interest <<- NULL - GSEA will be computed for all cell types

# path to custom metadata with additional groups used for DGE
# must contain 'dcc_filename' column to merge with geomx_obj metadata
# if more column names are identical to the existing ones, columns from the custom dt will be used
# if not needed, set to NULL
#custom_metadt_path <<- file.path(proj_dir, 'geomx-processing', 'data', 'b12_dcc_clinical_data.csv')
custom_metadt_path <<- file.path(output_dir, 'deconvolution/bayes_prism/bp_hitum_in_stroma.csv')

# DGE parameters
comparison_type <<- 'within' 
# 'within' when you compare different ROI types within sample
# between - comparisons between slides
main_var_name <<- 'Segment_hitumor06' # main variable to make comparison between
main_var_is_bin <<- FALSE # should variable be compared with all others at once (TRUE) or with each other separately
# if FALSE all labels in main_var_name will be compared as they are

#main_var_main_val <<- 'posCD8_posIBA1'
#main_var_main_val <<- 'CD8_.*Iba1' # if main_var_is_bin - TRUE - name of the main value (or regex - careful!)
main_var_main_val <<- NULL
#dge_categories <<- c('Segment', 'NACT_status') # categories to divide to when making DGE separately
dge_categories <<- c('NACT_status')

# don't change it - identifier of dge run
dge_name <<- paste0('dge_', comparison_type, '_slide_', main_var_name, 
                   '_bin_', main_var_is_bin, '_', gsub('\\*', '', main_var_main_val), '_',
                   paste0(dge_categories, collapse = '_'))

dge_logs_path <<- file.path(output_dir, 'dge', dge_name, 'dge_logs.txt')


run_unless_exists('Differential Gene Expression', dge_logs_path, 
                  file.path(proj_dir, 'geomx-processing', 'src', 'geomx_dge.R'))


################
# GSEA on DGE
#TODO make it loop over many msigdb subdbs
#TODO check if clustering is done separately for each data group
# 'HALLMARK', 'CP:KEGG_MEDICUS', 'CP:REACTOME', 'GO:BP'

signature_type <<- 'msigdb' # c('custom','msigdb') 
msigdb_subcat <<- 'GO:BP' # if signature type is msigdb, which subdatabase to use. one of: c('HALLMARK', 'CP:BIOCARTA', 'CP:KEGG_MEDICUS', 'CP:REACTOME', 'CP:PID', 'CP:WIKIPATHWAYS', 'GO:BP')
  # if signature_type == 'custom' set to NULL
# msigdb - on msigdb db specified in 'msigdb_subcat
# custom - on custom signatures list specified in custom_sign_path

# run DGE enrichment
# TODO make a proper pipeline step
source(file.path(proj_dir, 'geomx-processing', 'src','sidescripts', 'geomx_dge_enrichment.R'), local = TRUE)


# conditionally run L-R interactions analysis -----------------------------

# column name of cell type label in scRNAseq metadata
scrna_anno <<- 'mid_lvl_ct_updated' #either 'cell_type' / 'mid_lvl_ct' / 'mid_lvl_ct_updated' / 'low_lvl_ct'


