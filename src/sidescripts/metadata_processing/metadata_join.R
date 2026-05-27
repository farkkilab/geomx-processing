library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(tibble)


# TODO add:
# CNV
# LOH
# TCR/BCR diversity ?
# cycif whole slide metrics

# TODO do we need to recalculate fraction of immune per AOI?

# notes
# ct fractions are coming from sd (also used for clustering)
# roi clusters clusters_gmm_clustnr_5

#S131 info from iCAN: 79 coding mut, 5 tier3, 74 tier4
#TP53 TIER3 mut

# define paths ------------------------------------------------------------

proj_dir <<- '~/Documents/phd/st'
output_dir <<- file.path(proj_dir, 'geomx-processing', 'results', 'batch123-2808')
output_metadata_path <- file.path(output_dir, 'metadata_full_SENSITIVE.csv')

# define input files paths ------------------------------------------------

# main metadata
metadata_orig_path <- file.path(proj_dir, 'data/geomx/batch123/metadata/dcc_metadata_batch123_no_tls_cleaned.csv')

# clinical data
clinical_dt_path <<- file.path(proj_dir, 'data/geomx/clinical_data/9_eyemt_patient_clinical_data_SENSITIVE_upd_0426.csv')

# somatic mutations
snv_path <- file.path(proj_dir, 'data/geomx/clinical_data/9_eyemt_imtb_somatic_snv_summary.csv')

# ct fractions per aoi and roi from geomx_roi_hubs_integration
ct_frac_deconv_aoi_path <- file.path(output_dir, 'cycif_integration', 'b123_ct_frac_deconv_bcells.csv')
ct_frac_deconv_roi_path <- file.path(output_dir, 'cycif_integration', 'b123_ct_frac_deconv_roi_bcells.csv')
# TODO now SD fractions are merged - change if needed

# ct fractions of immune per roi (used for clustering) from geomx_relabel_roi_2nd_approach
ct_frac_deconv_roi_immunefrac_path <- file.path(output_dir, 'deconvolution', 'relabel-roi-deconv-dimred', 'sd_mye_lymph_b_ct_fractions_of_immune.csv')

# ROI clusters based on deconvolution ct fractions from geomx_relabel_roi_2nd_approach
# TODO change clust_type if needed - also , 'clusters_hclust_cut2' is ok
ct_frac_clust_path <- file.path(output_dir, 'deconvolution', 'relabel-roi-deconv-dimred', 'sd_mye_lymph_b_all_clustering_results.csv')
clust_type <- 'clusters_gmm_clustnr_5'
# descriptive labels for clusters - IN THIS CASE BOTH METHODS HAS THE SAME CLUSTERS DESCRIPTION
clust_labels <- list(CD8_Macro_domin = 1, mixed_w_CD4 = 2, mixed_w_others = 3, Macro_domin = 4, Bcell_domin = 5)

# clean data --------------------------------------------------------------

# main dcc metadata
metadt <- fread(metadata_orig_path) %>%
  mutate(sample_roi = paste0(Sample, '_', Roi_geomx))

# clinical data
clindt <- fread(clinical_dt_path, drop = seq(1, 6)) %>% 
  select(-tumor_purity, -`selected block`, -dataset, -`patient_in _ican`, -Site, -NACT_status)  %>% # removing sample related vars
  distinct()

# summary of somatic snv nr from iCAN reports
snv <- fread(snv_path) #TP53_mut - somatic SNV

# ct fractions from deconvolution per ROI and AOI  
ct_frac_aoi <- fread(ct_frac_deconv_aoi_path, select = c('dcc_filename', 'cell_type', 'ct_frac_sd')) %>%
  spread(key = 'cell_type', value = 'ct_frac_sd')
colnames(ct_frac_aoi) <- paste0('ct_frac_sd_aoi_', colnames(ct_frac_aoi))
colnames(ct_frac_aoi)[1] <- 'dcc_filename'

ct_frac_roi <- fread(ct_frac_deconv_roi_path, select = c('sample_roi', 'cell_type', 'ct_frac_sd')) %>%
  spread(key = 'cell_type', value = 'ct_frac_sd')
colnames(ct_frac_roi) <- paste0('ct_frac_sd_roi_', colnames(ct_frac_roi))
colnames(ct_frac_roi)[1] <- 'sample_roi'

ct_immunefrac_roi <- fread(ct_frac_deconv_roi_immunefrac_path)
colnames(ct_immunefrac_roi) <- paste0('ct_immunefrac_sd_roi_', colnames(ct_immunefrac_roi))
colnames(ct_immunefrac_roi)[1] <- 'sample_roi'

# clusters labels
roi_clust <- fread(ct_frac_clust_path, select = c('sample_roi', clust_type))
roi_clust$roi_cluster_label <- mapvalues(roi_clust[[clust_type]], 
                                      from=c(unname(unlist(clust_labels))),
                                      to=c(names(clust_labels)))


# merge data --------------------------------------------------------------

metadt <- left_join(metadt, clindt) %>%
  left_join(snv) %>%
  left_join(ct_frac_aoi) %>%
  left_join(ct_frac_roi) %>%
  left_join(ct_immunefrac_roi) %>%
  left_join(roi_clust[, c('sample_roi', 'roi_cluster_label')]) %>%
  as.data.frame()


# save --------------------------------------------------------------------

fwrite(metadt, output_metadata_path)
