# MHC-II in tumor vs T-cell abundance in corresponding stroma

# MHC-II in tumor vs IFNg in corresponding stroma

# MHC-II in tumor vs quirky epitopes in same AOI
library(data.table)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(reshape2)
library(tibble)
###################################
outp_dir <- '~/Documents/phd/st/geomx-processing/results/batch1-1903/downstream_for_fer'
dir.create(outp_dir)

#geomx_meta <- fread('~/Documents/phd/st/geomx-processing/results/batch1-1903/geomx_metadata.csv')
deconv_ct_frac <- data.frame(fread('~/Documents/phd/st/geomx-processing/results/batch1-1903/deconvolution/bayes_prism/bp_res_mid_lvl_ct_ct_fraction.csv'))
gsea_all <- data.frame(fread('~/Documents/phd/st/geomx-processing/results/batch1-1903/pathway_analysis/gsea/ssgsea_norm_harmony_batch_corr_all_custom_IFNg_pathways.csv.csv'))
gsea_deconv_tum <- data.frame(fread('~/Documents/phd/st/geomx-processing/results/batch1-1903/pathway_analysis/gsea/ssgsea_norm_harmony_batch_corr_deconv_tumor_custom_IFNg_pathways.csv.csv'))

deconv_tum_expr_mtx <- readRDS('~/Documents/phd/st/geomx-processing/results/batch1-1903/deconvolution/bayes_prism/bp_res_mid_lvl_ct_expr_mtx_cleaned_vst_harmony_batch_corr.RDS')$tumor

ct_of_interest <- c('Tcells', 'Bcells', 'DCs', 'Fibroblasts', 'NKcells', 'Macrophages')
meta_names <- c('Sample', 'Patient', 'NACT_status', 'PFS', 'PFS_months', 'Annotation_cell')

neoepitope_genes <- as.vector(fread('~/Documents/phd/st/geomx-processing/results/batch1-1903/downstream_for_fer/neoepitopes.csv'))
##########################
# prepare dfs
gsea_res <- gsea_all

gsea_res$Roi <- paste0(gsea_res$Sample, '_', gsea_res$Roi)
gsea_res <- left_join(gsea_res, deconv_ct_frac[, c('dcc_filename', ct_of_interest)]) 

####
# gsea_res_mhc <- gsea_res[gsea_res$pathway == 'FER_MHC_REVISITED', ] 
# 
# deconv_tum_expr_mtx <- data.frame(deconv_tum_expr_mtx[rownames(deconv_tum_expr_mtx) %in% neoepitope_genes[[1]], ])
# deconv_tum_expr_mtx <- rownames_to_column(deconv_tum_expr_mtx, 'gene_name')
# 
# fwrite(gsea_res_mhc, '~/Documents/phd/st/geomx-processing/results/batch1-1903/downstream_for_fer/metadata_with_tumor_deconv_mhc_gsea_and_ct_fractions.csv')
# fwrite(deconv_tum_expr_mtx, '~/Documents/phd/st/geomx-processing/results/batch1-1903/downstream_for_fer/deconv_tumor_expr_neoepitope_genes.csv')

####

gsea_res_str <- gsea_res[gsea_res$Segment == 'stroma', ]
gsea_res_tumor <- gsea_res[gsea_res$Segment == 'tumor', ]

# TODO instead of this make it wide before and join with paths of interest
tum_mhc <- gsea_res_tumor[gsea_res_tumor$pathway == 'FER_MHC_REVISITED', ]
tum_mhc <- rename(tum_mhc, ssgsea_mhc_tum = 'ssgsea_score')
tum_mhc <- rename_with(tum_mhc, ~paste0(., "_tum"), ct_of_interest)

stroma_ifn <- gsea_res_str[gsea_res_str$pathway == 'HALLMARK_INTERFERON_GAMMA_RESPONSE', ]
stroma_ifn <- rename(stroma_ifn, ssgsea_ifng_str = 'ssgsea_score')
stroma_ifn <- rename_with(stroma_ifn, ~paste0(., "_str"), ct_of_interest)

tum_str_dt <- left_join(tum_mhc[, c(meta_names, 'Roi', 'ssgsea_mhc_tum', paste0(ct_of_interest, '_tum'))], 
                        stroma_ifn[, c('Roi', 'ssgsea_ifng_str', paste0(ct_of_interest, '_str'))])


################################
# make scatters

inp_df <- tum_str_dt
xname <- 'ssgsea_mhc_tum'
ynames <- c(paste0(ct_of_interest, '_tum'), paste0(ct_of_interest, '_str'))

for(yname in ynames){
  for(color_colname in meta_names){
    print(color_colname)
    
    #rm NA if needed
    inp_df2 <- inp_df[!is.na(inp_df[[xname]]) & !is.na(inp_df[[yname]]), ]
    #manual_colours <- brewer.pal(length(unique(as.factor(inp_df[[color_colname]]))), 'Paired')
    M1 <- lm(get(yname) ~ exp(get(xname)) + Patient, data = inp_df2)
    
    
    ggplot(data = inp_df2, 
           aes(x = get(xname), y = get(yname))) +
      geom_point(aes(color = as.factor(get(color_colname))), alpha = 0.3) + 
      ggtitle(paste0(xname, ' vs ', yname),
              #subtitle = paste0('cor: ', round(cor(inp_df[[xname]], inp_df[[yname]], method = 'pearson'), 2))) +
              subtitle = paste('lm(y ~ exp(x) | Patient) : R2' , round(summary(M1)$r.squared, 2))) +
      xlab(xname) +
      ylab(yname) +
      guides(color=guide_legend(title=color_colname)) +
      #scale_color_manual(values=manual_colours) +
      #geom_smooth(method='glm', formula= y~x ) +
      stat_poly_line(method = 'lm', formula = y~exp(x)) +
      stat_poly_eq(use_label(c("R2")), formula = y~exp(x))
    
    ggsave(file.path(outp_dir, paste0('scatter_', xname, '_vs_', yname, '_by_', color_colname, '_tumor_deconvoluted_signal.svg')))
    
  }
}
