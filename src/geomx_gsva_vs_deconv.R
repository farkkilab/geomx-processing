library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(GSVA)
library(plyr)
library(dplyr)
library(data.table)
library(biomaRt)
library(DESeq2)
library(msigdbr)
library(tibble)
library(ggcorrplot)

# get variables -----------------------------------------------------------
data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'
output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2'

input_rds_path <- file.path(output_dir, 'geomx_qc_norm.RDS')

input_gsva_all_path <- file.path(output_dir, 'gsva', 'gsva_all_selected.csv')
input_gsva_deconv_macro_path <- file.path(output_dir, 'gsva', 'gsva_deconv_macro_mid_lvl_ct_selected.csv')
input_gsva_deconv_tcell_path <- file.path(output_dir, 'gsva', 'gsva_deconv_tcell_mid_lvl_ct_selected.csv')

input_bp_deconv_path <- file.path(output_dir, 'deconvolution', 'bp', 'bp_res_mid_lvl_ct_45.RDS')
input_sd_deconv_path <- file.path(output_dir, 'deconvolution', 'sd', 'sd_res_mid_lvl_ct_nofilt.rds')

# if not doing all hall_cp
selected_sig_path <- file.path('/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/immune_signatures_selected_forpaper_names.csv')

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')

dir.create(file.path(output_dir, 'gsva', 'lm_gsva_deconv'), showWarnings = T, recursive = T)


# read files --------------------------------------------------------------

selected_sig <- fread(selected_sig_path)

# ct_names <- c("Tcells", "Bcells", "Fibroblasts", "NKcells", "Macrophages", 
#               "DCs", "tumor")

ct_names <- c("Tcells", "Fibroblasts", "Macrophages", "tumor")

sd_deconv <- readRDS(input_sd_deconv_path)
sd_deconv <- data.frame(pData(sd_deconv)[, 'prop_of_all'])
sd_deconv <- sd_deconv[, c(ct_names, "Endothelial.cells")]

bp_deconv <- readRDS(input_bp_deconv_path)
bp_deconv <- BayesPrism::get.fraction(bp=bp_deconv,
                         which.theta="final",
                         state.or.type="type")
bp_deconv <- bp_deconv[, c(ct_names, "Endothelial cells")]

# iterate through 3 files
# iterate trough bp/sd deconv
gsva_type_list <- list('all' = input_gsva_all_path, 'deconv_macro' = input_gsva_deconv_macro_path,
                       'deconv_tcell' = input_gsva_deconv_tcell_path)

deconv_list <- list('bp' = bp_deconv, 'sd' = sd_deconv)

sapply(1:length(gsva_type_list), function(x){
  gsva <- fread(gsva_type_list[[x]])
  
  gsva <- gsva[gsva$pathway %in% selected_sig$pathway, ]
  
  gsva_type_name <- names(gsva_type_list)[x]
  print('XXXXXX')
  print(gsva_type_name)
  
  sapply(1:length(deconv_list), function(x){
    deconv <- deconv_list[[x]]
    deconv_name <- names(deconv_list)[x]
    print(deconv_name)
    
    dir.create(file.path(output_dir, 'gsva', 'lm_gsva_deconv', gsva_type_name, deconv_name), showWarnings = T, recursive = T)
    
    deconv <- as.data.frame(deconv)
    deconv <- tibble::rownames_to_column(deconv, 'dcc_filename')
    
    gsva_df <- left_join(gsva, deconv)
    # iterate trough all, pre, post, short, long, 
    
    gsva_list <- list('all' = gsva_df,
                      'pre' = gsva_df[gsva_df$`NACT status` == 'pre', ],
                      'post' = gsva_df[gsva_df$`NACT status` == 'post', ],
                      'long' = gsva_df[gsva_df$`NACT status` == 'post' & gsva_df$PFS == 'Long', ],
                      'short' = gsva_df[gsva_df$`NACT status` == 'post' & gsva_df$PFS == 'Short', ])
    
    
    sapply(1:length(gsva_list), function(x){
      gsva_inp <- gsva_list[[x]]
      gsva_name <- names(gsva_list)[x]
      print(gsva_name)
      
      # make scatterplot with lm 
      gsva_lm <- lapply(unique(as.vector(gsva_inp$pathway)), function(path_name){
        gsva_path <- gsva_inp[gsva_inp$pathway == path_name, ]
        
        lapply(colnames(deconv)[-1], function(ct){
          
          print(ct)
          # fit lm with ct fraction as explanatory var
          lm_res <- lm(gsva_score~get(ct),data=gsva_path)
          lm_coef <- summary(lm_res)$coefficients[2]
          lm_rsq <- summary(lm_res)$adj.r.squared
          
          cor_pe <- cor.test(gsva_path$gsva_score, gsva_path[[ct]], method = 'pearson')
          #cor_sp <- round(cor(gsva_path$gsva_score, gsva_path[[ct]], method = 'spearman'), 2)
          
          gsva_lm_res <- list('pathway' = path_name, 'deconv_ct' = ct,
                              lm_coef = lm_coef, lm_rsq = lm_rsq,
                              cor_pe_coeff = cor_pe$estimate, cor_pe_pval = cor_pe$p.value)
          
          png(file = file.path(output_dir, 'gsva', 'lm_gsva_deconv', gsva_type_name, deconv_name,
                               paste0('scatter_', path_name, '_', ct, '.png')))
          plot(gsva_path[[ct]], gsva_path$gsva_score, xlab = path_name, ylab = ct, sub = paste0('R2 = ', lm_rsq))
          abline(lm(gsva_score~get(ct),data=gsva_path),col='red')
          dev.off()
          
          return(gsva_lm_res)
        })
      })
      
      gsva_lm <- unlist(gsva_lm, recursive = F)
      gsva_lm_df <- rbindlist(gsva_lm, fill=TRUE)
      
      
      fwrite(gsva_lm_df, file.path(output_dir, 'gsva', 'lm_gsva_deconv', gsva_type_name, deconv_name,
                                   paste0('lm_gsva_', gsva_name, '.csv')))
      
      # make a corrplot
      
      gsva_corr_coeff <- as.data.frame(dcast(gsva_lm_df, pathway ~ deconv_ct,
                            value.var = 'cor_pe_coeff'))
      
      gsva_corr_pval <- as.data.frame(dcast(gsva_lm_df, pathway ~ deconv_ct,
                                             value.var = 'cor_pe_pval'))
      
      sapply(unique(selected_sig$path_type), function(path_type){
        gsva_corr_coeff_path_type <- dplyr::filter(gsva_corr_coeff, pathway %in% 
                                           selected_sig$pathway[selected_sig$path_type == path_type]) %>%
          tibble::column_to_rownames('pathway')
        
        gsva_corr_pval_path_type <- dplyr::filter(gsva_corr_pval, pathway %in% 
                                                     selected_sig$pathway[selected_sig$path_type == path_type]) %>%
          tibble::column_to_rownames('pathway')
        
        corrplot <- ggcorrplot(gsva_corr_coeff_path_type, hc.order = FALSE, outline.color = "white", # p.mat = gsva_corr_pval_path_type
                               title = path_type, tl.cex = 3, pch.cex = 2) +
          theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))
        
        pdf(file=file.path(output_dir, 'gsva', 'lm_gsva_deconv', gsva_type_name, deconv_name,
                           paste0('corrplot_gsva_', gsva_name,'_', path_type, '.pdf')),
            width=5, height=6) 
        
        plot(corrplot)
        dev.off()
        
      })
    })
  })
  
})




###################################33
#####################################

kk <- fread("/media/iganiemi/T7-iga/st/geomx-processing/results/nact2/gsva/lm_gsva_deconv/all/bp/lm_gsva_all.csv")






