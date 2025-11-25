# experimental: linear adjustment of GSEA scores calculated on the full signal 
# to the number of certain cell types calculated by spatialdecon

###############################################################################
###############################################################################
# adjusting gsva scores from full signal for spatialdecon cell freq -------
dir.create(file.path(output_dir, 'gsva', 'sd_lm_additional_msigdb_filt'), showWarnings = T, recursive = T)

sd_deconv <- readRDS(input_sd_deconv_path)
sd_deconv <- data.frame(pData(sd_deconv)[, 'prop_of_all'])
sd_deconv <- dplyr::select(sd_deconv, -Mast.cells, -other)
colnames(sd_deconv) <- paste0('deconv_', colnames(sd_deconv))
deconv_names <- colnames(sd_deconv)
sd_deconv <- tibble::rownames_to_column(sd_deconv, 'dcc_filename')

#gsva_all_long <- fread(file.path(output_dir, 'gsva', 'gsva_all_selected.csv'))
gsva_all_long <- gsva_list_long[[1]]
gsva_all_long <- left_join(gsva_all_long, sd_deconv)

selected_sig <- fread(selected_sig_path)
gsva_df <- filter(gsva_all_long, pathway %in% selected_sig$pathway)

##############
# calculate lm for each pathway vs deconvolution cell type
gsva_lm <- lapply(unique(as.vector(gsva_df$pathway)), function(path_name){
  gsva_path <- gsva_df[gsva_df$pathway == path_name, ]
  
  lapply(deconv_names, function(ct){
    
    # fit lm with ct fraction as explanatory var
    lm_res <- lm(gsva_score~get(ct),data=gsva_path)
    lm_coef <- summary(lm_res)$coefficients[2]
    lm_rsq <- summary(lm_res)$adj.r.squared
    
    gsva_lm_res <- list('pathway' = path_name, 'deconv_ct' = ct,
                        lm_coef = lm_coef, lm_rsq = lm_rsq)
    
    if(!do_gsva_hal_cp_all){
      png(file = file.path(output_dir, 'gsva', 'sd_lm_additional_msigdb_filt', paste0('scatter_', path_name, '_', ct, '.png')))
      plot(gsva_path[[ct]], gsva_path$gsva_score, xlab = path_name, ylab = ct)
      abline(lm(gsva_score~get(ct),data=gsva_path),col='red')
      dev.off()
    }
    return(gsva_lm_res)
  })
})

gsva_lm <- unlist(gsva_lm, recursive = F)
gsva_lm_df <- rbindlist(gsva_lm, fill=TRUE)

fwrite(gsva_lm_df, file.path(output_dir, 'gsva', paste0('sd_lm_gsva_', out_name, '.csv')))


##########################################
# correct gsva scores based on lm parameters
# adjusted for lin reg
# y = a + xb
# y = a // -xb

rsq_thr <- 0.4

gsva_lm_adj_tcell <- lapply(unique(gsva_df$pathway), function(path_name){
  gsva_df_path <- gsva_df[gsva_df$pathway == path_name, ]
  
  lm_tcell <- gsva_lm_df[gsva_lm_df$pathway == path_name & grepl(cd8_ct, gsva_lm_df$deconv_ct), ]
  
  if(lm_tcell$lm_rsq > rsq_thr){
    gsva_df_path$gsva_score <- gsva_df_path$gsva_score - (gsva_df_path[[paste0('deconv_', cd8_ct)]] * lm_tcell$lm_coef)
  }
  
  return(gsva_df_path)
})

gsva_lm_adj_macro <- lapply(unique(gsva_df$pathway), function(path_name){
  gsva_df_path <- gsva_df[gsva_df$pathway == path_name, ]
  
  lm_macro <- gsva_lm_df[gsva_lm_df$pathway == path_name & grepl(macro_ct, gsva_lm_df$deconv_ct), ]
  
  if(lm_macro$lm_rsq > rsq_thr){
    gsva_df_path$gsva_score <- gsva_df_path$gsva_score - (gsva_df_path[[paste0('deconv_', macro_ct)]] * lm_macro$lm_coef)
  }
  return(gsva_df_path)
})

gsva_lm_adj_tcell <- do.call(rbind, gsva_lm_adj_tcell)
gsva_lm_adj_macro <- do.call(rbind, gsva_lm_adj_macro)

fwrite(gsva_lm_adj_tcell, file.path(output_dir, 'gsva', 
                                    paste0('gsva_all_', out_name, '_sd_lm_adjusted_tcell_additional_msigdb_filt_', as.character(rsq_thr), '.csv')))

fwrite(gsva_lm_adj_macro, file.path(output_dir, 'gsva', 
                                    paste0('gsva_all_', out_name, '_sd_lm_adjusted_macro_additional_msigdb_filt_', as.character(rsq_thr), '.csv')))

