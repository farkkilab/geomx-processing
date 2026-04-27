library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggpmisc)

clin_path <- "~/Documents/phd/st/data/geomx/clinical_data/9_eyemt_patient_clinical_data_SENSITIVE_upd_0426.csv"

out_dir <- "~/Documents/phd/st/geomx-processing/results/clinical"
dir.create(out_dir)

source(file.path('~/Documents/phd/st/geomx-processing/src/geomx_utils.R'))

# load files and select vars ----------------------------------------------

clin <- as.data.frame(fread(clin_path))
colnames(clin)

vars_survival <- c('PFS_days', 'PFI_days', 'OS_days')
vars_labels <- c('stage', 'primary_surgery_residual', 'HRP_status', 'BRCA_status', 'primary_treatment_response')
vars_cont <- c('TMB', 'age_at_diagnosis', 'ovaHRDscar_score')

# change R0 vals into R=1 for both R<1 and r>0
clin$primary_surgery_residual <- gsub('<|>', '=', clin$primary_surgery_residual)

# filter to unique values per patient
clin <- clin[, c('Patient', 'deceased', vars_survival, vars_labels, vars_cont)] %>%
  distinct()


# PFS/PFI/OS + TMB/age correlations coloured by label variables ---------------------

vars_comb <- combn(c(vars_survival, vars_cont), 2, simplify = F)

#surv_vals <- surv_comb[[1]]
#col_var <- vars_labels[1]

for(vars in vars_comb){
  for(col_var in vars_labels){
    
    ggplot(data = clin, aes(x = get(vars[1]), get(vars[2]))) +
      geom_point(aes(color = get(col_var))) + 
      ggtitle(paste0(vars[1], ' vs ', vars[2])) + 
      xlab(vars[1]) +
      ylab(vars[2]) +
      labs(color=col_var) +
      geom_smooth(method='lm', formula= y~x, linewidth=0.5, color = 'black') +
      stat_poly_eq(use_label(c("R2")))
    
    ggsave(file.path(out_dir, paste0('scatter_corr_', vars[1], '_vs_', vars[2], '_', col_var, '.png')))
  } 
}


# PFS/PFI/OS + TMB/age boxplots with labels vars ------------------------------------

# surv_var <- vars_survival[1]
# col_var <- vars_labels[1]

for(surv_var in c(vars_survival, vars_cont)){
  for(col_var in vars_labels){
    
    ggplot(data = clin, aes(x = get(col_var), y = get(surv_var), fill = as.factor(get(col_var)))) +
      geom_boxplot() +
      geom_point(position= position_jitterdodge(dodge.width = 1, jitter.width= .3, jitter.height = 0),
                 size= 1, alpha = 0.6) +
      stat_summary(fun = "mean", geom = "point", colour = "red", position = position_dodge(0.9), size=0.3) +
      geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = FALSE, size = 0.2, label.size = 2.8) +
      theme(axis.text.x = element_text(angle=45, hjust=1, size = 5)) +
      ggtitle(paste0(surv_var, ' per ', col_var))+
      xlab(col_var) +
      ylab(surv_var) +
      labs(fill=col_var)
    
    ggsave(file.path(out_dir, paste0('boxpl_wilcox_', surv_var, '_', col_var, '.png')))
  }
}
