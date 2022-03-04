setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
source('utils.R')
source('load_data_cbp.R')
library(tidyverse)
config_name = 'ccle_inspect_only.yaml'
refresh_log = TRUE
save_intermediate = FALSE  # may spend extra time
use_cache = FALSE
use_cache_geno_pca = TRUE
source = TRUE

## Handle config
default_cfg = yaml.load_file(file.path('config', 'default.yaml'))
if (config_name != 'dafault.yaml'){
  config = yaml.load_file(file.path('config', config_name))
  config = merge_cfg(default_cfg, config)
}
config = pcs_cfg(config)
dt_cfg = config$dataset
eqtl_cfg = config$eQTL
preprocess_cfg = config$preprocess

# Load data
if (!use_cache){
  dt = load_data_cbp(dataset_home = dt_cfg$dataset_home,
                     exp_nm = dt_cfg$exp_nm, 
                     cna_nm = dt_cfg$cna_nm, 
                     mut_nm = dt_cfg$mut_nm,
                     sample_meta_nm = dt_cfg$sample_meta_nm,
                     patient_meta_nm = dt_cfg$patient_meta_nm,
                     case_complete_nm = dt_cfg$case_complete_nm,
                     case_list_dir_nm = dt_cfg$case_list_dir_nm,
                     na.str = dt_cfg$na.str,
                     gene_col_nm = dt_cfg$gene_col_nm,
                     normalize_genes = preprocess_cfg$norm_gene,
                     quantile_norm = preprocess_cfg$quantile,
                     z_score = preprocess_cfg$z_score,
                     has_log_ed = dt_cfg$has_log_ed,
                     rm_low_expr_gene = as.numeric(preprocess_cfg$rm_low_expr_gene),
                     diag_out = config$output,
                     filter_protein_coding = eqtl_cfg$genepos)
  gc()
  if (save_intermediate){
    save(dt, file = file.path(dirname(dt_cfg$dataset_home), 'clean_data_no_norm.RData'))
  }
} else {
  load(file = file.path(dirname(dt_cfg$dataset_home), 'clean_data_no_norm.RData'))
}
