library(yaml)
library(logging)
library(patchwork)
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
source('utils.R')
source('load_data_cbp.R')
gc()

HKGs = c('C1orf43', 'CHMP2A', 'EMC7', 'GPI', 'PSMB2', 
         'PSMB4', 'RAB7A', 'REEP5', 'SNRPD3','VCP', 'VPS29')

## Handle config
default_cfg_name = 'default_vsWT.yaml'
default_cfg = yaml.load_file(file.path('config', default_cfg_name))
config_name = 'tcga_brca.yaml'
if (config_name != default_cfg_name){
  config = yaml.load_file(file.path('config', config_name))
  config = merge_cfg(default_cfg, config)
}
config = pcs_cfg(config)
dt_cfg = config$dataset
eqtl_cfg = config$eQTL
preprocess_cfg = config$preprocess
noNorm.dt = load_data_cbp(dataset_home = dt_cfg$dataset_home,
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
                   diag_out = NULL,
                   filter_protein_coding = eqtl_cfg$genepos,
                   diploid_norm = preprocess_cfg$diploid_norm)

hasNorm.dt = load_data_cbp(dataset_home = dt_cfg$dataset_home,
                          exp_nm = dt_cfg$exp_nm, 
                          cna_nm = dt_cfg$cna_nm, 
                          mut_nm = dt_cfg$mut_nm,
                          sample_meta_nm = dt_cfg$sample_meta_nm,
                          patient_meta_nm = dt_cfg$patient_meta_nm,
                          case_complete_nm = dt_cfg$case_complete_nm,
                          case_list_dir_nm = dt_cfg$case_list_dir_nm,
                          na.str = dt_cfg$na.str,
                          gene_col_nm = dt_cfg$gene_col_nm,
                          normalize_genes = HKGs,
                          quantile_norm = preprocess_cfg$quantile,
                          z_score = preprocess_cfg$z_score,
                          has_log_ed = dt_cfg$has_log_ed,
                          rm_low_expr_gene = as.numeric(preprocess_cfg$rm_low_expr_gene),
                          diag_out = NULL,
                          filter_protein_coding = eqtl_cfg$genepos,
                          diploid_norm = preprocess_cfg$diploid_norm)

## Plot HKG in noNorm ===
library(pheatmap)
noNorm.dt[[1]]@colData$p53_state = 'Wildtype'
hasNorm.dt[[1]]@colData$p53_state = 'Wildtype'
ann = annotate_sample_mut(noNorm.dt[[2]]@data)
for (i in names(ann)){
  noNorm.dt[[1]]@colData[ann[[i]], 'p53_state'] = i
  hasNorm.dt[[1]]@colData[ann[[i]], 'p53_state'] = i
}
col_ann = noNorm.dt[[1]]@colData$p53_state
names(col_ann) = rownames(noNorm.dt[[1]]@colData)
col_ann = sort(col_ann)
col_ann = data.frame('p53_state'=col_ann, row.names = names(col_ann))
hkg_exp_noNorm = assay(noNorm.dt[[1]][HKGs,rownames(col_ann),'RNA'])
hkg_exp_hasNorm = assay(hasNorm.dt[[1]][HKGs,rownames(col_ann),'RNA'])
pheatmap(hkg_exp, annotation_col = col_ann, cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = F)
