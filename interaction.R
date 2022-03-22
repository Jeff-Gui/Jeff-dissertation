library(yaml)
library(logging)
library(patchwork)
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
source('utils.R')
source('load_data_cbp.R')
gc()

# config_name = 'metabric.yaml'
# tcga_stad.yaml, tcga_hnsc.yaml
# tcga_pan_nine.yaml
# ccle.yaml, tcga_lusc.yaml, tcga_blca.yaml, tcga_ov.yaml, tcga_lgg.yaml
config_name = 'tcga_pan_nine.yaml'  # tcga_luad.yaml, tcga_brca.yaml, metabric.yaml, tcga_coad.yaml
default_cfg_name = 'default.yaml'

# BATCH RUN
config_names = c('tcga_lusc.yaml', 'tcga_blca.yaml', 'tcga_ov.yaml', 'tcga_lgg.yaml',
                 'tcga_stad.yaml', 'tcga_luad.yaml', 'tcga_brca.yaml', 'tcga_coad.yaml',
                 'tcga_hnsc.yaml')

for (config_name in config_names) {
  # refresh_log = TRUE
  save_intermediate = FALSE  # may spend extra time
  use_cache = TRUE
  use_cache_geno_pca = TRUE
  source = TRUE
  
  ## Handle config
  default_cfg = yaml.load_file(file.path('config', default_cfg_name))
  if (config_name != default_cfg_name){
    config = yaml.load_file(file.path('config', config_name))
    config = merge_cfg(default_cfg, config)
  }
  config = pcs_cfg(config)
  dt_cfg = config$dataset
  eqtl_cfg = config$eQTL
  preprocess_cfg = config$preprocess
  # write_yaml(config, file.path(config$output, 'config.yaml'))
  
  ## Handle log
  # logpath = file.path(config$output, 'log.txt')
  # if (file.exists(logpath) & refresh_log){
  #   file.remove(logpath)
  # }
  # basicConfig(level = 'FINEST')
  # addHandler(writeToFile, file=logpath, level='DEBUG')
  # loginfo('Loading data...', logger = 'main')
  
  ## Load data
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
                       filter_protein_coding = eqtl_cfg$genepos,
                       diploid_norm = preprocess_cfg$diploid_norm)
    gc()
    if (save_intermediate){
      save(dt, file = file.path(dirname(dt_cfg$dataset_home), 'clean_data.RData'))
    }
  } else {
    load(file = file.path(dirname(dt_cfg$dataset_home), 'clean_data.RData'))
  }
  
  ## Prepare eQTL data
  eqtl_m = get_eQTL_m(dt, genes = strsplit(eqtl_cfg$genes, split = ',')[[1]],
                      sample_col_name = eqtl_cfg$maf_sample_col_nm,
                      min_sample_per_snp = eqtl_cfg$min_vaf,
                      meta_mut_fp = eqtl_cfg$meta_mut,
                      mode = eqtl_cfg$mode)
  snps = as.matrix(eqtl_m[[2]])
  snps[snps > 0] = 1 # binarize
  p53_exp = assay(dt[[1]]['TP53',,'RNA'])
  expr = assay(dt[[1]][,,'RNA'])
  
  nprobe = nrow(assay(dt[[1]][,,'RNA']))
  coll = list()
  hit_count = 0
  for (i in 1:nprobe){
    df_probe = as.data.frame(cbind(t(p53_exp), t(expr[i,])))
    colnames(df_probe) = c('p53', 'gene')
    df_probe = cbind(df_probe, t(snps))
    for (j in 3:ncol(df_probe)){
      snp_nm = colnames(df_probe)[j]
      form = as.formula(paste('gene~p53*', snp_nm))
      mod = summary(lm(form, data = df_probe))
      coe = as.data.frame(mod$coefficients)
      #if (sum(coe$`Pr(>|t|)`<0.05)>0){
      if (coe$`Pr(>|t|)`[4] < 0.05 & 
          coe$`Pr(>|t|)`[3] < 0.05){ # significant interaction
        int_nm = paste('p53:',snp_nm,sep='')
        coe['gene'] = rownames(expr)[i]
        coe['term'] = c('intercept', 'p53', snp_nm, int_nm)
        coll[[paste(rownames(expr)[i], int_nm, sep=':')]] = coe
        hit_count = hit_count + 1
      }
    }
    ggplot(df_probe) +
      geom_point(aes(x=p53, y=gene, color=as.factor(hot_spot)))
    if (i %% 1000 ==0){
      print(paste(i,hit_count,sep=', '))
    }
  }
  
}