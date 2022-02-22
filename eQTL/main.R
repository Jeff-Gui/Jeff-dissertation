library(MatrixEQTL)
library(yaml)
library(logging)
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
source('utils.R')
source('load_data_cbp.R')

# config_name = 'metabric.yaml'
config_name = 'tcga_brca.yaml'
refresh_log = TRUE
use_cache = FALSE

## Handle config
default_cfg = yaml.load_file(file.path('config', 'default.yaml'))
if (config_name != 'dafault.yaml'){
  config = yaml.load_file(file.path('config', config_name))
  config = merge_cfg(default_cfg, config)
}
dt_cfg = config$dataset
if (is.null(dt_cfg$na.str)){
  dt_cfg$na.str = ''
}
eqtl_cfg = config$eQTL
map = c(modelLINEAR, modelANOVA, modelLINEAR_CROSS)
names(map) = c('LINEAR', 'ANOVA', 'LINEAR_CROSS')
useModel = map[eqtl_cfg$model]
if (is.null(eqtl_cfg$output_file_name_tra)){
  eqtl_cfg$output_file_name_tra = tempfile()
}
if (is.null(eqtl_cfg$output_file_name_cis)){
  eqtl_cfg$output_file_name_cis = tempfile()
}
if (!dir.exists(config$output)){
  dir.create(config$output)
}
if (eqtl_cfg$meta_mut == '/'){
  eqtl_cfg$meta_mut = NULL
}
write_yaml(config, file.path(config$output, 'config.yaml'))

logpath = file.path(config$output, 'log.txt')
if (file.exists(logpath) & refresh_log){
  file.remove(logpath)
}

basicConfig(level = 'FINEST')
addHandler(writeToFile, file=logpath, level='DEBUG')
loginfo('Loading data...', logger = 'main')

preprocess_cfg = config$preprocess
if (!is.null(preprocess_cfg$norm_gene)){
  preprocess_cfg$norm_gene = strsplit(preprocess_cfg$norm_gene, split = ',')[[1]]
}

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
                     has_log_ed = dt_cfg$has_log_ed)
  gc()
} else {
  load(file = file.path(dirname(dt_cfg$dataset_home), 'clean_data.RData'))
}


## Quality control step


## Prepare eQTL data
eqtl_m = get_eQTL_m(dt, genes = strsplit(eqtl_cfg$genes, split = ',')[[1]],
                    sample_col_name = eqtl_cfg$maf_sample_col_nm,
                    min_sample_per_snp = eqtl_cfg$min_vaf,
                    meta_mut_fp = eqtl_cfg$meta_mut)
snps = eqtl_m$snp
gene = eqtl_m$gene
snpspos = eqtl_m$snp_lookup
genepos = read.table(eqtl_cfg$genepos, sep='\t', header = T)  # TODO
genepos = genepos[,c('geneid', 'chr', 'left', 'right')]

### Covariate
# TODO Maybe tailored to the analysis
if (is.null(eqtl_cfg$covariate_from_meta)){
  cvrt = SlicedData$new()
} else {
  cov_from_meta_col = strsplit(eqtl_cfg$covariate_from_meta, split = ',')[[1]]
  meta = data.frame(dt[[1]]@colData)
  cvrt = matrix(0, nrow = length(cov_from_meta_col), ncol = nrow(meta))
  colnames(cvrt) = rownames(meta)
  rownames(cvrt) = cov_from_meta_col
  for (cl in cov_from_meta_col){
    if (class(meta[,cl]) != 'numeric'){
      cvrt[cl,] = as.numeric(as.factor(meta[,cl])) - 1
    } else {
      cvrt[cl,] = meta[,cl]
    }
  }
  remove_cov = c()
  for (i in 1:nrow(cvrt)){
    if (length(unique(cvrt[i,]))==1){
      remove_cov = c(remove_cov, i)
      loginfo(logger = 'main', 'Covariate %s is not varying in the dataset.', rownames(cvrt)[i])
    }
  }
  if (length(remove_cov) > 0){
    cvrt = cvrt[-remove_cov,]
  }
  cvrt = SlicedData$new()$CreateFromMatrix(cvrt)
}


#### Error covariance matrix
#### Set to numeric() for identity.
errorCovariance = numeric()

## Run eQTL
gc()
me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name = file.path(config$output, eqtl_cfg$output_file_name_tra),
  pvOutputThreshold = as.numeric(eqtl_cfg$pvOutputThreshold_tra),
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = file.path(config$output, eqtl_cfg$output_file_name_cis),
  pvOutputThreshold.cis = as.numeric(eqtl_cfg$pvOutputThreshold_cis),
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = as.numeric(eqtl_cfg$cisDist),
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

# Checking the results
trans_eqtls = me$trans$eqtls
trans_eqtls = subset(trans_eqtls, trans_eqtls$FDR < 0.05)
snpids = unique(trans_eqtls$snps)
snpp = c()
for (i in snpids){
  if (length(grep('snp', i))>0){
    snpp = c(snpp, get_var_info_from_maf(dt[[2]]@data, i, eqtl_m[[3]]))
  } else {
    snpp = c(snpp, i)
  }
}
names(snpp) = snpids
trans_eqtls['protein_change'] = snpp[trans_eqtls$snps]
trans_eqtls = trans_eqtls[order(trans_eqtls$protein_change, trans_eqtls$FDR),]
write.table(trans_eqtls,file.path(config$output, 'trans_eqtl_fdr005.txt'),
            row.names = F, sep = '\t', quote = F)
gc()

## Visualise the results
snpid_to_plot = 'tissue_high'
gene_to_plot = 'IL28B'

if (length(grep('snp', snpid_to_plot))>0){
  plot_rs = get_single_eqtl_plot_dt(gene_to_plot, snpid_to_plot, dt,
                                    as.matrix(snps), lookup = eqtl_m[[3]])
  plot_title = strsplit(get_var_info_from_maf(dt[[2]]@data, snpid_to_plot, eqtl_m[[3]]),
                        split = '\\.')[[1]][2]
} else {
  plot_rs = get_single_eqtl_plot_dt(gene_to_plot, snpid_to_plot, dt,
                                    as.matrix(snps))
  plot_title = snpid_to_plot
}
plot_df = plot_rs[[1]]

ggplot(plot_df, aes(x=class, y=expression)) + theme_classic() +
  geom_violin() +
  geom_jitter(width = 0.2, alpha=0.5) +
  labs(x=plot_title, title = strsplit(config_name, split='\\.')[[1]][1], y=gene_to_plot)

# # positive controls
# pos_ctr_eqtl = subset(trans_eqtls, trans_eqtls$gene %in% 
#                         c('CCNA2', 'CHEK1', 'PSMA1', 'PSMC1', 'GMPS', 'DHFR', 'IMPDH'))
# 
# # Check distribution of expression
# expr_test = as.matrix(gene)
# expr_test = expr_test[,sample(colnames(expr_test), 10)]
# ggplot(melt(expr_test), aes(x=Var2, y=value)) + geom_boxplot()
# 
