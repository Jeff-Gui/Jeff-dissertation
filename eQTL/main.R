library(MatrixEQTL)
library(yaml)
library(logging)
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
source('utils.R')
source('load_data_cbp.R')

config_name = 'metabric.yaml'
config_name = 'tcga_brca.yaml'

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
write_yaml(config, file.path(config$output, 'config.yaml'))


basicConfig(level = 'FINEST')
addHandler(writeToFile, file=file.path(config$output,'log.txt'), level='DEBUG')
loginfo('Loading data...', logger = 'main')

## Load data
dt = load_data_cbp(dataset_home = dt_cfg$dataset_home,
                   exp_nm = dt_cfg$exp_nm, cna_nm = dt_cfg$cna_nm, mut_nm = dt_cfg$mut_nm,
                   sample_meta_nm = dt_cfg$sample_meta_nm,
                   case_complete_nm = dt_cfg$case_complete_nm,
                   case_list_dir_nm = dt_cfg$case_list_dir_nm,
                   na.str = dt_cfg$na.str,
                   gene_col_nm = dt_cfg$gene_col_nm)

## Quality control step

## Prepare eQTL data
eqtl_m = get_eQTL_m(dt, genes = strsplit(eqtl_cfg$genes, split = ',')[[1]],
                    sample_col_name = eqtl_cfg$maf_sample_col_nm,
                    min_sample_per_snp = eqtl_cfg$min_vaf)
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
  cvrt = SlicedData$new()$CreateFromMatrix(cvrt)
}


#### Error covariance matrix
#### Set to numeric() for identity.
errorCovariance = numeric()

## Run eQTL
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

trans_eqtls = me$trans$eqtls

plot_df = get_single_eqtl_plot_dt('CCNA2', 'snp7', dt, 
                                  as.matrix(snps), lookup = eqtl_m[[3]])[[1]]
ggplot(plot_df, aes(x=class, y=expression)) + theme_classic() +
  geom_violin() +
  geom_jitter(width = 0.2, alpha=0.5)
