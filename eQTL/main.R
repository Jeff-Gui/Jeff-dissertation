library(MatrixEQTL)
library(yaml)
setwd('/Users/jefft/Desktop/p53_project/eQTL')
source('utils.R')
source('load_data_cbp.R')

config_name = 'metabric.yaml'

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

## Load data
dt = load_data_cbp(dataset_home = dt_cfg$dataset_home,
                   exp_nm = dt_cfg$exp_nm, cna_nm = dt_cfg$cna_nm, mut_nm = dt_cfg$mut_nm,
                   sample_meta_nm = dt_cfg$sample_meta_nm,
                   case_complete_nm = dt_cfg$case_complete_nm,
                   case_list_dir_nm = dt_cfg$case_list_dir_nm,
                   na.str = dt_cfg$na.str,
                   gene_col_nm = dt_cfg$gene_col_nm)

## Quality control step
### VAF

## Prepare eQTL data
eqtl_m = get_eQTL_m(dt, genes = strsplit(eqtl_cfg$genes, split = ',')[[1]],
                    sample_col_name = eqtl_cfg$maf_sample_col_nm,
                    min_sample_per_snp = 20)
snps = eqtl_m$snp
gene = eqtl_m$gene
snpspos = eqtl_m$snp_lookup
genepos = read.table(eqtl_cfg$genepos, sep='\t', header = T)  # TODO
genepos = genepos[,c('geneid', 'chr', 'left', 'right')]

### Covariate
# TODO Maybe tailored to the analysis
cov_from_meta_col = c('ER_STATUS')
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

#### Error covariance matrix
#### Set to numeric() for identity.
errorCovariance = numeric()

## Run eQTL
me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name = eqtl_cfg$output_file_name_tra,
  pvOutputThreshold = as.numeric(eqtl_cfg$pvOutputThreshold_tra),
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = eqtl_cfg$output_file_name_cis,
  pvOutputThreshold.cis = as.numeric(eqtl_cfg$pvOutputThreshold_cis),
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = as.numeric(eqtl_cfg$cisDist),
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)


plot_df = get_single_eqtl_plot_dt('CHEK2', 'snp9', dt, 
                                  as.matrix(snps), lookup = lookup)[[1]]
ggplot(plot_df, aes(x=class, y=expression)) + theme_classic() +
  geom_boxplot()
