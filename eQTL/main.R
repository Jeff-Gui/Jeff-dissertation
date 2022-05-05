library(MatrixEQTL)
library(yaml)
library(logging)
library(patchwork)
library(getopt)
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
source('utils.R')
source('load_data_cbp.R')
gc()

# config_name = 'metabric.yaml'
# tcga_stad.yaml, tcga_hnsc.yaml
# tcga_pan_nine.yaml
# ccle.yaml, tcga_lusc.yaml, tcga_blca.yaml, tcga_ov.yaml, tcga_lgg.yaml
config_name = 'tcga_pan_nine.yaml'  # tcga_luad.yaml, tcga_brca.yaml, metabric.yaml, tcga_coad.yaml
# default_ult.yaml default_vsNull.yaml default_vsWT.yaml
default_cfg_name = 'default_vsWT.yaml'
load('/Users/jefft/Desktop/p53_project/datasets/UMD_TP53/SIFT_UMD.RData')

config_name = 'tcga_xena_brca.yaml'

# BATCH RUN
config_names = c('tcga_lusc.yaml', 'tcga_blca.yaml', 'tcga_ov.yaml', 'tcga_lgg.yaml',
                 'tcga_stad.yaml', 'tcga_luad.yaml', 'tcga_brca.yaml', 'tcga_coad.yaml',
                 'tcga_hnsc.yaml')

for (config_name in config_names) {
  print(config_name)
  refresh_log = TRUE
  save_intermediate = FALSE  # may spend extra time
  use_cache = TRUE # cached data (clean_data.RData)
  use_cache_geno_pca = TRUE # cached genotype PCs
  source = TRUE
  check_cna_wt = FALSE
  
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
  write_yaml(config, file.path(config$output, 'config.yaml'))
  
  ## Handle log
  logpath = file.path(config$output, 'log.txt')
  if (file.exists(logpath) & refresh_log){
    file.remove(logpath)
  }
  basicConfig(level = 'FINEST')
  addHandler(writeToFile, file=logpath, level='DEBUG')
  loginfo('Loading data...', logger = 'main')
  
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
                      mode = eqtl_cfg$mode,
                      SIFT = sift_pred)  # different mode may lead to subsampling
  snps = eqtl_m$snp
  dt[[1]] = dt[[1]][,colnames(as.matrix(snps)),]
  gene = eqtl_m$gene
  snpspos = eqtl_m$snp_lookup
  genepos = read.table(eqtl_cfg$genepos, sep='\t', header = T)  # TODO
  genepos = genepos[,c('geneid', 'chr', 'left', 'right')]
  
  ### Covariate
  # TODO Maybe tailored to the analysis
  remove_cov = c()
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
    for (i in 1:nrow(cvrt)){
      if (length(unique(cvrt[i,]))==1){
        remove_cov = c(remove_cov, i)
        loginfo(logger = 'main', 'Covariate %s is not varying in the dataset.', rownames(cvrt)[i])
      }
    }
  }
  if (!is.null(eqtl_cfg$covariate_code)){
    eqtl_cfg$covariate_code = strsplit(eqtl_cfg$covariate_code, split = ',')[[1]]
    if ('GENO' %in% eqtl_cfg$covariate_code){
      if (use_cache_geno_pca){
        gpca = read.table(file.path(dirname(dt_cfg$dataset_home), 'genotype_PC20.txt'), 
                          sep='\t', header=T)
      } else {
        gpca = genotype_pca(dt[[2]]@data, out_dic = dirname(dt_cfg$dataset_home), 
                            samples = rownames(dt[[1]]@colData), 
                            file_out = save_intermediate,
                            gene_col_name = dt_cfg$gene_col_nm, 
                            sample_col_name = eqtl_cfg$maf_sample_col_nm)
      }
      gpca = t(gpca[,1:5]) # Only using top five PCs.
      gpca = gpca[,colnames(cvrt)]
      cvrt = rbind(cvrt, gpca)
    }
  }
  if (length(remove_cov) > 0){
    cvrt = cvrt[-remove_cov,]
  }
  cvrt = SlicedData$new()$CreateFromMatrix(cvrt)
  
  #### Error covariance matrix
  #### Set to numeric() for identity.
  errorCovariance = numeric()
  
  ## Run eQTL
  map = c(modelLINEAR, modelANOVA, modelLINEAR_CROSS)
  names(map) = c('LINEAR', 'ANOVA', 'LINEAR_CROSS')
  useModel = map[eqtl_cfg$model]
  gc()
  
  ann = annotate_sample_mut(dt[[2]]@data)
  snps = as.matrix(snps)
  if (eqtl_cfg$mode=='mut-VS-null'){
    wt_spl = colnames(snps)[colnames(snps) %in% ann$nonsense]
  } else {
    wt_spl = colnames(snps)[!colnames(snps) %in% unlist(ann)]
    if (check_cna_wt){
      orln = length(wt_spl)
      cna = read.table(file.path(config$dataset$dataset_home,config$dataset$cna_nm), sep='\t', header = T)
      norm_wt = colnames(cna)[which(cna[cna$Hugo_Symbol=='TP53',]==0)]
      norm_wt = gsub('\\.', '-', norm_wt)
      wt_spl = intersect(wt_spl, norm_wt)
      loginfo(logger = 'main', 'Filtering %d copy-number normal wildtype p53 out of %d.', length(wt_spl), orln)
      if (length(wt_spl)<10){
        loginfo(logger = 'main', 'Too few WT, skipping this dataset.')
        next
      }
    }
  }
  gene = as.matrix(gene)
  cvrt = as.matrix(cvrt)
  loginfo(logger = 'main', paste('Using coveriates: ', rownames(cvrt), sep=''))
  
  cis_qtls_coll = list()
  trans_qtls_coll = list()
  for (i in 1:nrow(snps)){
    loginfo(logger = 'main', 'Processing mutant %s', rownames(snps)[i])
    mis_spl = colnames(snps)[which(snps[i,]!=0)]
    mis_spl = mis_spl[!mis_spl %in% ann$others]
    loginfo(logger = 'main', 'Mutant # %s, Control # %s', length(mis_spl), length(wt_spl))
    all_spl = c(wt_spl, mis_spl)
    sub_snps = t(as.matrix(snps[i,all_spl]))
    rownames(sub_snps) = rownames(snps)[i]
    
    me = Matrix_eQTL_main(
      snps = SlicedData$new()$CreateFromMatrix(sub_snps), 
      gene = SlicedData$new()$CreateFromMatrix(gene[,all_spl]), 
      cvrt = SlicedData$new()$CreateFromMatrix(cvrt[,all_spl]),
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
    cis_eqtls = me$cis$eqtls
    
    if (nrow(trans_eqtls) > 0){
      snpids = unique(trans_eqtls$snps)
      snpp = c()
      for (j in snpids){
        if (length(grep('snp', j))>0){
          snpp = c(snpp, get_var_info_from_maf(dt[[2]]@data, j, eqtl_m[[3]]))
        } else {
          snpp = c(snpp, j)
        }
      }
      names(snpp) = snpids
      trans_eqtls['protein_change'] = snpp[trans_eqtls$snps]
      
      before_flt = nrow(trans_eqtls)
      trans_eqtls = trans_eqtls[order(trans_eqtls$protein_change, trans_eqtls$FDR),]
      trans_qtls_coll[[i]] = trans_eqtls
      cis_qtls_coll[[i]] = cis_eqtls
      loginfo(logger = 'main', 'Identified %d genes, %d passed the FDR filter.', before_flt, sum(trans_eqtls$FDR < 0.05))
    }
  }
  trans_qtls_coll = Reduce(rbind, trans_qtls_coll)
  cis_qtls_coll = Reduce(rbind, cis_qtls_coll)
  if (nrow(trans_qtls_coll) > 0){
    write.table(trans_qtls_coll,file.path(config$output, 'trans_eqtl.txt'),
                row.names = F, sep = '\t', quote = F)
  }
  if (nrow(cis_qtls_coll) > 0){
    write.table(cis_qtls_coll,file.path(config$output, 'cis_eqtl.txt'),
                row.names = F, sep = '\t', quote = F)
  }
  loginfo(logger = 'main', 'In total, identified %d genes, %d passed the FDR filter.', nrow(trans_qtls_coll), sum(trans_qtls_coll$FDR < 0.05))
  trans_qtls_coll = subset(trans_qtls_coll, trans_qtls_coll$FDR < 0.05)
  if (nrow(trans_qtls_coll) > 0){
    write.table(trans_qtls_coll,file.path(config$output, 'trans_eqtl_fdr005.txt'),
                row.names = F, sep = '\t', quote = F)
  }
  gc()
}

if (!source){
  trans_eqtls = read.table('/Users/jefft/Desktop/p53_project/scripts/eQTL/outputs/tcga_brca_raw_seq/trans_eqtl_fdr005.txt', 
                           header=T, sep='\t')
  library(ggpubr)
  ## Visualise the results
  # overall
  ggplot(trans_eqtls, aes(x=beta, y=1/FDR)) + theme_classic() +
    geom_point(size=0.5) +
    scale_y_log10() +
    facet_wrap(~protein_change, ncol=2)
  
  # per-gene
  snpid_to_plot = 'breast_w2016'
  gene_to_plot = 'MEPE'
  
  plot_df = get_single_eqtl_plot_dt(gene_to_plot, snpid_to_plot, dt, as.matrix(snps))
  plot_title = unique(trans_eqtls$protein_change[which(trans_eqtls$snps == snpid_to_plot)])
  
  # CCLE only
  plot_df$site = dt[[1]]@colData[rownames(plot_df), 'PRIMARY_SITE']
  
  a = ggplot(plot_df, aes(x=as.logical(mut_state), y=expression)) + theme_classic() +
    geom_violin(width = 0.5) +
    # geom_boxplot(width = 0.05, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha=0.5) +
    labs(x=plot_title, title = strsplit(config_name, split='\\.')[[1]][1], y=gene_to_plot)
  
  ggplot(plot_df, aes(x=mut_cls, y=expression)) + theme_classic() +
    geom_violin(width = 0.5) +
    geom_jitter(width = 0.2, alpha=0.5) +
    labs(x=plot_title, title = strsplit(config_name, split='\\.')[[1]][1], y=gene_to_plot) +
    stat_compare_means(comparisons = list(c('eQTL mutant', 'frameshift'),
                                          c('eQTL mutant', 'nonsense'),
                                          c('eQTL mutant', 'other missense'),
                                          c('eQTL mutant', 'wildtype')),
                       aes(label=..p.signif..), method = 't.test')
  
  # CCLE only
  a + facet_wrap(~site, ncol=5)
  
  if (sum(plot_df$mut_state) != sum(as.logical(plot_df$mut_state))){
    b = ggplot(subset(plot_df, mut_state != 0), aes(x=mut_state, y=expression)) + 
      theme_classic() +
      geom_point() +
      labs(x='VAF', y=gene_to_plot)
    a + b
  } else {
    a
  }
  
  # ggplot(plot_df) + theme_classic() +
  #   geom_point(data=subset(plot_df, mut_state!=0), aes(x=mut_state, y=expression)) +
  #   geom_violin(data=subset(plot_df, mut_state==0), aes(x=0, y=expression), width = 0.1) +
  #   geom_boxplot(data=subset(plot_df, mut_state==0),aes(x=0, y=expression), width = 0.05, outlier.shape = NA) +
  #   geom_vline(xintercept = 0.07, linetype='dotted') +
  #   # geom_jitter(data=subset(plot_df, mut_state==0),aes(x=0, y=expression), width = 0.1, alpha=0.5) +
  #   labs(x=plot_title, title = strsplit(config_name, split='\\.')[[1]][1], y=gene_to_plot)
  
  # positive controls
  pos_ctr_eqtl = subset(trans_eqtls, trans_eqtls$gene %in%
                          c('CCNA2', 'CHEK1', 'PSMA1', 'PSMC1', 'GMPS', 'DHFR', 'IMPDH'))
  
  # Check distribution of expression
  expr_test = as.matrix(gene)
  expr_test = expr_test[,sample(colnames(expr_test), 10)]
  ggplot(melt(expr_test), aes(x=Var2, y=value)) + geom_boxplot()
}




