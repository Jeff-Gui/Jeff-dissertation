library(logging)
library(tidyverse)
library(maftools)
#library(limma)
library(data.table)
library(preprocessCore)
library(RNOmni)
# Ref: https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html
library(MultiAssayExperiment)
# Ref: https://bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/MultiAssayExperiment.html
# setwd('/Users/jefft/Desktop/p53_project')
# 
# dataset_home = '/Users/jefft/Desktop/p53_project/datasets/CCLE/ccle_broad_2019'
# dataset_home = '/Users/jefft/Desktop/p53_project/datasets/LUAD-TCGA/luad_tcga_pan_can_atlas_2018'
# case_complete_nm = 'cases_complete.txt'
# case_complete_nm = 'cases_all.txt'
# #exp_nm = 'data_mrna_seq_v2_rsem.txt'
# exp_nm = 'data_mrna_seq_rpkm.txt'
# exp_nm = 'data_mrna_seq_v2_rsem_zscores_ref_normal_samples.txt'
# cna_nm = 'data_cna.txt'
# mut_nm = 'data_mutations.txt'
# gene_col_nm = 'Hugo_Symbol'
# sample_meta_nm = 'data_clinical_sample.txt'
# patient_meta_nm = 'data_clinical_patient.txt'
# case_list_dir_nm = 'case_lists'
# na.str = ''
# has_log_ed = FALSE
# rm_low_expr_gene = NULL
# normalize_genes = NULL
# quantile_norm = TRUE
# z_score = FALSE
# diag_out = NULL
# filter_protein_coding = NULL

clean_matrix = function(mtx, complete_cases_dot, gene_col_nm){
  mtx = mtx[!is.na(mtx[gene_col_nm]),]
  # remove dup by taking the record with the largest mean
  dup_genes = unique(mtx[[gene_col_nm]][which(duplicated(mtx[[gene_col_nm]]))])
  remove_rows = c()
  for (i in dup_genes){
    sub = subset(mtx, mtx[[gene_col_nm]]==i)
    rms = rowMeans(sub[,complete_cases_dot])
    remove_rows = c(remove_rows, rownames(sub)[which(rms < max(rms))])
  }
  loginfo(logger = 'data.loader','Duplicated rows in matrix: %d', length(remove_rows))
  mtx = mtx[which(!rownames(mtx) %in% remove_rows),]
  remain_dup = which(duplicated(mtx[[gene_col_nm]]))
  if (length(remain_dup)>0){
    mtx = mtx[-remain_dup,]  # remove secondary dup, use first record.
  }
  rownames(mtx) = mtx[[gene_col_nm]]
  mtx = mtx[,complete_cases_dot]
  
  ## NA removal (?)
  mtx = na.omit(mtx)
  
  return(mtx)
}

load_data_cbp = function(dataset_home, case_complete_nm = 'cases_complete.txt',
                         exp_nm = 'data_mrna_seq_v2_rsem.txt',
                         cna_nm = 'data_cna.txt',
                         mut_nm = 'data_mutations.txt',
                         gene_col_nm = 'Hugo_Symbol',
                         sample_meta_nm = 'data_clinical_sample.txt',
                         patient_meta_nm = 'data_clinical_patient.txt',
                         case_list_dir_nm = 'case_lists',
                         na.str = '',
                         has_log_ed = TRUE,
                         rm_low_expr_gene = NULL,
                         normalize_genes = NULL,
                         quantile_norm = TRUE,
                         diploid_norm = FALSE,
                         z_score = FALSE, real_z_score = FALSE,
                         diag_out = NULL,
                         filter_protein_coding = NULL){
  if (diploid_norm){
    loginfo(logger = 'data.loader', 'Planning to normalize sample to diploid ones, loading CNA.')
  } else {
    loginfo(logger = 'data.loader', 'Not loading CNA data.')
  }
  patient_id_col_nm = 'PATIENT_ID'
  complete_cases = read.table(file.path(dataset_home, case_list_dir_nm, case_complete_nm), sep=':')
  complete_cases = trimws(strsplit(complete_cases[nrow(complete_cases),2], '\t')[[1]])
  
  exps = fread(file.path(dataset_home, exp_nm), sep='\t', header = T, na.strings = na.str)
  exps = as.data.frame(exps)
  loginfo(logger ='data.loader', 
          'Expression matrix dim before cleaning: %d genes X %d samples.', dim(exps)[1], dim(exps)[2])
  complete_cases = intersect(colnames(exps), complete_cases)  # sometimes, complete cases may not account for expression matrix. 
  
  sample_meta = read.table(file.path(dataset_home, sample_meta_nm), sep = '\t', comment.char = '#', header = T, quote = '\"', fill=T)
  patient_meta = read.table(file.path(dataset_home, patient_meta_nm), sep = '\t', comment.char = '#', header = T, quote = '\"', fill=T)
  sample_meta = left_join(sample_meta, patient_meta, by = patient_id_col_nm)
  
  exps = clean_matrix(exps, complete_cases_dot = complete_cases, gene_col_nm = gene_col_nm)
  # filter out protein coding
  if (!is.null(filter_protein_coding)){
    genepos = read.table(filter_protein_coding, sep='\t', header = T)
    exps = exps[which(rownames(exps) %in% genepos$geneid),]
  }
  
  if (diploid_norm){
    cnas = fread(file.path(dataset_home, cna_nm), sep='\t', header = T, na.strings = na.str)
    cnas = as.data.frame(cnas)
    complete_cases = intersect(complete_cases, colnames(cnas))
    cnas = clean_matrix(cnas, complete_cases_dot = complete_cases, gene_col_nm = gene_col_nm)
    co_genes = as.character(ordered(intersect(rownames(exps), rownames(cnas))))
    exps = exps[co_genes,complete_cases]
    cnas = cnas[co_genes,complete_cases]
  }
  
  loginfo(logger = 'data.loader',
          'Expression matrix dim after cleaning: %d genes X %d samples.', dim(exps)[1], dim(exps)[2])
  
  ## Normalization and filtering low-expressing genes
  # log2(x+1) if has not been log transformed. NA must has been addressed in above cleaning step.
  # However, if doing diploid norm, no log2.
  LOG_TRANS = has_log_ed
  if (!LOG_TRANS & !diploid_norm){
    loginfo('Log2 transforming expression data...', logger = 'data.loader')
    exps = log2(exps + 1)
  }
  
  # Remove low-expressing genes
    # for normalize by diploid samples, remove low-expressing at last
  exps = remove_low_genes(exps, rm_low_expr_gene = rm_low_expr_gene, diag_out = diag_out)
  
  gene_nm = rownames(exps)
  sample_nm = colnames(exps)
  # Normalize with house-keeping genes (divide)
  if (!is.null(normalize_genes)){
    normalize_genes = normalize_genes[which(normalize_genes %in% rownames(exps))]
    if (length(normalize_genes)==0){
      loginfo('No normalize gene record!', logger = 'data.loader')
    } else {
      loginfo(logger = 'data.loader', 'Normalize gene expression with %s.', normalize_genes)
      # because has been log transformed
      cm = log2(colMeans(2**exps[normalize_genes,]))
      for (j in 1:ncol(exps)){
        exps[,j] = exps[,j] - cm[j]
      }
    }
  }
  
  # Quantile normalization
  if (quantile_norm){
    loginfo('Quantile-normalizing across samples...', logger = 'data.loader')
    # exps = normalizeBetweenArrays(exps)
    exps = normalize.quantiles(as.matrix(exps))
  }
  
  # Rank-based normalization z-score
  if (z_score){
    loginfo('Z-scaling each gene based on rank...', logger = 'data.loader')
    # exps = t(scale(t(exps)))
    exps = t(apply(exps, 1, RankNorm))
  } else {
    if (real_z_score){
      loginfo('Z-scaling each gene', logger = 'data.loader')
      exps = t(scale(t(exps)))
    }
  }
  
  # Normalize to diploid
  if (diploid_norm){
    cnas = cnas[rownames(exps),]
    exps = norm_exp_to_diploid(cnas, exps)
  }
  
  rownames(exps) = gene_nm
  colnames(exps) = sample_nm
  exps = as.data.frame(exps)
  
  if (diploid_norm){
    # exps = remove_low_genes(exps, rm_low_expr_gene = rm_low_expr_gene, diag_out = diag_out)
  }
  
  ## Handling mutation data
  muts = read.maf(file.path(dataset_home, mut_nm))
  # Calculate VAF for filtering (not in this function)
  flag = FALSE
  if (('t_ref_count' %in% colnames(muts@data)) & ('t_alt_count' %in% colnames(muts@data))){
    if (!anyNA(muts@data$t_ref_count)){
      loginfo('Found t_ref_count and t_alt_count, calculating VAF.', logger = 'data.loader')
      muts@data[['VAF']] = muts@data$t_alt_count / (muts@data$t_alt_count + muts@data$t_ref_count)
      flag = TRUE
    }
  }
  if (!flag){
    loginfo('No VAF information in the MAF file.', logger = 'data.loader')
  }
  # muts = subsetMaf(muts, tsb=complete_cases, genes = co_genes)  # IF CNA here
  muts = subsetMaf(muts, tsb=complete_cases)
  sample_meta = subset(sample_meta, sample_meta$SAMPLE_ID %in% complete_cases)
  rownames(sample_meta) = sample_meta$SAMPLE_ID
  sample_meta = sample_meta[,-which(colnames(sample_meta)=='SAMPLE_ID')]
  sample_meta = sample_meta[colnames(exps),]
  
  g_s = muts@gene.summary
  mm_n = as.numeric(g_s[which(g_s[[gene_col_nm]]=='TP53'), 'Missense_Mutation'])
  t_n = as.numeric(g_s[which(g_s[[gene_col_nm]]=='TP53'), 'total'])
  ms_n = as.numeric(g_s[which(g_s[[gene_col_nm]]=='TP53'), 'MutatedSamples'])
  loginfo(logger = 'data.loader', 'p53 mutation: %d missense out of %d mutations. %d mutated samples.', mm_n, t_n, ms_n)
  
  ## Construct multi-assay-experiment
  exprdat = SummarizedExperiment(exps)
  colnames(exprdat) = complete_cases
  #cnadat = SummarizedExperiment(cnas)
  #colnames(cnadat) = complete_cases
  multiAssay = MultiAssayExperiment(
    #list('RNA'=exprdat, 'CNA'=cnadat),
    list('RNA'=exprdat),
    sample_meta)
  
  a = ls()
  rm(list=a[which(a != 'multiAssay' & a != 'muts')])
  return(list(multiAssay, muts))
}

integrate_dts = function(datasets, z_score_norm=FALSE){
  # z_score_norm: whether normalize RNA dataset with Z-score before merging
  #   if length is one, will apply to all datasets
  #   otherwise, each index match one dataset in the list
  if (length(z_score_norm)==1){
    z_score_norm = rep(z_score_norm, length(datasets))
  }
  if (length(datasets)==1){
    return(datasets)
  }
  # Merge genes detected in all datasets
  co_genes = c()
  co_exps = c()
  for (dt in datasets){
    co_exps = c(co_exps, names(dt@ExperimentList@listData))
  }
  co_exps = unique(co_exps)
  merged = list()
  for (cep in co_exps){
    co_rnms = c()
    ct = 0
    for (dt in datasets){
      rnm = rownames(dt@ExperimentList[[cep]])
      if (ct == 0){co_rnms = c(co_rnms, rnm)} else {
        co_rnms = intersect(co_rnms, rnm)
      }
      ct = ct + 1
    }
    merged_exp = data.frame(row.names = co_rnms)
    count = 0
    for (dt in datasets){
      count = count + 1
      dt_df = data.frame(dt@ExperimentList[[cep]]@assays@data)[co_rnms,]
      dt_df = dt_df[,-which(colnames(dt_df) %in% c('group', 'group_name'))]
      if (z_score_norm[count]){
        dt_df = t(scale(t(dt_df)))
      }
      merged_exp = cbind(merged_exp, dt_df)
    }
    
    # remove any row that contains NA (?)
    merged_exp = merged_exp[complete.cases(merged_exp),]
    
    merged[[cep]] = merged_exp
  }
  
  a = ls()
  rm(list=a[which(a != 'merged')])
  return(merged)
}


align_expression = function(merged_exp, sample_factor){
  # sample_factor: group information of the sample in merged df
  
  # Step 1 use limma package
  merged_exp = normalizeBetweenArrays(merged_exp)
}


# normalize to diploid samples
norm_exp_to_diploid = function(cnas, exps){
  col_spl = colnames(cnas)
  mask = as.matrix(cnas)
  exps = as.matrix(exps)
  mask[which(mask!=0)] = NA
  mask[which(mask==0)] = 1
  mask[which(is.na(mask))] = 0
  mask = exps * mask
  mask[which(mask==0)] = NA
  mus = rowMeans(mask, na.rm = T)
  stds = rowSds(mask, na.rm = T)
  for (i in 1:nrow(exps)){
    if (!is.na(mus[i])) {
      exps[i,] = (exps[i,] - mus[i]) / stds[i]
    } else {
      exps[i,] = NA
    }
  }
  return(exps)
}


# Remove low expressing genes
remove_low_genes = function(exps, rm_low_expr_gene, diag_out=NULL){
  expr_mean = rowMeans(exps)
  if (!is.null(rm_low_expr_gene)){
    rm_id = which(expr_mean < rm_low_expr_gene)
    names(rm_id) = NULL
    if (length(rm_id)>0){
      exps = exps[-rm_id,] 
    }
  }
  expr_mean_after = rowMeans(exps)
  pt = data.frame('value'=c(expr_mean, expr_mean_after))
  pt['cls'] = c(rep('before', length(expr_mean)), rep('after', length(expr_mean_after)))
  gp = ggplot(pt) + theme_classic() +
    geom_histogram(aes(x=value, fill=cls), bins=100) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    facet_wrap(~cls, ncol=1, scales = 'free') +
    labs(title = paste(length(expr_mean_after), length(expr_mean), 
                       'remain under threshold', rm_low_expr_gene, sep='-')) +
    theme(legend.position = 'none')
  if (!is.null(diag_out)){
    ggsave(file.path(diag_out, 'expression_hist.png'), gp, device = 'png')
  }
  return(exps)
}


# brca.tcga = load_data_cbp(dataset_home = '/Users/jefft/Desktop/p53_project/datasets/BRCA-TCGA/brca_tcga_pan_can_atlas_2018')
# 
# brca.cptac = load_data_cbp(dataset_home = '/Users/jefft/Desktop/p53_project/datasets/BRCA-CPTAC/brca_cptac_2020/',
#                            case_complete_nm = 'cases_3way_complete.txt',
#                            exp_nm = 'data_mrna_seq_fpkm.txt', na.str = 'NA')
# 
# brca.metabric = load_data_cbp(dataset_home = '/Users/jefft/Desktop/p53_project/datasets/METABRIC/brca_metabric/',
#                               exp_nm = 'data_mrna_agilent_microarray.txt', na.str = 'NA')
# 
# integrated = integrate_dts(list(brca.tcga[[1]], brca.cptac[[1]]), z_score_norm = T)
# sample_factor = factor(c(rep(1, dim(brca.tcga[[1]]@colData)[1]), 
#                          rep(2, dim(brca.cptac[[1]]@colData)[1])))
# integrated$RNA = align_expression(integrated$RNA, sample_factor)
# 
# a = prcomp(integrated$RNA)
# rot = a$rotation
# df = data.frame('PC1'=rot[,'PC1'], 'PC2'=rot[,'PC2'])
# df['dataset'] = F
# df$dataset[which(grep('TCGA', rownames(df)))] = T
# ggplot(df, aes(x=PC1,y=PC2)) + geom_point(aes(color=dataset))

# MAF can be converted to mae, but not summarizedexperiment.
