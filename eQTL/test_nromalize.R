library(data.table)
library(tidyverse)
source('load_data_cbp.R')

# AIM: test if normalizing RNA sequencing data matches pre-computed data
exp_target_nm = 'data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt'
exp_target_nm = 'data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt'

dataset_home = '/Users/jefft/Desktop/p53_project/datasets/BRCA-TCGA/brca_tcga_pan_can_atlas_2018'
case_complete_nm = 'cases_complete.txt'
exp_nm = 'data_mrna_seq_v2_rsem.txt'
cna_nm = 'data_cna.txt'
mut_nm = 'data_mutations.txt'
gene_col_nm = 'Hugo_Symbol'
sample_meta_nm = 'data_clinical_sample.txt'
patient_meta_nm = 'data_clinical_patient.txt'
case_list_dir_nm = 'case_lists'
na.str = ''
has_log_ed = FALSE
diag_out = NULL
filter_protein_coding = TRUE

patient_id_col_nm = 'PATIENT_ID'
complete_cases = read.table(file.path(dataset_home, case_list_dir_nm, case_complete_nm), sep=':')
complete_cases = trimws(strsplit(complete_cases[nrow(complete_cases),2], '\t')[[1]])

exps = fread(file.path(dataset_home, exp_nm), sep='\t', header = T, na.strings = na.str)
exps = as.data.frame(exps)
loginfo(logger ='data.loader',
        'Expression matrix dim before cleaning: %d genes X %d samples.', dim(exps)[1], dim(exps)[2])
cnas = fread(file.path(dataset_home, cna_nm), sep='\t', header = T, na.strings = na.str)
cnas = as.data.frame(cnas)

complete_cases = intersect(colnames(exps), complete_cases)  # sometimes, complete cases may not account for expression matrix. 
complete_cases = intersect(complete_cases, colnames(cnas))

sample_meta = read.table(file.path(dataset_home, sample_meta_nm), sep = '\t', comment.char = '#', header = T, quote = '\"', fill=T)
patient_meta = read.table(file.path(dataset_home, patient_meta_nm), sep = '\t', comment.char = '#', header = T, quote = '\"', fill=T)
sample_meta = left_join(sample_meta, patient_meta, by = patient_id_col_nm)

exps = clean_matrix(exps, complete_cases_dot = complete_cases, gene_col_nm = gene_col_nm)
cnas = clean_matrix(cnas, complete_cases_dot = complete_cases, gene_col_nm = gene_col_nm)
co_genes = as.character(ordered(intersect(rownames(exps), rownames(cnas))))
exps = exps[co_genes,]
cnas = cnas[co_genes,]
exps = as.matrix(exps)

# normalize to diploid samples
norm_exp_to_diploid = function(cnas, exps){
  col_spl = colnames(cnas)
  mask = as.matrix(cnas)
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

# normalize to all samples
mus = rowMeans(exps)
stds = rowSds(exps)
for (i in 1:nrow(exps)){
  exps[i,] = (exps[i,] - mus[i]) / stds[i]
}

# compare with pre-computed table
exps = assay(dt[[1]][,,'RNA'])
gene_col_nm = 'Hugo_Symbol'
exps_pre = fread('/Users/jefft/Desktop/p53_project/datasets/CCLE/ccle_broad_2019/data_mrna_seq_rpkm_zscores_ref_diploid_samples.txt',
                 sep='\t', header = T, na.strings = 'NA')
exps_pre = fread('/Users/jefft/Desktop/p53_project/datasets/BRCA-TCGA/brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt',
                 sep='\t', header = T, na.strings = 'NA')
exps_pre = as.data.frame(exps_pre)
exps_pre = clean_matrix(exps_pre, complete_cases_dot = colnames(exps), gene_col_nm = gene_col_nm)
exps_pre = exps_pre[rownames(exps),]
exps_pre = as.matrix(exps_pre)

pheatmap::pheatmap(cbind(exps[1:200,1:200], exps_pre[1:200,1:200]), 
                   cluster_rows = F, cluster_cols = F, labels_row = '', labels_col = '')
