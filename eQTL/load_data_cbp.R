library(tidyverse)
library(maftools)
# Ref: https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html
library(MultiAssayExperiment)
# Ref: https://bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/MultiAssayExperiment.html
# setwd('/Users/jefft/Desktop/p53_project')
# 
# dataset_home = '/Users/jefft/Desktop/p53_project/datasets/BRCA-TCGA/brca_tcga_pan_can_atlas_2018'
# case_complete_nm = 'cases_complete.txt'
# exp_nm = 'data_mrna_seq_v2_rsem.txt'
# cna_nm = 'data_cna.txt'
# mut_nm = 'data_mutations.txt'
# gene_col_nm = 'Hugo_Symbol'
# sample_meta_nm = 'data_clinical_sample.txt'

clean_matrix = function(mtx, complete_cases_dot, gene_col_nm){
  mtx = mtx[!is.na(mtx[gene_col_nm]),]
  # remove dup by taking the record with the largest mean
  dup_genes = mtx[[gene_col_nm]][which(duplicated(mtx[[gene_col_nm]]))]
  remove_rows = c()
  for (i in dup_genes){
    sub = subset(mtx, mtx[[gene_col_nm]]==i)
    rms = rowMeans(sub[,complete_cases_dot])
    remove_rows = c(remove_rows, rownames(sub)[which(rms < max(rms))])
  }
  mtx = mtx[which(!rownames(mtx) %in% remove_rows),]
  remain_dup = which(duplicated(mtx[[gene_col_nm]]))
  if (length(remain_dup)>0){
    mtx = mtx[-remain_dup,]  # remove secondary dup, use first record.
  }
  rownames(mtx) = mtx[[gene_col_nm]]
  mtx = mtx[,complete_cases_dot]
  return(mtx)
}

load_data_cbp = function(dataset_home, case_complete_nm = 'cases_complete.txt',
                         exp_nm = 'data_mrna_seq_v2_rsem.txt',
                         cna_nm = 'data_cna.txt',
                         mut_nm = 'data_mutations.txt',
                         gene_col_nm = 'Hugo_Symbol',
                         sample_meta_nm = 'data_clinical_sample.txt',
                         case_list_dir_nm = 'case_lists',
                         na.str = ''){
  
  complete_cases = read.table(file.path(dataset_home, case_list_dir_nm, case_complete_nm), sep=':')
  complete_cases = trimws(strsplit(complete_cases[nrow(complete_cases),2], '\t')[[1]])
  complete_cases_dot = gsub('-', '.', complete_cases)
  exps = read.table(file.path(dataset_home, exp_nm), sep='\t', header = T, na.strings = na.str)
  cnas = read.table(file.path(dataset_home, cna_nm), sep='\t', header = T, na.strings = na.str)
  muts = read.maf(file.path(dataset_home, mut_nm))
  sample_meta = read.table(file.path(dataset_home, sample_meta_nm), sep = '\t', comment.char = '#', header = T)
  
  exps = clean_matrix(exps, complete_cases_dot = complete_cases_dot, gene_col_nm = gene_col_nm)
  cnas = clean_matrix(cnas, complete_cases_dot = complete_cases_dot, gene_col_nm = gene_col_nm)
  co_genes = as.character(ordered(intersect(rownames(exps), rownames(cnas))))
  exps = exps[co_genes,]
  cnas = cnas[co_genes,]
  muts = subsetMaf(muts, tsb=complete_cases, genes = co_genes)
  sample_meta = subset(sample_meta, sample_meta$SAMPLE_ID %in% complete_cases)
  rownames(sample_meta) = sample_meta$SAMPLE_ID
  sample_meta = sample_meta[,-which(colnames(sample_meta)=='SAMPLE_ID')]
  rownames(sample_meta) = gsub('-', '.', rownames(sample_meta))
  sample_meta = sample_meta[colnames(exps),]
  rownames(sample_meta) = gsub('\\.', '-', rownames(sample_meta))
  gc()
  
  # construct multi-assay-experiment
  exprdat = SummarizedExperiment(exps)
  colnames(exprdat) = complete_cases
  cnadat = SummarizedExperiment(cnas)
  colnames(cnadat) = complete_cases
  multiAssay = MultiAssayExperiment(
    list('RNA'=exprdat, 'CNA'=cnadat),
    sample_meta)
  
  a = ls()
  rm(list=a[which(a != 'multiAssay' & a != 'muts')])
  return(list(multiAssay, muts))
}

integrate_dts = function(datasets){
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
    for (dt in datasets){
      dt_df = data.frame(dt@ExperimentList[[cep]]@assays@data)[co_rnms,]
      dt_df = dt_df[,-which(colnames(dt_df) %in% c('group', 'group_name'))]
      merged_exp = cbind(merged_exp, dt_df)
    }
    
    # remove any row that contains NA (?)
    merged_exp = merged_exp[complete.cases(merged_exp),]
    
    merged[[cep]] = merged_exp
  }
  
  a = ls()
  rm(list=a[which(a != 'merged')])
  gc()
  return(merged)
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
# integrated = integrate_dts(list(brca.tcga[[1]], brca.cptac[[1]]))
# # MAF can be converted to mae, but not summarizedexperiment. 


#=========================== Try TCGABiolink?

#===========================
### Below cbioportal api implementation. Too slow
# library(cBioPortalData)
# cbio = cBioPortal()
# studies = getStudies(api = cbio)
# # get all genes
# gt = geneTable(cbio, pageSize = 999999) %>% 
#   filter(type %in% c('protein-coding', 'miRNA', 'snoRNA', 'ncRNA', 'snRNA'))
# 
# # Cached in '/Users/jefft/Library/Caches/org.R-project.R/R/cBioPortalData'
# sid = 'brca_tcga'  # must in studyId in the table above
# tcga_brca_clin = clinicalData(cbio, studyId = sid)
# mprofile = molecularProfiles(api = cbio, studyId = sid)
# samples = sampleLists(cbio, studyId = sid)
# sl_id = samples$sampleListId[which(samples$name=='Complete samples')]
# dt = getDataByGenes(cbio, studyId = sid, genes = gt$entrezGeneId,
#               by = 'entrezGeneId', molecularProfileIds = 'brca_tcga_rna_seq_v2_mrna',
#               sampleListId = sl_id)
# 
# # Inspect below to select wanted profile
# mprofile_ftd = mprofile[which(mprofile$molecularAlterationType %in% 
#                                 c('MRNA_EXPRESSION', 'MUTATION_EXTENDED')),]
# mprofile_ftd = mprofile[which(mprofile$name %in% c('Mutations', 'mRNA expression (RNA Seq V2 RSEM)')),]
# 
# exps = molecularData(api = cbio, molecularProfileIds = mprofile_ftd$molecularProfileId[1],
#                      entrezGeneIds = 1:1000, sampleIds = c())
# muts = molecularData(api = cbio, molecularProfileIds = mprofile_ftd$molecularProfileId[2])