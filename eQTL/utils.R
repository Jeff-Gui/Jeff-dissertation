library(logging)
library(MatrixEQTL)
library(maftools)
library(tidyverse)
library(stringr)
library(MultiAssayExperiment)


get_eQTL_m = function(mae, genes='TP53', sample_col_name='Tumor_Sample_Barcode',
                      min_sample_per_snp=0.01, meta_mut_fp=NULL, mode='mut-VS-other',
                      use_p53_exp = FALSE, SIFT = NULL){
  # MAE: multiassay experiment object list, first object be gene expression and the second be MAF
  # Return Gene expression and SNP slicedData for eQTL mapping
  # mode: see get_SNP_m_from_maf
  expression = assay(mae[[1]][,,'RNA', drop=FALSE])
  if (use_p53_exp){
    p53_exp = as.numeric(expression['TP53',])
    p53_exp = (p53_exp - min(p53_exp)) / (max(p53_exp) - min(p53_exp))
    names(p53_exp) = colnames(expression)
  } else {
    p53_exp = NULL
  }
  
  rs = get_SNP_m_from_maf(mae[[2]]@data, genes=genes, samples=colnames(expression),
                          sample_col_name = sample_col_name,
                          min_sample_per_snp = min_sample_per_snp,
                          meta_mut_fp = meta_mut_fp, mode=mode,
                          p53_expr = p53_exp, SIFT = SIFT)
  gene = SlicedData$new()
  snp = SlicedData$new()
  snp = snp$CreateFromMatrix(rs[[2]])
  expression = as.matrix(expression)[,colnames(rs[[2]])]
  gene = gene$CreateFromMatrix(expression)
  loginfo(logger = 'eQTL_gen_matrix', 'eQTL mapping sample size: %d.', ncol(rs[[2]]))
  return(list('gene'=gene, 'snp'=snp, 'snp_lookup'=rs[[1]]))
}


get_null_samples = function(maf, sample_col_name='Tumor_Sample_Barcode', 
                            gene='TP53'){
  # get sample ID that has stop mutation for certain gene
  maf = subset(maf, maf$Hugo_Symbol==gene)
  rt = maf[[sample_col_name]][which(maf$Variant_Classification == 'Nonsense_Mutation')]
  return(as.character(rt))
}


get_SNP_m_from_maf = function(maf, sample_col_name,
                              genes=NULL, samples=NULL,
                              min_sample_per_snp=5,
                              meta_mut_fp=NULL, mode='mut-VS-other',
                              p53_expr=NULL, SIFT = NULL){
  # Get SNP matrix from MAF, also return SNP association table
  # Only consider missense mutation
  # genes: select specific gene list to use (Hugo symbol), if NULL, return all.
  # samples: tumour samples involved, if null will use all samples in maf.
  # mode: mut-VS-other: contrast mut+ to mut-; mut-VS-null: contrast mut+ to nonsense mutation.
  # SIFT: sift prediction score, a vector named by HGVSp_Short
  #   If consider SIFT: will consider Damaging mutant only
  
  # TODO: classify snp as 0, 1, or 2.
  
  meta_mut = load_meta_mut(meta_mut_fp)
  
  maf = as.data.frame(maf)
  null_spl = get_null_samples(maf, gene = 'TP53', sample_col_name = sample_col_name)
  maf = subset(maf, maf$Variant_Classification == 'Missense_Mutation')
  if (!is.null(SIFT)){
    maf = maf[which(SIFT[maf$HGVSp_Short] == 'Damaging'),]
  }
  sample_col_id = which(colnames(maf) == sample_col_name)
  maf[,sample_col_id] = as.character(maf[,sample_col_id])
  
  if (!is.null(genes)){
    maf = subset(maf, maf$Hugo_Symbol %in% genes)
  }
  if (!is.null(samples)){
    maf = subset(maf, maf[[sample_col_name]] %in% samples)
  } else {
    samples = unique(maf[[sample_col_name]])
  }
  # TODO: enfore filtering of VAF
  MIN_SAMPLE_PER_SNP = ceiling(min_sample_per_snp * length(samples))
  loginfo('Using VAF cutoff %f, at least %d / %d samples', min_sample_per_snp, MIN_SAMPLE_PER_SNP,
          length(samples), logger = 'eQTL_QC')
  
  hgvsc = maf[,c('HGVSc','Chromosome','Start_Position')]
  hgvsc = hgvsc[!duplicated(hgvsc),]
  lookup = hgvsc
  lookup['snpid'] = paste('snp', 1:nrow(hgvsc), sep='')
  lookup['Chromosome'] = paste('chr', lookup$Chromosome, sep='')
  colnames(lookup) = c('HGVSc', 'chr', 'pos', 'snpid')
  lookup_vector = lookup$snpid
  names(lookup_vector) = lookup$HGVSc
  
  snp_m = matrix(0, nrow = length(lookup_vector), ncol = length(samples))
  colnames(snp_m) = samples
  rownames(snp_m) = lookup_vector
  
  has_vaf = 'VAF' %in% colnames(maf)
  for (i in 1:length(lookup_vector)){
    nm = names(lookup_vector)[i]
    sp_idx = which(maf$HGVSc == nm)
    sp_name = maf[sp_idx, sample_col_id]
    # if (i %% 1000 == 0){
    #  print(i)
    #}
    for (sp in 1:length(sp_name)){
      if (!has_vaf){
        # If no VAF information, assume pure
        snp_m[lookup_vector[nm], sp_name[sp]] = 1
      } else {
        if (!is.null(p53_expr)){
          snp_m[lookup_vector[nm], sp_name[sp]] = maf[sp_idx[sp], 'VAF'] * p53_expr[sp_idx[sp]]
        } else {
          snp_m[lookup_vector[nm], sp_name[sp]] = maf[sp_idx[sp], 'VAF']
        }
      }
    }
  }
  rownames(lookup) = lookup$HGVSc
  lookup = lookup[,-1]
  lookup = lookup[,c(3,1,2)]
  
  # Resolve meta mutations
  mean_pos_snp = ceiling(mean(lookup$pos))
  if (!is.null(meta_mut_fp)){
    for (id in unique(meta_mut$meta_mut_id)){
      aa_pos = meta_mut$aa_pos[which(meta_mut$meta_mut_id==id)]
      if (length(grep('^p\\.', aa_pos) > 0)){
        maf_sub_meta_mut = subset(maf, maf$HGVSp_Short %in% aa_pos)
      } else {
        maf_sub_meta_mut = subset(maf, maf$Protein_position %in% aa_pos)
      }
      snp_m = rbind(snp_m, rep(0,ncol(snp_m)))
      rownames(snp_m)[nrow(snp_m)] = id
      if (!has_vaf){
        snp_m[id, maf_sub_meta_mut[[sample_col_name]]] = 1
      } else {
        # If VAF information available, average VAF of meta-mutants records for each sample
        for (barcode in unique(maf_sub_meta_mut[[sample_col_name]])){
          maf_sub_meta_mut_barcode = subset(maf_sub_meta_mut, maf_sub_meta_mut[[sample_col_name]] == barcode)
          snp_m[id, barcode] = mean(maf_sub_meta_mut_barcode[['VAF']])
        }
      }
      lookup = rbind(lookup, c(id, 'chr17', mean_pos_snp))
      rownames(lookup)[nrow(lookup)] = id
    }
  }
  lookup$pos = as.numeric(lookup$pos)
  
  # Filter low-frequency SNPs
  rm_snps = names(which(rowSums(snp_m) < MIN_SAMPLE_PER_SNP))
  loginfo('Filtering out %d / %d mutations. Remaining %d.', length(rm_snps), nrow(snp_m), 
          nrow(snp_m) - length(rm_snps),
          logger = 'eQTL_QC')
  
  if (length(rm_snps) > 0){
    snp_m = snp_m[which(!rownames(snp_m) %in% rm_snps),]
    lookup = subset(lookup, !lookup$snpid %in% rm_snps)
  }
  if (is.null(dim(snp_m))){
    snp_m = t(as.matrix(snp_m))
    rownames(snp_m) = lookup$snpid[1]
  }
  
  if (mode == 'mut-VS-null'){
    rm = which((colSums(snp_m) == 0) & (!colnames(snp_m) %in% null_spl))
    print(length(rm))
    if (length(rm)>0){
      snp_m = snp_m[,-rm]
    }
  }
  
  return(list('lookup'=lookup, 'SNP_m'=snp_m))
}


merge_cfg = function(cfg_a, cfg_b){
  # merge config b into config a, b will overlap with a
  # does not allow b having extra config than a
  for (i in names(cfg_a)){
    if (i %in% names(cfg_b)){
      if (class(cfg_a[[i]]) == 'list'){
        cfg_a[[i]] = merge_cfg(cfg_a[[i]], cfg_b[[i]])
      } else {
        cfg_a[[i]] = cfg_b[[i]]
      }
    }
  }
  return(cfg_a)
}


get_single_eqtl_plot_dt = function(gene_name, snpid, mae, snp_m,
                                   sample_col_name='Tumor_Sample_Barcode',
                                   gene_col_name='Hugo_Symbol'){
  # lookup: optional table for association between snp id and description
  # snp_m: snp matrix loaded to MatrixEQTL package
  rna = assay(mae[[1]][gene_name,,'RNA'])
  if (nrow(rna)==0){
    print('Gene not in expression matrix')
    return(NULL)
  }
  rna = t(rna)
  rna = data.frame(rna)
  rna['mut_state'] = snp_m[snpid,rownames(rna)]
  colnames(rna) = c('expression', 'mut_state')
  
  # annotate mut class
  rna[which(rna$mut_state==0), 'mut_cls'] = 'wildtype'
  mut_ann = annotate_sample_mut(mae[[2]]@data, sample_col_name, 'TP53', gene_col_name)
  # maf = subset(mae[[2]]@data, mae[[2]]@data[[gene_col_name]] == 'TP53')
  for (i in names(mut_ann)){
    rna[mut_ann[[i]], 'mut_cls'] = i
  }
  rna[which(rna$mut_state!=0), 'mut_cls'] = 'eQTL mutant'
  return(rna)
}


get_var_info_from_maf = function(maf, snpid, lookup){
  var_name = rownames(lookup)[which(lookup$snpid==snpid)]
  maf = subset(maf, maf$HGVSc == var_name)
  return(unique(maf$HGVSp_Short))
}


get_intersection_eqtl = function(eqtl_df, groups=NULL, group_col='snps', 
                                 doPlot=FALSE, check_beta=NULL){
  # Get eqtl that present in all groups specified
  # check_beta: either NULL, pos or neg, df must has beta column
  if (is.null(groups)){
    groups = unique(eqtl_df[[group_col]])
  }
  co_eqtl = list()
  for (snp in groups){
    sub_eqtl = subset(eqtl_df, eqtl_df[[group_col]] == snp)
    if (!is.null(check_beta)){
      if (check_beta == 'pos'){
        sub_eqtl = subset(sub_eqtl, sub_eqtl$beta > 0)
      } else {
        sub_eqtl = subset(sub_eqtl, sub_eqtl$beta < 0)
      }
    }
    co_eqtl[[snp]] = sub_eqtl$gene
  }
  co_eqtl_rd = Reduce(intersect, co_eqtl)
  eqtl_df = subset(eqtl_df, eqtl_df[[group_col]] %in% groups & eqtl_df$gene %in% co_eqtl_rd)
  eqtl_df = eqtl_df[,which(!colnames(eqtl_df) %in% c('statistic', 'FDR'))]
  eqtl_df = eqtl_df[order(eqtl_df$gene, eqtl_df$beta, decreasing = c(F,T)),]
  if (doPlot){
    plt = ggvenn::ggvenn(co_eqtl, stroke_color = 'white')
    return(list(eqtl_df, plt))
  } else {
    return(eqtl_df)
  }
}


genotype_pca = function(maf, out_dic=NULL, samples=NULL, file_out=FALSE,
                        gene_col_name='Hugo_Symbol', sample_col_name='Tumor_Sample_Barcode'){
  library(FactoMineR)
  # maf: maf object
  # out_dic: output directory
  maf = as.data.frame(maf)
  maf[[sample_col_name]] = as.character(maf[[sample_col_name]])
  gene_pool = unique(maf[[gene_col_name]])
  if (is.null(samples)){
    samples = unique(maf[[sample_col_name]])
  }
  gene_m = matrix(0, nrow=length(gene_pool), ncol=length(samples))
  rownames(gene_m) = gene_pool
  colnames(gene_m) = samples
  for (i in gene_pool){
    gene_m[i, maf[which(maf[[gene_col_name]] == i), sample_col_name]] = 1
  }
  if (!is.null(samples)){
    rest = samples[which(!samples %in% maf[[sample_col_name]])]
    if (length(rest) > 0){
      gene_m_aug = matrix(0, nrow=nrow(gene_m), ncol=length(rest))
      rownames(gene_m_aug) = rownames(gene_m)
      colnames(gene_m_aug) = rest
      gene_m = cbind(gene_m, gene_m_aug)
    }
  }
  
  # PCA
  library(factoextra)
  pca.res = PCA(t(gene_m), ncp=20, graph = F)
  ind = get_pca_ind(pca.res)
  scre_plot = fviz_eig(pca.res, addlabels = T)
  ggsave(file.path(out_dic, 'genotype_PCA.png'), scre_plot, device = 'png')
  #df = data.frame('PC1'=ind$coord[,1], 'PC3'=ind$coord[,20])
  #ggplot(df) +
  #  geom_point(aes(x=PC1, y=PC3))
  pca_coord = ind$coord
  rownames(pca_coord) = gsub('\\.', '-', rownames(pca_coord))
  if (file_out){
    write.table(pca_coord, file.path(out_dic, 'genotype_PC20.txt'), sep='\t', row.names = T, quote = F) 
  }
  return(pca_coord)
  
  # Simulate to test if do SNP PCA possible
  # dt = matrix(0, nrow=1000,ncol=10000)
  # mut_row = sample(1:1000,50000, replace=TRUE)
  # mut_col = sample(1:10000, 50000, replace=T)
  # for (i in 1:50000){
  #   dt[mut_row[i], mut_col[i]] = 1
  # }
  # pca.res = PCA(dt, ncp=20, graph = F)
  # fviz_eig(pca.res)  # if use SNP to do PCA, the explained varibility will be as low as 0.1%
}


annotate_sample_mut = function(maf, sample_col_name = 'Tumor_Sample_Barcode', gene='TP53', gene_col_name = 'Hugo_Symbol',
                               fp_SIFT = '/Users/jefft/Desktop/p53_project/datasets/UMD_TP53/SIFT_UMD.RData'){
  load(fp_SIFT)
  maf = subset(maf, maf[[gene_col_name]] %in% gene)
  maf[[sample_col_name]] = as.character(maf[[sample_col_name]])
  non_sample = unique(maf[[sample_col_name]][which(maf$Consequence=='stop_gained')])
  fsv_sample = unique(maf[[sample_col_name]][which(maf$Consequence=='frameshift_variant')])
  mis_sample = unique(maf[[sample_col_name]][which(maf$Consequence=='missense_variant' & 
                                                     sift_pred[maf$HGVSp_Short] == 'Damaging')])
  # Sample may have more than one kinds of mutations
  pool = c(mis_sample, non_sample, fsv_sample)
  others = pool[which(duplicated(pool))]
  big_pool = unique(maf[[sample_col_name]])
  dup_spl = maf[[sample_col_name]][which(duplicated(maf[[sample_col_name]]))]
  others = unique(c(others,  big_pool[which(!big_pool %in% pool)], dup_spl))
  return(list(
    'missense' = setdiff(mis_sample, others),
    'nonsense' = setdiff(non_sample, others),
    'frameshift' = setdiff(fsv_sample, others),
    'others' = others
  ))
}


pcs_cfg = function(config){
  config$output = file.path(config$experiment_home, config$output)
  dt_cfg = config$dataset
  eqtl_cfg = config$eQTL
  preprocess_cfg = config$preprocess
  
  if (is.null(dt_cfg$na.str)){
    dt_cfg$na.str = ''
  }
  if (is.null(eqtl_cfg$output_file_name_tra)){
    eqtl_cfg$output_file_name_tra = tempfile()
  }
  if (is.null(eqtl_cfg$output_file_name_cis)){
    eqtl_cfg$output_file_name_cis = tempfile()
  }
  if (!dir.exists(config$output)){
    dir.create(config$output)
  }
  if (eqtl_cfg$meta_mut[1] == '/'){
    eqtl_cfg$meta_mut = NULL
  }
  if (!is.null(preprocess_cfg$norm_gene)){
    preprocess_cfg$norm_gene = strsplit(preprocess_cfg$norm_gene, split = ',')[[1]]
  }
  
  config$dataset = dt_cfg
  config$eQTL = eqtl_cfg
  config$preprocess = preprocess_cfg
  return(config)
}


get_binary_SNP_m_from_maf = function(maf, snp_list, snp_col_nm=NULL,
                                     samples = NULL,
                                     gene='TP53', 
                                     sample_col_name='Tumor_Sample_Barcode',
                                     gene_col_name='Hugo_Symbol', mode='amino_acid',
                                     protein_change_col = 'Protein_Change',
                                     code_with_VAF = FALSE,
                                     fp_SIFT = '/Users/jefft/Desktop/p53_project/datasets/UMD_TP53/SIFT_UMD.RData'){
  maf = subset(maf, maf[[gene_col_name]] == gene)
  b4_length = nrow(maf)
  load(fp_SIFT)
  maf = maf[which(sift_pred[maf$HGVSp_Short] == 'Damaging'),]
  print(paste('SIFT filter:', nrow(maf),'out of', b4_length, sep=' '))
  # mode: if amino_acid, then should be specific mutation.
  # if position, then should be location only.
  if (mode=='position'){
    ext = sapply(maf[[protein_change_col]], function(x){strsplit(x, split = '\\.')[[1]][2]})
    maf[['position']] = as.numeric(str_extract(ext,"(?<=[A-Z])[0-9]+(?=[A-Z])"))
    snp_col_nm = 'position'
  } else {
    if (mode == 'position-tcga'){
      maf[['position']] = as.numeric(maf$Protein_position)
      snp_col_nm = 'position'
    }
  }
  
  snp_list_qc = list()
  for (i in 1:length(snp_list)){
    if (length(intersect(unlist(snp_list[[i]]), maf[[snp_col_nm]]))>=0){
      snp_list_qc[[as.character(i)]] = snp_list[[i]]
    }
  }
  if (length(snp_list_qc)<1){
    return(matrix())
  }
  if (is.null(samples)){
    samples = unique(maf[[sample_col_name]])
  }
  m = matrix(0, nrow=length(samples), ncol = length(snp_list_qc))
  if (!code_with_VAF){
    for (i in 1:nrow(m)){
      sub_maf = subset(maf, maf[[sample_col_name]] == samples[i])
      for (j in 1:ncol(m)){
        if (length(intersect(snp_list_qc[[j]], sub_maf[[snp_col_nm]]))>0){
          m[i,j] = 1
        }
      }
    }
  } else {
    for (i in 1:nrow(m)){
      sub_maf = subset(maf, maf[[sample_col_name]] == samples[i])
      for (j in 1:ncol(m)){
        its = intersect(snp_list_qc[[j]], sub_maf[[snp_col_nm]])
        if (length(its)>0){
          m[i,j] = mean(as.numeric(sub_maf[which(sub_maf[[snp_col_nm]] %in% its), 'VAF']))
        }
      }
    }
  }
  rownames(m) = samples
  colnames(m) = as.vector(t(lapply(snp_list_qc, function(x){
    return(paste(x, collapse = '_'))
  })))
  return(m)
}


load_meta_mut = function(meta_mut_fp){
  meta_mut = data.frame()
  if (!is.null(meta_mut_fp)){
    for (fp in meta_mut_fp){
      meta_mut_sub = read.table(fp, sep='\t', header = T) 
      meta_mut = rbind(meta_mut, meta_mut_sub)
    }
    colnames(meta_mut) = c('meta_mut_id', 'aa_pos')
  }
  if (sum(duplicated(meta_mut))>0){
    meta_mut = meta_mut[-duplicated(meta_mut),] 
  }
  return(meta_mut)
}


load_eQTL_output = function(dir_home, exclude=NULL, mode='fdr', beta=NULL, as_df=FALSE){
  beta_val = beta
  # beta: filter absolute effect size
  if (mode == 'fdr'){
    nm = 'trans_eqtl_fdr005.txt'
  } else {
    nm = 'trans_eqtl.txt'
  }
  output_dirs = list.files(dir_home)
  if (!is.null(exclude)){
    output_dirs = output_dirs[which(!output_dirs %in% exclude)]
  }
  coll = list()
  for (i in output_dirs){
    fp = file.path(dir_home, i, nm)
    if (file.exists(fp)){
      trans_eqtls = read.table(fp, sep='\t', header = T)
      if (nrow(trans_eqtls)==0){
        coll[[i]] = NA
      } else {
        trans_eqtls['experiment'] = i
        if (!is.null(beta_val)){
          trans_eqtls = subset(trans_eqtls, abs(trans_eqtls$beta) > beta_val)
        }
        if (nrow(trans_eqtls)==0){
          coll[[i]] = NA
        } else {
          coll[[i]] = trans_eqtls
        }
      }
    } else {
      coll[[i]] = NA
    }
  }
  if (!as_df){
    return(coll)
  } else {
    df_coll = data.frame()
    for (i in names(coll)){
      nm = toupper(strsplit(i, split='_')[[1]][2])
      df = coll[[i]]
      if (class(df)!='logical'){
        df = subset(df, abs(df$beta) > beta_cutoff)
        df_coll = rbind(df_coll, df)
      }
    }
    df_coll[['cancer']] = sapply(df_coll$experiment, function(x){
      return(strsplit(x, split='_')[[1]][2])
    })
    df_coll$cancer = toupper(df_coll$cancer)
    return(df_coll)
  }
}


load_clean_data = function(fp, ann_p53=TRUE, ann_bin_mut_list=NULL, mode=NULL, use_VAF_tcga=FALSE,
                meta_mut_home='/Users/jefft/Desktop/p53_project/datasets/meta_muts/'){
  # mode: tcga or ccle, if not specified, analyze from fp
  # load clean data while specifying the mutation state
  if (is.null(mode)){
    if (length(grep('tcga',fp, ignore.case = T))>0){
      mode = 'tcga'
    } else {
      mode = 'ccle'
    }
  }
  if (!mode %in% c('tcga', 'ccle')){
    print('Wrong mode!')
    return(NULL)
  }
  code_with_VAF = FALSE
  snp_col_nm = 'HGVSp_Short'
  if (mode == 'tcga'){
    protein_change_col = 'HGVSp_Short'
    if (use_VAF_tcga){
      code_with_VAF = TRUE
    }
  } else {
    protein_change_col = 'Protein_Change'
  }
  
  meta_mut = load_meta_mut(file.path(meta_mut_home, list.files(meta_mut_home)))
  load(fp)
  if (ann_p53){
    p53_ann = annotate_sample_mut(dt[[2]]@data)
    dt[[1]]@colData[['p53_state']] = 'Wildtype'
    for (i in names(p53_ann)){
      dt[[1]]@colData[p53_ann[[i]], 'p53_state'] = i
    }
  }
  sps = rownames(dt[[1]]@colData)
  if (!is.null(ann_bin_mut_list)){
    for (mut in ann_bin_mut_list){
      if (length(grep('p\\.', mut))>0){
        b_m = get_binary_SNP_m_from_maf(dt[[2]]@data, snp_list = list(mut),
                                        protein_change_col = protein_change_col,
                                        snp_col_nm = snp_col_nm,
                                        samples = sps, mode = 'amino_acid', code_with_VAF = code_with_VAF)
      } else {
        snp_list = unique(unlist(meta_mut$aa_pos[which(meta_mut$meta_mut_id==mut)]))
        if (length(snp_list)==0){
          print(paste(mut, 'not found in meta mutant directory.'))
          next
        }
        if (length(grep('^p\\.', snp_list))>0){
          b_m = get_binary_SNP_m_from_maf(dt[[2]]@data, snp_list = list(snp_list),
                                          protein_change_col = protein_change_col,
                                          snp_col_nm = snp_col_nm,
                                          samples = sps, mode = 'amino_acid', code_with_VAF = code_with_VAF)
        } else {
          b_m = get_binary_SNP_m_from_maf(dt[[2]]@data, snp_list = list(snp_list),
                                          protein_change_col = protein_change_col,
                                          snp_col_nm = snp_col_nm,
                                          samples = sps, mode = 'position', code_with_VAF = code_with_VAF)
        }
      }
      dt[[1]]@colData[[paste('has_', mut, sep='')]] = 0
      dt[[1]]@colData[rownames(b_m)[b_m==1],paste('has_', mut, sep='')] = 1
    }
  }
  return(dt)
}


ext_gene_GO = function(gene_lists, do_intersect=FALSE){
  gls = sapply(gene_lists, function(x){return(strsplit(x, split='/')[[1]])})
  if (length(gene_lists)==1){
    return(gls)
  }
  if (do_intersect){
    names(gls) = NULL
    return(Reduce(intersect, gls))
  } else {
    return(unique(unlist(gls)))
  }
}


annotate_mut_group = function(dt, check_muts, as.mtx=TRUE){
  # Wrapper of get binrary mutation result
  # get bulk mutation state at the same time.
  if (!'p53_state' %in% colnames(dt[[1]]@colData)){
    p53_ann = annotate_sample_mut(dt[[2]]@data)
    dt[[1]]@colData[['p53_state']] = 'Wildtype'
    for (i in names(p53_ann)){
      dt[[1]]@colData[p53_ann[[i]], 'p53_state'] = i
    }
  }
  
  # load meta mutant
  meta_mut_home = '/Users/jefft/Desktop/p53_project/datasets/meta_muts'
  meta_mut_fp = list.files(meta_mut_home)
  meta_mut_fp = sapply(meta_mut_fp, 
                       function(x){return(file.path(meta_mut_home, x))})
  meta_mut = load_meta_mut(meta_mut_fp)
  
  idx_ftd = rownames(dt[[1]]@colData)
  others = intersect(idx_ftd, rownames(dt[[1]]@colData)[which(dt[[1]]@colData$p53_state == 'others')])
  rt_list = list('others'=others)
  for (check_mut in check_muts){
    flag = FALSE
    if (length(grep('p\\.', check_mut))==0){
      mut_aa = unique(meta_mut$aa_pos[meta_mut$meta_mut_id == check_mut])
      if (length(grep('p\\.', mut_aa))==0){
        b_m = get_binary_SNP_m_from_maf(dt[[2]]@data, 
                                        snp_list = list(mut_aa),
                                        samples = idx_ftd, 
                                        mode = 'position')
      } else {
        b_m = get_binary_SNP_m_from_maf(dt[[2]]@data, 
                                        snp_list = list(mut_aa),
                                        snp_col_nm = 'Protein_Change',
                                        samples = idx_ftd, 
                                        mode = 'amino_acid')
      }
    } else {
      b_m = get_binary_SNP_m_from_maf(dt[[2]]@data, 
                                      snp_list = list(check_mut),
                                      snp_col_nm = 'Protein_Change',
                                      samples = idx_ftd, 
                                      mode = 'amino_acid')
    }
    mutant = rownames(b_m)[which(b_m==1)]
    mutant = mutant[!mutant %in% others]
    rt_list[[check_mut]] = mutant
  }
  rt_list[['wildtype']] = intersect(idx_ftd, rownames(dt[[1]]@colData)[which(dt[[1]]@colData$p53_state == 'Wildtype')])
  rt_list[['nonsense']] = intersect(idx_ftd, rownames(dt[[1]]@colData)[which(dt[[1]]@colData$p53_state == 'nonsense')])
  if (!as.mtx){
    return(rt_list)
  } else {
    mtx = matrix(0, nrow=length(unique(unlist(rt_list))), ncol=length(names(rt_list)))
    rownames(mtx) = unique(unlist(rt_list))
    colnames(mtx) = names(rt_list)
    mtx = as.data.frame(mtx)
    for (i in names(rt_list)){
      mtx[rt_list[[i]], i] = 1
    }
    return(mtx)
  }
}


gen_bed_GOI = function(GOI, chr_prefix=FALSE, output=NULL, expand=c(-1000,1000),
                       genepos='/Users/jefft/Desktop/p53_project/scripts/eQTL/hg38_gene_table_autosome.tsv'){
  if (length(expand)==1){
    expand = c(-expand[1], expand[1])
  }
  # generate bed file for gene of interest
  df = read.table(genepos, header = T)
  GOI = intersect(GOI, df$geneid)
  df = df[df$geneid %in% GOI,]
  rownames(df) = df$geneid
  df = df[GOI,]
  df$chr = gsub('chr','',df$chr)
  if (chr_prefix){
    df$chr = paste('chr', df$chr, sep='')
  }
  df = df[,c('chr','left','right')]
  df$left = df$left + expand[1]
  df$right = df$right + expand[2]
  for (i in 1:nrow(df)){
    if (df$left[i]<0){
      df$left[i] = 0
    }
  }
  rt_df = df
  rownames(df) = NULL
  if (!is.null(output)){
    write.table(df, file.path(output, 'GOI.bed'), sep='\t', quote = F, row.names = F, col.names = F)
    write.table(paste(rownames(rt_df), '   ', rt_df$chr, ':', rt_df$left,'-', rt_df$right,sep=''), 
                file.path(output, 'GOI_geneloc.txt'), quote = F, row.names = F, col.names = F)
  }
  return(rt_df)
}

