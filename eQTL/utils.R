library(logging)
library(MatrixEQTL)
library(maftools)
library(tidyverse)
library(MultiAssayExperiment)


get_eQTL_m = function(mae, genes='TP53', sample_col_name='Tumor_Sample_Barcode',
                      min_sample_per_snp=0.01, meta_mut_fp=NULL){
  # MAE: multiassay experiment object list, first object be gene expression and the second be MAF
  # Return Gene expression and SNP slicedData for eQTL mapping
  expression = assay(mae[[1]][,,'RNA', drop=FALSE])
  
  rs = get_SNP_m_from_maf(mae[[2]]@data, genes=genes, samples=colnames(expression),
                          sample_col_name = sample_col_name,
                          min_sample_per_snp = min_sample_per_snp,
                          meta_mut_fp = meta_mut_fp)
  gene = SlicedData$new()
  snp = SlicedData$new()
  snp = snp$CreateFromMatrix(rs[[2]])
  gene = gene$CreateFromMatrix(as.matrix(expression))
  return(list('gene'=gene, 'snp'=snp, 'snp_lookup'=rs[[1]]))
}


get_SNP_m_from_maf = function(maf, sample_col_name,
                              genes=NULL, samples=NULL,
                              min_sample_per_snp=5,
                              meta_mut_fp=NULL){
  # Get SNP matrix from MAF, also return SNP association table
  # genes: select specific gene list to use (Hugo symbol), if NULL, return all.
  # samples: tumour samples involved, if null will use all samples in maf.
  
  # TODO: classify snp as 0, 1, or 2.
  
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
  
  maf = as.data.frame(maf)
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
        snp_m[lookup_vector[nm], sp_name[sp]] = maf[sp_idx[sp], 'VAF']
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
      maf_sub_meta_mut = subset(maf, maf$Protein_position %in% aa_pos)
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
  maf = subset(mae[[2]]@data, mae[[2]]@data[[gene_col_name]] == 'TP53')
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


get_intersection_eqtl = function(eqtl_df, groups=NULL, group_col='snps'){
  # Get eqtl that present in all groups specified
  if (is.null(groups)){
    groups = unique(eqtl_df[[group_col]])
  }
  co_eqtl = list()
  for (snp in groups){
    sub_eqtl = subset(eqtl_df, eqtl_df[[group_col]] == snp)
    co_eqtl[[snp]] = sub_eqtl$gene
  }
  co_eqtl = Reduce(intersect, co_eqtl)
  eqtl_df = subset(eqtl_df, eqtl_df[[group_col]] %in% groups & eqtl_df$gene %in% co_eqtl)
  eqtl_df = eqtl_df[,which(!colnames(eqtl_df) %in% c('statistic', 'FDR'))]
  eqtl_df = eqtl_df[order(eqtl_df$gene, eqtl_df$beta, decreasing = c(F,T)),]
  return(eqtl_df)
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


annotate_sample_mut = function(maf, sample_col_name = 'Tumor_Sample_Barcode', gene='TP53', gene_col_name = 'Hugo_Symbol'){
  maf = subset(maf, maf[[gene_col_name]] %in% gene)
  maf[[sample_col_name]] = as.character(maf[[sample_col_name]])
  mis_sample = unique(maf[[sample_col_name]][which(maf$Consequence=='missense_variant')])
  non_sample = unique(maf[[sample_col_name]][which(maf$Consequence=='stop_gained')])
  fsv_sample = unique(maf[[sample_col_name]][which(maf$Consequence=='frameshift_variant')])
  # Sample may have more than one kinds of mutations
  pool = c(mis_sample, non_sample, fsv_sample)
  others = pool[which(duplicated(pool))]
  big_pool = unique(maf[[sample_col_name]])
  others = c(others,  big_pool[which(!big_pool %in% pool)])
  return(list(
    'missense' = mis_sample,
    'nonsense' = non_sample,
    'frameshift' = fsv_sample,
    'others' = others
  ))
}


pcs_cfg = function(config){
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

