library(logging)
library(MatrixEQTL)

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
  
  if (!is.null(meta_mut_fp)){
    meta_mut = read.table(meta_mut_fp, sep='\t', header = T)
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
  loginfo('Using VAF cutoff %f, at least %d / %d samples.', min_sample_per_snp, MIN_SAMPLE_PER_SNP,
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
  
  for (i in 1:length(lookup_vector)){
    nm = names(lookup_vector)[i]
    sp_idx = which(maf$HGVSc == nm)
    sp_name = maf[sp_idx, sample_col_id]
    for (sp in sp_name){
      snp_m[lookup_vector[nm], sp] = 1
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
      snp_m[id, maf_sub_meta_mut[[sample_col_name]]] = 1
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

get_single_eqtl_plot_dt = function(gene_name, snpid, mae, snp_m, lookup=NULL){
  # lookup: optional table for association between snp id and description
  # snp_m: snp matrix loaded to MatrixEQTL package
  rna = assay(mae[[1]][gene_name,,'RNA'])
  rna = t(rna)
  rna = data.frame(rna)
  rna['class'] = as.factor(snp_m[snpid,rownames(rna)])
  colnames(rna) = c('expression', 'class')
  if (!is.null(lookup)){
    var_name = rownames(lookup)[which(lookup$snpid==snpid)]
    return(list('plot_df'=rna, 'var_name'=var_name))
  }
  return(list('plot_df'=rna, 'var_name'=NULL))
}

get_var_info_from_maf = function(maf, snpid, lookup){
  var_name = rownames(lookup)[which(lookup$snpid==snpid)]
  maf = subset(maf, maf$HGVSc == var_name)
  return(unique(maf$HGVSp_Short))
}



