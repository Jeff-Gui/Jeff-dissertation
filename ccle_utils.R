library(tidyverse)
library(data.table)
library(ggpubr)

load_rnai = function(rnai_fp = '/Users/jefft/Desktop/p53_project/datasets/CCLE/extra_raw/RNAi_D2_combined_gene_dep_scores.csv',
                     hg19_ann = '/Users/jefft/Desktop/p53_project/scripts/eQTL/hg19_gene_table_autosome.tsv',
                     protein_coding_only = TRUE){
  tb = fread(rnai_fp, sep=',', header = T, na.strings = 'NA', quote = '\"')
  tb = as.data.frame(tb)
  rownames(tb) = sapply(tb$V1, function(x){
    return(strsplit(x, split = ' \\(')[[1]][1])
  })
  tb = tb[,-1]
  if (protein_coding_only){
    # CCLE is aligned to hg19
    genepos = read.table(hg19_ann, sep='\t', header = T)
    tb = tb[which(rownames(tb) %in% genepos$geneid),]
  }
  tb = tb[order(rownames(tb)),]
  return(tb)
}


load_protein_z_score = function(protein_fp='/Users/jefft/Desktop/p53_project/datasets/CCLE/ccle_broad_2019/data_protein_quantification_zscores.txt'){
  tb = fread(protein_fp, sep='\t', header = T, na.strings = 'NA', quote = '\"')
  tb = as.data.frame(tb)
  rownames(tb) = tb[,1]
  tb = tb[,-1]
  # collapse duplicated proteins to per-gene by taking mean
  gene_nms = sapply(rownames(tb), function(x){
    return(strsplit(x, split = '\\|')[[1]][1])
  })
  tb['gene_names'] = gene_nms
  dup_idx = which(duplicated(gene_nms))
  dup_nm = gene_nms[dup_idx]
  tb_non_dup = tb[-which(tb$gene_name %in% dup_nm),]
  tb_dup = tb[which(tb$gene_name %in% dup_nm),]
  tb_dedup = gather(tb_dup, key='cl', value='value', 1:(ncol(tb_dup)-1)) %>%
    group_by(gene_names, cl) %>% summarise(mean_value = mean(value, na.rm=T)) %>%
    spread(cl, mean_value)
  tb_dedup = as.data.frame(tb_dedup)
  rownames(tb_dedup) = tb_dedup$gene_names
  tb_dedup = tb_dedup[,-1]
  tb_dedup[is.nan(as.matrix(tb_dedup))] = NA
  tb_dedup['gene_names'] = rownames(tb_dedup)
  tb_dedup = tb_dedup[,colnames(tb)]
  
  # re-calculate the z-score before merging?
  tb_non_dup = rbind(tb_non_dup, tb_dedup)
  tb_non_dup = tb_non_dup[order(rownames(tb_non_dup)),]
  rownames(tb_non_dup) = tb_non_dup$gene_names
  tb_non_dup = tb_non_dup[,-which(colnames(tb_non_dup)=='gene_names')]
  return(tb_non_dup)
}


get_genes_plt = function(genes, ccle, tcga, mutation_groups, 
                         primary_site, rnai=NULL, comparison=NULL){
  source('/Users/jefft/Desktop/Manuscript/set_theme.R')
  # comparison: must be list('rna'=list(), 'rnai'=list())
  #   set to null if no any test
  if (is.null(rnai)){
    rnai = load_rnai()
  }
  if (!'p53_state' %in% colnames(tcga[[1]]@colData)){
    p53_ann = annotate_sample_mut(tcga[[2]]@data)
    tcga[[1]]@colData[['p53_state']] = 'Wildtype'
    for (i in names(p53_ann)){
      tcga[[1]]@colData[p53_ann[[i]], 'p53_state'] = i
    }
  }
  if (!'p53_state' %in% colnames(ccle[[1]]@colData)){
    p53_ann = annotate_sample_mut(ccle[[2]]@data)
    ccle[[1]]@colData[['p53_state']] = 'Wildtype'
    for (i in names(p53_ann)){
      ccle[[1]]@colData[p53_ann[[i]], 'p53_state'] = i
    }
  }
  for (i in 1:length(mutation_groups)){
    b_m = get_binary_SNP_m_from_maf(ccle[[2]]@data, 
                                    snp_list = list(mutation_groups[[i]]),
                                    samples = rownames(ccle[[1]]@colData), 
                                    mode = 'position')
    ccle[[1]]@colData[paste('mutation_binary_state',i, sep='.')] = b_m
  }
  for (i in 1:length(mutation_groups)){
    b_m = get_binary_SNP_m_from_maf(tcga[[2]]@data, 
                                    snp_list = list(mutation_groups[[i]]),
                                    samples = rownames(tcga[[1]]@colData),
                                    protein_change_col = 'HGVSp_Short',
                                    mode = 'position')
    tcga[[1]]@colData[paste('mutation_binary_state',i, sep='.')] = b_m
  }
  rt_list = list()
  for (gene in genes){
    ### Visualisation
    ccle[[1]]@colData['gene_expr'] = t(assay(ccle[[1]][gene,,'RNA'])[rownames(ccle[[1]]@colData)])
    tcga[[1]]@colData['gene_expr'] = t(assay(tcga[[1]][gene,,'RNA'])[rownames(tcga[[1]]@colData)])
    ccle[[1]]@colData['gene_rnai'] = NA
    tcga[[1]]@colData['gene_rnai'] = NA
    flag = FALSE
    if (gene %in% rownames(rnai)){
      flag = TRUE
      co_cell = intersect(rownames(ccle[[1]]@colData), colnames(rnai))
      ccle[[1]]@colData[co_cell, 'gene_rnai'] = as.numeric(t(rnai[gene,co_cell]))
    }
    
    if (is.null(primary_site)){
      idx = 1:nrow(ccle[[1]]@colData)
    } else {
      idx = which(ccle[[1]]@colData$PRIMARY_SITE %in% primary_site)
    }
    col_to_use = c('gene_expr', 'p53_state', 'gene_rnai',
                   paste('mutation_binary_state', 1:length(mutation_groups), sep='.'))
    df_plt = rbind(as.data.frame(ccle[[1]]@colData[idx,col_to_use]) %>%
                     mutate(db='ccle'),as.data.frame(tcga[[1]]@colData[,col_to_use]) %>%
                     mutate(db='tcga'))
    df_plt['itg_state'] = NA
    for (i in 1:length(mutation_groups)){
      df_plt[which(df_plt[[paste('mutation_binary_state', i, sep='.')]]==1),'itg_state'] = names(mutation_groups)[i]
    }
    
    df_plt[which(df_plt$p53_state=='Wildtype'),'itg_state'] = 'Wildtype'
    df_plt[which(df_plt$p53_state=='nonsense'), 'itg_state'] = 'Null'
    df_plt = subset(df_plt, !is.na(itg_state))
    # table(df_plt$itg_state)
    # table(df_plt$db)
    # summary(df_plt$gene_rnai)
    df_plt$db = toupper(df_plt$db)
    a = ggplot(df_plt, aes(x=itg_state, y=gene_expr)) +
      geom_violin() +
      geom_boxplot(width=0.3, outlier.shape = NA) +
      geom_jitter(width=0.1, alpha=0.7, size=0.5) +
      mytme +
      facet_wrap(~db) +
      labs(x='p53 state', y=paste(gene, 'expression')) +
      EnvStats::stat_n_text(fontface = "italic") + 
      theme(strip.background = element_rect(fill='transparent'),
            strip.text = element_text(size=12, face='bold'),
            axis.text.x = element_text(angle=45, hjust = 1))
    if (!is.null(comparison$rna)){
      a = a + stat_compare_means(comparisons = comparison$rna, method = 't.test')
    }
    
    if (flag){
      b = ggplot(subset(df_plt, df_plt$db=='CCLE' & !is.na(df_plt$gene_rnai)), 
                 aes(x=itg_state, y=gene_rnai)) +
        geom_violin() +
        geom_boxplot(width=0.3, outlier.shape = NA) +
        geom_jitter(width=0.1, alpha=0.7, size=0.5) +
        EnvStats::stat_n_text(fontface = "italic") + 
        mytme +
        labs(x='p53 state', y=paste(gene, 'RNAi dependency')) +
        theme(strip.background = element_rect(fill='transparent'),
              strip.text = element_text(size=16),
              axis.text.x = element_text(angle=45, hjust = 1))
      if (!is.null(comparison$rnai)){
        b = b + stat_compare_means(comparisons = comparison$rnai, method = 't.test')
      }
      rt_list[[paste(gene, 'rna', sep='_')]] = a
      rt_list[[paste(gene, 'rnai', sep='_')]] = b
    } else {
      rt_list[[paste(gene, 'rna', sep='_')]] = a
      rt_list[[paste(gene, 'rnai', sep='_')]] = ggplot()
    }
  }
  return(rt_list)
}


