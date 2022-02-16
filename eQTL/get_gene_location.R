### Get gene location table from GTF file for MatrixEQTL package
library(rtracklayer)
library(dplyr)

generate_gene_table = function(gtf_path, output_path, type_col='source', protein_coding_only=TRUE){
  gtf = import(gtf_path)
  gtf.df = as.data.frame(gtf)
  if ('transcript' %in% gtf.df$type){
    gtf.df = subset(gtf.df, gtf.df$type=='transcript')
  }
  gtf.df['source'] = gtf.df[[type_col]]
  gtf.df[which(is.na(gtf.df['source'])), 'source'] = 'NA'
  if (protein_coding_only){
    gtf.df = subset(gtf.df, gtf.df$source=='protein_coding')
  }
  smm = as.data.frame(gtf.df %>% group_by(transcript_id, gene_name) %>% 
                      summarise(tp_start=min(start), tp_end=max(end)))
  smm['tp_len'] = smm$tp_end - smm$tp_start
  rownames(smm) = as.character(1:nrow(smm))
  source_gene = as.data.frame(gtf.df %>% group_by(gene_name) %>%
    summarise(sc = unique(source)))
  sc_lookup = source_gene$sc
  names(sc_lookup) = source_gene$gene_name
  seq_gene = as.data.frame(gtf.df %>% group_by(gene_name) %>%
    summarise(seq = unique(seqnames)))
  seq_lookup = seq_gene$seq
  names(seq_lookup) = seq_gene$gene_name
  
  idx = 1:nrow(smm)
  lookup = idx
  names(lookup) = paste(smm$gene_name, smm$tp_len, sep='-')
  
  # set the longest transcript as selected
  smm_2 = as.data.frame(smm %>% group_by(gene_name) %>%
    summarise(gene_len=max(tp_len)))
  smm_ftd = smm[lookup[paste(smm_2$gene_name, smm_2$gene_len, sep='-')],]
  smm_ftd['source'] = sc_lookup[smm_ftd$gene_name]
  smm_ftd = smm_ftd[,which(colnames(smm_ftd) %in% c('gene_name', 'tp_start', 'tp_end', 'source'))]
  smm_ftd['chr'] = seq_lookup[smm_ftd$gene_name]
  # only autosome
  if (length(grep('chr', smm_ftd$chr))==0){
    smm_ftd = subset(smm_ftd, smm_ftd$chr %in% 1:22)
    smm_ftd$chr = paste('chr', smm_ftd$chr, sep='')
  }
  if (protein_coding_only){
    smm_ftd = smm_ftd[,c(1,5,2,3)]
    colnames(smm_ftd) = c('geneid', 'chr', 'left', 'right')
  } else {
    smm_ftd = smm_ftd[,c(1,5,2,3,4)]
    colnames(smm_ftd) = c('geneid', 'chr', 'left', 'right', 'source')
  }
  write.table(smm_ftd, output_path, 
              sep = '\t', quote = FALSE, row.names = FALSE)
}

# hg18
generate_gene_table(gtf_path = '/Users/jefft/Genome/NCBI36_annotation/Homo_sapiens.NCBI36.52.gtf.gz',
                    output_path = '/Users/jefft/Desktop/p53_project/scripts/eQTL/hg18_gene_table_autosome.tsv',
                    type_col = 'source', protein_coding_only = T)

# hg19 (GRCh37)
generate_gene_table(gtf_path = '/Users/jefft/Genome/GRCh37_annotation/gencode.v19.annotation.gtf.gz',
                    output_path = '/Users/jefft/Desktop/p53_project/scripts/eQTL/hg19_gene_table_autosome.tsv',
                    type_col = 'transcript_type', protein_coding_only = T)

# hg38 (GRCh38)
generate_gene_table(gtf_path = '/Users/jefft/Genome/GRCh38_annotation/Homo_sapiens.GRCh38.99.gtf',
                    output_path = '/Users/jefft/Desktop/p53_project/scripts/eQTL/hg38_gene_table_autosome.tsv',
                    type_col = 'transcript_biotype', protein_coding_only = T)

