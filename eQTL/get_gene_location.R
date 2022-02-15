### Get gene location table from GTF file for MatrixEQTL package
library(rtracklayer)
library(dplyr)

gtf = import('/Users/jefft/Genome/NCBI36/Homo_sapiens.NCBI36.52.gtf.gz')
gtf.df = as.data.frame(gtf)
smm = gtf.df %>% group_by(transcript_id, gene_name) %>% 
  summarise(tp_start=min(start), tp_end=max(end))
smm['tp_len'] = smm$tp_end - smm$tp_start
rownames(smm) = as.character(1:nrow(smm))
source_gene = gtf.df %>% group_by(gene_name) %>%
  summarise(sc = unique(source))
sc_lookup = source_gene$sc
names(sc_lookup) = source_gene$gene_name
seq_gene = gtf.df %>% group_by(gene_name) %>%
  summarise(seq = unique(seqnames))
seq_lookup = seq_gene$seq
names(seq_lookup) = seq_gene$gene_name

idx = 1:nrow(smm)
lookup = idx
names(lookup) = paste(smm$gene_name, smm$tp_len, sep='-')

# set the longest transcript as selected
smm_2 = smm %>% group_by(gene_name) %>%
  summarise(gene_len=max(tp_len))
smm_ftd = smm[lookup[paste(smm_2$gene_name, smm_2$gene_len, sep='-')],]
smm_ftd['source'] = sc_lookup[smm_ftd$gene_name]
smm_ftd = smm_ftd[,which(colnames(smm_ftd) %in% c('gene_name', 'tp_start', 'tp_end', 'source'))]
smm_ftd['chr'] = seq_lookup[smm_ftd$gene_name]
# only autosome
smm_ftd = subset(smm_ftd, smm_ftd$chr %in% 1:22)
smm_ftd$chr = paste('chr', smm_ftd$chr, sep='')
smm_ftd = smm_ftd[,c(1,5,2,3,4)]
colnames(smm_ftd) = c('geneid', 'chr', 'left', 'right', 'source')
write.table(smm_ftd, '/Users/jefft/Desktop/p53_project/eQTL/hg18_gene_table_autosome.tsv', 
            sep = '\t', quote = FALSE, row.names = FALSE)




