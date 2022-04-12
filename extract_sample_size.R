library(tidyverse)
dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt_SIFT_noHot'
dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt_noSIFT_noHot'
dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt'

### EXtract sample size ====
folder = list.files(file.path(dir_home, 'outputs'))
coll = data.frame()
for (fd in folder){
  log = read.table(file.path(dir_home, 'outputs', fd, 'log.txt'), sep=':', fill=T,
                   row.names = NULL, header = F, comment.char = '?')
  origin = grep('Processing',log$V5)
  qtl = read.table(file.path(dir_home, 'outputs', fd, 'trans_eqtl.txt'), header = T)
  if ('snps' %in% colnames(qtl)){
    df = data.frame('snp'=qtl$snps, 'pc'=qtl$protein_change)
  } else {
    df = data.frame('snp'=qtl$SNP, 'pc'=qtl$protein_change)
  }
  if (sum(duplicated(df))>0){
    df = df[-which(duplicated(df)),]
  }
  df = deframe(df)
  for (o in origin){
    id = strsplit(log[o,'V5'], split=' ')[[1]][3]
    sizes = trimws(strsplit(log[o+1,'V5'], split=',')[[1]])
    mut = strsplit(sizes[1], split = ' # ')[[1]][2]
    ctr = strsplit(sizes[2], split = ' # ')[[1]][2]
    nm = strsplit(log[o+2,'V5'], split=' ')[[1]][4]
    id = df[id]
    coll = rbind(coll, c(id,ctr,mut,nm,toupper(strsplit(fd, split = '_')[[1]][2])))
  }
}
colnames(coll) = c('Mutation', 'Control', 'Mutant','Gene_NM', 'Cancer')
write.table(coll, file.path(dir_home, 'sample_size.txt'), sep='\t', row.names = F, quote = F)

### compare ====
source('/Users/jefft/Desktop/Manuscript/set_theme.R')
plot_out = '/Users/jefft/Desktop/p53_project'
df_noSIFT_noHot = read.table('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt_noSIFT_noHot/sample_size.txt', header = T)
df_SIFT_noHot = read.table('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt_SIFT_noHot/sample_size.txt', header = T)
df_SIFT_hot = read.table('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt/sample_size.txt', header = T)

df = inner_join(df_SIFT_hot, df_SIFT_noHot, by=c('Mutation', 'Cancer'), suffix=c('SIFT_hot','SIFT_noHot'))
df = inner_join(df, df_noSIFT_noHot, by=c('Mutation', 'Cancer'))
ncol(df)
colnames(df)[9:11] = c('ControlnoSIFTnoHot', 'MutantnoSIFTnoHot','Gene_NMnoSIFTnoHot')
df = df[,-grep('Control', colnames(df))]

df_size = gather(df, key='gp', value='count', c(2,5,7))
df_size$gp = factor(df_size$gp, levels = c('MutantnoSIFTnoHot','MutantSIFT_noHot', 'MutantSIFT_hot'))
g = ggplot(df_size[df_size$Mutation %in% c('hot_spot', 'contact', 'sandwich', 'conformation'),]) +
  geom_bar(aes(x=Cancer, y=count, fill=gp), stat='identity', position = 'dodge') +
  facet_wrap(~Mutation, scales = 'free') +
  scale_y_continuous(expand = c(0,0)) + theme_classic() +
  mytme +
  theme(axis.text.x = element_text(angle=45, vjust = 0.5), legend.direction = 'horizontal',
        strip.text = element_text(size=14)) 
ggsave(file.path(plot_out, 'Mutant_sample_size.pdf'), bg = 'transparent',
       plot=g, height=8.27,width=11.69,units='in',device='pdf',dpi=300)

df_freq = df
df_freq$Gene_NMSIFT_hot = df_freq$Gene_NMSIFT_hot / df_freq$MutantSIFT_hot
df_freq$Gene_NMnoSIFTnoHot = df_freq$Gene_NMnoSIFTnoHot / df_freq$MutantnoSIFTnoHot
df_freq$Gene_NMSIFT_noHot = df_freq$Gene_NMSIFT_noHot / df_freq$MutantSIFT_noHot
df_freq = gather(df_freq[,c(1,4,3,6,8)], key='gp', value='count', 3:5)
df_freq$gp = factor(df_freq$gp, levels = c('Gene_NMnoSIFTnoHot','Gene_NMSIFT_noHot', 'Gene_NMSIFT_hot'))
df_freq = df_freq[df_freq$Mutation %in% c('hot_spot', 'contact', 'sandwich', 'conformation'),]
g = ggplot(df_freq) +
  geom_bar(aes(x=Cancer, y=count, fill=gp), stat='identity', position = 'dodge') +
  facet_wrap(~Mutation, scales = 'free') +
  scale_y_continuous(expand = c(0,0)) + theme_classic() +
  mytme +
  theme(axis.text.x = element_text(angle=45, vjust = 0.5), legend.direction = 'horizontal',
        strip.text = element_text(size=14)) 
ggsave(file.path(plot_out, 'Normalized_result_size.pdf'), bg = 'transparent',
       plot=g, height=8.27,width=11.69,units='in',device='pdf',dpi=300)
