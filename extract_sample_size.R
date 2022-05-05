library(tidyverse)
# dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt_SIFT_noHot'
# dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt_noSIFT_noHot'
dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt'

### EXtract sample size ====
folder = list.files(file.path(dir_home, 'outputs'))
coll = data.frame()
coll_sample = data.frame()
for (fd in folder){
  if (fd=='tcga_ov_raw_seq'){next}
  log = read.table(file.path(dir_home, 'outputs', fd, 'log.txt'), sep=':', fill=T,
                   row.names = NULL, header = F, comment.char = '?')
  
  mut_sample = strsplit(log$V6[grep('mutated sample', log$V6)], split='\\.')[[1]][2]
  mut_sample = as.numeric(strsplit(mut_sample, split=' ')[[1]][2])
  total_sample = as.numeric(log$V6[grep('mapping sample size', log$V5)])
  coll_sample = rbind(coll_sample, 
                      c(toupper(strsplit(fd, split='_')[[1]][2]),
                        mut_sample, total_sample, mut_sample/total_sample))
  
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
colnames(coll_sample) = c('Cancer', 'Mut_sample', 'Total_sample', 'Ratio')
write.table(coll_sample, file.path(dir_home, 'mut_size.txt'), sep='\t', row.names = F, quote = F)
write.table(coll, file.path(dir_home, 'sample_size.txt'), sep='\t', row.names = F, quote = F)

### plot mutant count ====
coll_sample = read.table(file.path(dir_home, 'mut_size.txt'), sep='\t', header = T)
coll_sample = coll_sample[order(coll_sample$Ratio, decreasing = T),]
coll_sample$Cancer = factor(coll_sample$Cancer, levels=coll_sample$Cancer)
g = ggplot(coll_sample, aes(y=Cancer)) +
  geom_point(aes(x=Ratio, size=Mut_sample)) +
  mytme +
  labs(y='', x='Frequency of p53 mutant samples') +
  scale_size_continuous(name = 'count of p53 mutant samples', limits = c(150,450), breaks=seq(200,400,100)) +
  theme(axis.text.y = element_text(angle=0, vjust=0.8),
        legend.direction = 'horizontal', panel.grid.major = element_line(color='grey', size=0.5),
        legend.key = element_rect(fill='transparent'))
ggsave(file.path(dir_home, 'p53_mut_rate.pdf'),
       plot=g, width=6,height=4,units='in',device='pdf',dpi=300)


   ### visualise ====
source('/Users/jefft/Desktop/Manuscript/set_theme.R')
coll_plt = coll
coll_plt$Gene_NM_perMut = as.numeric(coll_plt$Gene_NM) / as.numeric(coll_plt$Mutant)
coll_plt = coll_plt[coll_plt$Mutation %in% c('hot_spot', 'contact', 'sandwich', 'conformation'),]
g = ggplot(coll_plt) +
  geom_bar(aes(x=Cancer, y=Gene_NM_perMut), stat='identity', position = 'dodge') +
  facet_wrap(~Mutation, scales = 'free') +
  scale_y_continuous(expand = c(0,0)) + theme_classic() +
  mytme +
  theme(axis.text.x = element_text(angle=45, vjust = 0.5), legend.direction = 'horizontal',
        strip.text = element_text(size=14))
ggsave(file.path(dir_home,'plots', 'Normalized_result_size.pdf'), bg = 'transparent',
       plot=g, height=8.27,width=11.69,units='in',device='pdf',dpi=300)

### plot mutation rate


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
