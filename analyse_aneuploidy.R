setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
# TCGA-pan_VS-wt, TCGA-pan_VS-mutneg_ult
# ccle_home = '/Users/jefft/Desktop/p53_project/datasets/CCLE/processed_dt'
ccle_home = '/Users/jefft/Desktop/p53_project/datasets/CCLE_22Q1/pcd'
dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt_ploid'
eqtl_out = file.path(dir_home, 'outputs')
plot_out = file.path(dir_home, 'plots')
data_out = file.path(dir_home, 'data_out')
source('utils.R')
source('enrich_utils.R')
source('../ccle_utils.R')
source('../overlap_utils.R')
source('../revigo_utils.R')
source('../set_theme.R')
library(tidyverse)
library(stringr)
library(ggrepel)
library(ggsci)

# Count of controls ====
load('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt/data_out/control_summary.RData')
count_inA = sm_tb_stat_coll
load('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt_ploid/data_out/control_summary.RData')
count_exA = sm_tb_stat_coll
for (i in names(count_inA)){
  count_inA[[i]] = full_join(count_inA[[i]], count_exA[[i]], by='cancer', suffix=c('inA','exA'))
  count_inA[[i]]$group = i
}
count_coll = Reduce(rbind, count_inA)
count_coll = gather(count_coll, key='condition',value='count', c(2,5)) # two hit col
count_coll$cancer = factor(count_coll$cancer, levels = c('BRCA', 'LGG', 'BLCA', 'COAD', 'HNSC', 'STAD'))
g = ggplot(count_coll, aes(x=cancer, y=count)) +
  geom_bar(aes(fill=condition), position = 'dodge', stat = 'identity') +
  facet_wrap(~group, scale='free_y') +
  mytme +
  scale_y_continuous(expand=expansion(mult = c(0,0.1)), breaks=seq(0,25,5)) +
  scale_fill_manual(name='', values = c('#BCBD22FF','#E377C2FF'), labels=c('- ploidy', '+ ploidy')) +
  theme(strip.background = element_rect(fill='transparent'),
        legend.direction = 'horizontal',
        strip.text = element_text(face = 'bold', size=12),
        axis.text.x = element_text(angle=45, hjust = 1, vjust=1),
        axis.title.y = element_text(size=14)) +
  labs(x='', y='Count of\ndetected controls')
g
ggsave(file.path(plot_out,'control', 'compare_control_ploid_count.pdf'),
       plot=g, width=11.69*0.8, height=8.27*0.4,units='in',device='pdf',dpi=300)



# Load Gene ========
beta_cutoff = 0
coll = load_eQTL_output('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt/outputs', 
                        beta=beta_cutoff, exclude = c('tcga_nine_pool', 'metabric_raw_its'), as_df = T)
coll_ploid = load_eQTL_output('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt_ploid/outputs',
                              beta=beta_cutoff, exclude = c('tcga_nine_pool', 'metabric_raw_its'), as_df = T)
coll = coll[coll$protein_change=='hot_spot',]
coll_ploid = coll_ploid[coll_ploid$protein_change=='hot_spot',]
smm = full_join(coll %>% group_by(cancer) %>% summarise(count_pos=sum(beta>0), count_neg=sum(beta<0)),
          coll_ploid %>% group_by(cancer) %>% summarise(count_pos=sum(beta>0), count_neg=sum(beta<0)),
          by='cancer', suffix = c('.pldin','.pldout'))
smm[is.na(smm)] = 0
smm_plt = smm[smm$cancer!='LUSC',]
smm_plt = gather(smm_plt, key='group', value='count', 2:ncol(smm))
gps = as.data.frame(t(as.data.frame(strsplit(smm_plt$group, split='\\.'))))
colnames(gps) = c('sign', 'condition')
rownames(gps) = NULL
smm_plt = cbind(smm_plt, gps)
smm_plt$cancer = factor(smm_plt$cancer, levels=c('BRCA', 'LGG','BLCA','COAD', 'HNSC', 'STAD'))
g = ggplot(smm_plt) +
  geom_bar(aes(x=sign,y=count,fill=condition), stat = 'identity', position = 'dodge') +
  facet_wrap(~cancer, scale='free') + mytme +
  scale_fill_manual(name='', values=c('#E377C2FF','#BCBD22FF'), labels=c('+ ploidy', '- ploidy')) +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  theme(
        strip.background = element_rect(fill='transparent'),
        strip.text.x = element_text(face='bold', size=11)) +
  labs(x='',y='Gene count') +
  scale_x_discrete(labels=c('Negative', 'Positive'))
ggsave(file.path(plot_out, 'compare_ploid_count.pdf'),
       plot=g, width=11.69*0.7, height=8.27*0.8,units='in',device='pdf',dpi=300)

for (j in 1:nrow(smm)){
  print(smm$cancer[j])
  print(chisq.test(matrix(as.matrix(smm[j,2:5]),nrow=2))$p.value)
}

# test genes ====
coll$sign = 'neg'
coll$sign[coll$beta>0] = 'pos'
coll_ploid$sign = 'neg'
coll_ploid$sign[coll_ploid$beta>0] = 'pos'
use_col = c('cancer','gene','pvalue','FDR','beta','sign')
joinTb = full_join(coll[,use_col], coll_ploid[,use_col], 
        by=c('gene','cancer'), suffix = c('.pldin','.pldout'))
joinTb[is.na(joinTb)] = 0
joinTb = joinTb[joinTb$cancer!='LUSC',]
joinTb$missing = joinTb$beta.pldin * joinTb$beta.pldout==0
joinTb$cancer = factor(joinTb$cancer, levels = c('BRCA','LGG','BLCA','COAD','HNSC','STAD'))
joinTb$M = joinTb$beta.pldin - joinTb$beta.pldout
joinTb$A = joinTb$beta.pldin + joinTb$beta.pldout

for (check_cancer in unique(joinTb$cancer)){
  print(check_cancer)
  coad_join = joinTb[joinTb$cancer==check_cancer,]
  same_pos = sum(coad_join$beta.pldin>0 & coad_join$beta.pldout>0)
  same_neg = sum(coad_join$beta.pldin<0 & coad_join$beta.pldout<0)
  plod_pos = sum(coad_join$beta.pldin>0 & coad_join$beta.pldout<=0 |
                   coad_join$beta.pldin==0 & coad_join$beta.pldout<0)
  plod_neg = sum(coad_join$beta.pldin<0 & coad_join$beta.pldout>=0 |
                   coad_join$beta.pldin==0 & coad_join$beta.pldout>0)
  count_mtx = matrix(c(same_pos,same_neg,plod_pos,plod_neg),ncol=2)
  colnames(count_mtx) = c('sign_same', 'sign_dif')
  rownames(count_mtx) = c('pos','neg')
  print(count_mtx)
  print(chisq.test(count_mtx)$p.value) # is chi-square valid here?
}

# plotting
cel = ceiling(max(c(abs(joinTb$beta.pldin), abs(joinTb$beta.pldout))))
cel = ceiling(max(abs(joinTb$M), abs(joinTb$A)))
g = ggplot(joinTb[!joinTb$missing,], 
       aes(x=beta.pldin, y=beta.pldout)
       # aes(x=A, y=M)
       ) +
  # MA plot scale
  #scale_x_continuous(limits = c(-5,5), breaks = seq(-5,5,2.5)) +
  # scale_y_continuous(limits = c(-2,2), breaks = seq(-2,2,1)) +
  # Non-MA scale
  scale_x_continuous(limits = c(-3,3), breaks = seq(-3,3,2)) +
  scale_y_continuous(limits = c(-3,3), breaks = seq(-3,3,2)) +
  geom_vline(xintercept = 0, color='grey30') +
  geom_hline(yintercept = 0, color='grey30') +
  geom_abline(intercept = 0, slope = 1, color='grey30', linetype='dotted') +
  #geom_abline(intercept = 0, slope = -1, color='grey30', linetype='dotted') +
  geom_point(size=0.1, color='coral',alpha=0.5) + mytme +
  geom_jitter(data=joinTb[joinTb$missing & joinTb$beta.pldin==0,],color='royalblue', size=0.1, alpha=0.5,
              width = 0.2, height = 0) +
  geom_jitter(data=joinTb[joinTb$missing & joinTb$beta.pldout==0,],color='royalblue', size=0.1, alpha=0.5,
              width = 0, height = 0.2) +
  theme(strip.background = element_rect(fill=NA),
        strip.text = element_text(face='bold', size=12)) +
  labs(x='Beta not considering aneuploidy', y='Beta aneuploidy as a covariate') +
  facet_wrap(~cancer)
ggsave(file.path(plot_out, 'beta_overview.pdf'),
       plot=g, width=11.69*0.7, height=8.27*0.7,units='in',device='pdf',dpi=300)

missing_gene = joinTb[joinTb$missing,]
missing_gene$group = NA
for (i in 1:nrow(missing_gene)){
  if (missing_gene$beta.pldin[i]>0 | (missing_gene$beta.pldin[i]==0 & missing_gene$beta.pldout[i]<0)){
    missing_gene$group[i] = 'upAne'
  } else {
    missing_gene$group[i] = 'downAne'
  }
}

total_cancer = joinTb %>% group_by(cancer) %>% summarise(nm=length(gene))
total_cancer = deframe(total_cancer)
missing_smm = missing_gene %>% group_by(cancer,group) %>% summarise(pct=length(gene))
missing_smm$pct = missing_smm$pct / total_cancer[missing_smm$cancer]
missing_smm = missing_smm[order(missing_smm$pct),]
missing_smm$cancer = factor(missing_smm$cancer, levels = unique(missing_smm$cancer))
g = ggplot(missing_smm) + mytme +
  geom_bar(aes(x=cancer,y=pct*100,fill=group),
           stat = 'identity',position = 'stack') +
  labs(x='',y='% total genes') +
  scale_y_continuous(expand = c(0,0),breaks = seq(10,100,10)) +
  scale_fill_d3(palette = 'category20', labels = c('Negative interaction','Positive interaction')) +
  theme(legend.title = element_blank(), legend.direction = 'horizontal')
ggsave(file.path(plot_out, 'missing_gene.pdf'),
       plot=g, height=11.69*0.4, width=8.27*0.6,units='in',device='pdf',dpi=300)


# missing genes and their correlation to the aneuploidy ====
check_cancer = 'COAD'
coad_join = joinTb[joinTb$cancer==check_cancer,]
test = coad_join[coad_join$missing,]
## choose one !!!
test = test[order(test$beta.pldin, decreasing = T),]
test = test[order(test$beta.pldout, decreasing = T),]
## load data
coad = load_clean_data('/Users/jefft/Desktop/p53_project/datasets/4-COAD-TCGA/clean_data.RData',
                       ann_p53 = T, ann_bin_mut_list = c('hot_spot'))

## choose one !!!
goi = test$gene[nrow(test)] # negative interacttion
goi = test$gene[1] # positive interaction

for (goi in c('RAB26','NDN','TMEM170B','WWC1')){
  plt_dt = cbind(as.data.frame(coad[[1]]@colData), 
                 t(assay(coad[[1]][goi,,'RNA']))) # pick the first gene
  plt_dt = plt_dt[plt_dt$p53_state %in% c('missense', 'Wildtype'),]
  ggplot(plt_dt, aes(x=has_hot_spot,y=get(goi),color=ANEUPLOIDY_SCORE)) +
    geom_jitter(aes(group=has_hot_spot)) +
    labs(y=goi) +
    mytme + scale_color_gradientn(colours=colorRampPalette(c('red','blue'))(100))
  
  g = ggplot(plt_dt, aes(x=ANEUPLOIDY_SCORE,y=get(goi),color=as.factor(has_hot_spot),group=has_hot_spot)) +
    geom_jitter() + geom_smooth(method = 'lm') +
    labs(y=goi) +
    facet_wrap(~has_hot_spot, nrow=2) +
    scale_color_d3(palette = 'category20') +
    mytme +
    stat_cor()
  ggsave(file.path(plot_out, paste(goi, '_interact.pdf', sep='')),
         plot=g, height=11.69*0.4, width=8.27*0.6,units='in',device='pdf',dpi=300)
}

g = ggplot(joinTb[joinTb$gene %in% c('RAB26','NDN','TMEM170B','WWC1') & joinTb$cancer=='COAD',], 
           aes(x=beta.pldin, y=beta.pldout)) +
  geom_point(size=5, color='coral',alpha=0.5) + mytme +
  geom_text(aes(label=gene), nudge_x = 0.1,nudge_y=0.1) +
  labs(x='Beta not considering aneuploidy', y='Beta aneuploidy as a covariate')
ggsave(file.path(plot_out, 'exp_gene_loc.pdf'),
       plot=g, height=11.69*0.4, width=8.27*0.6,units='in',device='pdf',dpi=300)


## genome wide test ====
cor_dt = cbind(as.data.frame(coad[[1]]@colData)[,c('ANEUPLOIDY_SCORE','p53_state', 'has_hot_spot')], t(assay(coad[[1]][,,'RNA'])))
cor_dt = cor_dt[cor_dt$p53_state %in% c('Wildtype', 'missense'),]
aneu = cor_dt$ANEUPLOIDY_SCORE
ps = c()
cors = c()
for (i in 4:ncol(cor_dt)){
  tst = cor.test(aneu,cor_dt[,i])
  ps = c(ps, tst$p.value)
  cors = c(cors, tst$estimate)
}
aneu_corr = data.frame('gene'=colnames(cor_dt)[4:ncol(cor_dt)],'pvalue'=ps, 'cor'=cors)
rownames(aneu_corr) = aneu_corr$gene
coad_join$aneu_corr = aneu_corr[coad_join$gene,'cor']
coad_join$aneu_pvalue = aneu_corr[coad_join$gene,'pvalue']
g = ggplot(coad_join,aes(x=missing, y=abs(aneu_corr))) +
  geom_violin() + geom_boxplot(width=0.2) +
  labs(x='Interaction', y='Correlation to aneuploidy score', title='Gene expression in COAD') +
  stat_compare_means(method = 't.test') + mytme
ggsave(file.path(plot_out, 'COAD_correlation_aneu.pdf'),
       plot=g, height=11.69*0.4, width=8.27*0.5,units='in',device='pdf',dpi=300)


hist(coad_join[coad_join$missing,'aneu_corr'],breaks=100)
g = ggplot(coad_join[coad_join$missing,], 
       aes(x=beta.pldin, y=beta.pldout,color=aneu_corr)) +
  scale_color_gradient2(low='#17BECFFF',high='#D62728FF',mid='white') +
  geom_jitter(width = 0.2,height = 0.2,size=0.1) + mytme +
  theme(legend.direction = 'horizontal', legend.key.width = unit(0.4,'in'))
ggsave(file.path(plot_out, 'COAD_missing_overview.pdf'),
       plot=g, height=11.69*0.4, width=8.27*0.5,units='in',device='pdf',dpi=300)


