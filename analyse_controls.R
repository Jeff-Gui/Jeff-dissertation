setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
# TCGA-pan_VS-mutneg_ult, TCGA-pan_VS-wt, TCGA-pan_VS-null
dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-null'
eqtl_out = file.path(dir_home, 'outputs')
plot_out = file.path(dir_home, 'plots', 'control')
data_out = file.path(dir_home, 'data_out')
library(tidyverse)
library(pheatmap)
library(patchwork)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
library(ggsci)
source('utils.R')
source('../ccle_utils.R')
source('../set_theme.R')

gc()
# output_fp = 'outputs/TEST_BRCA'
output_fp = eqtl_out
coll = load_eQTL_output(output_fp, mode='fdr-no', exclude = 'tcga_nine_pool')
if (length(which(is.na(coll))>0)){
  coll = coll[-which(is.na(coll))]
}
coll = Reduce(rbind, coll)
# config which mutation group(s) to plot
table(coll$protein_change)
name_map = read.csv('/Users/jefft/Desktop/p53_project/datasets/mut_name_map.csv', header=T)
mutation_to_plt = name_map$label
names(mutation_to_plt) = name_map$code

mutation_to_plt = mutation_to_plt[c('hot_spot', 'contact', 'conformation', 'sandwich', 'p.R273H', 'p.R175H', 'breast_w2016')] # 
detected_exp = list()
for (ep in unique(coll$experiment)){
  detected_exp[[ep]] = intersect(mutation_to_plt, mutation_to_plt[unique(subset(coll, coll$experiment==ep)$protein_change)])
}


all_mut = names(table(coll$protein_change))
mut_exclude = all_mut[which(!all_mut %in% name_map$code)]
print(mut_exclude)  # this mutant will not be ploted

meta_mut_dir = '/Users/jefft/Desktop/p53_project/datasets/meta_muts'
meta_muts = load_meta_mut(file.path(meta_mut_dir, list.files(meta_mut_dir)))

# load control genes
library(readxl)
ctrs_raw = read_xlsx('/Users/jefft/Desktop/p53_project/Thesis/gene_signatures/collection.xlsx')
ctrs = na.omit(ctrs_raw)
ctrs = ctrs[order(ctrs$Gene_annotation),c('Gene', 'Gene_annotation')]
ctrs = ctrs[-which(duplicated(ctrs)),]
# check if control genes are expressed
load('/Users/jefft/Desktop/p53_project/datasets/TCGA-Pan-Nine/gene_matrix.RData')
mtx = as.data.frame(mtx)

# summarize discovered genes, for doing stat test
smm_cancer = data.frame()
for (i in unique(coll$experiment)){
  cr = toupper(strsplit(i, split='_')[[1]][2])
  all = sum(mtx[,cr]==1)
  pos = length(unique(coll[coll$experiment==i & coll$beta>0 & 
                           coll$protein_change=='hot_spot' & coll$FDR<0.05,'gene']))
  neg = length(unique(coll[coll$experiment==i & coll$beta<0 & 
                           coll$protein_change=='hot_spot' & coll$FDR<0.05,'gene']))
  smm_cancer = rbind(smm_cancer, c(cr,all,pos,neg))
}
colnames(smm_cancer) = c('Cancer', 'All', 'Positive', 'Negative')

no_exp = c()
for (gn in unique(ctrs$Gene)){
  if (!gn %in% rownames(mtx)){
    no_exp = c(no_exp, gn)
  }
}
# CXCL8,SUSD6,TIGAR,PTCHD4
print(paste('Genes not detected:', paste(no_exp, collapse = ',')))
print(length(no_exp))
ctrs = ctrs[which(!ctrs$Gene %in% no_exp),]
ctrs_clean = ctrs

ctrs = ctrs_clean
for (epr in unique(coll$experiment)){
  sub_epr = subset(coll, experiment == epr)
  for (snp in unique(names(mutation_to_plt))){
    nm = paste(epr, snp, sep='|')
    ctrs[[paste(nm, 'FDR', sep='|')]] = NA
    ctrs[[paste(nm, 'beta', sep='|')]] = NA
    sub_c = subset(sub_epr, protein_change==snp)
    rownames(sub_c) = sub_c$gene
    hit = intersect(ctrs$Gene, sub_c$gene)
    for (h in hit){
      ctrs[which(ctrs$Gene==h), paste(nm, 'FDR', sep='|')] = sub_c[h, 'FDR']
      ctrs[which(ctrs$Gene==h), paste(nm, 'beta', sep='|')] = sub_c[h, 'beta']
    }
  }
}

ctrs_long = gather(ctrs, key='group', value='value', 3:ncol(ctrs))
ctr_long = cbind(ctrs_long, t(as.data.frame(strsplit(ctrs_long$group, split = '\\|'))))
rownames(ctr_long) = NULL
colnames(ctr_long)[(ncol(ctr_long)-2):ncol(ctr_long)] = c('experiment', 'mutation', 'type')
ctr_long = ctr_long[,-which(colnames(ctr_long)=='group')]
ctr_long = spread(ctr_long, key=type, value=value)

gene_order = ctrs$Gene[order(ctrs$Gene_annotation, ctrs$Gene)]
ctr_long$Gene = factor(ctr_long$Gene, unique(gene_order))
ctr_long['wrong'] = NA
ctr_long[which(ctr_long$beta<0), 'wrong'] = T
ctr_long[which(ctr_long$beta>0), 'wrong'] = F

# ctr_long = subset(ctr_long, ctr_long$mutation %in% names(mutation_to_plt))
ctr_long$mutation = mutation_to_plt[ctr_long$mutation]
ctr_long$mutation = factor(ctr_long$mutation, levels = mutation_to_plt)

# ctr_long$beta[which(ctr_long$beta <= 1)] = NA

ctr_long$cancer = toupper(sapply(ctr_long$experiment, function(x){return(strsplit(x, split='_')[[1]][2])}))

# see how many genes cen be recovered in each cancer, hotspot only ====
sm_tb_stat_coll = list()
xtme = theme(axis.text.x = element_text(angle=45, size=12, vjust = 0.7),
             legend.text = element_text(size=12), legend.direction = 'horizontal')
ctr_long = ctr_long[ctr_long$cancer!='OV',]
#### wt control one
ctr1 = subset(ctr_long, ctr_long$Gene_annotation=='wt_control')
ctr1 = ctr1[ctr1$mutation=='Hotspots',]
ctr1[ctr1$FDR>0.05 & !is.na(ctr1$FDR), c('beta','FDR','wrong')] = NA
trh = ceiling(length(unique(ctr1$Gene))/2) # must get half of genes recovered
sm_tb = ctr1 %>% group_by(cancer) %>% summarise(hit=length(which(wrong==TRUE)),
                                        wrong_hit=length(which(wrong==FALSE)),
                                        not_detected=sum(is.na(beta)))
sm_tb = sm_tb[order(sm_tb$hit, decreasing = T),]
sm_tb_stat = sm_tb[!sm_tb$cancer %in% c('LUAD','LUSC','OV'),]
print(sd(sm_tb_stat$hit / rowSums(sm_tb_stat[,2:4])))
print(mean((sm_tb_stat$hit / rowSums(sm_tb_stat[,2:4]))))
sm_tb_stat_coll[['ctr1']] = sm_tb_stat
ftr_ctr1 = sm_tb$cancer
sm_tb = gather(sm_tb, key='gp', value='count', 2:4)
sm_tb$cancer = factor(sm_tb$cancer, levels = ftr_ctr1)
ctr1_count = ggplot(sm_tb, aes(x=cancer, y=count)) +
  geom_bar(aes(fill=gp), stat='identity',position = 'stack') +
  scale_y_continuous(breaks = trh, expand = c(0,0)) +
  geom_hline(yintercept = trh, linetype='dotted') +
  labs(y='Count', x='', title='Downregulated controls 1') +
  scale_fill_d3(palette = 'category20', name = '', labels = c('Recovered','Not recovered','Recovered with the wrong sign')) +
  mytme + xtme
sm_tb_ctr1 = sm_tb[sm_tb$gp=='hit',]

#### wt control 2
ctr1 = subset(ctr_long, ctr_long$Gene_annotation=='wt_control_2')
ctr1 = ctr1[ctr1$mutation=='Hotspots',]
ctr1[ctr1$FDR>0.05 & !is.na(ctr1$FDR), c('beta','FDR','wrong')] = NA
trh = ceiling(length(unique(ctr1$Gene))/2) # must get half of genes recovered
sm_tb = ctr1 %>% group_by(cancer) %>% summarise(hit=length(which(wrong==TRUE)),
                                                wrong_hit=length(which(wrong==FALSE)),
                                                not_detected=sum(is.na(beta)))
sm_tb = sm_tb[order(sm_tb$hit, decreasing = T),]
sm_tb_stat = sm_tb[!sm_tb$cancer %in% c('LUAD','LUSC','OV'),]
print(sd(sm_tb_stat$hit / rowSums(sm_tb_stat[,2:4])))
print(mean((sm_tb_stat$hit / rowSums(sm_tb_stat[,2:4]))))
sm_tb_stat_coll[['ctr2']] = sm_tb_stat
ftr_ctr2 = sm_tb$cancer
sm_tb = gather(sm_tb, key='gp', value='count', 2:4)
sm_tb$cancer = factor(sm_tb$cancer, levels = ftr_ctr2)
ctr2_count = ggplot(sm_tb, aes(x=cancer, y=count)) +
  geom_bar(aes(fill=gp), stat='identity',position = 'stack') +
  scale_y_continuous(breaks = trh, expand = c(0,0)) +
  geom_hline(yintercept = trh, linetype='dotted') +
  labs(y='Count', x='', title='Downregulated controls 2') +
  scale_fill_d3(palette = 'category20', name = '', labels = c('Hit','Not detected','Wrong hit')) +
  mytme + xtme + theme(legend.position = 'none')
sm_tb_ctr2 = sm_tb[sm_tb$gp=='hit',]

#### positive controls
ctr1 = subset(ctr_long, !ctr_long$Gene_annotation %in% c('wt_control', 'wt_control_2'))
ctr1 = ctr1[ctr1$mutation=='Hotspots',]
ctr1[ctr1$FDR>0.05 & !is.na(ctr1$FDR), c('beta','FDR','wrong')] = NA
trh = ceiling(length(unique(ctr1$Gene))/2) # must get half of genes recovered
sm_tb = ctr1 %>% group_by(cancer) %>% summarise(hit=length(which(wrong==FALSE)),
                                                wrong_hit=length(which(wrong==TRUE)),
                                                not_detected=sum(is.na(beta)))
sm_tb = sm_tb[order(sm_tb$hit, decreasing = T),]
sm_tb_stat = sm_tb[!sm_tb$cancer %in% c('LUAD','LUSC','OV'),]
print(sd(sm_tb_stat$hit / rowSums(sm_tb_stat[,2:4])))
print(mean((sm_tb_stat$hit / rowSums(sm_tb_stat[,2:4]))))
sm_tb_stat_coll[['pos']] = sm_tb_stat
ftr_pos = sm_tb$cancer
sm_tb = gather(sm_tb, key='gp', value='count', 2:4)
sm_tb$cancer = factor(sm_tb$cancer, levels = ftr_pos)
pos_count = ggplot(sm_tb, aes(x=cancer, y=count)) +
  geom_bar(aes(fill=gp), stat='identity',position = 'stack') +
  scale_y_continuous(breaks = trh, expand = c(0,0)) +
  geom_hline(yintercept = trh, linetype='dotted') +
  labs(y='Count', x='', title='Upregulated controls') +
  scale_fill_d3(palette = 'category20', name = '', labels = c('Hit','Not detected','Wrong hit')) +
  mytme + xtme + theme(legend.position = 'none')
sm_tb_pos = sm_tb[sm_tb$gp=='hit',]

### Rank cancers
mtx = matrix(0, nrow=length(unique(sm_tb$cancer)), ncol=3)
rownames(mtx) = unique(sm_tb$cancer)
colnames(mtx) = c('ctr1', 'ctr2', 'pos')
mtx = as.data.frame(mtx)
for (i in 1:nrow(sm_tb_ctr1)){
  mtx[as.character(sm_tb_ctr1$cancer[i]), 'ctr1'] = sm_tb_ctr1$count[i]
}
for (i in 1:nrow(sm_tb_ctr2)){
  mtx[as.character(sm_tb_ctr2$cancer[i]), 'ctr2'] = sm_tb_ctr2$count[i]
}
for (i in 1:nrow(sm_tb_pos)){
  mtx[as.character(sm_tb_pos$cancer[i]), 'pos'] = sm_tb_pos$count[i]
}
mtx$sumHit = rowSums(mtx)
mtx = mtx[order(mtx$sumHit, decreasing = T),]
mtx$cancer = factor(rownames(mtx), levels=rownames(mtx))
rank_pfl = ggplot(mtx) +
  geom_bar(aes(x=cancer, y=sumHit), stat = 'identity') +
  scale_y_continuous(expand = c(0,0)) +
  labs(x='', y='Sum of hit counts') +
  scale_fill_d3() +
  mytme +xtme

# Test if there is any enrichment in eQTL analysis
fishp = matrix(0, nrow=length(unique(sm_tb_stat$cancer)), ncol=length(unique(names(sm_tb_stat_coll))))
colnames(fishp) = names(sm_tb_stat_coll)
rownames(fishp) = unique(sm_tb_stat$cancer)
fishp = as.data.frame(fishp)
fishod = fishp
for (i in names(sm_tb_stat_coll)){
  fjn = full_join(sm_tb_stat_coll[[i]], smm_cancer, by=c('cancer'='Cancer'))
  for (j in 1:nrow(fjn)){
    if (i=='pos'){
      rvg = as.numeric(fjn$Positive[j])
    } else {
      rvg = as.numeric(fjn$Negative[j])
    }
    if (!is.na(fjn$hit[j])){
      mn = matrix(c(fjn$hit[j], fjn$not_detected[j]+fjn$wrong_hit[j], 
                    rvg-fjn$hit[j], 
                    as.numeric(fjn$All[j])-rvg-fjn$not_detected[j]-fjn$wrong_hit[j]),nrow=2)
      pv = fisher.test(mn)$p.value
      od = fisher.test(mn)$estimate
      fishp[fjn$cancer[j], i] = toSig(pv)
      fishod[fjn$cancer[j], i] = od
    } else {
      fishp[fjn$cancer[j], i] = toSig(NA)
      fishod[fjn$cancer[j], i] = NA
    }
  }
}

fishp[fishp==''] = 'na'
fishp$cancer = rownames(fishp)
fishp = fishp[fishp$cancer!='OV',]

ctr1_count = ctr1_count + geom_text(data=fishp, aes(x=cancer,y=28,label=ctr1), 
                                    color='white', fontface='bold', size=5, vjust=0)
ctr2_count = ctr2_count + geom_text(data=fishp, aes(x=cancer,y=12,label=ctr2), 
                                    color='white', fontface='bold',size=5, vjust=0)
pos_count = pos_count + geom_text(data=fishp, aes(x=cancer,y=33.7,label=pos), 
                                    color='white', fontface='bold',size=5, vjust=0)


marrangeGrob(grobs=list(rank_pfl, ctr1_count + theme(legend.position = 'none'), 
                        ctr2_count, pos_count),ncol=4,nrow=1, top='') %>% 
  ggsave(file.path(plot_out, 'TCGA-pan_VS-mutneg_ult-wFDR_hitCount.pdf'),
         plot=., width=11.69*1.5,height=8.27*0.5,units='in',device='pdf',dpi=300)
ggsave(file.path(plot_out, 'TCGA-pan_VS-mutneg_ult-wFDR_hitCount_legend.pdf'),
       plot=ctr1_count, width=11.69*0.7,height=8.27*0.5,units='in',device='pdf',dpi=300)
save(sm_tb_stat_coll, file = file.path(data_out, 'control_summary.RData'))

# marrangeGrob(grobs=list(rank_pfl, ctr1_count, pos_count, ctr2_count),ncol=2,nrow=2, top='') %>% 
#   ggsave(file.path(plot_out, 'TCGA-pan_VS-mutneg_ult-noFDR_hitCount.pdf'),
#          plot=., width=11.69*0.7,height=8.27*0.9,units='in',device='pdf',dpi=300)

# Annotate mutant passing QC or not
ctr_long$mut_ann = FALSE
for (i in 1:nrow(ctr_long)){
  if (ctr_long$mutation[i] %in% detected_exp[[ctr_long$experiment[i]]]){
    ctr_long$mut_ann[i] = TRUE
  }
}

ctr_long_breast = subset(ctr_long, ctr_long$cancer=='BRCA' & ctr_long$mutation %in% 
                           c('R175H', 'Walerych et al., 2016', 'Hotspots')) # for later
# WT control or Not !!! Choice
ctr_long = ctr_long[ctr_long$mutation != 'Walerych et al., 2016',]
ctr_long_plt_3 = subset(ctr_long, !ctr_long$Gene_annotation %in% c('wt_control', 'wt_control_2'))
ctr_long_plt_1 = subset(ctr_long, ctr_long$Gene_annotation=='wt_control')
ctr_long_plt_2 = subset(ctr_long, ctr_long$Gene_annotation=='wt_control_2')

# Plotting ====
plt.list = list()
myPalette = colorRampPalette(c("royalblue", '#D62728FF'))
up_p = ceiling(-log10(min(ctr_long$FDR, na.rm = T)))
print(up_p)
sc = scale_colour_gradientn(colours = myPalette(50), 
                            limits=c(0, up_p),breaks = c(0,8,16,32), trans = 'log1p')
ssp = scale_size_continuous(limits = c(0,ceiling(max(abs(ctr_long$beta), na.rm=T))))
count = 1
for (ctr_long_plt in list(ctr_long_plt_1, ctr_long_plt_2, ctr_long_plt_3)){
  for (i in unique(ctr_long_plt$experiment)){
    ctr_long_sub = subset(ctr_long_plt, ctr_long_plt$experiment==i)
    test = as.character(unique(ctr_long_sub$mutation[ctr_long_sub$mut_ann==FALSE]))
    nm_exp = toupper(strsplit(i, split = '_')[[1]][2])
    g = ggplot(ctr_long_sub) +
      geom_point(aes(x=Gene, y=mutation, size=abs(beta), 
                     color=-log10(FDR), shape=wrong), stroke=2) +
      scale_shape_manual(name='', label=c('Positive', 'Negative'), limits = c(FALSE, TRUE), values = c(1,2)) +  # 16 17
      # facet_rep_wrap(~experiment, ncol = 2, repeat.tick.labels = TRUE) +
      # scale_colour_gradient2(low = "blue", high = "purple", mid = "pink") +
      sc + ssp +
      mytme +
      geom_point(data = subset(ctr_long_sub, FDR < 0.05), 
                 aes(x=Gene, y=mutation), shape = '*', size=6, color='black') +
      # y = 'Mutation (group)'
      labs(x='', y='', title = nm_exp) +
      theme(strip.background = element_rect(fill=NA), 
            strip.text = element_text(face='bold', color='black', size=12),
            axis.text.x = element_text(hjust = 1, vjust=0.5, size=12, angle = 90, face='italic'),
            axis.text.y = element_text(hjust = 1, size=12),
            axis.title.y = element_text(size=18),
            panel.border = element_rect(color='black', size=1, fill=NA),
            panel.grid.major = element_line(color='grey',linetype='dotted'),
            legend.position = 'none'
      )
    if (length(test) > 0){
      wid = length(unique(ctr_long_sub$Gene)) / 2
      for (nd_gene in test){
        g = g + geom_hline(yintercept = nd_gene, color='black') +
          geom_label(x = wid, y=nd_gene, label.size = 0, label.r = unit(0,'lines'),
                     fill='white', label='Not passing VAF filter')
      }
    }
    plt.list[[paste(count, i, sep='_')]] = g
  }
  count = count + 1
}

marrangeGrob(grobs=plt.list,ncol=3,nrow=9, top = '') %>% 
  ggsave(file.path(plot_out, 'TCGA-pan_VS-mutneg_ult-noFDR_posCtr.pdf'),
       plot=., height=11.69*3.2, width=8.27*3.5,units='in',device='pdf',dpi=300)

ggsave(file.path(plot_out, 'TCGA-pan_VS-mutneg_ult-noFDR_legend.pdf'),
       plot=plt.list[[1]] + theme(legend.position = 'bottom', legend.direction = 'horizontal'), 
       width=11.69*2, height=8.27,units='in',device='pdf',dpi=300)


### Plot individual panel for breast cancer Fig 1B ====
plt.list = list()
myPalette = colorRampPalette(c("#9467BDFF", '#D62728FF'))
up_p = ceiling(-log10(min(ctr_long_breast$FDR, na.rm = T)))
print(up_p)
sc = scale_colour_gradientn(colours = myPalette(50), 
                            limits=c(0, up_p),breaks = c(0,2,4,8,16), trans = 'log1p')
sc = scale_colour_gradientn(colours = myPalette(50), 
                            limits=c(0, up_p),breaks = c(0,8,16))
ssp = scale_size_continuous(limits = c(0,ceiling(max(abs(ctr_long_breast$beta), na.rm=T))))
ctr_long_plt_3 = subset(ctr_long_breast, !ctr_long_breast$Gene_annotation %in% c('wt_control', 'wt_control_2'))
ctr_long_plt_1 = subset(ctr_long_breast, ctr_long_breast$Gene_annotation=='wt_control')
ctr_long_plt_2 = subset(ctr_long_breast, ctr_long_breast$Gene_annotation=='wt_control_2')
df_list = list(ctr_long_plt_1, ctr_long_plt_2, ctr_long_plt_3)
tle = c('BRCA downregulated controls 1', 'BRCA downregulated controls 2', 'BRCA upregulated controls')
for (j in 1:length(df_list)){
  ctr_long_plt = df_list[[j]]
  for (i in unique(ctr_long_plt$experiment)){
    ctr_long_sub = subset(ctr_long_plt, ctr_long_plt$experiment==i)
    ctr_long_sub$beta[ctr_long_sub$FDR>=0.05] = NA
    ctr_long_sub$FDR[ctr_long_sub$FDR>=0.05] = NA
    test = as.character(unique(ctr_long_sub$mutation[ctr_long_sub$mut_ann==FALSE]))
    g = ggplot(ctr_long_sub) +
      geom_point(aes(x=Gene, y=mutation, size=abs(beta), 
                     color=-log10(FDR), shape=wrong), stroke=2) +
      scale_shape_manual(name='', label=c('Positive', 'Negative'), limits = c(FALSE, TRUE), values = c(16,17)) +  # 16 17 or 1 2
      sc + ssp +
      mytme +
      # geom_point(data = subset(ctr_long_sub, FDR < 0.05), 
      #            aes(x=Gene, y=mutation), shape = '*', size=6, color='black') +
      labs(x='', y='', title = tle[j]) +
      scale_y_discrete(position = 'right') +
      theme(strip.background = element_rect(fill=NA), 
            strip.text = element_text(face='bold', color='black', size=12),
            axis.text.x = element_text(hjust = 1, vjust=0.5, size=12, angle = 90, face='italic'),
            axis.text.y = element_text(hjust = 1, size=12),
            axis.title.y = element_text(size=18),
            panel.border = element_rect(color='black', size=1, fill=NA),
            panel.grid.major = element_line(color='grey',linetype='dotted'),
            legend.position = 'none'
      )
    if (length(test) > 0){
      wid = length(unique(ctr_long_sub$Gene)) / 2
      for (nd_gene in test){
        g = g + geom_hline(yintercept = nd_gene, color='black') +
          geom_label(x = wid, y=nd_gene, label.size = 0, label.r = unit(0,'lines'),
                     fill='white', label='Not passing VAF filter')
      }
    }
    plt.list[[j]] = g
  }
}
marrangeGrob(grobs=plt.list,ncol=1,nrow=3, top = '') %>% 
  ggsave(file.path(plot_out, 'TCGA-pan_VS-mutneg_ult-noFDR_Ctr_BRCA.pdf'),
         plot=., height=11.69*0.8, width=8.27*1.1,units='in',device='pdf',dpi=300)

ggsave(file.path(plot_out, 'TCGA-pan_VS-mutneg_ult-noFDR_Ctr_BRCA_legend.pdf'),
       plot=plt.list[[3]] + theme(legend.position = 'bottom', legend.direction = 'horizontal'), 
       height=8, width=12,units='in',device='pdf',dpi=300)

### Conflicting result ===
tcga = load_clean_data('/Users/jefft/Desktop/p53_project/datasets/9-BRCA-TCGA/clean_data.RData',
                       check_cna_wt = T)

load('/Users/jefft/Desktop/p53_project/datasets/CCLE_22Q1/pcd/BRCA/clean_data.RData')
ccle = dt
remove(dt)
gc()
genes = c('BAG1', 'IGF1R', 'IGF2')

# mutation_groups = list(c(273,248,175,245,249,282))
mutation_groups = read.table('/Users/jefft/Desktop/p53_project/datasets/meta_muts/hot_spot.txt', header = T)
mutation_groups = list(mutation_groups$aa_pos)
names(mutation_groups) = c('Hotspots')

plt = get_genes_plt(genes=genes, ccle=ccle, tcga=tcga, mutation_groups, 
                    primary_site = 'Breast', rnai = NULL, 
                    comparison = list('rna'=list(c('Hotspots', 'Wildtype'))),
                    no_ccle = TRUE, plot_n = T, plot_nonsense = F)

for (i in 1:length(plt$plots)){
  plt$plots[[i]] = plt$plots[[i]] + labs(x='',y=toupper(strsplit(names(plt$plots)[i], split='_')[[1]][1]))
}

plt$plots %>% marrangeGrob(ncol=4, nrow=1, top = '',
                           layout_matrix = matrix(1:4,byrow = T, ncol=4)) %>%
  ggsave(file.path(plot_out, 'BRCA_conflict_controls.pdf'), bg = 'transparent',
         plot=., width=8.27*1.3,height=11.69*0.3,units='in',device='pdf',dpi=300)




### Contact VS Conform in Esposito, 2022 ====
source('/Users/jefft/Desktop/p53_project/scripts/ccle_utils.R')

# contact_cell = intersect(ccle[[1]]@colData$CELL_LINE_NAME, 
#                          c('MDA-MB-468','U373MG', 'U-251 MG', 'SF-295', 'HCC193', 'PC9'))
# conform_cell = intersect(ccle[[1]]@colData$CELL_LINE_NAME, 
#                          c('HCC1395', 'HCC1954', 'SK-MEL-2', 'SK-LMS-1'))
# genes = c('TEAD1', 'TEAD2', 'TEAD3', 'TEAD4')
# 
# comp = data.frame('name'=c(contact_cell, conform_cell), 'group'=c(rep('contact', 3), rep('conform', 3)))
# comp['ID'] = sapply(comp$name, function(x){ccle[[1]]@colData$PATIENT_ID[which(ccle[[1]]@colData$NAME==x)]})
# comp = cbind(comp, t(assay(ccle[[1]][genes,comp$ID,'RNA'])))
# 
# comp = gather(comp, key='gene', value='exp', 4:7)
# ggplot(comp, aes(x=group,y=exp)) +
#   geom_point() +
#   facet_wrap(~gene) +
#   geom_text(aes(label=name)) +
#   mytme

mutation_groups = list(c(273,248), c(175,245,249,282))
names(mutation_groups) = c('HS cont.', 'HS conf.')

genes = c('HMGCR', 'MVK', 'MVD', 'FDPS', 'SQLE', 'LSS', 'DHCR7')
genes = c('TEAD1', 'TEAD2', 'TEAD3', 'TEAD4')
# n=30 n=28 n=40 n=651; HS conf., HS cont., nonsense, wildtype
plt = get_genes_plt(genes=genes, ccle=ccle, tcga=tcga, mutation_groups, 
              primary_site = 'Breast', rnai = NULL, 
              comparison = list('rna'=list(c('HS cont.', 'Wildtype'),
                                           c('HS conf.', 'HS cont.'),
                                           c('HS conf.', 'Wildtype'))),
              no_ccle = TRUE, plot_n = TRUE, plot_nonsense = F,
              no_anova = T)

for (i in 1:length(plt$plots)){
  plt$plots[[i]] = plt$plots[[i]] + labs(x='',y=toupper(strsplit(names(plt$plots)[i], split='_')[[1]][1]))
}

plt$plots %>% marrangeGrob(ncol=4, nrow=2, top = '',
                           layout_matrix = matrix(1:8,byrow = T, ncol=4)) %>%
  ggsave(file.path(plot_out, 'MVA_controls.pdf'), bg = 'transparent',
         plot=., width=8.27*1.3,height=11.69*0.6,units='in',device='pdf',dpi=300)

plt$plots %>% marrangeGrob(ncol=4, nrow=1, top = '',
                     layout_matrix = matrix(1:4,byrow = T, ncol=4)) %>%
  ggsave(file.path(plot_out, 'TEAD_controls.pdf'), bg = 'transparent',
         plot=., width=8.27*1.3,height=11.69*0.3,units='in',device='pdf',dpi=300)

