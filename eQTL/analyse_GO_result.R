setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
exp_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-mutneg_ult'
plot_out = file.path(exp_home, 'plots')
source('utils.R')
source('/Users/jefft/Desktop/Manuscript/set_theme.R')
library(tidyverse)
library(stringr)

load(file.path(exp_home, 'GO_result_no_BG_filter.RData'))
go_df = data.frame()
for (i in 1:length(go_coll)){
  nm = names(go_coll)[i]
  sub_df = go_coll[[nm]]
  sub_df[['experiment']] = nm
  rownames(sub_df) = NULL
  go_df = rbind(go_df, sub_df)
}
go_df['cancer'] = sapply(go_df$experiment, 
              function(x){toupper(strsplit(x, split = '_')[[1]][2])})
go_df['sign'] = 'neg'
go_df[grep('pos', go_df$experiment), 'sign'] = 'pos'
go_df['mutation'] = sapply(go_df$experiment,
              function(x){strsplit(x, split='-')[[1]][2]})
go_df$qvalue[which(is.na(go_df$qvalue))] = go_df$p.adjust[which(is.na(go_df$qvalue))]

# QC gene ratio > 1%
go_df['gene_ratio_val'] = sapply(go_df$GeneRatio, function(x){eval(parse(text = x))})
hist(go_df$gene_ratio_val[grep('GO', go_df$ID)], breaks=100)
summary(go_df$gene_ratio_val[grep('GO', go_df$ID)])
go_df = go_df[which(go_df$gene_ratio_val > 0.01),]

# take intersect
library(ggvenn)
hot_spot_go = go_df[grep('hot_spot', go_df$experiment),]
hot_spot_go = hot_spot_go[grep('GO', hot_spot_go$ID),]
coll = list()
for (i in unique(hot_spot_go$cancer)){
  coll[[i]] = hot_spot_go$ID[which(hot_spot_go$cancer == i)]
}
g = ggvenn(coll, stroke_color = 'white')
g %>% ggsave(file.path(plot_out, 'pos_pan_GO_overlap.pdf'),
             plot=., width=5,height=5,units='in',device='pdf',dpi=300)
coll = coll[which(names(coll)!='STAD')]
co_go_four_cancer = Reduce(intersect, coll) # "GO:0006913" "GO:0051169"
test = subset(hot_spot_go, ID %in% co_go_four_cancer)
test %>% group_by(ID) %>% summarize(mn = mean(gene_ratio_val), sd=sd(gene_ratio_val))
co_genes = unlist(lapply(test$geneID, function(x){return(strsplit(x, split='/')[[1]])}))
intersect(co_genes, c('MDM2', ''))

### Conformational VS Contact mutant
ids = list()
for (expr in unique(hot_spot_go$experiment)){
  ids[[expr]] = subset(hot_spot_go, experiment == expr)$ID
}
g = ggvenn(ids[intersect(grep('pos', names(ids)),
                     grep('brca', names(ids)))]) +
  labs(title = 'BRCA')
g %>% ggsave(file.path('/Users/jefft/Desktop/p53_project/Plots/conformVScontact', 'BRCA_pos_GO.pdf'),
         plot=., width=5,height=5,units='in',device='pdf',dpi=300)
g = ggvenn(ids[intersect(grep('pos', names(ids)),
                     grep('coad', names(ids)))]) +
  labs(title = 'COAD')
g %>% ggsave(file.path('/Users/jefft/Desktop/p53_project/Plots/conformVScontact', 'COAD_pos_GO.pdf'),
             plot=., width=5,height=5,units='in',device='pdf',dpi=300)
g = ggvenn(ids[intersect(grep('pos', names(ids)),
                     grep('lgg', names(ids)))]) +
  labs(title = 'LGG')
g %>% ggsave(file.path('/Users/jefft/Desktop/p53_project/Plots/conformVScontact', 'LGG_pos_GO.pdf'),
             plot=., width=5,height=5,units='in',device='pdf',dpi=300)

#### dot plot of each top 5 terms
signn = 'pos'
top_go_plt = data.frame()
for (i in c('BLCA', 'COAD', 'LGG', 'BRCA')){
  for (j in c('hot_spot')){ # 'core', 'contact'
    sub_go_df = subset(go_df, (go_df$cancer == i) & (go_df$mutation == j) & (go_df$sign == signn))
    sub_go_df = sub_go_df[grep('GO', sub_go_df$ID),]
    if (nrow(sub_go_df)>0){
      top_go_plt = rbind(top_go_plt, sub_go_df[order(sub_go_df$p.adjust)[1:(min(5,nrow(sub_go_df)))],])
    }
  }
}

# top_go_plt$Description = factor(top_go_plt$Description, levels = unique(top_go_plt$Description))
g = ggplot(top_go_plt, aes(x=mutation)) +
  geom_point(aes(y=Description, color=-log10(p.adjust), size=gene_ratio_val, group=cancer)) +
  facet_wrap(~cancer, scale='free', nrow = 2) +
  scale_colour_gradient2(low = "blue", high = "purple", mid = "pink") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=25)) +
  theme(axis.text.x = element_text(angle=0, vjust = 1, size=16),
        strip.background = element_rect(fill = 'transparent'),
        strip.text = element_text(size=16, face='bold')) +
  mytme
g %>% ggsave(file.path(plot_out, paste(signn, '_pan_GO_dotPlot.pdf', sep='')),
             plot=., width=10,height=12,units='in',device='pdf',dpi=300)

### Get all genes positively regulated by conformational/contact changes
output_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-mutneg/outputs'
coll = list()
for (i in c('tcga_brca_raw_seq', 'tcga_coad_raw_seq', 'tcga_lgg_raw_seq')){
  fp = file.path(output_home, i, 'trans_eqtl_fdr005.txt')
  dtb = read.table(fp, sep='\t', header = T)
  dtb['cancer'] = toupper(strsplit(i, split='_')[[1]][2])
  coll[[i]] = dtb
}
coll = Reduce(rbind, coll)
coll = subset(coll, coll$protein_change %in% c('hot_spot_conform', 'hot_spot_contact'))
coll = subset(coll, coll$beta < 0)
gene_conform = list()
gene_contact = list()
for (i in c('BRCA', 'COAD', 'LGG')){
  dtb = subset(coll, coll$cancer ==i)
  gene_conform[[i]] = dtb$gene[which(dtb$protein_change == 'hot_spot_conform')]
  gene_contact[[i]] = dtb$gene[which(dtb$protein_change == 'hot_spot_contact')]
}
ggvenn(gene_conform) %>% 
  ggsave(file.path('/Users/jefft/Desktop/p53_project/Plots/conformVScontact', 'gene_neg_overlap_conform_BRCA_COAD_LGG.pdf'),
        plot=., width=5,height=5,units='in',device='pdf',dpi=300)
ggvenn(gene_contact) %>% 
  ggsave(file.path('/Users/jefft/Desktop/p53_project/Plots/conformVScontact', 'gene_neg_overlap_contact_BRCA_COAD_LGG.pdf'),
         plot=., width=5,height=5,units='in',device='pdf',dpi=300)

# There is < 5 overlappiong genes...
gene_conform = Reduce(intersect, gene_conform)
gene_contact = Reduce(intersect, gene_contact)

### G



