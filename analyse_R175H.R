library(tidyverse)
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-mutneg_ult'
eqtl_out = file.path(dir_home, 'outputs')
plot_out = file.path(dir_home, 'plots', 'R175H')
data_out = file.path(dir_home, 'data_out')
source('utils.R')
source('enrich_utils.R')
source('../ccle_utils.R')
source('/Users/jefft/Desktop/Manuscript/set_theme.R')
library(tidyverse)
library(stringr)
library(ggvenn)

# load GO
gl = load_GO_out(exp_home = dir_home, gene_ratio_min = 0)
go_df = gl$df
go_df = go_df[grep('GO', go_df$ID),]
go_df = go_df[grep('R175H', go_df$mutation),]
go_df_pos = subset(go_df, go_df$sign == 'pos')


# load gene
beta_cutoff=0
coll = load_eQTL_output(eqtl_out, beta=beta_cutoff, exclude = 'tcga_nine_pool')
df_coll = data.frame()
for (i in names(coll)){
  nm = toupper(strsplit(i, split='_')[[1]][2])
  df = coll[[i]]
  if (class(df)!='logical'){
    df = subset(df, abs(df$beta) > beta_cutoff)
    df_coll = rbind(df_coll, df)
  }
}
df_coll = subset(df_coll, df_coll$protein_change %in% c('p.R175H'))
pos_gene = subset(df_coll, df_coll$beta > 0)

# GO overlap
coll = list()
for (i in unique(go_df_pos$experiment)){
  coll[[toupper(strsplit(i, split='_')[[1]][2])]] = go_df_pos$ID[go_df_pos$experiment==i]
}
g = ggvenn(coll, stroke_color = 'white')
g %>% ggsave(file.path(plot_out, paste('R175H_BRCA-LGG-COAD_pos_co_GO.pdf', sep='')),
             plot=., width=4,height=4,units='in',device='pdf',dpi=300)

co_term = Reduce(intersect, coll)
top_go_plt = go_df_pos[which(go_df_pos$ID %in% co_term),]

g = ggplot(top_go_plt, aes(x=cancer)) +
  geom_point(aes(y=Description, color=-log10(p.adjust), size=Count, group=cancer)) +
  #facet_wrap(~cancer, scale='free', nrow = 2) +
  scale_colour_gradient2(low = "blue", high = "purple", mid = "pink") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=25)) +
  theme(axis.text.x = element_text(angle=0, vjust = 1, size=16),
        strip.background = element_rect(fill = 'transparent'),
        strip.text = element_text(size=16, face='bold')) +
  mytme
g
g %>% ggsave(file.path(plot_out, paste('R175H_BRCA-LGG-COAD_pos_co_GO_dotPlot.pdf', sep='')),
             plot=., width=10,height=5,units='in',device='pdf',dpi=300)

# Gene overlap
df_coll$cancer = toupper(sapply(df_coll$experiment, function(x){return(strsplit(x, split='_')[[1]][2])}))
df_coll_pos = subset(df_coll, df_coll$beta > 0 & df_coll$cancer != 'STAD') # too few in STAD
gene_count = table(df_coll_pos$gene)
print(table(gene_count)) # only 4 overlapping by three cancer
check = names(gene_count)[which(gene_count >= 2)] # at least overlap once.
names(gene_count)[which(gene_count >= 3)]
test = do_GO(check, background = unique(df_coll_pos$gene), ont='MF')
dotplot(test)
pcs_GO_out(test, dir = plot_out, filename = 'R175H_BRCA-LGG-COAD_pos_co_gene_GO_MF.pdf')

# Cancer specific signature
ont = 'BP'
coll = list()
go_list = list()
for (i in unique(df_coll_pos$cancer)){
  gs = subset(df_coll_pos, df_coll_pos$cancer==i)$gene
  coll[[i]] = gs
  test = do_GO(gs[which(!gs %in% check)], background = unique(df_coll_pos$gene), ont = ont)
  go_list[[i]] = pcs_GO_out(test, dir = plot_out, filename = paste('R175H_', i, '_', ont, '_specific.pdf', sep=''))
}
g = ggvenn(coll, stroke_color = 'white')
g %>% ggsave(file.path(plot_out, paste('R175H_BRCA-LGG-COAD_pos_co_gene.pdf', sep='')),
             plot=., width=4,height=4,units='in',device='pdf',dpi=300)
Reduce(intersect, coll)

# breast cancer and ion channel genes
ont = 'BP'
gs = subset(df_coll_pos, df_coll_pos$cancer=='BRCA')$gene
test = do_GO(gs[which(!gs %in% check)], background = unique(df_coll_pos$gene), ont = ont)
rs = test@result
rs = subset(rs, rs$p.adjust < 0.05)
tmport = strsplit(rs$geneID[1], split='/')[[1]]
print(length(grep('SLC', tmport)))
print(length(tmport))
print(length(grep('SLC', tmport)) / length(tmport))
## check expression in R175H, other mutant and WT
gene = tmport[grep('SLC', tmport)]
tcga = load_clean_data(fp = '/Users/jefft/Desktop/p53_project/datasets/9-BRCA-TCGA/clean_data.RData', 
                       ann_bin_mut_list = c('p.R175H'), mode='tcga')
table(tcga[[1]]@colData$has_p.R175H)
exp_SLC = assay(tcga[[1]][gene,,'RNA'])
exp_SLC = t(exp_SLC)
exp_SLC = cbind(exp_SLC, tcga[[1]]@colData[,c('p53_state', 'has_p.R175H')])
exp_SLC = as.data.frame(exp_SLC)
exp_SLC$gp = paste(exp_SLC$p53_state, exp_SLC$has_p.R175H, sep='-')
exp_SLC = subset(exp_SLC, !exp_SLC$p53_state %in% c('others', 'frameshift'))
table(exp_SLC$gp)
exp_SLC_plt = gather(exp_SLC, key='gene', value='exp', 1:19)
ggplot(exp_SLC_plt) +
  geom_tile(aes(x=gene, y=gp, fill=exp))

g = ggplot(exp_SLC_plt, aes(x=gp, y=exp)) +
  geom_boxplot() +
  facet_wrap(~gene, ncol=3) +
  stat_compare_means(comparisons = list(c('missense-0', 'missense-1'),
                                        c('missense-1', 'Wildtype-0'),
                                        c('missense-1', 'nonsense-0')),
                     method = 't.test') +
  theme(axis.text.x = element_text(angle=45)) +
  mytme
g %>% ggsave(file.path(plot_out, paste('BRCA-R175H-SLCgenes.pdf', sep='')),
             plot=., width=10,height=25,units='in',device='pdf',dpi=300)

exp_SLC = exp_SLC[order(exp_SLC$gp),]
ann_row = data.frame('gp'=exp_SLC[,22])
rownames(ann_row) = rownames(exp_SLC)
pheatmap(exp_SLC[,1:19], labels_row = F, labels_col = F, annotation_row = ann_row, cluster_rows = F)

