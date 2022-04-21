library(tidyverse)
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt'
eqtl_out = file.path(dir_home, 'outputs')
plot_out = file.path(dir_home, 'plots', 'COAD')
data_out = file.path(dir_home, 'data_out')
source('utils.R')
source('enrich_utils.R')
source('../ccle_utils.R')
source('../set_theme.R')
source('../overlap_utils.R')
source('../revigo_utils.R')
library(tidyverse)
library(stringr)
library(ggvenn)


# load gene
load('/Users/jefft/Desktop/p53_project/datasets/TCGA-Pan-Nine/gene_matrix.RData')
mtx = as.data.frame(mtx)
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
df_coll = df_coll[grep('p\\.',df_coll$protein_change),]
df_coll = df_coll[df_coll$experiment=='tcga_coad_raw_seq',]

# normalize number with sample size ====
ssize = read.table('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt/sample_size.txt', header = T)
ssize = ssize[ssize$Cancer=='COAD',]
ssize = ssize[grep('p\\.',ssize$Mutation),]
posnm = c()
negnm = c()
for (i in 1:nrow(ssize)){
  posnm = c(posnm, sum(df_coll$protein_change==ssize$Mutation[i] & df_coll$beta>0))
  negnm = c(negnm, sum(df_coll$protein_change==ssize$Mutation[i] & df_coll$beta<0))
}
ssize$posnm = posnm
ssize$negnm = negnm
ssize$Mutation = gsub('p\\.','',ssize$Mutation)
g = ggplot(ssize %>% gather(key='sign', value='Count', 6:7),
           aes(x=Mutation,y=Count/Mutant)) +
  geom_bar(aes(fill=sign),stat='identity', position = 'dodge') +
  scale_y_continuous(expand = expansion(mult = c(0,0.2))) +
  geom_text(aes(label=paste(Count, Mutant, sep='/'), group=sign), 
            position = position_dodge(width=0.9), vjust=-0.5, size=3) +
  labs(x='', y='Gene count / Mutant sample size') + scale_fill_d3(name='',labels = c('Negative', 'Positive')) +
  # geom_text(label=parse(text="'Not passing\nVAF filter'"), x='LUAD', y=2, size=2.5) +
  mytme + theme(legend.position = c(0.8,0.9))
ggsave(file.path(plot_out, 'count_overview.pdf'),
       plot=g, height=11.69*0.4,width=8.27*0.7,units='in',device='pdf',dpi=300)


# Gene overlap ====
df_coll$protein_change = gsub('p\\.','',df_coll$protein_change)
## !!! change sign !!!
df_coll_pos = subset(df_coll, df_coll$beta > 0 & df_coll$cancer != 'STAD') # too few in STAD
df_coll_pos = df_coll_pos[order(df_coll_pos$beta, decreasing = T),]
gene_coll = list()
for (i in unique(df_coll_pos$protein_change)){
  gene_coll[[i]] = df_coll_pos$gene[df_coll_pos$protein_change==i][1:min(200,sum(df_coll_pos$protein_change==i))]
}
# ggvenn(gene_coll)
um = gen_upSet_mtx(gene_coll)
# upset(um, intersect = colnames(um))

## run overlapping test
ovl_p_coll = list()
for (i in 2:5){
  cb = combn(1:5,i)
  for (j in 1:ncol(cb)){
    id = paste(colnames(um)[cb[,j]], collapse = ' X ')
    tcd = rep(NA,5)
    tcd[cb[,j]] = T
    ovl_p = run_overlap_test(um, tcd, n_loop = 1000)
    ovl_p_coll[[id]] = ovl_p
  }
}

uquery = list()
count = 1
for (j in names(ovl_p_coll)){
  if (ovl_p_coll[[j]]$p.left >= 0.95){
    ist = strsplit(j, split=' X ')[[1]]
    uquery[[count]] = upset_query(intersect = ist, color='#FF7F0EFF', fill='#FF7F0EFF')
    count = count + 1
  }
}
print(count-1)

g = upset(um, intersect = colnames(um), min_size=0,
          sort_intersections_by = 'degree', sort_intersections = 'ascending',
          mode = 'inclusive_intersection', width_ratio = 0.1, 
          min_degree=2, set_sizes = FALSE,
          queries = uquery,
          base_annotations=list(
            'Intersection size'=intersection_size(
              mode = 'inclusive_intersection',
              text=list(vjust=-0.2,hjust=0.5,angle=0,size=3),
              #mapping=aes(fill=mpaa),
              text_colors=c(on_background='black', on_bar='white'))
            + annotate(
              geom='text', x=Inf, y=Inf, color='#FF7F0EFF',
              label='Significant overlap', size=5,
              vjust=1, hjust=1
            ))
          # set_sizes = upset_set_size() + 
          #   scale_y_continuous(breaks=c(4000,2000,0))
) & theme(plot.background=element_rect(fill='transparent', color=NA))
g
ggsave(file.path(plot_out, 'panMut_pos_upset_GO.pdf'),
       plot=g, height=11.69*0.3,width=8.27*1,units='in',device='pdf',dpi=300)


### Get top five GO for each significant overlap ====
top5terms = list()
ids = c()
count = 0
for (i in uquery){
  its = i$intersect
  ids = c(ids, its)
  id = paste(its, collapse =' X ')
  genes = rownames(um)[which(rowSums(um[,its]) == length(its))]
  print(length(genes))
  if (length(genes)>=20){
    count = count + 1
    print(id)
    test = do_GO(genes, background = rownames(um))
    ps = pcs_GO_out(test, filename = paste('panMut_posOvl_GO',id,'.pdf',sep=''), 
                    dir = plot_out, cneplt = F,
                    size = c(8.27*0.7,11.69*0.5))
    
    rs = test@result
    rs = rs[rs$p.adjust < 0.05,]
    if (nrow(rs)>0){
      top5terms[[id]] = test@result$Description[1:min(5,nrow(test@result))]
    } else {
      top5terms[[id]] = NULL
    }
  } else {
    print('Fewer than 20')
    print(id)
  }
}
top5term_df = as.data.frame(sort(table(unlist(top5terms)), decreasing = F))
top5term_df$Var1 = factor(top5term_df$Var1, levels = top5term_df$Var1)
g = ggplot(top5term_df[(nrow(top5term_df)-5):nrow(top5term_df),]) + scale_x_continuous(expand = c(0,0)) +
  geom_bar(aes(y=Var1,x=Freq/count), stat = 'identity') + mytme + 
  labs(x='Frequency', y='') +
  geom_text(aes(label=Freq, y=Var1, x=Freq/count), nudge_x = -0.1, color='white') +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0), breaks=seq(0,1,0.5))
ggsave(file.path(plot_out, 'panMut_pos_upset_GO.pdf'),
       plot=g, height=11.69*0.3,width=8.27*0.5,units='in',device='pdf',dpi=300)

rownames(um)[get_idx_condt(um,c(T,T,T,T,T))]

# Cancer specific signature
ont = 'BP'
coll = list()
go_list = list()
for (i in unique(df_coll_pos$cancer)){
  gs = subset(df_coll_pos, df_coll_pos$cancer==i)$gene
  coll[[i]] = gs
  test = do_GO(gs[which(!gs %in% check)], background = unique(df_coll_pos$gene), ont = ont)
  rvg = run_revigo(test@result[,c('ID', 'pvalue')])
  test@result = subset(test@result, ID %in% rvg$`Term ID`)
  go_list[[i]] = pcs_GO_out(test, dir = plot_out, filename = paste('R273H_', i, '_', ont, '_specific.pdf', sep=''))
}
g = ggvenn(coll, stroke_color = 'white')
g %>% ggsave(file.path(plot_out, paste('R273H_BRCA-COAD_pos_co_gene.pdf', sep='')),
             plot=., width=4,height=4,units='in',device='pdf',dpi=300)

