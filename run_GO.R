library(tidyverse)
library(ggpubr)
library(clusterProfiler)
setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
source('utils.R')
source('enrich_utils.R')
source('/Users/jefft/Desktop/p53_project/scripts/ccle_utils.R')
source('/Users/jefft/Desktop/Manuscript/set_theme.R')
# TCGA-pan_VS-mutneg_ult TCGA-pan_VS-wt
dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt'
eqtl_out = file.path(dir_home, 'outputs')
plot_out = file.path(dir_home, 'plots', 'GO') # GO or GO_MF
data_out = file.path(dir_home, 'data_out')

flt_df = function(go_df, cutoff=0.01){
  bg = sapply(go_df$GeneRatio, function(x){eval(parse(text = x))})
  go_df = go_df[which(bg > cutoff),]
  return(go_df)
}

## !!! Config this!
mode = 'fdr' # fdr: collect FDR-filtered data
ont = 'BP'
go_rs_nm = 'GO_BP_result_no_BG_filter.RData'

## eQTL mapping outputs
### Experiment: contrast samples with p53 mutation to those without.
# does not filter beta here
beta_cutoff = 0 # 0
coll = load_eQTL_output(eqtl_out, beta=beta_cutoff, exclude = 'tcga_nine_pool', mode = mode)
study = c()
map_mut = c()
num_mut = c()
total_gene = c()
df_coll = data.frame()
for (i in names(coll)){
  nm = toupper(strsplit(i, split='_')[[1]][2])
  study = c(study, nm)
  df = coll[[i]]
  if (class(df)=='logical'){
    map_mut = c(map_mut, NA)
    num_mut = c(num_mut, NA)
    total_gene = c(total_gene, NA)
  } else {
    df = subset(df, abs(df$beta) > beta_cutoff)
    df_coll = rbind(df_coll, df)
    map_mut = c(map_mut, paste(unique(df$protein_change), collapse = '|'))
    num_mut = c(num_mut, length(unique(df$protein_change)))
    total_gene = c(total_gene, nrow(df))
  }
}
sry = data.frame('Study'=study, 'Mapped_mutant'=map_mut, 'Num_mutant'=num_mut, 
                 'Num_gene' = total_gene)
sry[['Avg_gene']] = sry$Num_gene / sry$Num_mutant
print(sry)


sry = sry[which(!is.na(sry$Study)),]
sry = sry[which(sry$Study != 'NINE'),]
df_coll = df_coll[grep('tcga', df_coll$experiment),]
if (length(grep('nine', df_coll$experiment))>0){
  df_coll = df_coll[-grep('nine', df_coll$experiment),]
}
write.table(sry, file.path(dir_home, 'output_summary.txt'), row.names = F, quote=F, sep='\t')


## do volcano plot
hist(df_coll$beta,breaks=100)
ggplot(df_coll) +
  geom_histogram(aes(x=beta), binwidth = 0.01)
g = ggplot(df_coll) +
  geom_point(aes(x=beta, y=-log10(FDR), color=experiment),size=0.5) +
  geom_vline(xintercept = 0.5, linetype='dotted', color='black') +
  mytme + theme(legend.text = element_text(size=12)) +
  geom_vline(xintercept = -0.5, linetype='dotted', color='black')
ggsave(file.path(dirname(plot_out), 'Volcano_overview.png'),
       plot=g, width=8.27*1.2,height=11.69*0.5,units='in',device='png',dpi=300)

# df_coll_hs = subset(df_coll, df_coll$protein_change=='hot_spot')
# ggplot(df_coll) +
#   # geom_hline(yintercept = -log10(0.05), linetype='dotted', color='grey') +
#   geom_point(aes(x=beta, y=-log10(FDR), color=protein_change)) +
#   geom_vline(xintercept = 0.5, linetype='dotted', color='black') +
#   geom_vline(xintercept = -0.5, linetype='dotted', color='black') +
#   facet_wrap(~experiment, scale='free') +
#   mytme + theme(
#     legend.text = element_text(size=10),
#     legend.title = element_blank(),
#     strip.background = element_blank(),
#     strip.text = element_text(size=14, color='black', face='bold'),
#     panel.border = element_rect(color='black', size=1, fill='transparent')
#     #legend.position = 'none'
#   )

## Gene ontology enrichment and also enrichment of transcription factors

### Run GO
cache = TRUE
if (cache){
  load(file.path(dir_home, 'GO_result.RData'))
  # load(file.path(dir_home, 'GO_result_no_BG_filter.RData'))
} else {
  load('/Users/jefft/Desktop/p53_project/datasets/TCGA-Pan-Nine/gene_matrix.RData')
  mtx = as.data.frame(mtx)
  term2gene = load_TRUST_term2gene()
  beta_cutoff = beta_cutoff # analyse genes with beta value above the cutoff
  min_gene = 20  # TODO?
  result = data.frame()
  go_coll = list()
  for (epr in unique(df_coll$experiment)){
    bg = rownames(mtx)[which(mtx[,toupper(strsplit(epr, split='_')[[1]][2])]==1)]
    df = subset(df_coll, df_coll$experiment == epr)
    for (mut in unique(df$protein_change)){
      print(paste(epr, mut, sep='_'))
      # num of pos/neg correlated genes and enriched GO terms
      msg = rep(0,6)
      df_sub = subset(df, df$protein_change == mut)
      df_sub = subset(df_sub, abs(df_sub$beta)> beta_cutoff)
      df_up = subset(df_sub, beta > beta_cutoff)
      df_down = subset(df_sub, beta < -beta_cutoff)
      msg[1] = nrow(df_up)
      msg[2] = nrow(df_down)
      cot = 3
      fsx = c('_pos_GO.pdf', '_pos_TF_GO.pdf', '_neg_GO.pdf', '_neg_TF_GO.pdf')
      fsx_shot = c('pos', 'pos_TF', 'neg', 'neg_TF')
      for (gp in list(df_up, df_down)){
        if (nrow(gp)>min_gene){
          up_GO = do_GO(gp, background = bg, ont = ont)
          if (is.null(up_GO)){
            msg[cot] = 0
          } else {
            po = pcs_GO_out(up_GO, pvaluecutoff = 0.05, 
                            filename = paste(epr, '_', mut, fsx[cot-2], sep=''), 
                            dir = plot_out)
            msg[cot] = po$n_sig_term
          }
          if (msg[cot]>0){
            go_coll[[paste(epr, mut, fsx_shot[cot-2], sep='-')]] = po$result
          }
        }
        if (nrow(gp)>min_gene){
          if (length(intersect(gp$gene, term2gene$gene))!=0){
            ern = enricher(gp$gene, TERM2GENE = term2gene,
                           universe = intersect(term2gene$gene, bg),
                           pvalueCutoff = 0.05)
            if (is.null(ern)){
              msg[cot+1] = 0
            } else {
              po = pcs_GO_out(ern, pvaluecutoff = 0.05, 
                              filename = paste(epr, '_', mut, fsx[cot-1], sep=''), 
                              dir = plot_out)
              msg[cot+1] = po$n_sig_term
            }
          } else {
            msg[cot+1] = 0
          }
          if (msg[cot+1]>0){
            go_coll[[paste(epr, mut, fsx_shot[cot-1], sep='-')]] = po$result
          }
        }
        cot = cot + 2
      }
      result = rbind(result, c(epr, mut, msg))
    }
  }
  colnames(result) = c('Experiment', 'Mutation', 'Pos_gene', 'Neg_gene', 'Pos_GO', 'Pos_TF_GO', 'Neg_GO', 'Neg_TF_GO')
  gc()
  result[,3:8] = as.numeric(as.matrix(result[,3:8]))
  save(go_coll, result, file = file.path(dir_home, go_rs_nm))
  
  ## QC gene ontology, require 5% gene ratio
  # to_remove = c()
  # sign_map = c('pos','neg','pos_TF','neg_TF')
  # names(sign_map) = c('Pos_GO', 'Neg_GO', 'Pos_TF_GO', 'Neg_TF_GO')
  # for (i in 1:nrow(result)){
  #   for (j in c('Pos_GO', 'Neg_GO', 'Pos_TF_GO', 'Neg_TF_GO')){
  #     if (result[[j]][i] > 0){
  #       sign = sign_map[j]
  #       cd = paste(result$Experiment[i], result$Mutation[i], sign, sep='-')
  #       go_coll[[cd]] = flt_df(go_coll[[cd]], cutoff = 0.05)
  #       result[[j]][i] = nrow(go_coll[[cd]])
  #       if (nrow(go_coll[[cd]])==0){
  #         to_remove = c(to_remove, cd)
  #       }
  #     }
  #   }
  # }
  # if (!is.null(to_remove)){
  #   go_coll = go_coll[which(!names(go_coll) %in% to_remove)]
  # }
  # save(go_coll, result, file = file.path(dir_home, 'GO_result.RData'))
}


load(file.path(dir_home, go_rs_nm))
# plot eQTL result
# print(result)
result['code'] = sapply(result$Experiment, function(x){
  return(toupper(strsplit(x, split = '_')[[1]][2]))
})
result_back = result
result = result_back
exclude = c()
# exclude = c('HNSC', 'STAD', 'LUAD') # too few mutations
name_map = read.csv('/Users/jefft/Desktop/p53_project/datasets/mut_name_map.csv', header=T)
mutation_to_plt = name_map$label
names(mutation_to_plt) = name_map$code
result = subset(result, !result$code %in% exclude)
result = subset(result, result$Mutation %in% names(mutation_to_plt))
result$Mutation = mutation_to_plt[result$Mutation]
result_long = gather(result, key='gp', value='count', 3:8)
result_long = result_long[-grep('TF', result_long$gp),]
result_long$count = as.numeric(result_long$count)
result_long$gp = factor(result_long$gp, levels = c('Pos_gene', 'Neg_gene', 'Pos_GO', 'Neg_GO'))
#to_plt = c('Pos_gene', 'Neg_gene')
#to_plt = c('Pos_GO', 'Neg_GO')
# strc_mut_plt = c('DNA contact', 'Core')
strc_mut_plt = c('DNA contact', 'Conformational', 'Sandwich')
other_mut_plt = mutation_to_plt[which(!mutation_to_plt %in% strc_mut_plt)]

group_mut = c('Hotspots','DNA contact', 'Conformational', 'Sandwich')


plt.list = list()
count = 1
for (to_plt in list(c('Pos_gene', 'Neg_gene'), c('Pos_GO', 'Neg_GO'))){
  if (length(grep('gene', to_plt))>0){
    lb = 'Gene'
  } else {
    lb = 'GO'
  }
  
  # add dummy rows
  rs_sub = subset(result_long, gp %in% to_plt & Mutation %in% group_mut)
  all_group = c()
  for (cnm in names(coll)){
    cnm = toupper(strsplit(cnm, split='_')[[1]][2])
    for (sn in to_plt){
      for (mut in group_mut){
        all_group = c(all_group, paste(cnm, sn, mut, sep='^'))
      }
    }
  }
  has_group = paste(rs_sub$code, rs_sub$gp, rs_sub$Mutation, sep='^')
  no_group = all_group[!all_group %in% has_group]
  for (ng in no_group){
    info = strsplit(ng, split='\\^')[[1]]
    rs_sub = rbind(rs_sub, c(NA, info[3], info[1], info[2], 0))
  }
  rs_sub$count = as.numeric(rs_sub$count)
  rs_sub$Mutation = factor(rs_sub$Mutation, levels = group_mut)
  rs_sub$code = factor(rs_sub$code, levels = c('BRCA', 'COAD', 'LGG', 'BLCA', 'STAD', 'LUSC', 'LUAD', 'HNSC', 'OV'))
  rs_sub$count_lb = rs_sub$count
  rs_sub$count_lb[rs_sub$count==0] = NA
  
  g = ggplot(rs_sub,
             aes(x=Mutation, y=count, fill=gp)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    facet_wrap(~code, scales = 'free_x', ncol=1, strip.position = 'right') +
    scale_fill_manual(values = c('coral', 'royalblue'), name = '', 
                      labels = c('Positive', 'Negative')) +
    labs(x='', y='Count', title = lb) +
    geom_vline(xintercept = seq(1.5,1.5*(length(group_mut)-1), 1), linetype='dotted') +
    geom_text(aes(label=count_lb), position = position_dodge(width=0.9), vjust=-0.5, size=3) +
    mytme +
    theme(strip.background = element_blank(),
          panel.spacing.y = unit(0, 'mm'),
          axis.text.x = element_blank(),
          axis.text.y = element_text(face='bold', size=10), 
          strip.text.y = element_text(face='bold', size=14, angle=0),
          ) +
    scale_y_sqrt(expand = expansion(mult = c(0,0.2)), breaks=c(100,1000,2000))
  plt.list[[count]] = g
  count = count + 1
}

plt.list %>% marrangeGrob(ncol=2, nrow=1, top = '') %>%
  ggsave(file.path(dirname(plot_out), 'Overview_mut_groups.pdf'),
         plot=., width=8.27*1.3,height=11.69,units='in',device='pdf',dpi=300)

plt.list = list()
single_mut = mutation_to_plt[which(!mutation_to_plt %in% group_mut)]
count = 1
for (to_plt in list(c('Pos_gene', 'Neg_gene'), c('Pos_GO', 'Neg_GO'))){
  if (length(grep('gene', to_plt))>0){
    lb = 'Gene'
  } else {
    lb = 'GO'
  }
  
  rs_sub = subset(result_long, gp %in% to_plt & Mutation %in% single_mut)
  rs_sub$count = as.numeric(rs_sub$count)
  rs_sub$Mutation = factor(rs_sub$Mutation)
  rs_sub$code = factor(rs_sub$code, levels = c('BRCA', 'COAD', 'LGG', 'BLCA', 'STAD', 'LUSC', 'LUAD', 'HNSC', 'OV'))
  rs_sub$count_lb = rs_sub$count
  rs_sub$count_lb[rs_sub$count==0] = NA
  
  g = ggplot(rs_sub,
             aes(x=Mutation, y=count, fill=gp)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    facet_wrap(~code, scales = 'free_x', ncol=2, strip.position = 'top') +
    scale_fill_manual(values = c('coral', 'royalblue'), name = '', 
                      labels = c('Positive', 'Negative')) +
    scale_y_sqrt(expand = expansion(mult = c(0,0.2)), breaks=c(10,100,1000)) +
    geom_text(aes(label=count_lb), position = position_dodge(width=0.9), vjust=-0.5, size=3) +
    labs(x='', y='Count', title = lb) +
    mytme +
    theme(strip.background = element_blank(),
          strip.text = element_text(face='bold', size=14, angle=0),
          axis.text.y = element_text(size=10, face='bold'),
          axis.text.x = element_text(angle=45, hjust = 1, size=10),
    )
  plt.list[[count]] = g
  count = count + 1
}

plt.list %>% marrangeGrob(ncol=2, nrow=1, top = '') %>%
  ggsave(file.path(dirname(plot_out), 'Overview_single_mut.pdf'),
         plot=., width=8.27*1.3,height=11.69,units='in',device='pdf',dpi=300)




