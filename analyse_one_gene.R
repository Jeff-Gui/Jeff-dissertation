setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
source('utils.R')
source('/Users/jefft/Desktop/Manuscript/set_theme.R')
source('/Users/jefft/Desktop/p53_project/scripts/ccle_utils.R')
library(tidyverse)
library(stringr)
library(ggpubr)
library(gridExtra)

# load meta mutant
meta_mut_home = '/Users/jefft/Desktop/p53_project/datasets/meta_muts'
meta_mut_fp = list.files(meta_mut_home)
meta_mut_fp = sapply(meta_mut_fp, 
                     function(x){return(file.path(meta_mut_home, x))})
meta_mut = load_meta_mut(meta_mut_fp)

gene = 'SNAI1'
primary_site = NULL
# Central_Nervous_System, Breast, Large intestine
primary_site = c('Breast') # set to NULL if do Pan
mut_contact = c(273,248)
mut_conform = c(175,282,249,245)
mut_contact = c(273,248,119,120,241,276,280)
source = FALSE
source = TRUE # set to TRUE if no loading

mut = c(273,248,245,282,175,249)
mutation_groups = list(mut)
names(mutation_groups) = c('Hotspots')

mutation_groups = list(mut_conform, mut_contact)
names(mutation_groups) = c('HS_conform', 'HS_contact')

mutation_groups = list(meta_mut$aa_pos[which(meta_mut$meta_mut_id=='conformation')],
                       meta_mut$aa_pos[which(meta_mut$meta_mut_id=='contact')],
                       meta_mut$aa_pos[which(meta_mut$meta_mut_id=='sandwich')])
names(mutation_groups) = c('Conformation', 'Contact', 'Sandwich')
for (i in 1:length(mutation_groups)){
  mutation_groups[[i]] = unique(mutation_groups[[i]])
}

tcga_fp = '/Users/jefft/Desktop/p53_project/datasets/9-BRCA-TCGA/clean_data.RData'
ccle_fp = '/Users/jefft/Desktop/p53_project/datasets/CCLE/clean_data_inspect.RData'
ccle_fp = '/Users/jefft/Desktop/p53_project/datasets/CCLE_22Q1/pcd/BRCA/clean_data.RData'
if (!source){
  # TCGA
  # load('/Users/jefft/Desktop/p53_project/datasets/9-BRCA-TCGA/clean_data_WGCNA.RData')
  load(tcga_fp)
  tcga = dt
  # CCLE
  load(ccle_fp)
  ccle = dt
  remove(dt)
  gc()
  
  ## Annotate p53 state
  p53_ann = annotate_sample_mut(tcga[[2]]@data)
  tcga[[1]]@colData[['p53_state']] = 'Wildtype'
  for (i in names(p53_ann)){
    tcga[[1]]@colData[p53_ann[[i]], 'p53_state'] = i
  }
  p53_ann = annotate_sample_mut(ccle[[2]]@data)
  ccle[[1]]@colData[['p53_state']] = 'Wildtype'
  for (i in names(p53_ann)){
    ccle[[1]]@colData[p53_ann[[i]], 'p53_state'] = i
  }
  
  ## Load rnai data
  rnai = load_rnai()
  ## Load crispr data
  crispr = load_crispr_score()
  crispr = crispr[,intersect(rownames(ccle[[1]]@colData),colnames(crispr))]
}

genes = c('HDGF','NOL10', 'EDA2R', 'MDM2', 'SPATA18')
genes = c('MYBL2', 'DPH2', 'DEPDC7', 'LFNG', 'BTG2')
genes = read.table('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-mutneg_ult/data_out/TCGA_BRCA-hot_spot-pos-CCLE_sig.txt', 
                   sep='\t', header = T)
genes = genes[which(genes$ccle.null.p < 0.05), ]
genes = genes$gene
genes = read.table('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt/data_out/TCGA_BRCA-hot_spot_CCLE_null_psig.txt', 
                   sep='\t', header = T)
genes = genes$gene

mycomp = list('rna'=list(c('Contact', 'Wildtype'),
                         c('Contact', 'Conformation'),
                         c('Contact', 'Nonsense'),
                         c('Contact', 'Sandwich')),
              'rnai'=list(c('Contact', 'Wildtype')))
mycomp = list('rna'=list(c('Hotspots', 'Wildtype')),
              'rnai' = list(c('Hotspots', 'Wildtype')))
names(table(ccle[[1]]@colData$PRIMARY_SITE))

plt = get_genes_plt(genes, ccle, tcga, mutation_groups, 
                    primary_site = 'Breast', rnai = rnai,
                    comparison = mycomp,
                    flt_other = T
                    )
plt$plots %>% marrangeGrob(ncol=2, nrow=3, top = '',
                     layout_matrix = matrix(1:6,byrow = T, ncol=2)) %>%
  ggsave('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt/plots/hotspot/BRCA_triple_sig.pdf',
         plot=., width=12,height=16,units='in',device='pdf',dpi=300)

#---------------- check plot data -----------------
plt_dt = plt$data
test = subset(plt_dt, plt_dt$db=='CCLE')
test = test[order(test$p53_state),]
write.table(test, '/Users/jefft/Desktop/21Q3.txt', sep='\t', quote=F, row.names = T)
summary(test$gene_rnai[test$p53_state!='Wildtype'], na.rm=T)

