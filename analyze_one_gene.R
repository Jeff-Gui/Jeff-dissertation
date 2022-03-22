setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
source('utils.R')
source('/Users/jefft/Desktop/Manuscript/set_theme.R')
source('/Users/jefft/Desktop/p53_project/scripts/ccle_utils.R')
library(tidyverse)
library(stringr)
library(ggpubr)
library(gridExtra)

gene = 'BCL2L2'
primary_site = NULL
# Central_Nervous_System, Breast, Large intestine
primary_site = c('Breast') # set to NULL if do Pan
mut = c(273,248,245,282,175,249)
mut_contact = c(273,248)
mut_conform = c(175,282,249,245)
source = FALSE
source = TRUE # set to TRUE if no loading

mutation_groups = list(mut)
names(mutation_groups) = c('Hotspots')

mutation_groups = list(mut_conform, mut_contact)
names(mutation_groups) = c('HS_conform', 'HS_contact')

if (!source){
  # TCGA
  # load('/Users/jefft/Desktop/p53_project/datasets/9-BRCA-TCGA/clean_data_WGCNA.RData')
  load('/Users/jefft/Desktop/p53_project/datasets/9-BRCA-TCGA/clean_data.RData')
  tcga = dt
  # CCLE
  load('/Users/jefft/Desktop/p53_project/datasets/CCLE/clean_data_inspect.RData')
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
}

genes = c('HDGF','NOL10', 'EDA2R', 'MDM2', 'SPATA18')
genes = c('MYBL2', 'DPH2', 'DEPDC7', 'LFNG', 'BTG2')
genes = read.table('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-mutneg_ult/data_out/TCGA_BRCA-hot_spot-pos-CCLE_sig.txt', 
                   sep='\t', header = T)
genes = genes[which(genes$ccle.null.p < 0.05), ]
genes = genes$gene
names(table(ccle[[1]]@colData$PRIMARY_SITE))
plt = get_genes_plt(genes, ccle, tcga, mutation_groups, 
                    primary_site = NULL, rnai = rnai, 
                    comparison = list('rna'=list(c('Hotspots', 'Wildtype')),
                                   'rnai'=list(c('Hotspots', 'Wildtype'))))
plt %>% marrangeGrob(ncol=2, nrow=3, top = 'Significant_TCGA-BRCA-hotspot',
                     layout_matrix = matrix(1:6,byrow = T, ncol=2)) %>%
  ggsave('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-mutneg_ult/plots/Pan_hot_spot_in_breast_CCLE.pdf',
         plot=., width=8,height=16,units='in',device='pdf',dpi=300)


