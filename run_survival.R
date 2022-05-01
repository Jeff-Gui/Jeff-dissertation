setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
# TCGA-pan_VS-mutneg_ult TCGA-pan_VS-wt
dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt'
ccle_home = '/Users/jefft/Desktop/p53_project/datasets/CCLE_22Q1/pcd'
eqtl_out = file.path(dir_home, 'outputs')
plot_out = file.path('/Users/jefft/Desktop/p53_project/eQTL_experiments/Meta')
data_out = file.path(dir_home, 'data_out')
source('utils.R')
source('enrich_utils.R')
source('../ccle_utils.R')
source('../set_theme.R')
library(survminer)
library(survival)
library(ggsci)

## load BRCA data ====
fp = '/Users/jefft/Desktop/p53_project/datasets/9-BRCA-TCGA/clean_data.RData'
dt = load_clean_data(fp = fp, ann_bin_mut_list = c('conformation', 'contact', 'sandwich'), mode='tcga')

## add stratify info in meta ====
meta = as.data.frame(dt[[1]]@colData)
meta = subset(meta, meta$p53_state %in% c('missense', 'Wildtype') 
              # & meta$has_contact==1
              )
meta = subset(meta, meta$SUBTYPE %in% c('BRCA_LumA', 'BRCA_LumB'))

goi = 'ARHGEF2'
expr = as.numeric(assay(dt[[1]][goi,rownames(meta),'RNA']))
expr_bin = rep('low', length(expr))
expr_bin[which(expr>median(expr))] = 'high'
meta[paste(goi,'_group',sep='')] = expr_bin

## plot survival ====
meta$os_status_num = as.numeric(sapply(meta$OS_STATUS, function(x){return(strsplit(x, split=':')[[1]][1])}))
meta$pfs_status_num = as.numeric(sapply(meta$PFS_STATUS, function(x){return(strsplit(x, split=':')[[1]][1])}))

# stratify by p53 mutation
fit = survfit(Surv(meta$OS_MONTHS, meta$os_status_num)~p53_state, data=meta)
fit = survfit(Surv(meta$PFS_MONTHS, meta$pfs_status_num)~p53_state, data=meta)


# stratify by gene expression
gnm = paste(goi,'_group',sep='')
fit = survfit(Surv(meta$OS_MONTHS, meta$os_status_num)~get(gnm), data=meta)

# plot
fit
# summary(fit)
gdt = ggsurvplot(fit, data = meta, conf.int = F, pval = T, 
                #legend.labs=c(paste(goi,'high'),paste(goi,'low'))
                )
g = gdt$plot
g + scale_color_d3(palette = 'category20')

# Mutation frequency in different subtype ====
ggplot(meta[meta$SUBTYPE!='',], aes(x=SUBTYPE, fill=p53_state)) +
  geom_bar(position = 'stack') +
  coord_polar(theta = 'y')





