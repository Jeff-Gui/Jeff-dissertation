setwd('/Users/jefft/Desktop/p53_project/scripts/eQTL')
# TCGA-pan_VS-mutneg_ult TCGA-pan_VS-wt
dir_home = '/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt'
ccle_home = '/Users/jefft/Desktop/p53_project/datasets/CCLE_22Q1/pcd'
eqtl_out = file.path(dir_home, 'outputs')
plot_out = file.path(dir_home, 'plots', 'coreVScontact')
data_out = file.path(dir_home, 'data_out')
source('utils.R')
source('enrich_utils.R')
source('../ccle_utils.R')
source('/Users/jefft/Desktop/Manuscript/set_theme.R')
library(tidyverse)
library(stringr)
library(ggvenn)
library(enrichplot)
library(ggsci)
library(ComplexUpset)
source('../overlap_utils.R')
exp_plt_out = file.path(dir_home, 'plots', 'coreVScontact')

# load genes and upset matrix ====
load('/Users/jefft/Desktop/p53_project/datasets/TCGA-Pan-Nine/gene_matrix.RData')
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
df_coll = subset(df_coll, df_coll$protein_change %in% c('contact', 'conformation', 'sandwich'))

load(file.path(data_out, 'UpsetMtx_strc_pos.RData'))
hs_mtx_coll_pos = hs_mtx_coll
load(file.path(data_out, 'UpsetMtx_strc_neg.RData'))
hs_mtx_coll_neg = hs_mtx_coll

coad_pos = hs_mtx_coll_pos$BRCA
test =  do_GO(rownames(coad_pos)[get_idx_condt(coad_pos, test_condition = c(F,T,F))],
              background = rownames(coad_pos))
dotplot(test)
goi = ext_gene_GO(test@result$geneID[c(3,4,6,7)])
goi = rownames(coad_pos)[1:40]

# STRINGdb ====
library(STRINGdb)
stringdb = STRINGdb$new(species=9606, score_threshold=400, version='11.5') # human = 9606
# map gene to protein
string_dt = stringdb$map(my_data_frame=data.frame('SYMBOL'=goi), 
                         my_data_frame_id_col_names = "SYMBOL", 
                         removeUnmappedRows = TRUE)
print(paste(nrow(string_dt),'/',length(goi),' matched.', sep=''))
# stringdb$plot_network(string_dt$STRING_id) # default vis
data_links = string_dt$STRING_id %>% stringdb$get_interactions()
data_links = data_links[-which(duplicated(data_links)),]
gc()

# generate net from STRING ====
library(igraph)
library(ggraph)
library(tidyverse)

links = data_links %>%
  mutate(from = string_dt[match(from, string_dt$STRING_id), "SYMBOL"]) %>% 
  mutate(to = string_dt[match(to, string_dt$STRING_id), "SYMBOL"]) %>%  
  dplyr::select(from, to , last_col()) %>% 
  dplyr::rename(weight = combined_score)
# remove single link (one from, one to)
links_2 = links %>% mutate(from_c = dplyr::count(., from)$n[match(from, dplyr::count(., from)$from)]) %>%
  mutate(to_c = dplyr::count(., to)$n[match(to, dplyr::count(., to)$to)]) %>%
  filter(!(from_c == 1 & to_c == 1)) %>%
  dplyr::select(1,2,3)

nodes = links_2 %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
net = graph_from_data_frame(d=links_2,vertices=nodes,directed = F)
mem = components(net)$membership
net_ids = unique(components(net)$membership)  # there may be multiple connected nets

for (i in net_ids){
  sub_net = subgraph(net, names(mem[mem==i]))
  V(sub_net)$deg = degree(sub_net)
  V(sub_net)$size = degree(sub_net)/5
  E(sub_net)$width = E(sub_net)$weight/10
  ggraph(sub_net,layout = "centrality", cent = deg)+
    geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
    geom_node_point(aes(size=size), color="orange", alpha=0.7)+
    geom_node_text(aes(filter=deg>3, label=name), size = 3, repel = T)+
    scale_edge_width(range = c(0.2,1))+
    scale_size_continuous(range = c(1,10) )+
    guides(size='none')+
    theme_graph()
}

# ggraph(net,layout = "linear", circular = TRUE)+
#   geom_edge_arc(aes(edge_width=width), color = "lightblue", show.legend = F)+
#   geom_node_point(aes(size=size), color="orange", alpha=0.7)+
#   geom_node_text(aes(filter=deg>5, label=name), size = 5, repel = F)+
#   scale_edge_width(range = c(0.2,1))+
#   scale_size_continuous(range = c(1,10) )+
#   guides(size='none')+
#   theme_graph()



