library(httr)
library(stringi)
library(rvest)

## Revigo: filtering GO
run_revigo = function(GO_pvalue_df, measure='SIMREL', cutoff='0.7'){
  # input: GO dataframe with at least two columns: GO term ID, pvalue (no need to be padjust)
  ustring = c()
  for (i in 1:nrow(GO_pvalue_df)){
    ustring = c(ustring, paste(GO_pvalue_df[i,1], GO_pvalue_df[i,2], sep=' '))
  }
  userData = paste(ustring, collapse = '\n')
  
  # Submit job to Revigo
  POST(
    url = "http://revigo.irb.hr/Revigo.aspx",
    body = list(
      cutoff = cutoff,
      valueType = "pvalue",
      speciesTaxon = "9606", # human species ID
      measure = measure,
      goList = userData
    ),
    encode = "form"
  ) -> res
  
  dat = content(res, encoding = "UTF-8")
  a = html_table(dat)[[1]]
  print(paste('Remaining GO terms:',sum(a$Eliminated=='False'),'/', nrow(a)))
  a = a[a$Eliminated=='False',]
  # the value column is log10(pvalue)
  return(a)
}

# Testing ====
# load('/Users/jefft/Desktop/p53_project/eQTL_experiments/TCGA-pan_VS-wt/GO_BP_result_no_BG_filter.RData')
# userData = go_coll$`tcga_brca_raw_seq-contact-pos`[,c('ID', 'pvalue')]



