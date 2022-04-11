library(reticulate)
setwd('/Users/jefft/Desktop/p53_project/scripts/partition_test')
use_condaenv('/opt/anaconda3/envs/bio')
py_config()
source_python('/Users/jefft/Desktop/p53_project/scripts/partition_test/main.py')

# R call python ref: https://blog.csdn.net/qq_31342997/article/details/89433255

# initialize random matrix
mtx = matrix(sample(rep(sample(0:4, 100, replace = T),3)), ncol=3)
test = do_random_partition(n_time=2000, mtx=mtx, MIN_SAMPLE_PER_LEAF=10)

# New version ====
library(utils)

lappend = function (lst, ...){
  lst = c(lst, list(...))
  return(lst)
}

pmt = function(vt){
  # vt must be unique values
  if (length(vt)==1){
    return(list(vt, NULL))
  } else {
    rest = pmt(vt[2:length(vt)])
    return(append(rest, sapply(rest, function(x){return(c(vt[1],x))})))
  }
}

pmt_curtree = function(ct, min_sample=20){
  # ct must be named numerical vector
  if (length(ct)==1){
    return(list(ct))
  } else {
    above_trh = names(ct)[ct>=min_sample]
    below_trh = names(ct)[ct<min_sample]
    
    if (sum(ct[2:length(ct)])<min_sample){
      
    }
    rest = pmt(ct[2:length(ct)])

    return(append(rest, sapply(rest, function(x){return(c(ct[1],x))})))
  }
}

get_possible_group = function(vt, min_sample){
  # vt: numeric vector of discrete group ID
  ct = table(vt)
  gp_rt = list()
  all_gp = pmt(names(ct))
  for (gp in all_gp){
    if (sum(ct[gp])>=min_sample){
      gp_rt = lappend(gp_rt, gp)
    }
  }
  return(gp_rt)
}

ruc_cho_layer = function(mtx,inh_gp=NULL){
  MIN_SAMPLE = 20
  if (nrow(mtx) < MIN_SAMPLE){
    return()
  } else {
    all_gp = get_possible_group(mtx[,1], min_sample = MIN_SAMPLE)
    for (gp in all_gp){
      idx = which(mtx[,1] %in% gp)
      ## run objective, update lookup
      obj = rnorm(1,0,1)  # TODO EDIT objective function here!!!
      count <<- count + 1
      inh_gp_now = lappend(inh_gp, gp)
      if (obj > cur_best){
        cur_best <<- obj  # global variables
        cur_best_gp <<- inh_gp
      }
      if (ncol(mtx)>=2){
        ruc_cho_layer(as.matrix(mtx[idx,2:ncol(mtx)],nrow=length(idx)), inh_gp=inh_gp_now)
      }
    }
  }
}


# TESTING ====
cur_best <<- -Inf
cur_best_gp <<- NA
mtx = matrix(sample(1:10,300,replace = T),ncol=3)
# config bar
bar = txtProgressBar(style = 3)
est_max = 1
for (i in 1:ncol(mtx)){
  est_max = est_max * (2**length(unique(mtx[,i])))
}
count <<- 0
ruc_cho_layer(mtx)




