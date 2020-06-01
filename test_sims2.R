library(tidyverse)
library(trendsegmentR)
library(purrr)
library(intervals)
# library(clue)
library(mclust)  # for adjusted Rand index
library(extraDistr)
# library(stm) 
library(lpSolve)
library(hdp) # Nicola Robert's hdp package
source("net_common.R")
source('hdp_module.R')
source('hdp_inference_simple.R')
source('docgen.R')
source("helpers.R")
source("blei_module.R")
# reticulate::source_python("dai/dai_hdp.py")  # this might cause problem if run twice.

gam_mean = 1
gam_var = 1

nicrob_wrapper = function(curr_y, max_iter) {
  y_table <- create_docword_freq(curr_y)
  # require(hdp)
  hdp_chain = hdp_quick_init(y_table, alphaa=gam_mean^2/gam_var, alphab=gam_mean/gam_var) 
  hdp_post = hdp_posterior(hdp_chain, burnin=1, n=max_iter, space=1, verbosity=0)
  zh = lapply(hdp_post@clust_dp_counts, function(y) label_mat2vec(y[-1,]))  # label_mat2vec() is in net_common.R
  return(zh)
}

slice_wrapper = function(curr_y, max_iter) {
  zh <- hdp_slice_samplerC(curr_y, beta0=0.5, gam0=1, ITRmax=max_iter, Kcap=20, Tcap=20) 
  temp <- lapply(1:max_iter, function(itr) compute_doc_labels(zh[[itr]]) + 1)
  return( temp )
}

sm_wrapper = function(curr_y, max_iter) {
  blei_splitmerge_wrapper(curr_y, max_iter, remove_dir = T)
}

nosm_wrapper = function(curr_y, max_iter) {
  blei_splitmerge_wrapper(curr_y, max_iter, split_merge = F, remove_dir = T)
}


# j = 2
# table(zh[[20]][j])
# which_most_frequent(zh[[20]][j])
# temp[[20]][j]
# temp
# compute_doc_labels(zh[[3]])+1

# seed <- .Random.seed
set.seed(1234)

J <- 50;
W <- 15;
n <- 100;

# data generation
# out <- sample_hdp(n=rep(n,J), J=J, gam0=1, beta0=3, categorical = T, W=W)
# curr_y <- out$y
# curr_z <- out$z
out = gen_doc_data(J=J, K=3, n=n, W=W)
curr_y = out$corpus1
curr_z = out$word_label

y_table <- create_docword_freq(curr_y) # do.call(rbind, lapply(1:J, function(j) tabulate(curr_y[[j]], nbins = W)))
doc_labels <- compute_doc_labels(curr_z)
# tabulate(doc_labels)
# compute_aggregate_nmi(doc_labels, tru_doc_label)

# perf_meas = "cRand"
arand = function(z1,z2) mclust::adjustedRandIndex(unlist(z1), unlist(z2))
nmi = function(z1,z2) compute_mutual_info(unlist(z1), unlist(z2))
nrep = 5

methods = list(slice = slice_wrapper, nicrob = nicrob_wrapper, sm=sm_wrapper, nosm = nosm_wrapper)
method_names = names(methods)
nmethods = length(methods)
ITRmax <- 50L
result = NULL
for (r in 1:nrep) {
  for (m in 1:nmethods) {
    id = (r-1)*nmethods + m
    print(id)
    method = methods[[m]]
    zh <- method(curr_y, ITRmax)  
    result[[id]] <- tibble(arand = sapply(1:ITRmax, function(itr) arand(zh[[itr]], doc_labels)), 
                          nmi = sapply(1:ITRmax, function(itr) nmi(zh[[itr]], doc_labels)), 
                          itr = 1:ITRmax,
                          method = method_names[m], rep=r)
   
  }
  
  # zh <- hdp_slice_samplerC(curr_y, beta0=0.5, gam0=1, ITRmax=ITRmax, Kcap=20, Tcap=20) 
  # 
  # result[[r]] <- tibble(perf=temp, itr=1:ITRmax, method="slice", rep=r)
  # 
  # zh <- apply_nicrob_hdp(y_table, ITRmax=ITRmax, burnin = 1, alphaa=gam_mean^2/gam_var, alphab= gam_mean/gam_var)
  # temp <- sapply(1:ITRmax, function(itr) perf_meas(zh[[itr]], doc_labels))
  # result[[r]] <- result[[r]] %>% add_row(tibble(perf=temp, itr=1:ITRmax, method="NicRob", rep=r))
  # 
  # # dtof = dai_warpper(curr_y, ITRmax) # dtof : doc topic f
  # # zh = lapply(dtof, function(freq_table) max.col(freq_table))
  # # temp <- sapply(1:ITRmax, function(itr) perf_meas(zh[[itr]], doc_labels))
  # # result[[r]] <- result[[r]] %>% add_row(tibble(perf=temp, itr=1:ITRmax, method="DAI", rep=r))
  # zh = blei_splitmerge_wrapper(curr_y, max_iter=ITRmax, remove_dir = T)
  # temp <- sapply(1:ITRmax, function(itr) perf_meas(zh[[itr]], doc_labels))
  # result[[r]] <- result[[r]] %>% add_row(tibble(perf=temp, itr=1:ITRmax, method="Blei (sp-mr)", rep=r))
}
res = bind_rows(result)
saveRDS(res, "test2.rds")

mean_res = res %>% group_by(itr,method) %>% summarise(arand = mean(arand), nmi = mean(nmi))


ggplot(mean_res, aes(x=itr, y=arand, color=method)) + geom_point() + 
  theme_bw() +
  xlab("iteration") + ylab("Adjusted Rand") +
  theme(legend.background=element_blank(), 
        legend.title=element_blank(), 
        legend.position = c(0.8, 0.3),
        legend.text = element_text(size=10),
        text = element_text(size=12)) +
  guides(fill=guide_legend(keywidth=0.25,keyheight=0.25,default.unit="inch"))
# ggsave("test.png",width=4,height=5) 


ari = mean_res %>% dplyr::filter(method =="slice") %>% pull(mean_perf)
detect_change(ari)
# ggplot(res, aes(x=factor(itr), y=perf, color=method)) + geom_boxplot()