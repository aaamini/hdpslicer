library(tidyverse)
library(trendsegmentR)
library(purrr)
library(intervals)
# library(clue)
library(mclust)  # for adjusted Rand index
library(extraDistr)
library(stm) 
library(hdp) # Nicola Robert's hdp package
source("net_common.R")
source('hdp_module.R')
source('hdp_inference_simple.R')
source('docgen.R')
source("helpers.R")
reticulate::source_python("dai/dai_hdp.py")  # this might cause problem if run twice.

# Try it with: 
# breaks = c(0,1,2,69,100)
intervals_from_breaks = function(breaks) {
  m = length(breaks)
  idx = c(1, rep(2:(m-1), each=2), m)
  Intervals( matrix(breaks[idx], ncol=2, byrow = T) )
}

detect_change = function(y) {
  # require(intervals)
  # require(purrr)
  n = length(y)
  # interp = approx(1:n, y, n=5*length(y))
  out = trendsegment(y)
  
  breaks = c(0, out$cpt, n)
  ints = intervals_from_breaks(breaks) # split into intervals according to break points
  # new_breaks = c(0, as.matrix(contract(ints, delta=1))[,1], n) # merge short intervals
  new_breaks = c(as.matrix(contract(ints, delta=1))[,1]-1, n) # merge short intervals
  # breaks = breaks[c(T, diff(breaks) > 2)]  # remove short intervals
  
  interval_sds = unlist( purrr::map( split(y, cut(seq_along(y), new_breaks)), sd) )
  return( new_breaks[ which.max(interval_sds) + 1 ] )
}

apply_nicrob_hdp = function(y_table, ITRmax = 50, burnin = 1, alphaa=1, alphab=1) {
  # require(hdp)
  hdp_chain = hdp_quick_init(y_table, alphaa=alphaa, alphab=alphab) 
  hdp_post = hdp_posterior(hdp_chain, burnin=burnin, n=ITRmax, space=1, verbosity=0)
  zh = lapply(hdp_post@clust_dp_counts, function(y) label_mat2vec(y[-1,]))  # label_mat2vec() is in net_common.R
  return(zh)
}

seed <- .Random.seed
#set.seed(14586)
J <- 50;
W <- 15;
n <- 100;

# data generation
# out <- sample_hdp(n=rep(n,J), J=J, gam0=1, beta0=3, categorical = T, W=W)
# curr_y <- out$y
# curr_z <- out$z
out = gen_doc_data(J=50, K=3, n=n, W=W)
curr_y = out$corpus1
curr_z = out$word_label

y_table <- create_docword_freq(curr_y) # do.call(rbind, lapply(1:J, function(j) tabulate(curr_y[[j]], nbins = W)))
doc_labels <- compute_doc_labels(curr_z)
# tabulate(doc_labels)
# compute_aggregate_nmi(doc_labels, tru_doc_label)


gam_mean = 1
gam_var = 1
# perf_meas = "cRand"
perf_meas = function(z1,z2) mclust::adjustedRandIndex(unlist(z1), unlist(z2))
nrep = 10


ITRmax <- 50L
result = NULL
for (r in 1:nrep) {
  zh <- hdp_slice_samplerC(curr_y, beta0=0.5, gam0=1, ITRmax=ITRmax, Kcap=20, Tcap=20) 
  temp <- sapply(1:ITRmax, function(itr) perf_meas(compute_doc_labels(zh[[itr]]), doc_labels))
  result[[r]] <- tibble(perf=temp, itr=1:ITRmax, method="slice", rep=r)
  
  zh <- apply_nicrob_hdp(y_table, ITRmax=ITRmax, burnin = 1, alphaa=gam_mean^2/gam_var, alphab= gam_mean/gam_var)
  temp <- sapply(1:ITRmax, function(itr) perf_meas(zh[[itr]], doc_labels))
  result[[r]] <- result[[r]] %>% add_row(tibble(perf=temp, itr=1:ITRmax, method="NicRob", rep=r))
  
  dtof = dai_warpper(curr_y, ITRmax) # dtof : doc topic f
  zh = lapply(dtof, function(freq_table) max.col(freq_table))
  temp <- sapply(1:ITRmax, function(itr) perf_meas(zh[[itr]], doc_labels))
  result[[r]] <- result[[r]] %>% add_row(tibble(perf=temp, itr=1:ITRmax, method="DAI", rep=r))
}
res = bind_rows(result)
mean_res = res %>% group_by(itr,method) %>% summarise(mean_perf = mean(perf))

ggplot(mean_res, aes(x=itr, y=mean_perf, color=method)) + geom_point() + 
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