library(trendsegmentR)
library(purrr)
library(intervals)
library(clue)
library(extraDistr)
library(hdp) # Nicola Robert's hdp package
source("net_common.R")
source('hdp_module.R')
source('hdp_inference_simple.R')

compute_aggregate_nmi <- function (z_estim, z, method = "NMI"){
  require(clue)
  tru_labs <- as.cl_hard_partition( unlist(z) )
  estim_labs <- as.cl_hard_partition( unlist(z_estim) )
  
  temp = cl_agreement(estim_labs, tru_labs, method = method)
  temp[is.nan(temp)] = 0
  temp
}

compute_doc_labels = function(curr_z) {
  K0 <- max(unlist(curr_z))
  return( sapply(1:length(curr_z), function(j) which.max(tabulate(curr_z[[j]], K0))) )
}

create_docword_freq = function(curr_y) {
  W = max(unlist(curr_y))
  do.call(rbind, lapply(1:length(curr_y), function(j) tabulate(curr_y[[j]], nbins = W)))
}

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
  new_breaks = c(0, as.matrix(contract(ints, delta=1))[,1], n) # merge short intervals
  # breaks = breaks[c(T, diff(breaks) > 2)]  # remove short intervals
  
  interval_sds = unlist( map( split(y, cut(seq_along(y), new_breaks)), sd) )
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

library(Matrix)
K = 3 # number of clusters
tru_doc_label = sample(K, J, replace=T) # true document labels
word_label_dist = 0.7*label_vec2mat(tru_doc_label) + 0.3*Matrix(rdirichlet(J,rep(1,K))) # sample word distributions among labels in each document
image(Matrix(word_label_dist))
# z = lapply(1:J, function(j) label_mat2vec(t(rmultinom(n, 1, prob=word_label_dist[j,]))))
z = lapply(1:J, function(j) sample(K, n, replace=T,  prob=word_label_dist[j,])) # sample word label in each document

shifts = seq(1,0.9*W,length.out = K)
shapes = purrr::map(1:K, ~dnorm(1:W,shifts[.],W/6))
shapes = purrr::map(1:K, ~(shapes[[.]])/sum(shapes[[.]]))
shape_mat = do.call(rbind, shapes)
plot(shape_mat[3,])
# phi = rdirichlet(K, rep(1/W,W))
alpha0 = 3
phi = rdirichlet(K, alpha0*shape_mat)
image(Matrix(phi))
y = lapply(1:J, function(j) rcat(n, phi[z[[j]],]))
curr_z = z
curr_y = y

y_table <- create_docword_freq(curr_y) # do.call(rbind, lapply(1:J, function(j) tabulate(curr_y[[j]], nbins = W)))
doc_labels <- compute_doc_labels(curr_z)
tabulate(doc_labels)
compute_aggregate_nmi(doc_labels, tru_doc_label)

## methods
ITRmax <- 100
zh <- hdp_slice_samplerC(curr_y, beta0=.5, gam0=1, ITRmax=ITRmax, Kcap=30, Tcap=30, W=W)
# nmi <- sapply(1:ITRmax, function(itr) compute_aggregate_nmi(zh[[itr]], curr_z))
nmi <- sapply(1:ITRmax, function(itr) compute_aggregate_nmi(compute_doc_labels(zh[[itr]]), doc_labels, method="cRand"))
plot(nmi)
(cpt = detect_change(nmi))
lines(rep(cpt,2), range(nmi), col="red", lty=20)
table(compute_doc_labels(zh[[ITRmax]]))
cbind(compute_doc_labels(zh[[ITRmax]]), doc_labels)
compute_confusion_matrix(compute_doc_labels(zh[[ITRmax]]), doc_labels)

ITRmax <- 1000
gam_mean = 5
gam_var = .1
zhnr <- apply_nicrob_hdp(y_table, ITRmax=ITRmax, burnin = 1, alphaa=gam_mean^2/gam_var, alphab= gam_mean/gam_var)
nmi <- sapply(1:ITRmax, function(itr) compute_aggregate_nmi(zhnr[[itr]], doc_labels, method="cRand"))
plot(nmi)
cbind(compute_doc_labels(zhnr[[ITRmax]]), doc_labels)
compute_confusion_matrix(zhnr[[ITRmax]], doc_labels)

# nmi <- sapply(1:ITRmax, function(itr) compute_mutual_info(unlist(zh[[itr]]), unlist(curr_z)))






# plot(out$x, out$est)


# library(hdp)
# example_data_hdp
# 
# set.seed(10)
# quick_chain <- hdp_posterior(out$y, burnin=1, n=50, space=1, seed=1234)
# temp <- hdp_extract_components(quick_chain)
# str(temp)
# plot_comp_size(temp, bty="L", lab=c(3, 5, 7))
# base(quick_chain)
# 
