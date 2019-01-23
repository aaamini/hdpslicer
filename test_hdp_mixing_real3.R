library(plotly)
plot1 <-function(y) plot_ly(x=seq(1,length(y)),y=~y)

source('hdp_module.R')
source('hdp_inference_simple.R')

library(MASS)
library(mvtnorm)
library(clue)

compute_aggregate_nmi <- function (z_estim, z){
  require(clue)
  tru_labs <- as.cl_hard_partition( unlist(z) )
  estim_labs <- as.cl_hard_partition( unlist(z_estim) )
  
  cl_agreement(estim_labs, tru_labs, method = 'NMI') 
}

load(file = "HDP_data/paper_labels_HDP.RData")
z <- as.integer(factor(papers_label))
z
load(file = "HDP_data/paper_words_HDP.RData")
library(Matrix)
image(Matrix(tdm.mat5))
(W <- dim(tdm.mat5)[1])
(D <- dim(tdm.mat5)[2])
y <- list()
for (j in 1:D) {
  idx <- which(tdm.mat5[,j] > 0)
  y[[j]] <- rep(idx, tdm.mat5[idx,j] )
}
lapply(y,length)

y

J <- 100
idx <-sample(1:length(y),J)
#y[[idx[10]]]
curr_y <- y[idx]
curr_z <- z[idx]
length(curr_z)

ITRmax <- 50
zh <- hdp_slice_sampler_simple(curr_y, beta0=3, gam0=1, ITRmax=ITRmax, W=W, randinit=F) 
 
res <- hdp_infer(curr_y, beta0=3, gam0=1, 
                 ITRmax=ITRmax, Kcap=20, Tcap=20, W=W)
n <- sapply(curr_y, length)
zh <- lapply(1:J, function(j) res$zb[j,1:n[j]])
       
# zh <- hdp_slice_sampler(y, beta0=5, gam0=1, 
#                         ITRmax=ITRmax, Kcap=10, Tcap=10,
#                         categorical=T, W=W, randinit=T)
  
nmi <- sapply(3:ITRmax, function(itr) compute_aggregate_nmi(zh[[itr]], zh[itr-1]))

compute_doc_labels <- function(zh) {
  sapply(1:length(zh), function(j) as.integer( names(which.max(table(zh[j]))) ) )
}


zh_flattend <- lapply(1:ITRmax, function(itr) compute_doc_labels(zh[[itr]]) )
nmi_tru <- sapply(3:ITRmax, function(itr) compute_aggregate_nmi(zh_flattend[[itr]], curr_z))
  
p <- plot_ly(x=3:ITRmax, width = 600, height = 400) %>% 
  #add_markers(y = nmi[[1]], name = 'n = 30', marker=list(symbol="diamond-open-dot")) %>%
  add_markers(y = nmi_tru, name = 'tru_nmi') %>%
  #add_markers(y = nmi, name = 'consec. nmi',marker=list(symbol="x")) %>%
  #add_trace(y = nmi[[1]], name = 'trace 1', mode = 'lines+markers') %>%
  layout(xaxis = list(title="iteration"), yaxis=list(title="NMI")) %>%
  layout(legend = list(x = 0.7, y = 0.1)) 

p


export(p, file=paste("nmi", fname, ".pdf"))
export(p, file=paste("nmi", fname, ".png"))

save(out, file="real_niloo1.RData")

p <- plot_ly(nmi_df, x=2:ITRmax, y=nmi) %>% layout(xaxis = list(title="iteration"), yaxis=list(title="NMI"))
p

#pdf("nmi_1.png")
#print(p)
#dev.off()

library(hdp)
example_data_hdp

set.seed(10)
quick_chain <- hdp_posterior(out$y, burnin=1, n=50, space=1, seed=1234)
temp <- hdp_extract_components(quick_chain)
str(temp)
plot_comp_size(temp, bty="L", lab=c(3, 5, 7))
base(quick_chain)

#plot_lik(quick_chain, bty="L")
# cRand <- cl_agreement(estim_labs, tru_labs, method = "cRand")
# acc <- cl_agreement(estim_labs, tru_labs, method = "diag")
# cat("NMI     = ", round(nmi,3), "\nacc.    = ", round(acc,3), "\nadjRand = ", round(cRand,3), "\n")

