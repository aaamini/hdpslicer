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

# p <- plot_ly(out$Y, alpha=0.75) %>%
#   add_trace(x=~x1, y=~x2, color=~factor(id), marker=list(size=5)) %>% #, marker=list(size=~2*id+3)) %>% #
#   add_trace(x=~out$phi[,1],y=~out$phi[,2], marker=list(symbol='cross',size=10,color='black'))
# p


# out$Kmax

# export(p, file=paste("true_clusters", fname, ".pdf"))
# export(p, file=paste("true_clusters", fname, ".png"))
# #attach(out)

seed <- .Random.seed
#set.seed(14586)
nvec <- c(30,100,300)
nmax <- max(nvec)
J <- 20;
W <- 10;
out <- sample_hdp(n=rep(nmax,J), J=J, gam0=1, beta0=3, categorical = T, W=W)
fname <- paste("n_", nmax, "_K_", out$Kmax,"_J_",J,"_W_",W,sep="")

# ## create feq. table 8/23/19
# y <- out$y
# z <- out$z
# y_table <- do.call(rbind, lapply(1:J, function(j) tabulate(y[[j]])))
# K0 <- max(max(unlist(z)))
# doclevel_labels <- sapply(1:J, function(j) which.max(tabulate(z[[j]], K0)))
# doclevel_labels

##
ITRmax <- 60
nmi <- list()
for (k in 1:length(nvec)) {
  cat(sprintf('--- n =  %3d---\n', nvec[k]))
  #curr_y <- lapply(out$y, function(x) x[1:nvec[k],])
  curr_y <- lapply(out$y, function(x) x[1:nvec[k]])
  curr_z <- lapply(out$z, function(x) x[1:nvec[k]])
  
  #zh <- hdp_slice_sampler(curr_y, beta0=3, gam0=1, ITRmax=ITRmax,categorical=T, W=W) 
   zh <- hdp_slice_samplerC(curr_y, beta0=3, gam0=1, ITRmax=ITRmax, Kcap=30, Tcap=30, W=W)
  
  # zh <- hdp_slice_sampler(curr_y, beta0=3, gam0=1, ITRmax=ITRmax, categorical=T, W=W)#, cat_prior_alpha=rep(1/10,10))
  #itrMAX <- length(zh)
  nmi[[k]] <- sapply(2:ITRmax, function(itr) compute_aggregate_nmi(zh[[itr]], curr_z))
  
}

p <- plot_ly(x=2:ITRmax, width = 600, height = 400) %>% 
  add_markers(y = nmi[[1]], name = 'n = 30', marker=list(symbol="diamond-open-dot")) %>%
  add_markers(y = nmi[[2]], name = 'n = 100') %>%
  add_markers(y = nmi[[3]], name = 'n = 300',marker=list(symbol="x")) %>%
  #add_trace(y = nmi[[1]], name = 'trace 1', mode = 'lines+markers') %>%
  layout(xaxis = list(title="iteration"), yaxis=list(title="NMI")) %>%
  layout(legend = list(x = 0.7, y = 0.1)) 

p

save(out, file="example_cat9c.RData")

export(p, file=paste("nmi", fname, ".pdf"))
export(p, file=paste("nmi", fname, ".png"))

p <- plot_ly(nmi_df, x=2:ITRmax, y=nmi) %>% layout(xaxis = list(title="iteration"), yaxis=list(title="NMI"))
p

save(out, file="example_cat8.RData")

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

