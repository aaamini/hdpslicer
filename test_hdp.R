library(plotly)
plot1 <-function(y) plot_ly(x=seq(1,length(y)),y=~y)

source('hdp_module.R')

library(MASS)
library(mvtnorm)
library(clue)

nmax <- 100
out <- sample_hdp(n=rep(nmax,5), J=5, gam0=1, beta0=3)


p <- plot_ly(out$Y, alpha=0.75) %>%
  add_trace(x=~x1, y=~x2, color=~factor(id), marker=list(size=5)) %>% #, marker=list(size=~2*id+3)) %>% #
  add_trace(x=~out$phi[,1],y=~out$phi[,2], marker=list(symbol='cross',size=10,color='black'))
p

# save(out, file="example1.RData")
out$Kmax
fname <- paste("n_", nmax, "_K_", out$Kmax,sep="")
fname

export(p, file=paste("true_clusters", fname, ".pdf"))
export(p, file=paste("true_clusters", fname, ".png"))
#attach(out)


zh <- hdp_slice_sampler(out$y, beta0=3, gam0=1, ITRmax=75)

# cbind(unlist(zh), unlist(out$z))

compute_aggregate_nmi <- function (z_estim, z){
  require(clue)
  tru_labs <- as.cl_hard_partition( unlist(z) )
  estim_labs <- as.cl_hard_partition( unlist(z_estim) )
  
  cl_agreement(estim_labs, tru_labs, method = 'NMI') 
}

itrMAX <- length(zh)
nmi <- sapply(2:itrMAX, function(itr) compute_aggregate_nmi(zh[[itr]], out$z))

p <- plot_ly(x=2:itrMAX, y=nmi) %>% layout(xaxis = list(title="iteration"), yaxis=list(title="NMI"))
p

export(p, file=paste("nmi", fname, ".pdf"))
export(p, file=paste("nmi", fname, ".png"))


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

