library(plotly)
plot1 <-function(y) plot_ly(x=seq(1,length(y)),y=~y)

source('hdp_module.R')

library(MASS)
library(mvtnorm)

out <- sample_hdp(n=rep(30,5), J=5, gam0=1, beta0=3)

plot_ly(out$Y, alpha=0.75) %>%
  add_trace(x=~x1, y=~x2, color=~factor(id), marker=list(size=5)) %>% #, marker=list(size=~2*id+3)) %>% #
  add_trace(x=~out$phi[,1],y=~out$phi[,2], marker=list(symbol='cross',size=10,color='black'))

out$Kmax
#attach(out)

library(lineprof)

l <- lineprof(hdp_slice_sampler(out$y, beta0=3, gam0=1, ITRmax=50))

shine(l)
zh <- hdp_slice_sampler(out$y, beta0=3, gam0=1, ITRmax=50)


cbind(unlist(zh), unlist(out$z))

library(clue)
tru_labs <- as.cl_hard_partition( unlist(out$z) )
estim_labs <- as.cl_hard_partition( unlist(zh) )

nmi <- cl_agreement(estim_labs, tru_labs, method = 'NMI')
cRand <- cl_agreement(estim_labs, tru_labs, method = "cRand")
acc <- cl_agreement(estim_labs, tru_labs, method = "diag")
cat("NMI     = ", round(nmi,3), "\nacc.    = ", round(acc,3), "\nadjRand = ", round(cRand,3), "\n")


# cc <- runif(10000, 0, u_upper)
# temp.df <- data.frame(u = cc, gam= u_upper)
# aggregate(.~gam, temp.df, function(x) 2*mean(x))



#k <- 9
#phi_vec[k,]
# temp_y <-  matrix(2*rnorm(200),ncol=2)
# temp_z <- f_vec[[k]](temp_y)
#plot_ly(x=temp_y[,1], y=temp_y[,2], z=temp_z, marker=list(size=2))
# plot_ly(x=temp_y[,1], y=temp_y[,2], z=temp_z, type="contour")#, autocontour = F, contours = list(start = .5, end=2, size=.5))

# f_vec <- list()
# for (k in 1:K) {
#   f_vec[[k]] <-  function(y) Ker(y, phi_vec[k,]) 
# }
# f_vec
