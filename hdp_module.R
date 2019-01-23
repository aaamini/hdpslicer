### HDP inference -- slice sampler ###
source('hdp_inference.R')

### Sampling from HDP ###
sample_gem <- function(alpha, tol=1e-8) {
  
  n <- 100*max(alpha,1)
  enough_atoms <- FALSE
  while (!enough_atoms) {
    bet <- rbeta(n, 1, alpha)
    temp <- c(1,cumprod(1-bet))
    w <- temp[1:n] * bet
    Ks <- max(which(1 - cumsum(w) > tol))
    if (Ks < n) {
      enough_atoms <- TRUE
    } else {
      n <- n*1.5
      warning('doubling \n')
    }
  }
  
  w[1:Ks]
}

hdp_gauss_likelihood <- function(nj) {
  matrix(rnorm(2*nj, mean=0, sd=.2), ncol=2)
}

sample_hdp <- function(n, J, gam0=2, beta0=2, categorical=F, W=10, cat_prior_alpha=1/W) {
  # J : number of rest. 
  # n : a vector of length J
  
  gam <- list()
  t <- list()
  k <- list()
  z <- list()
  y <- list()
  
  beta <- sample_gem(beta0)
  # plot_ly(x=seq(1,length(beta)),y=~beta)
  # hist(samplePmf(1000,beta))
  for (j in 1:J) {
    gam[[j]] <- sample_gem(gam0)
    t[[j]]   <- samplePmf(n[j],gam[[j]])
    
  }
  
  Tmax <- max(sapply(t, max))
  
  for (j in 1:J) {
    k[[j]] <- samplePmf(Tmax,beta)
    z[[j]] <- sapply(seq(1,n[j]), function(i) k[[j]][ t[[j]][i] ])
  }
  
  
  if (categorical) {
    # multinomial
    require(extraDistr)
    phi_samp <- function(K) rdirichlet(K, rep(cat_prior_alpha,W))
    likelihood <- function(phi, zvec, nj)  rcat(nj, phi[zvec,])
    
  } else { # Gaussian
    W <- 2
    phi_samp <- function(K) matrix(rnorm(2*K,mean=0,sd=1.5),ncol=2) 
    likelihood <- function(phi, zvec, nj)  {
      temp <- phi[zvec,] + matrix(rnorm(2*nj, mean=0, sd=.2), ncol=2)
      colnames(temp) <- c('x1','x2')
      temp
    }
  }
  
  Kmax <- max(sapply(k, max))
  phi <-phi_samp(Kmax)
  
  Y <- data.frame()
  for (j in 1:J) {
    #y[[j]] <- phi[z[[j]],] + matrix(rnorm(2*n[j], mean=0, sd=.2), ncol=2)
    y[[j]] <- likelihood(phi, z[[j]], n[j])
    # colnames(y[[j]]) <- c('x1','x2')
    
    Y <- rbind(Y, data.frame(id=j, y[[j]]))
    # Y <- rbind(Y, cbind(y[[j]], data.frame(id=rep(j,n[[j]]))) )
  }
  
  
  list(Y=Y, y=y, z=z, kb=k, tb=t, beta=beta, gam=gam, phi=phi, Tmax=Tmax, Kmax=Kmax)
}


