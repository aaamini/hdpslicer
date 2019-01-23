
### HDP inference -- slice sampler ###
library(Rcpp)
Rcpp::sourceCpp('hdp_inference.cpp')

hdp_slice_samplerC <- function(y, beta0=3, gam0=1, 
                               ITRmax=50, Kcap=20, Tcap=20, W=10) {
  
  res <- hdp_infer_C_base(y, beta0=beta0, gam0=gam0, 
                          ITRmax=ITRmax, Kcap=Kcap, Tcap=Tcap, W=W)
  n <- sapply(y, length)
  
  J <- length(y)
  lapply(1:ITRmax,  function(itr) lapply(1:J, function(j) res$zb_list[[itr]][j,1:n[j]]))
}


hdp_slice_sampler_simple <- function(y, beta0=3, gam0=1, ITRmax=50, Kcap=20, Tcap=20, W=10, cat_prior_alpha=NA, randinit=F) {
  # y should be list of length J, where each y[[j]] contains the observations for j^th restaurant, i.e., y[[j]] is a matrix of size  n[j] x d where n[j] is the # of observation or customers in res. j
  # for categorical data W should be given and correct!
  J <- length(y)
  
  n <- sapply(y, length)
  
  
  if (randinit) {
    tb <- lapply(1:J, function(j) sample(1:Tcap, size=n[j], replace=T))
    kb <- lapply(1:J, function(j) sample(1:Kcap, size=Tcap, replace=T))
    z <-  lapply(1:J, function(j) sapply(1:n[j], function(i) kb[[j]][ tb[[j]][i] ] ))
                 
  } else {
    tb <- lapply(1:J, function(j) rep(1,n[j]))
    kb <- lapply(1:J, function(j) rep(1,Tcap))
    z  <- lapply(n,   function(nj) rep(1,nj))
    
  }
   
  u <- lapply(1:J, function(j) runif(n[j]))
  v <- lapply(1:J, function(j) runif(Tcap))
  
 
    require(extraDistr)
    #if (any(is.na(cat_prior_alpha))) cat_prior_alpha <- rep(1/W,W)
    if (any(is.na(cat_prior_alpha))) cat_prior_alpha <- rep(1/W,W)
    
    update_phi <- function(y, z, K) {
      tab <- table( factor(unlist(z), levels=1:K), factor(unlist(y), levels=1:W) )
      cat_post_alpha <- sweep(tab, 2, cat_prior_alpha, '+')
      phi <- rdirichlet(K, cat_post_alpha)
     
      list(phi=phi, cat_post_alpha=cat_post_alpha) 
    }
    
    # Ker <- function(y, phi) {
    #   if (length(y) == 0) return(1)
    #   phi[y]
    # }
    # 
  
  z_hist <- list()
  converged <- FALSE  # not used now
  itr <- 1
  while (itr < ITRmax && !converged) {
   
    # undpate gamma
    gamp <- list()
    gam <- list()
    T_all <- list()
    Tv <- rep(0,J)
    Tj_overflow <- F
    for (j in 1:J) {
      # update gamma
      g_counts <- beta_counts( tb[[j]], Tcap )
      gamp[[j]] <- rbeta( Tcap, 1 + g_counts[1,], gam0 + g_counts[2,] )
      gam[[j]] <- stick_break_func(gamp[[j]])
     
      T_all[[j]] <- sapply( 1:n[j], function(i) find_tunc_idx(gam[[j]], u[[j]][i]) )
    }
    
    
    k_counts <- beta_counts( unlist(kb), Kcap )
    betap <- rbeta(Kcap, 1 + k_counts[1,], beta0 + k_counts[2,])
    beta <- stick_break_func( betap )
    K_all <- list()
    Kv <- rep(0,J)
    for (j in 1:J) {
      K_all[[j]] <- sapply( 1:Tcap, function(t) find_tunc_idx(beta, v[[j]][t]) )
     
      
      Kv[j] <- max( K_all[[j]] )
    }
    K <- max(Kv)
    
    phi_vec <- update_phi(y, z, Kcap)$phi
    f_vec <- sapply(1:Kcap, function(k) { function(y) Ker(y, phi_vec[k,]) } )
    
    for (j in 1:J) {      
      # update k
      for (t in 1:Tcap) {
        alphap <- table( factor(y[[j]], levels=1:W), factor(tb[[j]], levels=1:Tcap) )
        
        temp1 <- log(phi_vec+1e-10) %*% alphap[,t] 
        
        temp <- temp1[ 1:K_all[[j]][t] ]
        temp <- temp - max(temp)
       
        prob_vec <- exp(temp)+1e-11
        
        if (any(is.na(prob_vec))) {
          cat('phi_vec',phi_vec,'\n')
          cat('K_all',K_all[[j]][t],'\n')
          cat('yj',y[[j]],'\n')
          cat('tbj',tb[[j]],'\n')
          cat('alphap',alphap,'\n')
          cat('temp',temp,'\n')
          cat('temp1',temp1,'\n')
          cat('prob_vec',prob_vec,'\n')
          return(list(temp=temp,alphap=alphap))
        }
        
        # prob_vec <- exp( -log(phi_vec) %*% alphap[,t] )[ 1:K_all[[j]][t] ]
        
        #prod_list <- lapply( 1:K_all[[j]][t], function(k) f_vec[[k]]( y[[j]][tb[[j]] == t, ] ) )
        #prob_vec <- safe_list_prod(prod_list)
        #prob_vec <- sapply( 1:K_all[[j]][t], function(k) prod( f_vec[[k]]( y[[j]][tb[[j]] == t, ] ) ) )
        kb[[j]][t] <- samplePmf( 1, prob_vec)
      }
      
      # update t
      for (i in 1:n[j]){
        prob_vec <- sapply(1:T_all[[j]][i], function(t) phi_vec[ kb[[j]][t], y[[j]][i] ] ) + 1e-11
        # prob_vec <- sapply(1:T_all[[j]][i], function(t) f_vec[[  kb[[j]][t] ]]( y[[j]][i,] ) )
        tb[[j]][i] <- samplePmf( 1, prob_vec)
      }
      
      # update u
      u_upper <- sapply(seq(1,n[j]), function(i) gam[[j]][ tb[[j]][i] ])
      u[[j]] <- runif(n[j], 0, u_upper)
      
      # update v
      v_upper <- sapply(seq(1,Tcap), function(t) beta[ kb[[j]][t] ])
      v[[j]] <- runif(Tcap, 0, v_upper)
      
      # update z
      z[[j]] <-  sapply(1:n[j], function(i) kb[[j]][ tb[[j]][i] ] )
      
    } # endfor j
    
    itr <- itr + 1
    if (itr %% 5 == 0) {
      cat(sprintf("%6d: ",itr),'\n')
      # cat(table(round(beta,2)),"\n")
      beta_hist <- hist(beta)
      cat(beta_hist$mids,"\n")
      cat(beta_hist$counts,"\n")
    }
    
    z_hist[[itr]] <- z
      
  } # end while
  
  #z
  z_hist
} # hdp_slice_sampler


safe_list_prod <- function(prod_list){
  log_sums <- sapply(prod_list, function(x) sum(log(x)))
  sapply(log_sums, function(x) exp(x-max(log_sums)))
}


find_tunc_idx <- function(beta, threshold) {
  # # which(cumsum(beta) > 1-threshold)[1]
  # # we return the index after the length of beta if there is insufficient atoms
  # temp <- which(c(cumsum(beta),1) > 1-threshold)[1]
  # if (any(is.na(temp))){
  #   print('---')
  #   print(beta)
  #   print(threshold)
  # }
  #   
  # temp
  
  temp <- rev( which(c(0,cumsum(beta)) < 1-threshold) )[1]
  
  if (temp > length(beta)) temp <- temp-1
  
  temp
}

### Auxiliary functions ###
# samplePmf <- function(n, pmf, normalize=F){
#   if (normalize) { 
#     pmf <- pmf / sum(pmf) 
#   }
#   sample(x = seq(1,length(pmf)), n, replace = T, prob=pmf)  
# }
samplePmf <- function(n, pmf){
  # automatic normalization of pmf by "sample" (?)
  # samples from a single pmf
  sample(x = seq(1,length(pmf)), n, replace = T, prob=pmf)  
}


stick_break_func <- function(x) {
  temp <- c(1,cumprod(1-x))
  
  temp[1:length(x)] * x
}

beta_counts <- function(z,K) {
  # z is a vector of labels, 
  # K is the maximum label number to be counted in z
  # output[1,j]: hom many labels == j are in z
  # output[2,j]: hom many labels > j are in z
  zcounts <- tabulate(z,K)
  
  rbind(zcounts, c(rev(cumsum(rev(zcounts)))[-1],0))    
}