### HDP inference -- slice sampler ###

hdp_slice_sampler <- function(y, beta0=3, gam0=1, ITRmax=50, doubling_factor=1.5) {
  # y should be list of length J, where each y[[j]] contains the observations for j^th restaurant, i.e., y[[j]] is a matrix of size  n[j] x d where n[j] is the # of observation or customers in res. j
  J <- length(y)
  n <- sapply(y, function(x) dim(x)[1])
  
  #Tv <- rep(1,J)
  #Kv <- rep(1,J)
  #K <- max(Kv)
  
  Kcap <- 20  # hard-coded for now
  Tjcap <- rep(20,J) # hard-coded for now
  #Tcap <- Tv
  #Kcap <- K
  tb <- lapply(1:J, function(j) rep(1,n[j]))
  kb <- lapply(1:J, function(j) rep(1,Tjcap[j]))
  z  <- lapply(n,   function(nj) rep(1,nj))
   
  u <- lapply(1:J, function(j) runif(n[j]))
  v <- lapply(1:J, function(j) runif(Tjcap[j]))
  
  kb_old <- kb
  tb_old <- tb
  u_old <- u
  v_old <- v
  z_old <- z
  
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
      g_counts <- beta_counts( tb[[j]], Tjcap[j] )
      gamp[[j]] <- rbeta( Tjcap[j], 1 + g_counts[1,], gam0 + g_counts[2,] )
      gam[[j]] <- stick_break_func(gamp[[j]])
      # plot1(gam[[j]])
      #print(u[[j]])
      #print('-')
      T_all[[j]] <- sapply( 1:n[j], function(i) find_tunc_idx(gam[[j]], u[[j]][i]) )
      Tv[j] <- max(T_all[[j]])
      #print(Tv)
      
      if ( Tv[j] > Tjcap[j] ) {
        Tj_overflow <- T
        Tjcap_old <- Tjcap[j]
        Tjcap[j] <- round( doubling_factor*Tjcap[j] )
        # pad v with enough atoms
        v_old[[j]] <- c( v_old[[j]], runif(Tjcap[j]-Tjcap_old) )
        break
      }
    }
    if (Tj_overflow)  {
      cat('Doubling Tjcap\n')
      kb <- kb_old
      tb <- tb_old
      u <- u_old
      v <- v_old
      z <- z_old
      next
    }
    # round(gam[[j]],2)
    
    k_counts <- beta_counts( unlist(kb), Kcap )
    betap <- rbeta(Kcap, 1 + k_counts[1,], beta0 + k_counts[2,])
    beta <- stick_break_func( betap )
    K_all <- list()
    Kv <- rep(0,J)
    for (j in 1:J) {
      K_all[[j]] <- sapply( 1:Tjcap[j], function(t) find_tunc_idx(beta, v[[j]][t]) )
      Kv[j] <- max( K_all[[j]] )
    }
    K <- max(Kv)
    if (K > Kcap)  {
      cat('Doubling Kcap\n')
      Kcap <- round(doubling_factor*Kcap)
      kb <- kb_old
      tb <- tb_old
      u <- u_old
      v <- v_old
      z <- z_old
      next
    }
    
    # if we got here, it is safe to save current state as old state
    kb_old <- kb
    tb_old <- tb
    u_old <- u
    v_old <- v
    z_old <- z
    
    phi_vec <- update_phi(y, z, Kcap)$phi
    f_vec <- sapply(1:Kcap, function(k) { function(y) Ker(y, phi_vec[k,]) } )
    
    for (j in 1:J) {      
      # update t
      for (i in 1:n[j]){
        prob_vec <- sapply(1:T_all[[j]][i], function(t) f_vec[[  kb[[j]][t] ]]( y[[j]][i,] ) )
        tb[[j]][i] <- samplePmf( 1, prob_vec)
      }
      
      # update k
      for (t in 1:Tjcap[j]) {
        prod_list <- lapply( 1:K_all[[j]][t], function(k) f_vec[[k]]( y[[j]][tb[[j]] == t, ] ) )
        prob_vec <- safe_list_prod(prod_list)
        #prob_vec <- sapply( 1:K_all[[j]][t], function(k) prod( f_vec[[k]]( y[[j]][tb[[j]] == t, ] ) ) )
        kb[[j]][t] <- samplePmf( 1, prob_vec)
      }
      
      # update u
      u_upper <- sapply(seq(1,n[j]), function(i) gam[[j]][ tb[[j]][i] ])
      u[[j]] <- runif(n[j], 0, u_upper)
      
      # if (any(is.na(u[[j]]))){
      #   print(">>>>>>")
      #   print('tb[j]')
      #   print(tb[[j]])
      #   print('gam[j]')
      #   print(gam[[j]])
      #   print('u[[j]]')
      #   print( u[[j]] )
      #   print("<<<<<<")
      # }
      
      # update v
      v_upper <- sapply(seq(1,Tjcap[j]), function(t) beta[ kb[[j]][t] ])
      v[[j]] <- runif(Tjcap[j], 0, v_upper)
      
      # update z
      z[[j]] <-  sapply(1:n[j], function(i) kb[[j]][ tb[[j]][i] ] )
      
    } # endfor j
    
    itr <- itr + 1
    if (itr %% 5 == 0) {
      cat(sprintf("%3d: ",itr))
      #cat(round(beta,2),"\n\n")
      cat(table(round(beta,2)),"\n\n")
      #cat('.')
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

Ker <- function(y, phi, prec2_y=5^2) {
  #require(mvtnorm)
  dmvnorm(y, mean = phi, sigma=diag(1,2)/prec2_y)
}

update_phi <- function(y, z, K, prec2_y = 5^2, prec2_phi = 1/(1.5^2)) {
  require(MASS)
  #require(mvtnorm)
  
  z_flat <- unlist(z)
  y_flat <- do.call(rbind,y)
  z_freq <- tabulate(z_flat,K)
  # sapply(1:K, function (k) sum(z_flat==k))
  
  yz_df <- data.frame(y_flat, z=z_flat) 
  temp <- which(z_freq == 0)
  
  # adding observation (y) zero for missing labels, so that we effectively sample phi from the prior for those with missing labels. 
  yz_df <- rbind(data.frame(x1=0,x2=0, z=temp), yz_df) # 
  
  phi_mean <- aggregate(.~z, yz_df, function(x) sum(x)*(prec2_y/(prec2_phi + length(x)*prec2_y)) )  # mean of phi-posterior 
  phi_mean <- as.matrix(phi_mean[,-1])
  prec2_post <- prec2_phi + z_freq*prec2_y
  
  phi_new <- do.call(rbind, lapply(1:K, function(k) mvrnorm(1, mu = phi_mean[k,], Sigma = diag(1,2)/prec2_post[k]) ) )
  
  list(phi=phi_new, prec2=prec2_post)
}

find_tunc_idx <- function(beta, threshold) {
  # which(cumsum(beta) > 1-threshold)[1]
  # we return the index after the length of beta if there is insufficient atoms
  temp <- which(c(cumsum(beta),1) > 1-threshold)[1]
  if (any(is.na(temp))){
    print('---')
    print(beta)
    print(threshold)
  }
    
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