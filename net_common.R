# require(igraph)

label_vec2mat <- function(z, K=NULL){
  if (is.null(K)) {
    K <- max(unique(z))
  }
  # diag(K)[z,]
  Diagonal(K)[z,]
}

label_mat2vec <- function(Z){
  max.col(Z)
}

map_to_seq <- function(z, nmax=NA) {
  as.integer(factor(z, nmax=nmax))
}

compute_confusion_matrix <- function (z, y) {
  # Compute the confusion matrix between labels "y" and "z"
  # z,y Two sets of labels
  # K   number of labels in both "c" and "e"
  z = map_to_seq(z)
  y = map_to_seq(y)
  # if (is.null(K)) K = max(c(z,y))
  # t(label_vec2mat(z,K)) %*% label_vec2mat(y,K)
  t(label_vec2mat(z)) %*% label_vec2mat(y)
}


# M = label_vec2mat(c,K)'*label_vec2mat(e,K);

label_vec2mat <- function(z, K=NULL, sparse=F) {
  if (is.null(K)) K <- max(z)
  
  if (K==1) 
    return( as.matrix( rep(1,length(z)) , ncol=1) )
  else {
    if (sparse) {
      return( Diagonal(K)[z,] )  
    } else {
      return( diag(K)[z,] )
    }
  }
}


compute_mutual_info  <- function(z,y) {
  # normMUI Computes the normalized mutual information between two clusters
  #  Labels should be either vectors or n x k matrices
  
  # c = turn_into_column_ifvec(c);
  # e = turn_into_column_ifvec(e);
  
  # if ( !is.null(dim(z)) ) z = label_mat2vec(z)
  # if ( !is.null(dim(y)) ) z = label_mat2vec(y)
  
  CM = compute_confusion_matrix(z,y)
  normCM = CM / sum(CM); # normalized confusion matrix
  IDX = CM == 0 # index of zero elements of CM so that we can avoid them
  
  jointEnt = - sum( (normCM[!IDX])*log(normCM[!IDX]) )
  indpt = matrix(rowSums(normCM),ncol=1) %*% matrix(colSums(normCM),nrow=1)
  muInfo = sum(normCM[!IDX] * log(normCM[!IDX] / indpt[!IDX]) )
  
  muInfo / jointEnt
}

