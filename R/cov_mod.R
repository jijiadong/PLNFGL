#' Function for calculating the covariance matrix by moment method
#' @param Y the list of matrices of scRNA-seq count data
#' @param Smin a small positive number added to prevent the estimation of library size from being zero
#' @return a list of covariance matrices by moment method
#' @export
cov_q <- function(Y,Smin=0.1){

  K <- length(Y)

  # estimate the library size Si by total sum scaling
  Si <- list()
  for (k in 1:K){
    Si[[k]] <- rowSums(Y[[k]]) + Smin/ncol(Y[[k]])
  }

  # estimate the alpha[j] and sample covariance matrix
  alpha <- list()
  sigma = list()
  for (k in 1:K){
    n <- dim(Y[[k]])[1]
    YS <- Y[[k]]/Si[[k]]
    alpha[[k]] <- colMeans(YS)
    sigma[[k]] <- log(t(YS) %*% YS / n) - log(outer(alpha[[k]],alpha[[k]],"*"))
    diag(sigma[[k]]) <- log(colMeans(YS^2) - colMeans(YS/Si[[k]])) -log(alpha[[k]]^2)
  }
  sigma
}

#' Function for modifying the covariance matrix to positive definiteness
#' @param Y the list of matrices of scRNA-seq count data
#' @param smin a small positive number added to prevent the estimation of library size from being zero
#' @param m a replacement value for -Inf when all elements of the covariance matrix for moment estimation are -Inf
#' @return a list of covariance matrices by moment method and positive definite modification
#' @export
cov_mod <- function(Y,smin=0.1,m=NULL){
  p <- dim(Y[[1]])[2]
  K <- length(Y)
  n <- rep(0,K)
  for( k in 1:K ){ n[k] <- dim(Y[[k]])[1] }
  S0 <- cov_q(Y, Smin=smin)
  a <- list()
  for (k in 1:K){
    a[[k]] <- S0[[k]][upper.tri(S0[[k]])]
    if ( length(which(!is.infinite(a[[k]])))==0 ){
      S0[[k]][is.infinite(S0[[k]])] <- m[[k]]
    }else{
      S0[[k]][is.infinite(S0[[k]])] <- -max(abs(a[[k]][!is.infinite(a[[k]])]))
    }
  }
  #modified sample covariance matrix
  sigma_hat1 <- list()
  for (k in 1:K){
    obj <- list()
    prob <- list()
    A <- Variable(p, p, PSD=TRUE )
    obj[[k]] <- max(abs(A - S0[[k]]))
    prob[[k]] <- Problem(Minimize(obj[[k]]))
    sigma_hat1[[k]] <- CVXR::solve(prob[[k]], solver = 'SCS')
  }
  sigma_hat2 <- list()
  for (k in 1:K){
    I <- diag(rep(1, p))
    sigma_hat2[[k]] <- sigma_hat1[[k]][[1]] + max(abs(sigma_hat1[[k]][[1]] - S0[[k]]))*I
  }
  S <- sigma_hat2
  return(S)
}

