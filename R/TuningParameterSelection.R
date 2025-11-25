#' Tuning parameter selection by AIC
#'
#' @param mat.list.t.gau the list of matrices, which contains the Gaussian distributed counts
#' @param lam1 numeric, Tuning parameter for the sparsity
#' @param lam2 numeric, Tuning parameter for the similarities across the conditions
#' @param returnJGL logical, whether return the JGL model results. default is F.
#' @return a list of 1) a vector of tuning parameter values and the AIC value; 2) the JGL result object
#' @export
AIC_select <- function(mat.list.t.gau, res.S, m, lam1 = 0.1, lam2 = 0.1, returnJGL = F){
  res.jgl <- PNJGL(mat.list.t.gau, S=res.S, lambda1 = lam1, lambda2 = lam2, return.whole.theta = T, m=m)
  AIC <- 0; bic <- 0; ebic <- 0
  p <- ncol(res.jgl$theta[[1]])
  nk <- list()
  for (k in 1:length(mat.list.t.gau)){
    nk[[k]] <- nrow(mat.list.t.gau[[k]])
    traceST <- sum(diag(res.S[[k]] %*% res.jgl$theta[[k]]))
    Ek <- sum(! res.jgl$theta[[k]] == 0)
    detT <- det(res.jgl$theta[[k]])
    AICk <- (traceST*nk[[k]] - log(detT)*nk[[k]] + 2*Ek)/1e4
    AIC <- AIC + AICk
    # bick <- (traceST*nk - log(detT)*nk + log(nk)*Ek)/1e4
    # bic <- bic + bick
    # ebick <- (traceST*nk - log(detT)*nk + log(nk)*Ek + 4/2*Ek*log(p))/1e4
    # ebic <- ebic + ebick
  }
  res <- c(lam1, lam2, AIC)
  if (returnJGL){
    return(list(res,
                res.jgl))
  }else{
    return(res)
  }
}

#' Derive the Joint Graphical Lasso estimation result by AIC tuning parameter selection
#'
#' @param GauList a list of gaussian distributed matrices
#' @param l1vec a vector of candidate values for lambda1, if NULL, the values will be 1:5/20
#' @param l2vec a vector of candidate values for lambda2, if NULL, the values will be 1:5/50
#' @return AIC.vec: a table of AIC values ; jgl.res: the JGL result object
#' @export
getJGLTuningParamResult <- function(GauList, res.S, m, l1vec = NULL, l2vec = NULL){
  # Gaulist: samples by genes
  if (is.null(l1vec)){
    l1vec = c(.05,.1,.15,.2,.25)
  }
  if (is.null(l2vec)){
    l2vec = c(.02,.04,.06,.08,.1)
  }
  AIC.vec = NULL
  jgl.res = NULL
  for (l1 in l1vec){
    for (l2 in l2vec){
      #cat("Tuning parameter l1=", l1,", l2=", l2, "\n")
      tps <- AIC_select(mat.list.t.gau = GauList, res.S=res.S, m=m, lam1 = l1, lam2 = l2, returnJGL = T)
      #print(tps)

      if (is.null(AIC.vec)){
        jgl.res <- tps[[2]]
      } else if (tps[[1]][3] < min(AIC.vec[,3])){
        jgl.res <- tps[[2]]
      }
      AIC.vec <- rbind(AIC.vec, tps[[1]])
    }
  }
  return(list(AIC.vec = AIC.vec,
              jgl.res = jgl.res))
}

