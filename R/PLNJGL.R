#' The main function that jointly estimates multiple gene regulatory networks based on scRNA-seq data
#'
#' @param Y A list of n×p data matrices, one matrix per class.
#' @param S A list of empirical covariance matrices corresponding to \code{Y}.
#' @param penalty Character string specifying the penalty type. Either 
#'   \code{"fused"} or \code{"group"}.
#' @param lambda1 The tuning parameter for the graphical lasso penalty.
#' @param lambda2 The tuning parameter for the fused or group lasso penalty.
#' @param rho ADMM step-size parameter.
#' @param weights Class weights. Either a numeric vector of length K or 
#'   \code{"equal"} or \code{"sample.size"}.
#' @param penalize.diagonal Logical; whether to penalize diagonal entries.
#' @param maxiter Maximum number of ADMM iterations.
#' @param tol Convergence tolerance for ADMM.
#' @param warm Optional warm-start list of K matrices.
#' @param return.whole.theta Logical; if TRUE, returns full p×p estimates;
#'   otherwise returns only connected components.
#' @param screening Screening rule: \code{"fast"} or \code{"memory.efficient"}.
#' @param truncate Threshold below which estimated coefficients are set to zero.
#'
#' @return A list containing:
#'   \item{theta}{Estimated inverse covariance matrices.}
#'   \item{diag.theta.unconnected}{Diagonal estimates for unconnected nodes 
#'     (returned only if return.whole.theta = FALSE).}
#'   \item{connected}{Logical vector indicating connected nodes.}
#'
#' @export
PNJGL <-
  function(Y,S,penalty="fused",lambda1,lambda2,rho=1,weights="equal",penalize.diagonal=FALSE,maxiter=500,tol=1e-5,warm=NULL,return.whole.theta=FALSE, screening="fast",truncate=1e-5)
  {
    ## initialize:
    p = dim(Y[[1]])[2]
    K = length(Y)
    n = rep(0,K)
    for(k in 1:K) {n[k] = dim(Y[[k]])[1]}
    
    # assign feature names if none exist:
    if(length(dimnames(Y[[1]])[[2]])==0)
    {
      for(k in 1:K)
      {
        dimnames(Y[[k]])[[2]]=paste("V",1:p,sep="")
      }
    }
    
    
    # set weights:
    if(length(weights)==1){if(weights == "equal"){
      weights = rep(1,K)
    }}
    if(length(weights)==1){if(weights == "sample.size"){
      weights = n/sum(n)
    }}
    
    ## if p is very high, identify connected nodes without storing full S matrices
    if( screening == "memory.efficient" ) 
    { 
      if(penalty=="fused"){connected = screen.fgl(Y,lambda1,lambda2,weights) }	
      if(penalty=="group"){connected = screen.ggl(Y,lambda1,lambda2,weights) }	
    }
    if( screening == "fast" ) { connected = rep(TRUE,p) }
    if(!((screening == "memory.efficient")|(screening == "fast"))){stop("screening must equal \"fast\" or \"memory.efficient\".")}
    
    ### now get criterion over connected S:
    
    
    # if a penalty matrix is entered, only take its appropriate rows:
    lam1 = lambda1
    lam2 = lambda2
    if(length(lam1)>1) {lam1 = lam1[connected,connected]}
    if(length(lam2)>1) {lam2 = lam2[connected,connected]}
    
    ## examine criteria:  (value 0 where S allows theta=0 to satisfy KKT; value 1 where theta must be connected)
    if(penalty=="fused")
    {
      if(K==2)  #use bi-conditional screening rule to identify block structure exactly
      {
        crit1 = list()
        for(k in 1:K) { crit1[[k]] =  abs(S[[k]])*weights[k] > lam1 + lam2 }  
        S.sum = matrix(0,sum(connected),sum(connected))
        for(k in 1:K) {S.sum = S.sum + weights[k]*S[[k]]}
        S.sum = abs(S.sum)
        crit2 = S.sum > 2*lam1
      }
      
      if(K>2)  #use sufficient screening rule to identify larger-grained block structure
      {
        crit1 = list()
        for(k in 1:K) { crit1[[k]] =  abs(S[[k]])*weights[k] > lam1 }  
        crit2 = matrix(0,sum(connected),sum(connected))
      }
      
      # are both criteria met?
      critboth = crit2
      for(k in 1:K) {critboth = critboth + crit1[[k]]}
      critboth = (critboth!=0)				
      diag(critboth) = 1
    }
    
    if(penalty=="group")
    {
      ## examine criteria:  (value 0 where S allows theta=0 to satisfy KKT; value 1 where theta must be connected)
      tempsum = matrix(0,sum(connected),sum(connected))
      for(k in 1:K) {tempsum = tempsum + (pmax(weights[k]*abs(S[[k]]) - lam1,0))^2 }    
      critboth = tempsum > lam2^2
      diag(critboth) = 1
    }
    
    ## now identify block structure using igraph:
    g1 <- graph.adjacency(critboth)	
    cout = clusters(g1)
    blocklist = list()
    # identify unconnected elements, and get blocks:
    unconnected = c()
    #	for(i in 2:(cout$no+1))
    #	{
    #		if(sum(cout$membership==(i-1))==1) { unconnected <- c(unconnected,which(cout$membership==(i-1))) }
    #		if(sum(cout$membership==(i-1))>1) { blocklist[[length(blocklist)+1]] <- which(cout$membership==(i-1)) }
    #	}
    
    # adapt cout$membership to start with index 1:
    if(min(cout$membership)==0){cout$membership=cout$membership+1}
    for(i in 1:(cout$no))
    {
      if(sum(cout$membership==i)==1) { unconnected <- c(unconnected,which(cout$membership==i)) }
      if(sum(cout$membership==i)>1) { blocklist[[length(blocklist)+1]] <- which(cout$membership==i) }
    }
    
    # final set of connected nodes
    connected[unconnected] = FALSE
    # connected indices of connected nodes:  0 for unconnected nodes, and 1:length(connected) for the rest.  
    # maps features 1:p to their order in the connected features
    connected.index = rep(0,p)
    connected.index[connected] = 1:sum(connected)
    # regular indices of connected nodes: map connected nodes onto 1:p indexing:
    
    # redefine unconnected as !connected (up until now it's been extra nodes caught as unconnected)
    unconnected=!connected
    
    ## define theta on all connected:   (so theta is really theta.connected).
    theta = list()
    for(k in 1:K) 
    {
      theta[[k]] = matrix(0,sum(connected),sum(connected))
      if(sum(connected)>0)
      {
        dimnames(theta[[k]])[[1]]=dimnames(theta[[k]])[[2]]=dimnames(Y[[k]])[[2]][connected]	
      }
    }
    
    ## get solution on unconnected nodes
    # data:
    Yu = list()
    for(k in 1:K){Yu[[k]] = Y[[k]][,unconnected]}
    # penalty vectors:
    # note: for admm.iters.unconnected, we use the penalize.diagonal argument before calling the function.  for admm.iters, we use it IN the function.  
    if(length(lambda1)==1) { lam1.unconnected = lambda1 }
    if(length(lambda1)>1) { lam1.unconnected = diag(lambda1)[unconnected] }
    if(length(lambda2)==1) { lam2.unconnected = lambda2 }
    if(length(lambda2)>1) { lam2.unconnected = diag(lambda2)[unconnected] }
    # if penalize.diagonal==FALSE, then set the appropriate penalty vectors to zero:
    if(!penalize.diagonal)
    { 
      lam1.unconnected = lam1.unconnected * 0
      if(penalty=="group") {lam2.unconnected = lam2.unconnected * 0}
    }
    # get the unconnected portion of theta:
    if(sum(unconnected)>0) 
    {
      theta.unconnected = PLNadmm.iters.unconnected(Yu,unconnected,lambda1=lam1.unconnected,lambda2=lam2.unconnected,penalty=penalty,rho=rho,weights=weights,maxiter=maxiter,tol=tol)$Z
      for(k in 1:K) { names(theta.unconnected[[k]])=dimnames(Y[[k]])[[2]][!connected] }
    }
    if(sum(unconnected)==0) {theta.unconnected = NULL}
    
    ## now run JGL on each block of the connected nodes to fill in theta:
    if(length(blocklist)>0){
      for(i in 1:length(blocklist)){
        # the variables in the block
        bl <- blocklist[[i]] 
        Ybl = list()
        # get the data on only those variables
        for(k in 1:K) 
        {
          Ybl[[k]] = Y[[k]][,bl]
        }  
        # penalty matrices:
        if(length(lambda1)==1) { lam1.bl = lambda1 }
        if(length(lambda1)>1) { lam1.bl = lambda1[bl,bl] }
        if(length(lambda2)==1) { lam2.bl = lambda2 }
        if(length(lambda2)>1) { lam2.bl = lambda2[bl,bl] }
        # initialize lambdas:
        lam1.bl = penalty.as.matrix(lam1.bl,dim(Ybl[[1]])[2],penalize.diagonal=penalize.diagonal)
        if(penalty=="fused") {lam2.bl = penalty.as.matrix(lam2.bl,dim(Ybl[[1]])[2],penalize.diagonal=TRUE)}
        if(penalty=="group") {lam2.bl = penalty.as.matrix(lam2.bl,dim(Ybl[[1]])[2],penalize.diagonal=penalize.diagonal)}
        
        # implement warm start if desired
        if(length(warm)==0) {warm.bl = NULL}
        if(length(warm)>0)
        {
          warm.bl = list()
          for(k in 1:K) { warm.bl[[k]] = warm[[k]][bl,bl] }
        }
        # run JGL on the block:
        Thetabl = PLNadmm.iters(Ybl,lam1.bl,lam2.bl,penalty=penalty,rho=rho,weights=weights,penalize.diagonal=TRUE,maxiter=maxiter,tol=tol,warm=warm.bl)
        # update Theta with Thetabl's results:
        for(k in 1:K) {theta[[k]][connected.index[bl],connected.index[bl]] = Thetabl$Z[[k]]}   
      }}
    
    # round very small theta entries down to zero:
    if(dim(theta[[1]])[1]>0)
    {
      for(k in 1:K)
      {
        rounddown = abs(theta[[k]])<truncate; diag(rounddown)=FALSE
        theta[[k]]=theta[[k]]*(1-rounddown)
      }}
    
    # update connected:
    #stillconnected = rep(FALSE,length(connected))
    #for(k in 1:K)
    #{
    #	stillconnected = stillconnected + (rowSums(theta[[k]]!=0)>1)
    #}
    
    # return output: theta on connected nodes, diagonal theta on unconnected nodes, and the identities of the connected nodes
    if(!return.whole.theta) 
    {
      out = list(theta=theta,theta.unconnected=theta.unconnected,connected=connected)
    }
    if(return.whole.theta) 
    {
      whole.theta = list()
      for(k in 1:K) 
      {
        whole.theta[[k]] = matrix(0,p,p)
        diag(whole.theta[[k]])[unconnected] = theta.unconnected[[k]]
        whole.theta[[k]][connected,connected] = theta[[k]]
        dimnames(whole.theta[[k]])[[1]] = dimnames(whole.theta[[k]])[[2]] = dimnames(Y[[k]])[[2]]
      }
      out = list(theta=whole.theta,connected=connected)
    }
    class(out)="jgl"
    return(out)
  }


### ADMM for FGL when Theta is known to be diagonal (unconnected nodes):
#	call:	admm.iters.unconnected(Yu,lambda1=lam1.unconnected,lambda2=lam2.unconnected,rho,weights,
#   penalize.diagonal=FALSE,maxiter,tol)
# lam1.unconnected, lam2.u...  vectors of penalties.
PLNadmm.iters.unconnected = function(Y,unconnected,lambda1,lambda2,penalty="fused",rho=1,rho.increment=1,weights,maxiter = 1000,tol=1e-5,warm=NULL)
{
  K = length(Y)
  for(k in 1:K){Y[[k]]=as.matrix(Y[[k]])}
  p = dim(Y[[1]])[2]
  n=weights
  
  # define S vectors: empirical covs of each element
  S00 = list()
  S = cov_mod(Y)
  for(k in 1:K)
  {
    S00[[k]]=rep(0,p)
    for(j in 1:p){ S00[[k]][j] = S[[k]][j,j] }
  }
  
  # initialize theta:
  theta = list()
  for(k in 1:K){theta[[k]] = 1/S00[[k]]}
  # initialize Z:
  Z = list(); for(k in 1:K){Z[[k]] = rep(0,p)}
  # initialize W:
  W = list();	for(k in 1:K){W[[k]] = rep(0,p)}
  
  # initialize lambdas:
  lam1 = lambda1
  lam2 = lambda2
  
  # iterations:
  iter=0
  diff_value = 10
  while((iter==0) || (iter<maxiter && diff_value > tol))
  {
    # update theta to minimize -logdet(theta) + <S,theta> + rho/2||theta - Z + W ||^2_2:
    theta.prev = theta
    for(k in 1:K)
    {
      B = n[k]*S00[[k]] - rho*(Z[[k]] - W[[k]])  
      theta[[k]] = 1/(2*rho) * ( -B + sqrt(B^2 + 4*rho*n[k]) )  
      #		B = S00[[k]] - rho/n[k]*(Z[[k]] - W[[k]])  
      #		theta[[k]] = n[k]/(2*rho) * ( -B + sqrt(B^2 + 4*rho/n[k]) )  # written as in our paper
    }
    
    # update Z:
    # define A matrices:
    A = list()
    for(k in 1:K){ A[[k]] = theta[[k]] + W[[k]] }
    if(penalty=="fused")
    {
      # use flsa to minimize rho/2 ||Z-A||_F^2 + P(Z):
      if(K==2){Z = flsa2(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
      if(K>2){Z = flsa.general(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
    }
    if(penalty=="group")
    {
      # minimize rho/2 ||Z-A||_F^2 + P(Z):
      Z = dsgl(A,rho,lam1,lam2,penalize.diagonal=TRUE)
    }
    
    # update the dual variable W:
    for(k in 1:K){W[[k]] = W[[k]] + (theta[[k]]-Z[[k]])}
    
    # bookkeeping:
    iter = iter+1
    diff_value = 0
    for(k in 1:K) {diff_value = diff_value + sum(abs(theta[[k]] - theta.prev[[k]])) / sum(abs(theta.prev[[k]]))}
    # increment rho by a constant factor:
    rho = rho * rho.increment
  }
  
  diff = 0; for(k in 1:K){diff = diff + sum(abs(theta[[k]]-Z[[k]]))}
  out = list(theta=theta,Z=Z,diff=diff,iters=iter)
  return(out)
}

### ADMM for FGL:
PLNadmm.iters = function(Y,lambda1,lambda2,penalty="fused",rho=1,rho.increment=1,weights,penalize.diagonal,maxiter = 1000,tol=1e-5,warm=NULL)
{
  K = length(Y)
  p = dim(Y[[1]])[2]
  n=weights
  
  ns = c(); for(k in 1:K){ns[k] = dim(Y[[k]])[1]}
  S = cov_mod(Y)
  #S = list(); for(k in 1:K){S[[k]] = cov(Y[[k]])*(ns[k]-1)/ns[k]}
  # initialize theta:
  theta = list()
  for(k in 1:K){theta[[k]] = diag(1/diag(S[[k]]))}
  # initialize Z:
  Z = list(); for(k in 1:K){Z[[k]]=matrix(0,p,p)}
  # initialize W:
  W = list();	for(k in 1:K) {W[[k]] = matrix(0,p,p) }
  
  # initialize lambdas:  (shouldn't need to do this if the function is called from the main wrapper function, JGL)
  lam1 = penalty.as.matrix(lambda1,p,penalize.diagonal=penalize.diagonal)
  if(penalty=="fused") {lam2 = penalty.as.matrix(lambda2,p,penalize.diagonal=TRUE)}
  if(penalty=="group") {lam2 = penalty.as.matrix(lambda2,p,penalize.diagonal=penalize.diagonal)}
  # iterations:
  iter=0
  diff_value = 10
  while((iter==0) || (iter<maxiter && diff_value > tol))
  {
    # reporting
    #	if(iter%%10==0)
    if(FALSE)
    {
      print(paste("iter=",iter))
      if(penalty=="fused")
      {
        print(paste("crit=",crit(theta,S,n=rep(1,K),lam1,lam2,penalize.diagonal=penalize.diagonal)))
        print(paste("crit=",crit(Z,S,n=rep(1,K),lam1,lam2,penalize.diagonal=penalize.diagonal)))
      }
      if(penalty=="group"){print(paste("crit=",gcrit(theta,S,n=rep(1,K),lam1,lam2,penalize.diagonal=penalize.diagonal)))}
    }
    
    # update theta:
    theta.prev = theta
    for(k in 1:K){
      edecomp = eigen(S[[k]] - rho*Z[[k]]/n[k] + rho*W[[k]]/n[k])
      D = edecomp$values
      V = edecomp$vectors
      D2 = n[k]/(2*rho) * ( -D + sqrt(D^2 + 4*rho/n[k]) )
      theta[[k]] = V %*% diag(D2) %*% t(V)
    }
    
    # update Z:
    # define A matrices:
    A = list()
    for(k in 1:K){ A[[k]] = theta[[k]] + W[[k]] }
    if(penalty=="fused")
    {
      # use flsa to minimize rho/2 ||Z-A||_F^2 + P(Z):
      if(K==2){Z = flsa2(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
      if(K>2){Z = flsa.general(A,rho,lam1,lam2,penalize.diagonal=TRUE)}  # the option to not penalize the diagonal is exercised when we initialize the lambda matrices
    }
    if(penalty=="group")
    {
      #  minimize rho/2 ||Z-A||_F^2 + P(Z):
      Z = dsgl(A,rho,lam1,lam2,penalize.diagonal=TRUE)
    }
    
    # update the dual variable W:
    for(k in 1:K){W[[k]] = W[[k]] + (theta[[k]]-Z[[k]])}
    
    # bookkeeping:
    iter = iter+1
    diff_value = 0
    for(k in 1:K) {diff_value = diff_value + sum(abs(theta[[k]] - theta.prev[[k]])) / sum(abs(theta.prev[[k]]))}
    # increment rho by a constant factor:
    rho = rho*rho.increment
  }
  diff = 0; for(k in 1:K){diff = diff + sum(abs(theta[[k]]-Z[[k]]))}
  out = list(theta=theta,Z=Z,diff=diff,iters=iter)
  return(out)
}

flsa2 <-
  function(A,L,lam1,lam2,penalize.diagonal)  #A is a list of 2 matrices from which we apply an L2 penalty to departures
  {
    # 1st apply fused penalty:
    S1 = abs(A[[1]]-A[[2]])<=2*lam2/L
    X1 = (A[[1]]+A[[2]])/2
    Y1 = X1
    
    S2 = (A[[1]] > A[[2]]+2*lam2/L)
    X2 = A[[1]] - lam2/L
    Y2 = A[[2]] + lam2/L
    
    S3 = (A[[2]] > A[[1]]+2*lam2/L)
    X3 = A[[1]] + lam2/L
    Y3 = A[[2]] - lam2/L
    
    X = soft(a = S1*X1 + S2*X2 + S3*X3, lam = lam1/L, penalize.diagonal=penalize.diagonal)
    Y = soft(a = S1*Y1 + S2*Y2 + S3*Y3, lam = lam1/L, penalize.diagonal=penalize.diagonal)
    
    return(list(X,Y))
  }

make.adj.matrix <-
  function(theta, separate=FALSE)
  {
    K = length(theta)
    adj = list()
    if(separate)
    {
      for(k in 1:K)
      {
        adj[[k]] = (abs(theta[[k]])>1e-5)*1
      }
    }
    if(!separate)
    {
      adj = 0*theta[[1]]
      for(k in 1:K)
      {
        adj = adj+(abs(theta[[k]])>1e-5)*2^(k-1)
      }
    }
    return(adj)
  }

penalty.as.matrix <-
  function(lambda,p,penalize.diagonal)
  {
    # for matrix penalties:  check dim and symmetry:
    if(is.matrix(lambda))
    {
      if(sum(lambda!= t(lambda))>0) {stop("error: penalty matrix is not symmetric")}
      if(sum(abs(dim(lambda)-p))!=0 ) {stop("error: penalty matrix has wrong dimension")}
    }
    # for scalar penalties: convert to matrix form:
    if(length(lambda)==1) {lambda=matrix(lambda,p,p)}
    # apply the penalize.diagonal argument:
    if(!penalize.diagonal) {diag(lambda)=0}
    return(lambda)
  }

soft <-
  function(a,lam,penalize.diagonal){ # if last argument is FALSE, soft-threshold a matrix but don't penalize the diagonal
    out <- sign(a)*pmax(0, abs(a)-lam)
    if(!penalize.diagonal) diag(out) <- diag(a)
    return(out)
  }

dsgl <-
  function(A,L,lam1,lam2,penalize.diagonal)
  {
    lam1 = lam1*1/L
    lam2 = lam2*1/L
    
    if(is.matrix(A[[1]])) {p=dim(A[[1]])[1]}
    if(is.vector(A[[1]])) {p=length(A[[1]])}
    K=length(A)
    softA = A
    for(k in 1:K) {softA[[k]] = soft(A[[k]],lam1,penalize.diagonal=penalize.diagonal) }   #if penalize.diagonal=FALSE was used in ggl(), then this will not penalize the diagonal.
    normsoftA = A[[1]]*0
    for(k in 1:K) {normsoftA = normsoftA + (softA[[k]])^2}
    
    normsoftA = sqrt(normsoftA)
    
    notshrunk = (normsoftA>lam2)*1
    # reset 0 elements of normsoftA to 1 so we don't get NAs later. 
    normsoftA = normsoftA + (1-notshrunk)
    
    out = A
    for(k in 1:K)
    {
      out[[k]] = softA[[k]]*(1-lam2/normsoftA)
      out[[k]] = out[[k]]*notshrunk
    }
    return(out)
  }

flsa.general <-
  function(A,L,lam1,lam2,penalize.diagonal)
  {
    trueA = A
    if(is.matrix(A[[1]])) {p=dim(A[[1]])[1]}
    if(is.vector(A[[1]])) {p=length(A[[1]])}
    K = length(A)
    # results matrices:
    X = list()
    #for(k in 1:K) {X[[k]] = matrix(NA,p,p)} 
    for(k in 1:K) {X[[k]] = A[[1]]*NA} 
    if(is.matrix(A[[1]])) {fusions = array(FALSE,dim=c(K,K,p,p))}  
    if(is.vector(A[[1]])) {fusions = array(FALSE,dim=c(K,K,p,1))}
    
    # get starting newc: list of matrices.  newc[[k]][i,j] gives the (2*ordermats[[k]]-K-1) adjustment for lam2 at that element.  Basically, how many rank higher vs how many rank lower?
    newc = list()
    for(k in 1:K)
    {
      others = setdiff(1:K,k)
      others.smaller.k = 1:(k-1)
      newc[[k]] = A[[1]]*0               
      for(o in others) {newc[[k]] = newc[[k]] + (A[[o]]-A[[k]]< -1e-4) - (A[[o]]-A[[k]]>1e-4)} 
    }
    
    ######### start the loop here:
    for(iter in 1:(K-1))
    {
      
      # create order matrices:
      ordermats = list()
      for(k in 1:K)
      {
        others = setdiff(1:K,k)
        others.smaller.k = 1:(k-1)
        ordermats[[k]] = A[[1]]*0   
        for(o in others) {ordermats[[k]] = ordermats[[k]] + (A[[k]]-A[[o]]>1e-4)} 
        # to deal with ties, also add a unit to ordermat[[k]] if a class with a lower k has a tie at element i,j:
        if(k>1)
        {
          for(o in others.smaller.k) {ordermats[[k]] = ordermats[[k]] + (abs(A[[o]]-A[[k]])<1e-4)} 
        }
        ordermats[[k]] = ordermats[[k]] + 1
      }
      
      # create beta.g matrices, holding the solution to Holger's "unconstrained problem" 
      #  (prending we're not constraining the order of the solution to match the order of the A matrices)
      betas.g = list()
      for(k in 1:K)
      {
        betas.g[[k]] = A[[k]] - lam2/L*newc[[k]]
      }
      
      # identify and fuse all elements for which the betas.g are out of order:
      new.ordermats = list()
      for(k in 1:K)
      {
        others = setdiff(1:K,k)
        others.smaller.k = 1:(k-1)
        new.ordermats[[k]] = A[[1]]*0   
        for(o in others) {new.ordermats[[k]] = new.ordermats[[k]] + (betas.g[[k]]-betas.g[[o]]>1e-4)}  
        # to deal with ties, also add a unit to ordermat[[k]] if a class with a lower k has a tie at element i,j:
        if(k>1)
        {
          for(o in others.smaller.k) {new.ordermats[[k]] = new.ordermats[[k]] + (abs(betas.g[[o]]-betas.g[[k]])<1e-4)} 
        }
        new.ordermats[[k]] = new.ordermats[[k]] + 1
      }
      
      # identify neighboring fusions:  "fusions": K x K x p x p array: K x K matrices, T/F for fusions 
      for(k in 1:K){
        for(kp in 1:K){
          #given k,kp, declare a fusion when their ordermats entries are adjacent, and their new.ordermats entries have the opposite direction:
          fusions[k,kp,,] = fusions[k,kp,,]+           
            ((ordermats[[k]]-1==ordermats[[kp]])&(new.ordermats[[k]]<new.ordermats[[kp]]))+
            ((ordermats[[k]]+1==ordermats[[kp]])&(new.ordermats[[k]]>new.ordermats[[kp]]))+
            (abs(A[[k]]-A[[kp]])<1e-4)
          #(existing fusions, neighboring fusions, and ties)
          fusions = (fusions>0)*1
        }}
      
      
      # now we've noted fusions between all entries with adjacent ordermats entries and reversed new.ordermats entries
      # next: extend fusions to non-adjecent entries: if a-b and b-c, then connect a-c:
      for(k in 1:K){
        for(kp in 1:K){
          others = setdiff(1:K,c(k,kp))
          for(o in others)
          {
            #identify elements in o which are fused with the same element in both k and kp, then add them to the list of k-kp fusions:
            bothfused = fusions[k,o,,] & fusions[kp,o,,]    
            fusions[k,kp,,] = fusions[k,kp,,] | bothfused       
          }
        }}
      
      # now recalculate A with the new fused entries:
      # to recalculate A, for each non-zero entry, identify the classes k with which it must be fused, and get their average:
      for(k in 1:K)
      {
        others = setdiff(1:K,k)
        #fusemean and denom: the mean value of all the trueA to be fused, and the number of values to be fused:
        fusemean = trueA[[k]]
        denom = A[[1]]*0+1   
        for(o in others)
        {
          fusemean = fusemean+fusions[k,o,,]*trueA[[o]]  #add the values of the elements which must be fused to fusemean   
          denom = denom+fusions[k,o,,]     
        }	
        # now redefine A[[k]]: unchanged from trueA if there's no fusion, and the mean of the fused elements when there is fusion:
        A[[k]] = fusemean/denom
      }
      
      #newc: list of matrices.  newc[[k]][i,j] gives the (2*ordermats[[k]]-K-1) adjustment for lam2 at that element.  Basically, how many rank higher vs how many rank lower?
      newc = list()
      for(k in 1:K)
      {
        others = setdiff(1:K,k)
        others.smaller.k = 1:(k-1)
        newc[[k]] = A[[1]]*0   
        for(o in others) {newc[[k]] = newc[[k]] + (A[[o]]-A[[k]]< -1e-4) - (A[[o]]-A[[k]]>1e-4)}  
      }
      
    } #end loop here
    
    # final version of betas.g:
    for(k in 1:K)
    {
      betas.g[[k]] = A[[k]] - lam2/L*newc[[k]]
    }
    # now soft-threshold the solution matices:
    for(k in 1:K)
    {
      X[[k]] = soft(betas.g[[k]],lam=lam1/L,penalize.diagonal=penalize.diagonal)
    }
    return(X)
  }

crit <-
  function(theta, S,n, lam1, lam2, penalize.diagonal)  # theta = list of pXp matrices, length k
  {
    p = dim(S[[1]])[1]
    K = length(S)
    lam1 = penalty.as.matrix(lam1,p,penalize.diagonal=penalize.diagonal)
    lam2 = penalty.as.matrix(lam2,p,penalize.diagonal=TRUE)	
    crit = 0
    for(k in 1:length(theta))
    {
      # add log det that was entered as an argument, or else calculate it
      crit = crit+n[k]*log(det(theta[[k]]))-n[k]*sum(S[[k]]*theta[[k]])-sum(lam1*abs(theta[[k]])) 
      for(kp in k:length(theta))
      {
        crit = crit - sum(lam2*abs(theta[[k]]-theta[[kp]]))
      }
    }
    return(crit)
  }
gcrit <-
  function(theta, S,n, lam1, lam2, penalize.diagonal)  # theta = list of pXp matrices, length k
  {
    p = dim(S[[1]])[1]
    K = length(S)
    lam1 = penalty.as.matrix(lam1,p,penalize.diagonal=penalize.diagonal)
    lam2 = penalty.as.matrix(lam2,p,penalize.diagonal=penalize.diagonal)	
    crit = 0
    sumsquared = matrix(0,p,p)
    for(k in 1:length(theta))
    {
      # add log det that was entered as an argument, or else calculate it
      crit = crit+n[k]*log(det(theta[[k]]))-n[k]*sum(S[[k]]*theta[[k]])-sum(lam1*abs(theta[[k]])) 
      # add to term that will become group L2 norm:
      sumsquared = sumsquared + theta[[k]]^2
    }
    crit = crit - sum(lam2*sqrt(sumsquared))
    return(crit)
  }

plot.jgl <-
  function(x,...)
  {
    .env = "environment: namespace:JGL"
    #UseMethod("plot")
    theta=x$theta
    K=length(theta)
    adj = make.adj.matrix(theta)
    diag(adj)=0
    gadj = graph.adjacency(adj,mode="upper",weighted=TRUE)
    #weight the edges according to the classes they belong to
    E(gadj)$color = 2^(K)-get.edge.attribute(gadj,"weight")
    #plot the net using igraph
    plot(gadj, vertex.frame.color="white",layout=layout.fruchterman.reingold,
         vertex.label=NA, vertex.label.cex=3, vertex.size=1)
  }


print.jgl <-
  function(x, ...)
  {
    #	.env = "environment: namespace:JGL"
    #	UseMethod("print")
    p = dim(x$theta[[1]])[1]
    K = length(x$theta)
    #cat("Call: ")
    #dput(x$call)
    cat("\n")

    # return number of connected nodes:
    cat("Number of connected nodes: ",sum(x$connected),"\n")
    subnetlengths = nedges = c()
    for(k in 1:K)
    {
      subnetlengths[k] = length(subnetworks(x$theta)[[k]])
      nedges[k] = (sum(x$theta[[k]]!=0)-p)/2
    }
    # number of subgraphs:
    cat("Number of subnetworks in each class: ", subnetlengths, "\n")
    # number of edges per class:
    cat("Number of edges in each class: ", nedges, "\n")

    # number of shared edges:
    sharededges = matrix(0,p,p)
    for(k in 1:K){sharededges = sharededges + (x$theta[[k]]!=0)}
    nsharededges = (sum(sharededges==K)-p)/2
    cat("Number of edges shared by all classes: ", nsharededges, "\n")

    # L1 norm of classes:
    norm = c()
    for(k in 1:K) {norm[k] = sum(abs(x$theta[[k]]))-sum(abs(diag(x$theta[[k]])))}
    cat("L1 norm of off-diagonal elements of classes' Thetas: ", norm, "\n")
  }




