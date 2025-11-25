# PLNFGL: Poisson Log-Normal Joint Graphical Lasso

The PLNFGL R package provides a framework for joint inference of gene regulatory networks from single-cell RNA-seq (scRNA-seq) data across multiple conditions.
It addresses the challenges posed by dropout events in scRNA-seq data by modeling gene expression with a multivariate Poisson log-normal distribution.

## Installation
```
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("jijiadong/PLNFGL")
library(PLNFGL)
```
## Simulating Multi-Class Poisson Log-Normal Data
### 1. Generate True Precision Matrices
We first generate scale-free networks for three classes:
```
library(igraph)

p <- 100                 # number of nodes
K <- 3                   # number of classes

# scale-free graph
Siginv0 <- graph.generate(p, graph = "scale-free", m = 10)
Siginv0 <- Siginv0 * runif(p*p, 0.2, 0.4)
Siginv0[lower.tri(Siginv0)] <- t(Siginv0)[lower.tri(Siginv0)]
diag(Siginv0) <- 1

# create class-specific network differences
Siginv1 <- Siginv0
Siginv2 <- Siginv0; Siginv2[(0.9*p):p, (0.9*p):p] <- 0; diag(Siginv2) <- 1
Siginv3 <- Siginv0; Siginv3[(0.8*p):p, (0.8*p):p] <- 0; diag(Siginv3) <- 1

Siginv <- list(Siginv1, Siginv2, Siginv3)
```

### 2. Generate Poisson Log-Normal Data
```
n <- c(200, 200, 200)

X <- list()
Y <- list()
Sig <- list()

for(k in 1:K){
  Sig[[k]] <- cov2cor(solve(Siginv[[k]]))  # convert precision to covariance
  logX <- MASS::mvrnorm(n[k], rep(-2.3, p), Sig[[k]])  # Gaussian layer
  X[[k]] <- exp(logX)

  # library size
  size_k <- exp(rnorm(n[k], mean = log(10), sd = 0.1))
  
  # Poisson observations
  Y[[k]] <- matrix(0, n[k], p)
  for(i in 1:n[k]){
    for(j in 1:p){
      Y[[k]][i,j] <- rpois(1, lambda = size_k[i] * X[[k]][i,j])
    }
  }
}
```
What the simulated data look like:
- Y[[k]] is an n<sub>k</sub> × p matrix of count data
- Overdispersed counts due to Poisson-lognormal mixing
- Three related but not identical networks

## Compute Sample Covariance Matrices
PLNFGL uses a modified covariance estimator for count data
```
library(CVXR)

S0 <- cov_mod(Y)   # returns a list of K covariance matrices
S1 = cov2cor(S0[[1]]);S2 = cov2cor(S0[[2]]);S3 = cov2cor(S0[[3]])
S = list(S1, S2, S3)
```
## Joint Estimation of Multiple Networks
You can estimate all K class-specific networks using:
```
fit_fgl <- PNJGL(
  Y = Y,
  S = S,
  penalty = "fused",
  lambda1 = 0.6,      # sparsity
  lambda2 = 0.2,      # fusion across classes
  rho = 1,
  weights = "equal",
  maxiter = 200,
  tol = 1e-5,
  return.whole.theta = TRUE
)

fit_fgl$theta   # list of K estimated precision matrices
```
which is a p × p precision matrix for class k.
- Non-zero off-diagonal entries represent edges in the inferred network
- Shared edges across classes indicate conserved regulatory relationships
- Differences between classes reveal class-specific network rewiring
