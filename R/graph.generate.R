#' Generate adjacency matrices under various network structures
#'
#' @param p Integer. The dimension of the adjacency matrix; must be at least 2.
#' @param graph Character string specifying the type of network structure. 
#'   Available options include \code{"star"}, \code{"hub"}, \code{"banded"}, 
#'   \code{"scale-free"}, \code{"Barabasi-Albert"}, \code{"small-world"}, 
#'   \code{"Watts-Strogatz"}, \code{"random"}, and \code{"Erdos-Renyi"}.
#' @param m Integer. The number of blocks into which the \code{p} nodes are split
#'   when constructing a block-structured network. Defaults to 1.
#' @param neighb Integer. The neighborhood size used when generating small-world 
#'   or Watts–Strogatz networks, or the number of off-diagonal bands to retain when 
#'   \code{graph = "banded"}. Must not exceed the block size.
#'
#' @return A \code{p × p} adjacency matrix representing the generated network 
#'   structure, with block-wise construction if \code{m > 1}.
#'
#' @export
graph.generate <- function(p, graph = "hub", m = 1, neighb = 2) {
  graph <- match.arg(graph,c("star", "hub", "banded", "scale-free", "Barabasi-Albert",
                             "small-world", "Watts-Strogatz", "random","Erdos-Renyi"))

  if (graph == "hub" | graph == "star") {

    generator <- function(p) {
      adj <- diag(p)
      adj[1, ] <- adj[, 1] <- 1
      adj
    }

  } else if (graph == "banded") {

    generator <- function(p) {
      if (neighb > p) {
        stop("The number of bands cannot exceed the dimension of each block")
      }
      # adj <- 1/(abs(outer(1:p,1:p,"-"))+1)
      adj <- matrix(1, nrow=p, ncol=p)
      adj[abs(row(adj)-col(adj)) > neighb] <- 0
      adj
    }

  } else if (graph == "scale-free" | graph == "Barabasi-Albert") {

    generator <- function(p) {
      G <- sample_pa(p, power = 1, directed = FALSE)
      adj <- get.adjacency(G, sparse = FALSE) + diag(p)
      adj
    }

  } else if (graph == "small-world" | graph == "Watts-Strogatz") {

    generator <- function(p) {
      G <- sample_smallworld(1, p, neighb, 0.05)
      adj <- get.adjacency(G, sparse = FALSE) + diag(p)
      adj
    }

  } else if (graph == "random" | graph == "Erdos-Renyi") {

    generator <- function(p) {
      G <- sample_gnp(p, 1/p)
      adj <- get.adjacency(G, sparse = FALSE) + diag(p)
      adj
    }

  }

  # Construct the block split
  ind <- split(seq(p), ceiling(m*seq(p)/p))

  # Fill in blocks to construct full network structure
  net <- diag(p)
  for (i in ind) {
    net[i, i] <- generator(length(i))
  }
  net
}

#' Generate a nearest-neighbor network
#'
#' This function generates a symmetric adjacency (precision-like) matrix for a 
#' \code{p}-node network by placing nodes at random coordinates in a 2D plane and 
#' connecting each node to its \code{k} nearest neighbors. Edge weights are drawn 
#' uniformly from \code{[0.2, 0.4]}, and diagonal entries are set to enforce 
#' strict diagonal dominance.
#'
#' @param p Integer. The number of nodes in the network.
#' @param k Integer. The number of nearest neighbors each node is connected to. 
#'   Must satisfy \code{k < p}. Default is 5.
#'
#' @return A \code{p × p} symmetric matrix representing the weighted adjacency 
#'   structure (or precision-like matrix). Off-diagonal entries correspond to 
#'   edge weights, and diagonal entries are set to the sum of absolute weights of 
#'   each row plus 0.1 to ensure diagonal dominance.
#'
#' @export
networkNN <- function(p, k=5){
  x <- runif(p); y <- runif(p)
  coord <- cbind(x,y)
  dist <- matrix(0, p, p)
  for (i in 2:p){
    for (j in 1:(i-1)){
      dist[i,j] <- sqrt((coord[i,1]-coord[j,1])^2 + (coord[i,2]-coord[j,2])^2)
    }
  }
  dist <- dist + t(dist)
  B <- matrix(0, p, p)
  for (i in 1:p){
    ind <- order(dist[i,])
    B[i, ind[1:k]] <- 1
  }
  A <- matrix(0, p, p)
  for (i in 2:p){
    for (j in 1:(i-1)) A[i,j] = A[j,i] <- ifelse(B[i,j] == 1 | B[j,i] == 1, runif(1, min = 0.2, max = 0.4), 0)
  }
  for (i in 1:p) A[i,i] <- sum(abs(A[i,]))+0.1

  return(A)
}



