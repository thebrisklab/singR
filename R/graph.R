### Functions for graph.

# function to create network matrices from vectorized lower diagonals:

#' Create network matrices from vectorized lower diagonals
#' \code{vec2net} transfer the matrix vectorized lower diagonals into net to show the component image.
#' @param invector vectorized lower diagonals.
#' @param make.diag default value = 1.
#'
#' @return a net matrx
#' @export
#'
#' @examples net = vec2net(1:10)
vec2net = function(invector,make.diag=1) {
  #invector: choose(p,2) x 1, where p is the number of nodes
  nNode = (1 + sqrt(1+8*length(invector)))/2
  outNet = matrix(0,nNode,nNode)
  outNet[lower.tri(outNet)] = invector
  dim(outNet) = c(nNode,nNode)
  outNet = outNet + t(outNet)
  diag(outNet) = make.diag
  outNet
}


############
#‘ Function for plotting networks with ggplot
#' create graph dataset with netmat and mmp_order
#' a data.frame called with vectorization of reordered netmat by mmp_order.
#' @param gmatrix netmat
#' @param sort_indices mmp_order
#'
#' @return a data.frame with vectors:
#'  ## X1: vector of numerics.
#'  ## X2: vector of numerics.
#'  ## value: vectorization of reordered netmat by mmp_order.
#'
#' @export
#'
create.graph.long = function(gmatrix,sort_indices=NULL) {
  nnode = nrow(gmatrix)
  X1 = c(1:nnode)%x%rep(1,nnode) # row number comes from Kronecker product
  X2 =  rep(1,nnode)%x%c(1:nnode) # column number
  if (!is.null(sort_indices)) {
    gmatrix = gmatrix[sort_indices,sort_indices]
  }
  value = as.vector(as.matrix(gmatrix))
  data.frame(X1,X2,value)
}





