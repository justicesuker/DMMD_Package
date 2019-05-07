projection <- function(X){
  return(X %*% solve(crossprod(X),t(X)))
}

#' Function that calculates the optima matrix 
#'  
#' @param X The original matrix
#' @param M The matrix that contains 
#' @param r The rank of output matrix
#' 
#' @details The resulting matrix is the solution the following problem:
#' \deqn{min_{Y}{||X - Y||^2_F}}
#' subject to \deqn{M is a subset of column space of Y}
#' 
#' @return A matrix that is the solution to the optimization problem.
#'
#' @examples 
#' X = matrix(c(2,1,1,3,2,2),nrow = 3)
#' M = matrix(c(1,1,0),nrow = 3)
#' r = 2
#' FindOptMatrix(X,M,r)
FindOptMatrix <- function(X, M, r){
  s = dim(M)[2]
  n = dim(M)[1]
  if (n != dim(X)[1]){stop("The dimension of M and X does not match.")}
  if (s > r){stop("The number of columns in M should not exceed r.")}
  temp = X - projection(M) %*% X
  svd_result = svd(temp)$u
  Mtilde = cbind(M,svd_result[,1:(r-s)])
  result = projection(Mtilde) %*% X 
  return(result)
}
