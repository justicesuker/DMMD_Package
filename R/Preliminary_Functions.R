# Center and Scale for each row in a matrix.
MatscaleRow <- function(X, center = TRUE, scale = TRUE){
  result = X
  if (center){
    Xmean = rowMeans(X)
    result = X - Xmean
  }
  n = dim(X)[1]
  p = dim(X)[2]
  sdvec = rep(0,n)
  if (scale){
    for (i in 1:n){
      tempstd = sd(X[i,])
      if (tempstd != 0 && !is.na(tempstd)){
        result[i,] = result[i,]/(tempstd * sqrt((p - 1) / p))
      }
      else{
        stop("The matrix either contains NA or there is no variation in some rows.")
        break
      }
    }
  }
  return(result)
}

# General function that does center and scale for a matrix, either row-wise or column-wise.
Matscale = function(X, center = TRUE, scale = TRUE, att = 'row'){
  if (att == 'row'){
    result = MatscaleRow(X, center = center, scale = scale)
  }
  if (att == 'col'){
    temp = MatscaleRow(t(X), center = center, scale = scale)
    result = t(temp)
  }
  return(result)
}

#' Function that double standardizes a matrix
#'
#' @param X Matrix to be standardized.
#' @param tol Tolerance. Default is square root of machine precision.
#' @param maxIter Maximum iteration. Default is 500.
#' @details After double standardization, the matrix will have mean zero of each row and column. Also, the l_2 norm of each row is p; l_2 norm of each column is n, where n is the number of rows and p is the number of columns of the original matrix.
#' 
#' @return A list that contains:
#' \item{Result}{The matrix after double standardization.}
#' \item{Iter}{Number of iterations the function actually runs.}
#'
#' @examples
#' X = matrix(c(1,0,3,1,-1,4,5,0,6), nrow = 3, ncol = 3)
#' DoubleStandardize(X)
#' 
DoubleStandardize <- function(X, tol = .Machine$double.eps^0.5, maxIter = 500){
  # Initialize 
  # Xaft stands for the matrix after one iteration
  Xaft = X
  # Xprev stands for the matrix before doing the trasformation.
  # Initialize with X-1 to make sure it doesn't converge in the first step
  Xprev = X - 1
  count = 1
  while (sum((Xprev - Xaft)^2) > tol){
    if (count == maxIter){
      print("Warning : function does not converge")
      break
    }
    else{
      # Do a row-wise center and scale
      tempX = Matscale(Xaft, att = 'row')
      Xprev = Xaft
      # Do a column-wise center and scale
      Xaft = Matscale(tempX, att = 'col')
      count = count + 1
    }
  }
  rownames(Xaft) = rownames(X)
  colnames(Xaft) = colnames(X)
  return(list("Result" = Xaft, "Iter" = count))
}

# A function that calculates the projection matrix of specified matrix X
projection <- function(X, ortho = FALSE){
  if (!ortho){
    return(X %*% solve(crossprod(X),t(X)))
  }
  else{
    return(tcrossprod(X))
  }
}

# Generate pseudo random orthogonal matrix (n by p (n >= p))
GenOrthoMatrix <- function(n, p = n, tol = 1e-1){
  if (n < p){stop("n must be greater than p")}
  # initialize the rank to be 0.
  rank = 0
  while (rank < p){
    # Generate a random Gaussian matrix
    temp = matrix(rnorm(n^2), nrow = n, ncol = n)
    # SVD on the random matrix
    svd_temp = svd(temp)
    # Calculate the rank of the random matrix
    rank = sum(svd_temp$d > tol)
    # If the rank is equal to p, stops.
  }
  # The output is the full rank u matrix of the random Gaussian matrix
  result = svd_temp$u[,1:p]
  return(result)
}

# A function that does rank-r svd approximation of a specified matrix X.
svd_recover <- function(X, svd_result = NULL, r){
  if (r < 1){
    stop("The rank must be a positive integer!")
  }
  if (is.null(svd_result)){
    svd_result = svd(X)
  }
  # When rank equals 1, perform a scalar multiplication
  if (r == 1){
    result = svd_result$d[1] * svd_result$u[,1] %*% t(svd_result$v[,1])
  }
  # Normal case. We use matrix multiplication
  if (r != 1){
    result = svd_result$u[,1:r] %*% diag(svd_result$d[1:r]) %*% t(svd_result$v[,1:r])
  }
  return(result)
}