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
    }
  }
  return(result)
}

Matscale = function(X, center = TRUE, scale = TRUE, att = 'row'){
  if (att == 'row'){
    result = MatscaleRow(X, center = TRUE, scale = TRUE)
  }
  if (att == 'col'){
    temp = MatscaleRow(t(X), center = TRUE, scale = TRUE)
    result = t(temp)
  }
  return(result)
}

#' Double standardize a matrix
#'
#' @param X Matrix to be standardized.
#' @param tol Tolerance. Default is square root of machine precision.
#' @param maxIter Maximum iteration. Default is 100.
#'
#' @return A list that contains:
#' \item{Result}{The matrix after double standardization.}
#' \item{Iter}{Number of iterations the function actually runs.}
#'
#' @examples
DoubleStandardize <- function(X, tol = .Machine$double.eps^0.5, maxIter = 100){
  Xaft = X
  Xprev = X - 1
  count = 1
  while (sum((Xprev - Xaft)^2) > tol){
    if (count == maxIter){print("Warning : function does not converge")}
    else{
      tempX = Matscale(Xaft, att = 'row')
      Xprev = Xaft
      Xaft = Matscale(tempX, att = 'col')
    }
    count = count + 1
  }
  rownames(Xaft) = rownames(X)
  colnames(Xaft) = colnames(X)
  return(list("Result" = Xaft, "Iter" = count))
}

