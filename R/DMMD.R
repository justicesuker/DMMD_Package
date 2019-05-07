# This algorithm is still in development. 
# Do not use it until another update. Commented on May 6, 2019.

DMMD <- function(X1, X2, tol = .Machine$double.eps^0.5){
  result_col = SMMD(X1, X2, tol = tol)
  result_row = SMMD(X1, X2, tol = tol)
  return(list(result_col, result_row))
}