#' Function that calculates the optimal matrix in single matched case
#'  
#' @param X The original matrix
#' @param M The column space that the output matrix should have 
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
#' FindOpt_SM(X,M,r)
FindOpt_SM <- function(X, M, r){
  M = cbind(M)
  s = dim(M)[2]
  n = dim(M)[1]
  if (n != dim(X)[1]){stop("The dimension of M and X does not match.")}
  if (s > r){stop("The number of columns in M should not exceed r.")}
  temp = X - projection(M) %*% X
  svd_result = svd(temp)$u
  if (r == s){
    Mtilde = M
  }
  else{
    Mtilde = cbind(M,svd_result[,1:(r-s)])
  }
  result = projection(Mtilde) %*% X 
  return(result)
}

#' This function only works when either the rank of M or N equal to the specified rank r (\deqn{r = max(r_c, r_r)}). It's a prerequisite of FindOpt_DM_Iterative
#' @param X The original matrix
#' @param M The column space that the output matrix should have 
#' @param N The row space that the output matrix should have
#' @param r The rank of output matrix
#' 
#' @details The resulting matrix is the solution the following problem:
#' \deqn{min_{Y}{||X - Y||^2_F}}
#' subject to \deqn{M is a subset of column space of Y, N is a subset of row space of Y}
#' 
#' @return A matrix that is the solution to the optimization problem.
FindOpt_SimplifiedCase<- function(X, M, N, r){
  r_m = dim(M)[2]
  r_n = dim(N)[2]
  m = dim(M)[1]
  n = dim(N)[1]
  projM = projection(M)
  projN = projection(N)
  if (m != dim(X)[1]){stop("The dimension of M and X does not match.")}
  if (n != dim(X)[2]){stop("The dimension of N and X does not match.")}
  if (r != max(r_m,r_n)){stop("This is not a simplied case, please use FindOpt_DM_Iterative function")}
  # Simple case
  if (r_m == r_n){
    return(projM %*% X %*% projN)
  }
  # If r = r_m > r_n, we need to find the rest of the row space and form them into a new matrix that contains the full row space and then calculate projM %*% X %*% proj(Nnew)
  else if (r_m > r_n){
    svd_result = svd(projM %*% X %*% (diag(n) - projN))
    S = svd_result$v[,1:(r_m - r_n), drop = FALSE]
    sol = projM %*% X %*% (projN + projection(S,ortho = TRUE))
    return(sol)
  }
  # If r_m < r_n = r, we need to find the rest of the column space and form them into a new matrix that contains the full column space and then calculate proj(Mnew) %*% X %*% projN
  else{
    svd_result = svd((diag(m) - projM) %*% X %*% projN)
    R = svd_result$u[,1:(r_n - r_m), drop = FALSE]
    sol = (projM + projection(R,ortho = TRUE)) %*% X %*% projN
    return(sol)
  }
}

# An iterative function that solves the optimization problem of (2) in the manuscript.
# M is the column space that the solution should have
# N is the row space that the solution should have
# r is the rank
FindOpt_DM_Iterative <- function(X, M, N, r, maxiter = 1e4, tol = .Machine$double.eps^0.5){
  M = as.matrix(M)
  N = as.matrix(N)
  X = as.matrix(X)
  r_m = dim(M)[2]
  r_n = dim(N)[2]
  m = dim(M)[1]
  n = dim(N)[1]
  # Sanity check
  if (m != dim(X)[1]){stop("The dimension of M and X does not match.")}
  if (n != dim(X)[2]){stop("The dimension of N and X does not match.")}
  if (r_m > r){stop("The number of columns in M should not exceed r.")}
  if (r_n > r){stop("The number of columns in N should not exceed r.")}
  # If this is a simplified case, call the FindOpt_SimplifiedCase function
  if (r == max(r_m,r_n)){
    return(list(result = FindOpt_SimplifiedCase(X, M, N, r)))
  }
  # Main body of the function
  else{
    projN = projection(N)
    projM = projection(M)
    # Warm start. Initialize
    temp_init = (diag(m) - projM) %*% X
    # Initialize the rest of column space.
    R = svd(temp_init)$u[,1:(r - r_m), drop = FALSE]
    # Initialize the objective function difference after one iteration as +inf to make sure it does not converge at the first step. 
    opt_diff = +Inf
    # Record the number of iterations
    num_iter = 1
    # Initialize the values of the objective function
    opt_function_val = c()
    while(opt_diff > tol){
      # Stop the loop if the current number of iteration exceeds the specified maximum iteration number 
      if (num_iter > maxiter){
        warning('Algorithm falis to converge. Consider increasing maximum iteration numbers.')
        break
      }
      # Calculate the best possible rest of row space for the specified column space 
      temp1 = (projM + projection(R, ortho = TRUE)) %*% X %*% (diag(n) - projN)
      Snew = svd(temp1)$v[,1:(r - r_n), drop = FALSE]
      # Calculate the best possible rest of column space for the calculated row space 
      temp2 = (diag(m) - projM) %*% X %*% (projN + projection(Snew, ortho = TRUE))
      Rnew = svd(temp2)$u[,1:(r - r_m), drop = FALSE]
      
      # After one iteration
      S = Snew
      R = Rnew
      result = (projM + projection(R, ortho = TRUE)) %*% X %*% (projN + projection(S, ortho = TRUE))
      # Calculate the objective value
      opt_val = sum((result - X)^2)
      opt_function_val = append(opt_function_val,opt_val)
      # After one iteration, record objective function difference for the convergence criteria
      # If the number of iteration is 1, force the opt_diff = +inf for the second iteration.
      if (num_iter == 1){
        opt_diff = +Inf
      }
      else{
        opt_diff = abs(opt_val - opt_function_val[num_iter - 1])
      }
      # iteration number + 1
      num_iter = num_iter + 1
    }
    return(list(result = result, R = R, S = S, num_iter = num_iter, opt_function_val = opt_function_val))
  }
}
