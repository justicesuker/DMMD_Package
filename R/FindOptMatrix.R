#' Function that calculates the optimal low-rank signal matrix given partial column space basis
#'  
#' @param X The noisy data matrix
#' @param M The column space that the output signal matrix contains
#' @param r The rank of output signal matrix
#' 
#' @details The resulting matrix is the solution the following problem:
#' \deqn{min_{Y}{||X - Y||^2_F}}
#' subject to \deqn{M is a subset of column space of Y}
#' 
#' @return A matrix that is the solution to the optimization problem.
#' @export
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

#' Function that calculates the optimal low-rank signal matrix given partial column and row space basis.
#' @details This function only works when either the rank of M (column space basis) or N (row space basis) equal to the specified rank r (\deqn{r = max(r_c, r_r)}). It's a prerequisite of 'FindOpt_DM_Iterative' function.
#' The resulting matrix is the solution to the following problem:
#' \deqn{min_{Y}{||X - Y||^2_F}}
#' subject to \deqn{M is a subset of column space of Y, N is a subset of row space of Y}
#' @param X The noisy matrix (n by p)
#' @param M The column space that the output signal matrix contains (n by \deqn{r_c})
#' @param N The row space that the output signal matrix contains (p by \deqn{r_r})
#' @param r The rank of output signal matrix (\deqn{r = max(r_c, r_r)})
#' @return A list with the following elements:
#' \item{result}{A matrix that is the solution to the optimization problem}
#' \item{R}{The remaining \deqn{r - r_c} number of column bases}
#' \item{S}{The remaining \deqn{r - r_r} column bases}

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
    result = projM %*% X %*% projN
    S = matrix(rep(0, n), nrow = n)
    R = matrix(rep(0, m), nrow = m)
    return(list(result = result, R = R, S = S))
  }
  # If r = r_m > r_n, we need to find the rest of the row space and form them into a new matrix that contains the full row space and then calculate projM %*% X %*% proj(Nnew)
  else if (r_m > r_n){
    svd_result = svd(projM %*% X %*% (diag(n) - projN))
    S = svd_result$v[,1:(r_m - r_n), drop = FALSE]
    sol = projM %*% X %*% (projN + projection(S,ortho = TRUE))
    R = matrix(rep(0, m), nrow = m)
    return(list(result = sol, R = R, S = S))
  }
  # If r_m < r_n = r, we need to find the rest of the column space and form them into a new matrix that contains the full column space and then calculate proj(Mnew) %*% X %*% projN
  else{
    svd_result = svd((diag(m) - projM) %*% X %*% projN)
    R = svd_result$u[,1:(r_n - r_m), drop = FALSE]
    sol = (projM + projection(R,ortho = TRUE)) %*% X %*% projN
    S = matrix(rep(0, n), nrow = n)
    return(list(result = sol, R = R, S = S))
  }
}

#' Function that iteratively calculates the optimal low-rank signal matrix given partial column and row space basis.
#' @details This function solves the optimization problem of (2) in Dongbang Yuan & Irina Gaynanova (2022) Double-Matched Matrix Decomposition for Multi-View Data, Journal of Computational and Graphical Statistics, DOI: 10.1080/10618600.2022.2067860, that is:
#' \deqn{min_{A}{||X - A||^2_F}}
#' subject to \deqn{M is a subset of column space of A, N is a subset of row space of A}
#' @param X The noisy matrix. The dimension is n by p
#' @param M The column space that the signal matrix contains. The dimension is n by \deqn{r_c}
#' @param N The row space that the signal matrix contains. The dimension is p by \deqn{r_r}
#' @param r The rank of output signal matrix
#' @param maxiter The maximum number of iterations allowed in the calculation. Default is 1e4
#' @param tol The tolerance used to monitor the convergence. Default is the square root of machine precision
#' @return A list with the following elements:
#' \item{result}{A matrix that is the estimated solution to the optimization problem}
#' \item{R}{The remaining \deqn{r - r_c} column bases}
#' \item{S}{The remaining \deqn{r - r_r} row bases}
#' \item{num_iter}{The actual number of iteration run}
#' \item{opt_function_val}{A vector of objective values calculated at each iteration. It is supposed to be monotonically decreasing}
#' @examples
#' X = matrix(c(2,1,1,3,2,2),nrow = 3)
#' M = matrix(c(1,1,0), nrow = 3)
#' N = matrix(c(0,1), nrow = 2)
#' r = 2
#' FindOpt_DM_Iterative(X, M, N, r)
#' @export
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
    temp = FindOpt_SimplifiedCase(X, M, N, r)
    return(list(result = temp$result, R = temp$R, S = temp$S, num_iter = 0, opt_function_val = NA))
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