# Function for updating individual column space
updateR <- function(X, P_rall, P_c, dimR){
  n = nrow(X)
  if (dimR == 0){
    R = matrix(rep(0, n), nrow = n)
  }
  # Otherwise, use the usual update.
  else{
    temp_init = (diag(n) - P_c) %*% X %*% P_rall
    R = svd(temp_init)$u[, 1:dimR, drop = F]
  }
  return(R)
}

# Function for updating individual row space
updateS <- function(X, P_call, P_r, dimS){
  p = ncol(X)
  # Edge case when there is no individual.
  if (dimS == 0){
    S = matrix(rep(0, p), nrow = p)
  }
  # Otherwise, use the usual update.
  else{
    temp_init = P_call %*% X %*% (diag(p) - P_r)
    S = svd(temp_init)$v[, 1:dimS, drop = F]
  }
  return(S)
}

# Function for updating joint column structure M given two individuals.
updateM <- function(X1, X2, R1, R2, rc){
  n = nrow(X1)
  if (rc == 0){
    M = matrix(rep(0, n), nrow = n)
    return(M)
  }
  else{
    PR1 = tcrossprod(R1)
    PR2 = tcrossprod(R2)
    M1 = (diag(n) - PR1) %*% X1
    M2 = (diag(n) - PR2) %*% X2
    # Be careful of these edge cases since we need to calculate projection of (R1, R2):
    # If R1 is zero vector
    if (sum(R1^2) < 1e-8){
      # If R2 is zero vector
      if (sum(R2^2) < 1e-8){
        # we don't need to care about the orthogonality constraint of M and R1 or R2.
        PR = diag(n)
      }
      else{
        # R is only R2
        PR = diag(n) - projection(R2, ortho = TRUE)
      }
    }
    # If R1 is not zero vector
    else{
      # If R2 is simply zero vector
      if (sum(R2^2) < 1e-8){
        # We don't need to care about the orthogonality constraint of M and R2.
        PR = diag(n) - projection(R1, ortho = TRUE)
      }
      else{
        # Usual case
        R = cbind(R1, R2)
        PR = diag(n) - projection(R)
      } 
    }
    M = svd(PR %*% cbind(M1, M2))$u[, 1:rc, drop = F]
    return(M)
  }
}

# Function for updating joint row structure N given two individuals.
updateN <- function(X1, X2, S1, S2, rr){
  p = ncol(X1)
  if (rr == 0){
    N = matrix(rep(0, p), nrow = p)
    return(N)
  }
  else{
    PS1 = tcrossprod(S1)
    PS2 = tcrossprod(S2)
    N1 = X1 %*% (diag(p) - PS1) 
    N2 = X2 %*% (diag(p) - PS2) 
    # Be careful of these edge cases since we need to calculate projection of (S1, S2):
    # If S1 is zero vector
    if (sum(S1^2) < 1e-8){
      # If S2 is zero vector
      if (sum(S2^2) < 1e-8){
        # we don't need to care about the orthogonality constraint of N and S1 or S2.
        PS = diag(p)
      }
      else{
        # S is only S2
        PS = diag(p) - projection(S2, ortho = TRUE)
      }
    }
    # If S1 is not zero vector
    else{
      # If S2 is simply zero vector
      if (sum(S2^2) < 1e-8){
        # S is only S1
        PS = diag(p) - projection(S1, ortho = TRUE)
      }
      else{
        # Usual case
        S = cbind(S1, S2)
        PS = diag(p) - projection(S)
      } 
    }
    N = svd(rbind(N1, N2) %*% PS)$v[, 1:rr, drop = F]
    return(N)
  }
}

#' Main function of iterative DMMD algorithm
#' @details This function decomposed double-matched matrices according to Lemma 1 in Dongbang Yuan & Irina Gaynanova (2022) Double-Matched Matrix Decomposition for Multi-View Data, Journal of Computational and Graphical Statistics, DOI: 10.1080/10618600.2022.2067860. This function also updates the joint structures M and N compared to DMMD function.
#' @param X1 The first noisy data matrix
#' @param X2 The second noisy data matrix
#' @param r1 Total rank for X1. Default is NULL, meaning unknown, which will be estimated by rank estimation procedure determined by 'method'
#' @param r2 Total rank for X2. Default is NULL, meaning unknown, which will be estimated by rank estimation procedure determined by 'method'
#' @param rc Joint column rank. Default is NULL, meaning unknown, which will be estimated by profile likelihood method
#' @param rr Joint row rank. Default is NULL, meaning unknown, which will be estimated by profile likelihood method
#' @param angle_threshold The threshold angle for principal angles. Principal angles greater than the threshold are not considered as joint signal. Default is 90 degree
#' @param variance1 Either "equal" or "unequal". Default is "equal". This argument is the variance assumption used in the profile likelihood method for determining the total rank
#' @param variance2 Either "equal" or "unequal". Default is "equal". This argument is the variance assumption used in the profile likelihood method for determining the joint rank
#' @param method The method used for determining the total ranks r1 and r2. Either "PL" (profile likelihood) or "ED" (edge distribution). Default is "PL"
#' @param tol The tolerance used to determine convergence. Default is the square root of the machine precision 
#' @param maxiter Maximum number of iterations allowed in the iterative algorithm. Default is 1000
#' @param verbose Do you want to see the calculating progress of the function? Default is FALSE, which means the function stays silent.

#' @return A list with the following elements:
#' \item{r1}{The estimated total rank for X1, if not specified}
#' \item{r2}{The estimated total rank for X2, if not specified}
#' \item{rc}{The estimated joint column rank for X1, if not specified}
#' \item{rr}{The estimated joint row rank for X1, if not specified}
#' \item{A1}{The estimated low-rank signal matrix of X1}
#' \item{A2}{The estimated low-rank signal matrix of X2}
#' \item{E1}{The noise matrix of X1, \deqn{X_1 = A_1 + E_1}}
#' \item{E2}{The noise matrix of X2, \deqn{X_2 = A_2 + E_2}}
#' \item{Jc1}{The estimated low-rank joint column signal for X1}
#' \item{Jc2}{The estimated low-rank joint column signal for X2}
#' \item{Jr1}{The estimated low-rank joint row signal for X1}
#' \item{Jr2}{The estimated low-rank joint row signal for X2}
#' \item{Ic1}{The estimated low-rank individual column signal for X1}
#' \item{Ic2}{The estimated low-rank individual column signal for X2}
#' \item{Ir1}{The estimated low-rank individual row signal for X1}
#' \item{Ir2}{The estimated low-rank individual row signal for X2} 
#' \item{obj_vec}{A vector of objective values calculated at each iteration}
#' @export

#' @examples
#' data = DoubleDataGen(n = 20, p = 16, rank = c(4, 3), rc = 2, rr = 1, nrep = 1)
#' X1 = data$X1[[1]]
#' X2 = data$X2[[1]]
#' result = DMMD_i(X1,X2, verbose = TRUE)
DMMD_i <- function(X1, X2, r1 = NULL, r2 = NULL, rc = NULL, rr = NULL, angle_threshold = 90 * pi/180, variance1 = c("equal", "unequal"), variance2 = c("equal", "unequal"), method = c("PL", "ED"), tol = .Machine$double.eps^0.5, maxiter = 1000, verbose = FALSE){
  variance1 = match.arg(variance1)
  variance2 = match.arg(variance2)
  method = match.arg(method)
  # Check if the column names are equal
  if (!identical(colnames(X1), colnames(X2))){
    warning("This is an algorithm for double matched matrices. The column names of given matrices do not match")
  }
  # Check if the row names are equal
  if (!identical(rownames(X1), rownames(X2))){
    warning("This is an algorithm for double matched matrices. The row names of given matrices do not match")
  }
  n = nrow(X1)
  p = ncol(X1)
  # Check if the specified total ranks are legal
  if (!is.null(r1) | !is.null(r2)){
    if (max(r1,r2) > min(n, p) | min(r1, r2) <= 0){
      stop("The specified rank is not legal, please check.")
    }
  }
  # Check if the specified joint ranks are legal
  if (!is.null(rc) | !is.null(rr)){
    if (min(rc, rr) <= 0){
      stop("The specified joint rank is not legal, please check.")
    }
    if (!is.null(r1) | !is.null(r2)){
      if (max(rc, rr) > min(r1, r2)){
        stop("The specified joint rank is not legal, please check.")
      }
    }
  }
  # Save the svd result of the original matrices
  svd_x1 = svd(X1)
  svd_x2 = svd(X2)

  # Get the estimated total rank of X1 and X2. Store it as r1 and r2.
  if (is.null(r1)){
    if (method == "PL"){
      r1 = ProfileLikCluster(svd_x1$d, variance = variance1)$index
    }
    if (method == "ED"){
      r1 = Select_ED_Rank(svd_x1$d, maxiter = maxiter)
    }
  }
  if (is.null(r2)){
    if (method == "PL"){
      r2 = ProfileLikCluster(svd_x2$d, variance = variance1)$index
    }
    if (method == "ED"){
      r2 = Select_ED_Rank(svd_x2$d, maxiter = maxiter)
    }
  }
  # Check if the specified joint rank is legal
  if (!is.null(rc) | !is.null(rr)){
    if (max(rc,rr) > min(r1, r2)){
      stop("The specified joint rank is not legal, please check.")
    }
  }
  # Get the estimated column/row space of X1 and X2
  X1_est_c = as.matrix(svd_x1$u[,1:r1])
  X2_est_c = as.matrix(svd_x2$u[,1:r2])
  X1_est_r = as.matrix(svd_x1$v[,1:r1])
  X2_est_r = as.matrix(svd_x2$v[,1:r2])
  
  # Calculate joint column space
  # Get the principal angles
  angle_result_c = angle_cal(X1_est_c, X2_est_c)
  # Get the principal vectors
  pv1_c = angle_result_c$principal_vector1
  pv2_c = angle_result_c$principal_vector2
  # If the specified joint column rank is NULL. Calculate it using the PL method specified.
  if (is.null(rc)){
    principal_angle_c = angle_result_c$angle
    rc = joint_angle_cluster(principal_angle_c, angle_threshold = angle_threshold, variance = variance2)$joint_rank
  }
  
  # Calculate joint row space
  angle_result_r = angle_cal(X1_est_r, X2_est_r)
  # Get the principal vectors
  pv1_r = angle_result_r$principal_vector1
  pv2_r = angle_result_r$principal_vector2
  # If the specified joint row rank is NULL. Calculate it using the PL method specified.
  if (is.null(rr)){
    principal_angle_r = angle_result_r$angle
    rr = joint_angle_cluster(principal_angle_r, angle_threshold = angle_threshold, variance = variance2)$joint_rank
  }
  # Get initial estimates for M and N by averaging
  # Calculate joint column space projection matrix (this is MM')
  if (rc == 0){
    joint_space_c = matrix(rep(0, n), nrow = n)
    P_c = projection(joint_space_c, ortho = TRUE)
  }
  else{
    joint_space_c = (pv1_c[,1:rc] + pv2_c[,1:rc])/2
    P_c = projection(joint_space_c)
  }
  
  # Calculate joint row space projection matrix (this is NN')
  if (rr == 0){
    joint_space_r = matrix(rep(0, p), nrow = p)
    P_r = projection(joint_space_r, ortho = TRUE)
  }
  else{
    joint_space_r = (pv1_r[,1:rr] + pv2_r[,1:rr])/2
    P_r = projection(joint_space_r)
  }
  
  # Initialize R1, R2 and complete full column space
  R1 = updateR(X1, P_rall = diag(p), P_c = P_c, dimR = r1 - rc)
  P_call1 = P_c + projection(R1, ortho = TRUE)

  R2 = updateR(X2, P_rall = diag(p), P_c = P_c, dimR = r2 - rc)
  P_call2 = P_c + projection(R2, ortho = TRUE)

  # Initialize S1, S2 and complete full row space
  S1 = updateS(X1, P_call = P_call1, P_r = P_r, dimS = r1 - rr)
  P_rall1 = P_r + projection(S1, ortho = TRUE)
  
  S2 = updateS(X2, P_call = P_call2, P_r = P_r, dimS = r2 - rr)
  P_rall2 = P_r + projection(S2, ortho = TRUE)

  # Calculate estimated signals and objective value
  A1 = P_call1 %*% X1 %*% P_rall1
  A2 = P_call2 %*% X2 %*% P_rall2
  obj_old = sum((X1-A1)^2) + sum((X2-A2)^2)
  
  obj_vec = c(obj_old) 
  k = 0
  error = 1000
  # Main while loop
  while((error > tol)&(k < maxiter)){
    # Count iterations
    k = k + 1
    if (verbose){
      print(paste('Iteration', k, sep = ' '))
      print(obj_old)
      print("Update M")
    }
    
    ### Update joint M ###
    M = updateM(X1 = X1 %*% P_rall1, X2 = X2 %*% P_rall2, R1 = R1, R2 = R2, rc = rc)
    P_c = projection(M, ortho = TRUE)
    P_call1 = P_c + projection(R1, ortho = TRUE)
    P_call2 = P_c + projection(R2, ortho = TRUE)
    
    A1 = P_call1 %*% X1 %*% P_rall1
    A2 = P_call2 %*% X2 %*% P_rall2
    obj_new = sum((X1-A1)^2) + sum((X2-A2)^2)
    if (verbose){
      print(obj_new)
      print("Update N")
    }
    
    ### Update joint N ###
    N = updateN(X1 = P_call1 %*% X1, X2 = P_call2 %*% X2, S1 = S1, S2 = S2, rr = rr)
    P_r = projection(N, ortho = TRUE)
    P_rall1 = P_r + projection(S1, ortho = TRUE)
    P_rall2 = P_r + projection(S2, ortho = TRUE)

    A1 = P_call1 %*% X1 %*% P_rall1
    A2 = P_call2 %*% X2 %*% P_rall2
    obj_new = sum((X1-A1)^2) + sum((X2-A2)^2)
    
    if (verbose){
      print(obj_new)
      print("Update R1, R2")
    }
    
    # update R and full column space
    R1 = updateR(X1, P_rall = P_rall1, P_c = P_c, dimR = r1 - rc)
    P_call1 = P_c + projection(R1, ortho = TRUE)
    R2 = updateR(X2, P_rall = P_rall2, P_c = P_c, dimR = r2 - rc)
    P_call2 = P_c + projection(R2, ortho = TRUE)
    A1 = P_call1 %*% X1 %*% P_rall1
    A2 = P_call2 %*% X2 %*% P_rall2
    obj_new = sum((X1-A1)^2) + sum((X2-A2)^2)
    if (verbose){
      print(obj_new)
      print("Update S1, S2")
    }
    
    # update S and full row space
    S1 = updateS(X1, P_call = P_call1, P_r = P_r, dimS = r1 - rr)
    P_rall1 = P_r + projection(S1, ortho = TRUE)
    S2 = updateS(X2, P_call = P_call2, P_r = P_r, dimS = r2 - rr)
    P_rall2 = P_r + projection(S2, ortho = TRUE)

    # calculate objective function and prepare for next iteration
    A1 = P_call1 %*% X1 %*% P_rall1
    A2 = P_call2 %*% X2 %*% P_rall2
    obj_new = sum((X1-A1)^2) + sum((X2-A2)^2)
    if (verbose){
      print(obj_new)
    }
    error = abs(obj_new - obj_old)

    obj_old = obj_new
    obj_vec = c(obj_vec,obj_new)
  }
  
  # Get the decomposition and return results
  J1_c = P_c %*% A1
  J2_c = P_c %*% A2
  J1_r = A1 %*% P_r
  J2_r = A2 %*% P_r
  I1_c = A1 - J1_c
  I2_c = A2 - J2_c
  I1_r = A1 - J1_r
  I2_r = A2 - J2_r 
  E1 = X1 - A1
  E2 = X2 - A2

  return(list(r1 = r1, r2 = r2, rc = rc, rr = rr,
              A1 = A1, A2 = A2, E1 = E1, E2 = E2, 
              Jc1 = J1_c, Jc2 = J2_c, Jr1 = J1_r, Jr2 = J2_r,
              Ic1 = I1_c, Ic2 = I2_c, Ir1 = I1_r, Ir2 = I2_r, obj_vec = obj_vec))
  
}