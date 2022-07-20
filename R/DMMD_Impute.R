# Function to impute missing value with column mean
NaiveImpute <- function(data){
  for(i in 1:ncol(data)){
    temp_na_row = is.na(data[,i])
    data[temp_na_row, i] <- mean(data[,i], na.rm = TRUE)
  }
  return(data)
}

# Function to impute missing value with specific imputation method, that is "mean","unif" or "normal". 
RandomImpute <- function(data, method = c('mean', 'unif', 'normal')){
  if (method == 'mean'){
    return(NaiveImpute(data))
  }
  # Impute with uniform distributed random variables in the same range
  if (method == 'unif'){
    for(i in 1:ncol(data)){
      temp_na_row = is.na(data[,i])
      temp_min = min(data[,i], na.rm = TRUE)
      temp_max = max(data[,i], na.rm = TRUE)
      # The upperbound of the uniform distribution is the maximum value. The lower bound is the minimum value observed
      data[temp_na_row, i] <- runif(sum(temp_na_row), min = temp_min, max = temp_max)
    }
    return(data)
  }
  
  # Impute with normal distributed random variables with mean and variance to be the same as observed data
  if (method == 'normal'){
    for(i in 1:ncol(data)){
      temp_na_row = is.na(data[,i])
      temp_mean = mean(data[,i], na.rm = TRUE)
      temp_sd = sd(data[,i], na.rm = TRUE)
      # The mean of the normal distribution is the sample mean. The standard deviation is the sample standard deviation.
      data[temp_na_row, i] <- rnorm(sum(temp_na_row), mean = temp_mean, sd = temp_sd)
    }
    return(data)
  }
}

# A function to impute missing value with known matrix entries
# 'data' is the matrix with NAs. 'newdata' is the complete matrix
ImputeKnown <- function(data, newdata){
  data_imputed = data
  if (!identical(dim(data),dim(newdata))){
    stop("The dimension does not match.")
  }
  for(i in 1:ncol(data_imputed)){
    temp_na_row = is.na(data_imputed[,i])
    data_imputed[temp_na_row, i] <- newdata[temp_na_row, i]
  }
  return(data_imputed)
}

#' An iterative algorithm to do imputation based on iterative DMMD method.
#'
#' @param X1 The first data matrix with NAs
#' @param X2 The second data matrix with NAs
#' @param init The imputation method for initializing X1 and X2. Could be 'mean', 'unif' or 'normal'. Default 'mean'
#' @param r1 Total rank for X1. Default is NULL, meaning unknown, which will be estimated by rank estimation procedure determined by 'method'. Passed to DMMD_i function
#' @param r2 Total rank for X2. Default is NULL, meaning unknown, which will be estimated by rank estimation procedure determined by 'method'. Passed to DMMD_i function
#' @param rc Joint column rank. Default is NULL, meaning unknown, which will be estimated by profile likelihood method. Passed to DMMD_i function
#' @param rr Joint row rank. Default is NULL, meaning unknown, which will be estimated by profile likelihood method. Passed to DMMD_i function
#' @param angle_threshold The threshold angle for principal angles. Principal angles greater than the threshold are not considered as joint signal. Default is 90 degree. Passed to DMMD_i function
#' @param variance1 Either "equal" or "unequal". Default is "equal". This argument is the variance assumption used in the profile likelihood method for determining the total rank. Passed to DMMD_i function
#' @param variance2 Either "equal" or "unequal". Default is "equal". This argument is the variance assumption used in the profile likelihood method for determining the joint rank. Passed to DMMD_i function
#' @param method The method used for determining the total ranks r1 and r2. Either "PL" (profile likelihood) or "ED" (edge distribution). Default is "PL". Passed to DMMD_i function
#' @param tol The tolerance used to determine convergence. Default is the square root of the machine precision. Passed to DMMD_i function
#' @param maxiter1 Maximum number of iterations allowed in the iterative imputation. Default is 1000
#' @param maxiter2 Maximum number of iterations allowed in each called DMMD_i function. Default is 1000
#' @return A list with the following elements:
#' \item{Imputed_X1}{The final imputed and complete X1}
#' \item{Imputed_X2}{The final imputed and complete X2}
#' @export
#'
#' @examples
#' data = DoubleDataGen(n = 20, p = 16, rank = c(4, 3), rc = 2, rr = 1, nrep = 1)
#' X1_missing = data$X1[[1]]
#' X2_missing = data$X2[[1]]
#' X1_missing[1,1] = NA
#' X2_missing[3,2] = NA
#' result = DMMD_Impute(X1_missing, X2_missing)
#' result$Imputed_X1 - data$X1[[1]]
DMMD_Impute <- function(X1, X2, init = c('mean', 'unif', 'normal'), r1 = NULL, r2 = NULL, rc = NULL, rr = NULL, angle_threshold = 90 * pi/180, variance1 = c("equal", "unequal"), variance2 = c("equal", "unequal"), method = c("PL", "ED"), tol = .Machine$double.eps^0.5, maxiter1 = 1e3, maxiter2 = 1e3){
  init = match.arg(init)
  variance1 = match.arg(variance1)
  variance2 = match.arg(variance2)
  method = match.arg(method)
  if (!identical(colnames(X1), colnames(X2))){
    warning("This is an algorithm for double matched matrices. The column names of given matrices do not match")
  }
  if (!identical(rownames(X1), rownames(X2))){
    warning("This is an algorithm for double matched matrices. The row names of given matrices do not match")
  }
  n = dim(X1)[1]
  p = dim(X1)[2]
  # Initialize with the method specified 
  X1_old = RandomImpute(X1, method = init)
  X2_old = RandomImpute(X2, method = init)
  err = +Inf
  iter = 0
  while (iter <= maxiter1 & err > tol){
    # Use DMMD algorithm to get the signal matrix
    temp_result = DMMD_i(X1 = X1_old, X2 = X2_old, r1 = r1, r2 = r2, rc = rc, rr = rr, 
                         angle_threshold = angle_threshold, variance1 = variance1, variance2 = variance2, method = method, tol = tol, maxiter = maxiter2)
    signal1_new = temp_result$A1
    signal2_new = temp_result$A2
    # Impute the missing values of old matrix with the entries of signal matrices calculated 
    X1_new = ImputeKnown(X1, signal1_new)
    X2_new = ImputeKnown(X2, signal2_new)
    # record the error and see if it converges
    err = max(sum(X1_new - X1_old)^2, sum(X2_new - X2_old)^2)
    # Plus one iteration
    iter = iter + 1
    X1_old = X1_new
    X2_old = X2_new
  }
  if (iter > maxiter1){
    warning('DMMD impute fails to converge.')
  }
  return(list("Imputed_X1" = X1_new, "Imputed_X2" = X2_new))
}