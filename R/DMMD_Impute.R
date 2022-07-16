# A simple function to impute missing value with column mean
NaiveImpute <- function(data){
  for(i in 1:ncol(data)){
    temp_na_row = is.na(data[,i])
    data[temp_na_row, i] <- mean(data[,i], na.rm = TRUE)
  }
  return(data)
}

# A function to impute missing value with "mean","unif" or "normal" 
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
# data is the matrix with NAs. newdata is the complete matrix
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

# An iterative algorithm to do imputation based on DMMD method
# init is "mean", "unif" or "normal"
DMMD_Impute <- function(X1, X2, init = c('mean', 'unif', 'normal'), r1 = NULL, r2 = NULL, joint_rank_c = NULL, joint_rank_r = NULL, angle_threshold = 90 * pi/180, variance1 = "equal", variance2 = "equal", throw = FALSE, method = method, tol = .Machine$double.eps^0.5, maxiter1 = 1e3, maxiter2 = 1e3){
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
    temp_result = DMMD_i(X1 = X1_old, X2 = X2_old, r1 = r1, r2 = r2, rc = joint_rank_c, rr = joint_rank_r, 
                         angle_threshold = angle_threshold, variance1 = variance1, variance2 = variance2, throw = throw, method = method, eps = tol, kmax = maxiter2)
    signal1_new = X1_old - temp_result$Error$Error1
    signal2_new = X2_old - temp_result$Error$Error2
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
