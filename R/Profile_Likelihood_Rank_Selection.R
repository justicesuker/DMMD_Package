#' Calculate the MLE of variance of two vectors.
#' @importFrom stats var
#' @param x The first vector
#' @param y The second vector
#' @param variance Either "equal" or "unequal", i.e., whether the assumption is equal variance or unequal variance. If no input is given, equal variance is assumed.
#' 
#' @return A vector of length 2 that contains the estimated variance for the two vectors
VarianceMLE <- function(x, y, variance = c("equal", "unequal")){
  variance = match.arg(variance)
  m = length(x)
  n = length(y)
  a = var(x)
  b = var(y)
  if (any(is.na(x)) || any(is.na(y))) {
    stop("There is NA in input vectors.")
  }
  else{
    if (variance == "equal"){
      if (m + n < 3){stop("Input vector is of length 2 or less.")}
      if (m == 1){
        result = (n - 1) * b / (1 + n)
        resultvec = c(result,result)
        return(resultvec)
      }
      if (n == 1){
        result = (m - 1) * a / (m + 1)
        resultvec = c(result,result)
        return(resultvec)
      }
      if (m != 1 && n != 1){
        result = ((m - 1) * a + (n - 1) * b) / (m + n)
        resultvec = c(result,result)
        return(resultvec)
      }
    }
    if (variance == "unequal"){
      if (m < 2 || n < 2){
        stop("Length of vector must be at least 2 for unequal variance assumption.")
      }
      result1 = (m - 1) * a / m
      result2 = (n - 1) * b / n
      resultvec = c(result1,result2)
      return(resultvec)
    }
  }
}

#' Calculate the profile log-likelihood of given two data vectors assuming normal distributions.
#' @importFrom stats dnorm
#' @param x The first data vector
#' @param y The second data vector
#' @param variance Either "equal" or "unequal", i.e., whether the assumption is equal variance or unequal variance. If no input is given, equal variance is assumed.
#'
#' @return The calculated profile log-likelihood
ProfileLoglik <- function(x, y, variance = c("equal", "unequal")){
  variance = match.arg(variance)
  mu1 = mean(x)
  mu2 = mean(y)
  sigmaMLE = sqrt(VarianceMLE(x,y,variance))
  part1 = sum(log(dnorm(x, mean = mu1, sd = sigmaMLE[1])))
  part2 = sum(log(dnorm(y, mean = mu2, sd = sigmaMLE[2])))
  result = part1 + part2
  return(result)
}

#' Preliminary function that separates a vector into two groups using profile likelihood assuming equal variance.
#'
#' @param x An ordered (max to min) vector to be clustered into two groups
#'
#' @return A list that contains the following elements:
#' \item{index}{The index that separates the vector \code{x}. The element at the index is considered in the first cluster}
#' \item{profileloglikvec}{A vector containing profile log-likelihood for each index}
ProfileLikClusterEqual <- function(x){
  l = length(x)
  if (l == 1){
    warning("Length of vector is 1. There is only one cluster.")
    return(list(index = 1, profileloglikvec = NA))
  }
  else if (l == 2){
    warning("Length of vector is 2. Each cluster contains one point without running the algorithm.")
    return(list(index = 1, profileloglikvec = NA))
  }
  # Normal case
  else{
    maxindex = l - 1
    # Initialize the profile likelihood vector
    profileloglikvec = rep(NA,maxindex)
    # Consider each index that is less than or equal to maxindex
    for (p in 1:maxindex){
      # Separate the vector into two using the index p
      data1 = x[1:p]
      data2 = x[(p+1):l]
      # Calculate the corresponding profile log-likelihood
      profileloglikvec[p] = ProfileLoglik(data1,data2,"equal")
    }
    # Find the maximizer which will be the output
    index = which.max(profileloglikvec)
    return(list(index = index, profileloglikvec = profileloglikvec))
  }
}

#' Function that separates a vector into two groups using profile likelihood
#'
#' @param x An ordered (max to min) vector to be clustered into two groups
#' @param variance Either "equal" or "unequal", i.e., whether the assumption is equal variance or unequal variance. If no input is given, equal variance is assumed.
#' @return A list that contains the following elements:
#' \item{index}{The index where separates the vector \code{x}. The element at the index is considered in the first cluster}
#' \item{profileloglikvec}{A vector containing profile log-likelihood for each index}

#' @export
#' @examples
#' x = c(20,9.5,9,8,6,5,4.5,3)
#' ProfileLikCluster(x)
ProfileLikCluster <- function(x, variance = c("equal", "unequal")){
  variance = match.arg(variance)
  l = length(x)
  # if variance = "equal", simply call the basic ProfileLikClusterEqual function 
  if (!variance %in% c('equal', 'unequal')){
    stop('The input of variance must be either equal or unequal.')
  }
  if (variance == "equal"){
    return(ProfileLikClusterEqual(x))
  }
  # Assuming unequal variance
  else{
    # Extreme cases
    if (l == 1){
      warning("Length of vector is 1. There is only one cluster.")
      return(list(index = 1, profileloglikvec = NA))
    }
    else if (l == 2){
      warning("Length of vector is 2. Each cluster contains one point without running the algorithm.")
      return(list(index = 1, profileloglikvec = NA))
    }
    else if (l == 3){
      warning("Length of vector is 3. The result is based on equal variance assumption.")
      return(ProfileLikClusterEqual(x))
    }
    # Normal cases
    else{
      maxindex = l - 1
      # Initialize the profile likelihood vector
      profileloglikvec = rep(NA,maxindex)
      for (p in 2:(maxindex - 1)){
        # Separate the vector into two using the index p
        data1 = x[1:p]
        data2 = x[(p+1):l]
        # Calculate the corresponding profile log-likelihood
        profileloglikvec[p] = ProfileLoglik(data1, data2, variance = "unequal")
      }
      # Calculate the profile log-likelihood assuming equal variance when index = 1. Because for only one point we cannot assume unequal variance.
      profileloglikvec[1] = ProfileLoglik(x[1],x[2:l], variance = "equal")
      # Calculate the profile log-likelihood assuming equal variance when index = l - 1. Because for only one point we cannot assume unequal variance.
      profileloglikvec[maxindex] = ProfileLoglik(x[1:(maxindex)],x[(maxindex + 1):l],variance = "equal")
      # Combine all the profile log-likelihood and find the maximizer which will be the output
      index = which.max(profileloglikvec)
    }
    if (index == 1 || index == l - 1){
      warning("Profile likelihood assuming unequal variance is not greater than likelihood assuming equal variance. Printed result is based on equal variance assumption.")
      return(ProfileLikClusterEqual(x))
    }
    else{return(list(index = index, profileloglikvec = profileloglikvec))}
  }
}