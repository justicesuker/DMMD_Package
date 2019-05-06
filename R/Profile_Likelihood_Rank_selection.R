#' Function that calculates MLE of variance of two vectors.
#' 
#' @param x The first vector
#' @param y The second vector
#' @param variance Either "equal" or "unequal", 
#' corresponding to equal or unequal variance assumption.
#' Default is "equal".
#'
#' @return A vector of length two. If assumes equal variance, the result is the MLE version of pooled standard deviation repeating twice. If assumes unequal variance, the result is two MLE version of sample standard deviation.
#'
#' @examples
#' x = c(1,2,3)
#' y = c(1,3,4)
#' VarianceMLE(x,y,variance = "equal")
#' VarianceMLE(x,y,variance = "unequal")
#' 
VarianceMLE <- function(x, y, variance = "equal"){
  m = length(x)
  n = length(y)
  if (variance == "equal"){
    if (m + n < 3){stop("Input vector is of length 2 or less")}
    if (m == 1){
      result = (n - 1) * var(y) / (1 + n)
      resultvec = c(result,result)
      return(resultvec)
    }
    if (n == 1){
      result = (m - 1) * var(x) / (m + 1)
      resultvec = c(result,result)
      return(resultvec)
    }
    if (m != 1 && n != 1){
      result = ((m - 1) * var(x) + (n - 1) * var(y)) / (m + n)
      resultvec = c(result,result)
      return(resultvec)
    }
  }
  if (variance == "unequal"){
    if (m < 2 || n < 2){
      stop("Input vector is too short for unequal variance calculation")
    }
    result1 = (m - 1) * var(x) / m
    result2 = (n - 1) * var(y) / n
    resultvec = c(result1,result2)
    return(resultvec)
  }
}

ProfileLoglik <- function(x,y,variance = "equal"){
  mu1 = mean(x)
  mu2 = mean(y)
  sigmaMLE = sqrt(VarianceMLE(x,y,variance))
  if (sum(sigmaMLE == 0) > 0){
    stop("There is no variance in data, please try another method")
  }
  part1 = sum(log(dnorm(x, mean = mu1, sd = sigmaMLE[1])))
  part2 = sum(log(dnorm(y, mean = mu2, sd = sigmaMLE[2])))
  result = part1 + part2
  return(result)
}

# Input : An ordered (max to min) vector that is going to be clustered into two groups
# Basic function that assumes equal variance.
ProfileLikClusterEqual <- function(x){
  l = length(x)
  if (l < 3){
    stop("Length of vector must be at least 3")
  }
  profileloglikvec = rep(NA,l-1)
  for (p in 1:(l-1)){
    data1 = x[1:p]
    data2 = x[(p+1):l]
    profileloglikvec[p] = ProfileLoglik(data1,data2,"equal")
  }
  index = which.max(profileloglikvec)
  return(list(index = index, profileloglikvec = profileloglikvec))
}

#' Function that separates a vector into two groups using profile likelihood
#'
#' @param x Vector to be separated
#' @param variance Either "equal" or "unequal", 
#' corresponding to equal or unequal variance assumption.
#' Default is "equal".
#'
#' @return A list that contains the following elements:
#' \item{index}{The index where separates the vector \code{x}}
#' \item{profileloglikvec}{A vector that contains profile log-likelihood for each index.}
#'
#' @examples
#' x = c(20,9.5,9,8,6,5,4.5,3)
#' ProfileLikCluster(x)

ProfileLikCluster <- function(x, variance = "equal"){
  l = length(x)
  if (variance == "equal"){
    return(ProfileLikClusterEqual(x))
  }
  if (variance == "unequal"){
    if (l < 4){
      stop("Sample size must be at least 4 to use unequal variance assumption")
    }
    profileloglikvec = rep(NA,l-1)
    for (p in 2:(l-2)){
      data1 = x[1:p]
      data2 = x[(p+1):l]
      profileloglikvec[p] = ProfileLoglik(data1,data2,variance = "unequal")
    }
    profileloglikvec[1] = ProfileLoglik(x[1],x[2:l],variance = "equal")
    profileloglikvec[l-1] = ProfileLoglik(x[1:(l-1)],x[l],variance = "equal")
    index = which.max(profileloglikvec)
    if (index == 1 || index == l){
      warning("Profile likelihood assuming unequal variance is not greater than likelihood assuming equal variance. Printed result is based on equal variance assumption.")
      return(ProfileLikClusterEqual(x))
    }
    else{return(list(index = index, profileloglikvec = profileloglikvec))}
  }
}

