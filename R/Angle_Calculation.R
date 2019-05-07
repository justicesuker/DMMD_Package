require(matlib)

#' Function that calculates principal angles between column spaces of two matrices
#'
#' @param X The first matrix.
#' @param Y The second matrix.
#' @param tol Tolerence, default is the square root of machine precision.
#'
#' @return A list that contains the following:
#' \item{angle}{A vector of principal angles with increasing order.}
#' \item{cos_angle}{A vector of cosine principal angles}
#' \item{principal_vector1}{Principal vectors of matrix \code{X}}
#' \item{principal_vector2}{Principal vectors of matrix \code{Y}}
#'
#' @examples
#' X = matrix(c(1,1,1,1,1,0),nrow = 3, ncol = 2)
#' Y = matrix(c(1,1,1,2,1,0),nrow = 3, ncol = 2)
#' angle_cal(X,Y)
#' 
angle_cal <- function(X, Y, tol = .Machine$double.eps^0.5){
  X = as.matrix(X)
  Y = as.matrix(Y)
  svd_x = svd(X)
  svd_y = svd(Y)
  rank1 = sum(svd_x$d > tol)
  rank2 = sum(svd_y$d > tol)
  Xnorm = svd_x$u[,1:rank1]
  Ynorm = svd_y$u[,1:rank2]
  M = crossprod(Xnorm,Ynorm)
  svd_result = svd(M)
  cos_angle = svd_result$d
  l = length(cos_angle)
  principal_angle = rep(NA, l)
  for (i in 1:l){
    if (cos_angle[i] >= 1){principal_angle[i] = 0}
    if (cos_angle[i] <= 0){principal_angle[i] = pi/2}
    if (cos_angle[i] > 0 && cos_angle[i] < 1){
      principal_angle[i] = acos(cos_angle[i])
    }
  }
  principal_mat1 = Xnorm %*% svd_result$u
  principal_mat2 = Ynorm %*% svd_result$v
  return(list("angle" = principal_angle, "cos_angle" = cos_angle,
              "principal_vector1" = principal_mat1, 
              "principal_vector2" = principal_mat2))
}

#' Function that estimates joint rank by method of profile likelihood.
#' @param angle_vec A vector of principal angles
#' @param angle_threshold Optional Threshold (radians) that is used to truncate principal angles. Default is 90 degrees, i.e., no threshold is used.
#' @param variance Either "equal" or "unequal", i.e whether the assumption is equal variance or unequal variance. Default is unequal.
#' @param cos_method TRUE or FALSE, i.e., whether use angles or cosine angles as the data to cluster. Default is FALSE.
#'
#' @return A list with the following elements:
#' \item{joint_rank}{Estimated joint rank.}
#' \item{profileloglikvec}{Profile log likelihood calculated at each index.
#' The function will return NA with a warning message, if less or equal to 2 principal angles are smaller than the threshold.}
#'
#' @examples
#' 
#' data = datagen(n = 100, p = 20, joint_rank = 5, individual_rank = c(5,7), nrep = 1)
#' X1 = data$X1_list[[1]]
#' X2 = data$X2_list[[1]]
#' angle_vec = angle_cal(X1,X2)$angle
#' joint_angle_cluster(angle_vec)
#' 
joint_angle_cluster <- function(angle_vec, angle_threshold = 90 * pi/180, variance = "unequal", cos_method = FALSE){
  angle_vec = angle_vec[angle_vec < angle_threshold]
  l = length(angle_vec)
  if (l <= 2){
    if (l == 0){
      stop("There is no angle that is smaller than threshold, the joint rank is 0.")
    }
    if (l == 1){
      warning("There is only one angle that is smaller than threshold, the joint rank is recommended as 1 without using profile likelihood method.")
      return(list("joint_rank" = 1, "profileloglikvec" = NA))
    }
    if (l == 2){
      warning("There are only two angles that are smaller than threshold, the joint rank is recommended to be 2 without using profile likelihood method.")
      return(list("joint_rank" = 2, "profileloglikvec" = NA))
    }
  }
  else{
    if (var(angle_vec) < .Machine$double.eps^0.5){
      warning("The variance of angles that are smaller than threshold is too small. It is treated that all of them are joint signals. Please check manually.")
      return(list("joint_rank" = l, "profileloglikvec" = NA))
    }
    else{
      # Add 0 and angle_threshold artificially to include 0 joint rank case and full joint rank case
      angle_vec = c(0,angle_vec,angle_threshold)
      if (cos_method){angle_vec = cos(angle_vec)}
      result = ProfileLikCluster(angle_vec, variance)
      index = result$index - 1
      profileloglikvec = result$profileloglikvec
      profileloglikvec = profileloglikvec[-c(1,(l+2))]
    }
    return(list("joint_rank" = index, "profileloglikvec" = profileloglikvec))
  }
}
