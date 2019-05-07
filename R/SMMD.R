#' Single Matched Matrix Decomposition
#' @description Function that does single-matched matrix factorization. The matching is in the number of rows. 
#' @details The decomposition is based on the theory in the paper of
#' Feng, Qing, et al. "Angle-based joint and individual variation explained." Journal of multivariate analysis 166 (2018): 241-265.
#' One can also refer to 
#' Lock, Eric F., et al. "Joint and individual variation explained (JIVE) for integrated analysis of multiple data types." The annals of applied statistics 7.1 (2013): 523.
#' for details.
#' @param X1 The first matrix
#' @param X2 The second matrix
#' @param angle_threshold Optional threshold (radians) feeded to \code{joint_angle_cluster} function. Default is 90 degrees, i.e., no threshold is used.
#' @param variance1 Either "equal" or "unequal". The assumption of equal or unequal variance used to estimate the total rank.
#' @param variance2 Either "equal" or "unequal". The assumption of equal or unequal variance used to estimate the joint rank.
#' @param tol Tolerence, default is the square root of machine precision.
#' 
#' @return A list that contains the following:
#' \item{r1}{The estimated total rank of \code{X1}}
#' \item{r2}{The estimated total rank of \code{X2}}
#' \item{joint_rank}{The estimated rank of joint signals between \code{X1} and \code{X2}}
#' \item{individual_rank}{A vector of length 2. The vector contains the estimated ranks of individual signals of \code{X1} and \code{X2}}
#' \item{J1}{The joint signal of \code{X1}}
#' \item{J2}{The joint signal of \code{X2}}
#' \item{I1}{The individual signal of \code{X1}}
#' \item{I2}{The individual signal of \code{X2}}
#' 
#' @examples
#' data = datagen(n = 100, p = 20, joint_rank = 5, individual_rank = c(5,7), nrep = 1)
#' X1 = data$X1_list[[1]]
#' X2 = data$X2_list[[1]]
#' result_SMMD = SMMD(X1,X2)
#' result_SMMD$joint_rank
#' result_SMMD$individual_rank
#' 
SMMD <- function(X1, X2, angle_threshold = 90 * pi/180, variance1 = "unequal", variance2 = "unequal", tol = .Machine$double.eps^0.5){
  svd_x1 = svd(X1)
  svd_x2 = svd(X2)
  r1 = ProfileLikCluster(svd_x1$d, variance = variance1)$index
  r2 = ProfileLikCluster(svd_x2$d, variance = variance1)$index
  angle_result = angle_cal(X1, X2, tol = tol)
  principal_angle = angle_result$angle
  pv1 = angle_result$principal_vector1
  pv2 = angle_result$principal_vector2
  joint_rank = joint_angle_cluster(
    principal_angle, angle_threshold = angle_threshold, variance = variance2)$joint_rank
  joint_space = (pv1[,1:joint_rank] + pv2[,1:joint_rank])/2
  J1 = projection(joint_space) %*% X1
  J2 = projection(joint_space) %*% X2
  X1_per =  X1 - projection(joint_space) %*% X1
  X2_per =  X2 - projection(joint_space) %*% X2
  svd_x1_per = svd(X1_per)
  svd_x2_per = svd(X2_per)
  individual_rank = c(r1 - joint_rank, r2 - joint_rank)
  individual_space1 = svd_x1_per$u[,1:individual_rank[1]]
  individual_space2 = svd_x2_per$u[,1:individual_rank[2]]
  I1 = projection(individual_space1) %*% X1
  I2 = projection(individual_space2) %*% X2 
  return(list("r1" = r1, "r2" = r2, "joint_rank" = joint_rank, "individual_rank" = individual_rank, "J1" = J1, 
              "J2" = J2, "I1" = I1, "I2" = I2))
}